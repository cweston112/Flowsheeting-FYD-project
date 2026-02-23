from __future__ import annotations

from typing import Dict, List, Optional, Tuple
import numpy as np

from ..flowsheet_tools import UnitOp, EPS


class MultiStageCentrifugalContactorTBP(UnitOp):
    """
    Multi-stage centrifugal contactor (TBP extraction), counter-current stage model.

    D is defined as (org/aq). This model treats partitioning as TRUE / BIDIRECTIONAL:
      - If a species enters in the aqueous feed, some transfers to organic.
      - If a species enters in the organic feed, some transfers to aqueous.

    Stage-count performance uses constant-flow counter-current "xN/x0" with extraction factor E:
      - For aqueous -> organic extraction:   E_ex = D * (O/A) = D * OA
        aqueous-out / aqueous-in = xN_over_x0(E_ex, N)

      - For organic -> aqueous stripping:    E_st = (A/O) / D = AO / D = (1/OA) / D
        organic-out / organic-in = xN_over_x0(E_st, N)

    DESIGN UPDATE:
      We choose (N, OA_design) that minimizes a combined objective:
          J = w_N * N + w_OA * OA
      subject to:
        (1) required_recovery constraints are met
        (2) if tbp_stoich provided, enough TBP is present in organic to support extracted metals
            (including optional tbp_free_min).

    Notes:
      - TBP limiting is handled by sizing OA (not by post-hoc composition hacking).
      - This unit does not internally scale inlet flows during apply(); the flowsheet should size
        the organic feed upstream (as your flowsheet already does).
    """

    def __init__(
        self,
        name: str,
        D: Dict[str, float],
        required_recovery: Dict[str, float],
        nonextract_to_aq: Optional[List[str]] = None,
        nonextract_to_org: Optional[List[str]] = None,
        N_max: int = 50,
        OA_max: float = 1e3,
        OA_min: float = 0.0,
        N_balance_cap: int = 20,
        # --- TBP stoichiometry (optional) ---
        tbp_species: str = "TBP",
        tbp_stoich: Optional[Dict[str, float]] = None,  # solute -> mol TBP / mol solute in organic
        tbp_free_min: float = 2.0,  # optional minimum "free TBP" required in organic out
        # --- NEW: combined objective weights ---
        w_N: float = 1.0,
        w_OA: float = 1.0,
    ):
        super().__init__(name)
        self.D = {k: float(v) for k, v in D.items()}
        self.required_recovery = {k: float(v) for k, v in required_recovery.items()}
        self.nonextract_to_aq = set(nonextract_to_aq or [])
        self.nonextract_to_org = set(nonextract_to_org or [])
        self.N_max = int(N_max)
        self.OA_max = float(OA_max)
        self.OA_min = float(OA_min)
        self.N_balance_cap = int(N_balance_cap)

        self.tbp_species = str(tbp_species)
        self.tbp_stoich = {k: float(v) for k, v in (tbp_stoich or {}).items()}
        self.tbp_free_min = float(tbp_free_min)

        # objective weights
        self.w_N = float(w_N)
        self.w_OA = float(w_OA)

        # reported fields
        self.N_used: int = 1
        self.OA_design: float = 0.0
        self.OA_used: float = 0.0          # as provided by org_in
        self.OA_effective: float = 0.0     # equals min(OA_design, OA_used)

        # additional reporting
        self.org_scale_factor: float = 1.0
        self.O_tot_required: float = 0.0

        self.tbp_required: float = 0.0
        self.tbp_available: float = 0.0
        self.tbp_limited: bool = False

    @staticmethod
    def _xN_over_x0(E: float, N: int) -> float:
        if N <= 0:
            return 1.0
        if E <= 0.0:
            return 1.0
        if abs(E - 1.0) < 1e-12:
            return 1.0 / (N + 1.0)
        s = (E ** (N + 1) - 1.0) / (E - 1.0)
        return 1.0 / s

    @classmethod
    def _recovery_aq_to_org(cls, D: float, OA: float, N: int) -> float:
        return 1.0 - cls._xN_over_x0(D * OA, N)

    @classmethod
    def _required_OA_for_recovery(cls, D: float, N: int, recovery: float, OA_max: float) -> float:
        lo, hi = 0.0, 1.0
        while cls._recovery_aq_to_org(D, hi, N) < recovery:
            hi *= 2.0
            if hi > OA_max:
                raise ValueError("Required O/A exceeds OA_max.")
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if cls._recovery_aq_to_org(D, mid, N) >= recovery:
                hi = mid
            else:
                lo = mid
        return hi

    @staticmethod
    def _as_fraction(v: float) -> float:
        """
        Accept either:
          - fraction in [0,1]
          - percent in (1,100]
        """
        v = float(v)
        return v / 100.0 if v > 1.0 else v

    def _simulate_split_for_OA(
        self,
        reg,
        idx: Dict[str, int],
        A_base: np.ndarray,
        O_base: np.ndarray,
        OA: float,
        N: int,
        A_tot: float,
        O_tot_base: float,
    ) -> Tuple[np.ndarray, np.ndarray, float]:
        """
        Simulate the cascade at a target OA by scaling the organic inlet flow (and composition)
        proportionally. Returns (Aq, Org, org_scale_factor).
        """
        O_tot_req = OA * A_tot
        s = O_tot_req / max(O_tot_base, EPS)
        O = s * O_base
        A = A_base  # aqueous fixed

        AO = (1.0 / OA) if OA > EPS else float("inf")

        Aq = np.zeros(reg.n(), dtype=float)
        Org = np.zeros(reg.n(), dtype=float)

        for sp in reg.species:
            j = idx[sp]
            nA_in = float(A[j])
            nO_in = float(O[j])
            nT = nA_in + nO_in
            if nT <= EPS:
                continue

            # Forced routing first
            if sp in self.nonextract_to_aq:
                Aq[j] = nT
                continue
            if sp in self.nonextract_to_org:
                Org[j] = nT
                continue

            if sp in self.D:
                Di = float(self.D[sp])
                if Di <= 0.0:
                    Aq[j] = nT
                    continue

                # aqueous-fed contribution (extraction)
                E_ex = Di * OA
                xratio_aq = self._xN_over_x0(E_ex, N)
                nA_out_from_A = xratio_aq * nA_in
                n_to_org_from_A = max(nA_in - nA_out_from_A, 0.0)

                # organic-fed contribution (stripping)
                if AO == float("inf"):
                    xratio_org = 1.0
                else:
                    E_st = (AO / Di)
                    xratio_org = self._xN_over_x0(E_st, N)

                nO_out_from_O = xratio_org * nO_in
                n_to_aq_from_O = max(nO_in - nO_out_from_O, 0.0)

                Aq[j] = max(nA_out_from_A + n_to_aq_from_O, 0.0)
                Org[j] = max(nO_out_from_O + n_to_org_from_A, 0.0)
            else:
                Aq[j] = nA_in
                Org[j] = nO_in

        return Aq, Org, float(s)

    def _tbp_balance_at_OA(
        self,
        reg,
        idx: Dict[str, int],
        A_base: np.ndarray,
        O_base: np.ndarray,
        OA: float,
        N: int,
        A_tot: float,
        O_tot_base: float,
    ) -> Tuple[float, float, float, float]:
        """
        Returns (tbp_req, tbp_avail, org_scale_factor, O_tot_required) for a candidate OA.
        """
        Aq, Org, s = self._simulate_split_for_OA(reg, idx, A_base, O_base, OA, N, A_tot, O_tot_base)

        if self.tbp_species not in idx:
            raise ValueError(f"tbp_species '{self.tbp_species}' not found in registry.")
        j_tbp = idx[self.tbp_species]

        tbp_in_org = float(Org[j_tbp])
        tbp_avail = max(tbp_in_org - self.tbp_free_min, 0.0)

        tbp_req = 0.0
        for sp, nu in self.tbp_stoich.items():
            if sp not in idx:
                continue
            j = idx[sp]
            n_org = float(Org[j])
            if n_org > EPS and nu > 0:
                tbp_req += float(nu) * n_org

        O_tot_required = OA * A_tot
        return float(tbp_req), float(tbp_avail), float(s), float(O_tot_required)

    def design_N_and_OA(self) -> Tuple[int, float]:
        """
        Choose (N, OA_design) that minimizes:
            J = w_N * N + w_OA * OA
        subject to:
          - required_recovery constraints (mass transfer)
          - AND, if tbp_stoich provided, TBP stoichiometry feasibility
            at that OA when organic inlet is scaled proportionally to meet OA.
        """
        aq_in = self.inlets["aq_in"]
        org_in = self.inlets["org_in"]

        reg = aq_in.reg
        idx = reg.index()

        A_tot = aq_in.total_molar_flow()
        O_tot_base = org_in.total_molar_flow()

        A_base = aq_in.to_dense()
        O_base = org_in.to_dense()

        N_cap = min(self.N_max, self.N_balance_cap)

        best_N: Optional[int] = None
        best_OA: Optional[float] = None
        best_J: Optional[float] = None

        for N in range(1, N_cap + 1):
            # --- mass-transfer required OA to meet recoveries ---
            OA_req_mt = 0.0
            for sp, rec_in in self.required_recovery.items():
                if sp not in self.D:
                    raise ValueError(f"required_recovery species '{sp}' missing from D.")
                rec = self._as_fraction(float(rec_in))
                rec = min(max(rec, 0.0), 1.0)
                OA_req_mt = max(
                    OA_req_mt,
                    self._required_OA_for_recovery(self.D[sp], N, rec, self.OA_max),
                )

            OA_req = max(OA_req_mt, self.OA_min)

            # --- TBP stoichiometry feasibility may push OA higher ---
            if self.tbp_stoich:
                if self.tbp_species not in idx:
                    raise ValueError(f"tbp_species '{self.tbp_species}' not found in registry.")

                def ok(OA_test: float) -> bool:
                    tbp_req, tbp_avail, _, _ = self._tbp_balance_at_OA(
                        reg, idx, A_base, O_base, OA_test, N, A_tot, O_tot_base
                    )
                    return tbp_req <= tbp_avail + 1e-12

                if not ok(OA_req):
                    lo = OA_req
                    hi = max(1.0, OA_req)
                    while not ok(hi):
                        hi *= 2.0
                        if hi > self.OA_max:
                            raise ValueError("TBP stoichiometry infeasible within OA_max.")
                    for _ in range(80):
                        mid = 0.5 * (lo + hi)
                        if ok(mid):
                            hi = mid
                        else:
                            lo = mid
                    OA_req = hi

            # --- combined objective ---
            J = self.w_N * float(N) + self.w_OA * float(OA_req)

            if best_J is None or J < best_J - 1e-12:
                best_J = float(J)
                best_N = int(N)
                best_OA = float(OA_req)
            # deterministic tie-breakers (optional but helpful)
            elif abs(J - best_J) <= 1e-12:
                # prefer fewer stages, then lower OA
                if best_N is None or N < best_N or (N == best_N and OA_req < (best_OA or float("inf"))):
                    best_N = int(N)
                    best_OA = float(OA_req)

        if best_N is None or best_OA is None:
            raise ValueError("No feasible (N, O/A) found.")
        return int(best_N), float(best_OA)

    def apply(self) -> None:
        aq_in = self.inlets["aq_in"]
        org_in = self.inlets["org_in"]
        aq_out = self.outlets["aq_out"]
        org_out = self.outlets["org_out"]

        reg = aq_in.reg
        idx = reg.index()

        # ----- Design (includes TBP feasibility if enabled) -----
        N, OA_design = self.design_N_and_OA()
        self.N_used = int(N)
        self.OA_design = float(OA_design)

        # ----- Provided availability (for reporting only) -----
        A_tot = aq_in.total_molar_flow()
        O_tot_base = org_in.total_molar_flow()
        self.OA_used = O_tot_base / max(A_tot, EPS)

        # Use design OA for xratio, but cap by available organic at the inlet
        self.OA_effective = min(float(self.OA_design), float(self.OA_used))

        # Reporting only: what would be required vs what is provided
        self.O_tot_required = float(self.OA_design) * A_tot
        self.org_scale_factor = self.O_tot_required / max(O_tot_base, EPS)

        A = aq_in.to_dense()
        O = org_in.to_dense()

        OA = float(self.OA_effective)
        AO = (1.0 / OA) if OA > EPS else float("inf")

        Aq = np.zeros(reg.n(), dtype=float)
        Org = np.zeros(reg.n(), dtype=float)

        # ----- Split each species -----
        for sp in reg.species:
            j = idx[sp]
            nA_in = float(A[j])
            nO_in = float(O[j])
            nT = nA_in + nO_in
            if nT <= EPS:
                continue

            if sp in self.nonextract_to_aq:
                Aq[j] = nT
                continue
            if sp in self.nonextract_to_org:
                Org[j] = nT
                continue

            if sp in self.D:
                Di = float(self.D[sp])
                if Di <= 0.0:
                    Aq[j] = nT
                    continue

                # aqueous-fed contribution (extraction)
                E_ex = Di * OA
                xratio_aq = self._xN_over_x0(E_ex, self.N_used)
                nA_out_from_A = xratio_aq * nA_in
                n_to_org_from_A = max(nA_in - nA_out_from_A, 0.0)

                # organic-fed contribution (stripping)
                if AO == float("inf"):
                    xratio_org = 1.0
                else:
                    E_st = (AO / Di)
                    xratio_org = self._xN_over_x0(E_st, self.N_used)

                nO_out_from_O = xratio_org * nO_in
                n_to_aq_from_O = max(nO_in - nO_out_from_O, 0.0)

                Aq[j] = max(nA_out_from_A + n_to_aq_from_O, 0.0)
                Org[j] = max(nO_out_from_O + n_to_org_from_A, 0.0)
            else:
                Aq[j] = nA_in
                Org[j] = nO_in

        # ----- TBP feasibility reporting (design should make this feasible; apply caps OA by availability) -----
        self.tbp_limited = False
        self.tbp_required = 0.0
        self.tbp_available = 0.0
        if self.tbp_stoich:
            if self.tbp_species not in idx:
                raise ValueError(f"tbp_species '{self.tbp_species}' not found in registry.")
            j_tbp = idx[self.tbp_species]
            tbp_in_org = float(Org[j_tbp])
            tbp_avail = max(tbp_in_org - self.tbp_free_min, 0.0)

            tbp_req = 0.0
            for sp, nu in self.tbp_stoich.items():
                if sp not in idx:
                    continue
                j = idx[sp]
                n_org = float(Org[j])
                if n_org > EPS and nu > 0:
                    tbp_req += float(nu) * n_org

            self.tbp_required = float(tbp_req)
            self.tbp_available = float(tbp_avail)
            self.tbp_limited = tbp_req > tbp_avail + 1e-9

        aq_out.from_dense(Aq)
        org_out.from_dense(Org)