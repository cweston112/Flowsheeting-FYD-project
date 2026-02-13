from __future__ import annotations

from typing import Dict, List, Optional, Tuple
import numpy as np

from ..flowsheet_tools import UnitOp, EPS


class MultiStageCentrifugalContactorTBP(UnitOp):
    """
    Multi-stage centrifugal contactor (TBP extraction), counter-current stage model.

    D is defined as (org/aq). This updated model treats partitioning as TRUE / BIDIRECTIONAL:
      - If a species enters in the aqueous feed, some transfers to organic.
      - If a species enters in the organic feed, some transfers to aqueous.

    Stage-count performance is computed with the standard constant-flow counter-current
    "xN/x0" factor using an extraction factor E:
      - For aqueous -> organic extraction:   E_ex = D * (O/A) = D * OA
        aqueous-out / aqueous-in = xN_over_x0(E_ex, N)

      - For organic -> aqueous stripping:    E_st = (A/O) / D = AO / D = (1/OA) / D
        organic-out / organic-in = xN_over_x0(E_st, N)

    IMPORTANT: We make xratio follow the DESIGN O/A, but we cannot exceed available organic:
      OA_effective = min(OA_design, OA_used)

    If you include TBP in self.D (e.g. D["TBP"] = ...), TBP will partition accordingly.
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
        N_balance_cap: int = 20,
    ):
        super().__init__(name)
        self.D = {k: float(v) for k, v in D.items()}
        self.required_recovery = {k: float(v) for k, v in required_recovery.items()}
        self.nonextract_to_aq = set(nonextract_to_aq or [])
        self.nonextract_to_org = set(nonextract_to_org or [])
        self.N_max = int(N_max)
        self.OA_max = float(OA_max)
        self.N_balance_cap = int(N_balance_cap)

        # reported fields
        self.N_used: int = 1
        self.OA_design: float = 0.0
        self.OA_used: float = 0.0
        self.OA_effective: float = 0.0

    @staticmethod
    def _xN_over_x0(E: float, N: int) -> float:
        """
        x_N/x_0 for a counter-current cascade with constant phase flows and extraction factor E.
        """
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
        # recovery to organic from aqueous feed
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

    def design_N_and_OA(self) -> Tuple[int, float]:
        """
        Choose (N, OA_design) that minimizes J = OA*N subject to required_recovery (to organic)
        for the species listed in required_recovery.
        """
        N_cap = min(self.N_max, self.N_balance_cap)
        best_N, best_OA, best_J = None, None, float("inf")

        for N in range(1, N_cap + 1):
            OA_req = 0.0
            for sp, rec in self.required_recovery.items():
                if sp not in self.D:
                    raise ValueError(f"required_recovery species '{sp}' missing from D.")
                OA_req = max(OA_req, self._required_OA_for_recovery(self.D[sp], N, float(rec), self.OA_max))
            J = OA_req * N
            if J < best_J:
                best_J, best_N, best_OA = J, N, OA_req

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

        # ----- Design -----
        N, OA_design = self.design_N_and_OA()
        self.N_used = int(N)
        self.OA_design = float(OA_design)

        # ----- Actual availability -----
        A_tot = aq_in.total_molar_flow()
        O_tot = org_in.total_molar_flow()
        self.OA_used = O_tot / max(A_tot, EPS)

        # Use design OA for xratio, but cap by available organic
        self.OA_effective = min(self.OA_design, self.OA_used)

        # Convenience
        OA = float(self.OA_effective)
        AO = (1.0 / OA) if OA > EPS else float("inf")

        # Dense vectors
        A = aq_in.to_dense()
        O = org_in.to_dense()

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

            # Forced routing first
            if sp in self.nonextract_to_aq:
                Aq[j] = nT
                continue
            if sp in self.nonextract_to_org:
                Org[j] = nT
                continue

            # True partitioning if D provided
            if sp in self.D:
                Di = float(self.D[sp])

                # Handle degenerate D
                if Di <= 0.0:
                    # D -> 0 means overwhelmingly aqueous: everything ends aqueous
                    Aq[j] = nT
                    continue

                # (1) Contribution from aqueous feed (aq -> org extraction)
                E_ex = Di * OA
                xratio_aq = self._xN_over_x0(E_ex, self.N_used)  # aq_out / aq_in for aq-fed solute
                nA_out_from_A = xratio_aq * nA_in
                n_to_org_from_A = max(nA_in - nA_out_from_A, 0.0)

                # (2) Contribution from organic feed (org -> aq stripping)
                # extraction factor for stripping: E_st = AO / D
                if AO == float("inf"):
                    # no organic present -> organic feed cannot be stripped meaningfully; everything stays aqueous anyway
                    xratio_org = 1.0
                else:
                    E_st = (AO / Di)
                    xratio_org = self._xN_over_x0(E_st, self.N_used)  # org_out / org_in for org-fed solute

                nO_out_from_O = xratio_org * nO_in
                n_to_aq_from_O = max(nO_in - nO_out_from_O, 0.0)

                # Totals
                Aq[j] = max(nA_out_from_A + n_to_aq_from_O, 0.0)
                Org[j] = max(nO_out_from_O + n_to_org_from_A, 0.0)
            else:
                # Not modelled as transferable: stays where it is
                Aq[j] = nA_in
                Org[j] = nO_in

        aq_out.from_dense(Aq)
        org_out.from_dense(Org)
