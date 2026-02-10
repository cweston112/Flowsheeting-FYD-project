from __future__ import annotations

from typing import Dict, List, Optional, Tuple
import numpy as np

from ..flowsheet_tools import UnitOp, EPS


class MultiStageCentrifugalContactorAHA_Strip(UnitOp):
    """
    X2 AHA STRIP contactor (countercurrent, ideal stages).

    Convention:
      - D is provided as org/aq distribution coefficient: D = C_org / C_aq.
      - Stripping removes solute from organic -> aqueous, so we target recovery-to-aqueous.

    Extraction factor for stripping:
      E_strip = (A/O) / D

    Ports:
      inlets:  "org_in" (loaded organic), "aq_in" (strip solution)
      outlets: "org_out" (stripped organic), "aq_out" (loaded aqueous)
    """

    def __init__(
        self,
        name: str,
        D_org_over_aq: Dict[str, float],
        required_recovery_to_aq: Dict[str, float],  # fraction or percent
        nontransfer_keep_in_org: Optional[List[str]] = None,
        nontransfer_keep_in_aq: Optional[List[str]] = None,
        N_max: int = 50,
        AO_max: float = 1e3,          # A/O max used for design feasibility checks
        N_balance_cap: int = 20,
    ):
        super().__init__(name)
        self.D = {k: float(v) for k, v in D_org_over_aq.items()}

        # accept either fractions (<=1.0) or percent (>1.0)
        self.required_recovery_to_aq = {
            k: (float(v) / 100.0 if float(v) > 1.0 else float(v))
            for k, v in required_recovery_to_aq.items()
        }

        self.nontransfer_keep_in_org = set(nontransfer_keep_in_org or [])
        self.nontransfer_keep_in_aq = set(nontransfer_keep_in_aq or [])

        self.N_max = int(N_max)
        self.AO_max = float(AO_max)
        self.N_balance_cap = int(N_balance_cap)

        # reported / computed
        self.N_used: int = 1
        self.AO_design: float = 0.0
        self.AO_used: float = 0.0

    @staticmethod
    def _xN_over_x0(E: float, N: int) -> float:
        # Ideal countercurrent stage relation
        if abs(E - 1.0) < 1e-12:
            return 1.0 / (N + 1.0)
        s = (E ** (N + 1) - 1.0) / (E - 1.0)
        return 1.0 / s

    @classmethod
    def _recovery_to_aq(cls, D_org_over_aq: float, AO: float, N: int) -> float:
        # E_strip = (A/O)/D
        if D_org_over_aq <= 0.0:
            # D ~ 0 => overwhelmingly aqueous; stripping is essentially complete
            return 1.0
        E = (AO / D_org_over_aq)
        return 1.0 - cls._xN_over_x0(E, N)

    @classmethod
    def _required_AO_for_recovery_to_aq(cls, D: float, N: int, recovery: float, AO_max: float) -> float:
        # Find AO such that recovery_to_aq >= target
        lo, hi = 0.0, 1.0
        while cls._recovery_to_aq(D, hi, N) < recovery:
            hi *= 2.0
            if hi > AO_max:
                raise ValueError("Required A/O exceeds AO_max.")
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if cls._recovery_to_aq(D, mid, N) >= recovery:
                hi = mid
            else:
                lo = mid
        return hi

    def design_N_and_AO(self) -> Tuple[int, float]:
        """
        Choose (N, AO_design) meeting all required recovery-to-aqueous targets,
        minimizing J = AO_design * N.
        """
        N_cap = min(self.N_max, self.N_balance_cap)
        best_N, best_AO, best_J = None, None, float("inf")

        for N in range(1, N_cap + 1):
            AO_req = 0.0
            for sp, rec in self.required_recovery_to_aq.items():
                if sp not in self.D:
                    raise KeyError(f"{self.name}: required recovery species '{sp}' missing from D.")
                AO_req = max(AO_req, self._required_AO_for_recovery_to_aq(self.D[sp], N, rec, self.AO_max))
            J = AO_req * N
            if J < best_J:
                best_J, best_N, best_AO = J, N, AO_req

        if best_N is None or best_AO is None:
            raise ValueError("No feasible (N, A/O) found.")
        return int(best_N), float(best_AO)

    def apply(self) -> None:
        org_in = self.inlets["org_in"]
        aq_in = self.inlets["aq_in"]
        org_out = self.outlets["org_out"]
        aq_out = self.outlets["aq_out"]

        reg = org_in.reg
        idx = reg.index()

        # Design (for reporting / optional sizing checks)
        N, AO_design = self.design_N_and_AO()
        self.N_used = int(N)
        self.AO_design = float(AO_design)

        # Actual A/O used by this unit
        O_tot = org_in.total_molar_flow()
        A_tot = aq_in.total_molar_flow()
        self.AO_used = float(A_tot / max(O_tot, EPS))
        AO_eff = min(max(self.AO_used, 0.0), self.AO_max)

        O = org_in.to_dense()
        A = aq_in.to_dense()

        Org = np.zeros(reg.n(), dtype=float)
        Aq = np.zeros(reg.n(), dtype=float)

        for sp in reg.species:
            j = idx[sp]
            nO = float(O[j])
            nA = float(A[j])
            nT = nO + nA
            if nT <= 0.0:
                continue

            if sp in self.nontransfer_keep_in_org:
                Org[j] = nT
                continue
            if sp in self.nontransfer_keep_in_aq:
                Aq[j] = nT
                continue

            if sp in self.D:
                D = self.D[sp]
                if D <= 0.0:
                    # D ~ 0 => goes to aqueous
                    Org[j] = 0.0
                    Aq[j] = nT
                    continue

                # Stripping removes from organic:
                # organic raffinate fraction = xN/x0 where E = (A/O)/D
                E = AO_eff / D
                xratio = self._xN_over_x0(E, N)

                nO_out = xratio * nO
                nA_out = (nO - nO_out) + nA

                Org[j] = max(nO_out, 0.0)
                Aq[j] = max(nA_out, 0.0)
            else:
                # Default: keep each phase's own inventory
                Org[j] = nO
                Aq[j] = nA

        org_out.from_dense(Org)
        aq_out.from_dense(Aq)
