from __future__ import annotations

from typing import Dict, List, Optional, Tuple
import numpy as np

from ..flowsheet_tools import UnitOp, EPS


class MultiStageCentrifugalContactorTBP(UnitOp):
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
        self.N_used: int = 1
        self.OA_design: float = 0.0
        self.OA_used: float = 0.0

    @staticmethod
    def _xN_over_x0(E: float, N: int) -> float:
        if abs(E - 1.0) < 1e-12:
            return 1.0 / (N + 1.0)
        s = (E ** (N + 1) - 1.0) / (E - 1.0)
        return 1.0 / s

    @classmethod
    def _recovery(cls, D: float, OA: float, N: int) -> float:
        return 1.0 - cls._xN_over_x0(D * OA, N)

    @classmethod
    def _required_OA_for_recovery(cls, D: float, N: int, recovery: float, OA_max: float) -> float:
        lo, hi = 0.0, 1.0
        while cls._recovery(D, hi, N) < recovery:
            hi *= 2.0
            if hi > OA_max:
                raise ValueError("Required O/A exceeds OA_max.")
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if cls._recovery(D, mid, N) >= recovery:
                hi = mid
            else:
                lo = mid
        return hi

    def design_N_and_OA(self) -> Tuple[int, float]:
        N_cap = min(self.N_max, self.N_balance_cap)
        best_N, best_OA, best_J = None, None, float("inf")
        for N in range(1, N_cap + 1):
            OA_req = 0.0
            for sp, rec in self.required_recovery.items():
                OA_req = max(OA_req, self._required_OA_for_recovery(self.D[sp], N, rec, self.OA_max))
            J = OA_req * N
            if J < best_J:
                best_J, best_N, best_OA = J, N, OA_req
        if best_N is None or best_OA is None:
            raise ValueError("No feasible (N, O/A) found.")
        return best_N, best_OA

    def apply(self) -> None:
        aq_in = self.inlets["aq_in"]
        org_in = self.inlets["org_in"]
        aq_out = self.outlets["aq_out"]
        org_out = self.outlets["org_out"]

        reg = aq_in.reg
        idx = reg.index()

        N, OA = self.design_N_and_OA()
        self.N_used = N
        self.OA_design = OA

        A_tot = aq_in.total_molar_flow()
        O_tot = org_in.total_molar_flow()
        self.OA_used = (O_tot / max(A_tot, EPS))

        A = aq_in.to_dense()
        O = org_in.to_dense()

        Aq = np.zeros(reg.n(), dtype=float)
        Org = np.zeros(reg.n(), dtype=float)

        for sp in reg.species:
            j = idx[sp]
            nA, nO = A[j], O[j]
            nT = nA + nO
            if nT == 0.0:
                continue
            if sp in self.nonextract_to_aq:
                Aq[j] = nT
                continue
            if sp in self.nonextract_to_org:
                Org[j] = nT
                continue
            if sp in self.D:
                xratio = self._xN_over_x0(self.D[sp] * self.OA_used, N)
                nA_out = xratio * nA
                nOrg_out = (nA - nA_out) + nO
                Aq[j] = max(nA_out, 0.0)
                Org[j] = max(nOrg_out, 0.0)
            else:
                Aq[j] = nA
                Org[j] = nO

        aq_out.from_dense(Aq)
        org_out.from_dense(Org)
