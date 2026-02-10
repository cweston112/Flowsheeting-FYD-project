from __future__ import annotations

from typing import Dict
import numpy as np

from ..flowsheet_tools import UnitOp


class SpeciesSplitter(UnitOp):
    """Split selected species to outlet A by specified fractions. Unspecified species default to 0 (stay in B)."""

    def __init__(self, name: str, frac_to_A: Dict[str, float], default_to_A: float = 0.0):
        super().__init__(name)
        self.frac_to_A = {k: float(v) for k, v in frac_to_A.items()}
        self.default_to_A = float(default_to_A)
        for sp, f in self.frac_to_A.items():
            if not (0.0 <= f <= 1.0):
                raise ValueError(f"{name}: frac_to_A[{sp}] must be in [0,1]")

    def apply(self) -> None:
        Sin = self.inlets["in"]
        SA = self.outlets["A"]
        SB = self.outlets["B"]

        reg = Sin.reg
        idx = reg.index()
        Fin = Sin.to_dense()

        FA = np.zeros(reg.n(), dtype=float)
        FB = np.zeros(reg.n(), dtype=float)

        for sp in reg.species:
            j = idx[sp]
            n = Fin[j]
            if n == 0.0:
                continue
            fA = self.frac_to_A.get(sp, self.default_to_A)
            FA[j] = fA * n
            FB[j] = (1.0 - fA) * n

        SA.from_dense(FA)
        SB.from_dense(FB)
