from __future__ import annotations

from ..flowsheet_tools import UnitOp


class Splitter(UnitOp):
    def __init__(self, name: str, frac_to_A: float):
        super().__init__(name)
        if not (0.0 <= frac_to_A <= 1.0):
            raise ValueError("frac_to_A must be in [0,1]")
        self.frac_to_A = float(frac_to_A)

    def apply(self) -> None:
        Sin = self.inlets["in"]
        SA = self.outlets["A"]
        SB = self.outlets["B"]

        Fin = Sin.to_dense()
        FA = self.frac_to_A * Fin
        FB = (1.0 - self.frac_to_A) * Fin

        SA.from_dense(FA)
        SB.from_dense(FB)
