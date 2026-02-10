from __future__ import annotations
from ..flowsheet_tools import UnitOp, Stream

class HeaterCooler(UnitOp):
    """1-in/1-out: copy composition/p/phase, set T."""
    def __init__(self, name: str, *, set_T: float):
        super().__init__(name)
        if set_T <= 0:
            raise ValueError("set_T must be > 0 K")
        self.set_T = float(set_T)

    def apply(self) -> None:
        Sin: Stream = self.inlets["in"]
        Sout: Stream = self.outlets["out"]
        Sout.mol = dict(Sin.mol)
        Sout.p = Sin.p
        Sout.phase = Sin.phase
        Sout.T = self.set_T
