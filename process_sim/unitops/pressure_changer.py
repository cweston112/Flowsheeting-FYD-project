from __future__ import annotations
from ..flowsheet_tools import UnitOp, Stream

class PressureChanger(UnitOp):
    """1-in/1-out: copy composition/T/phase, set p."""
    def __init__(self, name: str, *, set_p: float):
        super().__init__(name)
        if set_p <= 0:
            raise ValueError("set_p must be > 0 Pa")
        self.set_p = float(set_p)

    def apply(self) -> None:
        Sin: Stream = self.inlets["in"]
        Sout: Stream = self.outlets["out"]
        Sout.mol = dict(Sin.mol)
        Sout.T = Sin.T
        Sout.phase = Sin.phase
        Sout.p = self.set_p
