from __future__ import annotations
from ..flowsheet_tools import UnitOp

class BufferVessel(UnitOp):
    def __init__(self, name: str, *, set_T: float | None = None, set_p: float | None = None, set_phase: str | None = None):
        super().__init__(name)
        self.set_T = set_T
        self.set_p = set_p
        self.set_phase = set_phase

    def apply(self) -> None:
        fin = self.inlets["in"]
        fout = self.outlets["out"]
        fout.mol = dict(fin.mol)
        fout.T = fin.T if self.set_T is None else self.set_T
        fout.p = fin.p if self.set_p is None else self.set_p
        fout.phase = fin.phase if self.set_phase is None else self.set_phase
