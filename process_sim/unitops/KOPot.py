from __future__ import annotations

from dataclasses import dataclass

from ..flowsheet_tools import UnitOp, EPS

class KnockOutPotWater(UnitOp):
    """
    Simple knock-out pot that removes all water from a gas stream.
      - Inlet: mixed gas (may contain H2O)
      - Outlet "liq": all H2O
      - Outlet "gas": everything else

    Assumes complete condensation/knockout of water at operating conditions.
    """
    def __init__(self, name: str, *, water_name: str = "H2O", print_diagnostics: bool = False):
        super().__init__(name)
        self.water_name = str(water_name)
        self.print_diagnostics = bool(print_diagnostics)

    def apply(self) -> None:
        inp = self.inlets["in"]
        gas = self.outlets["gas"]
        liq = self.outlets["liq"]

        # clear
        gas.mol = {}
        liq.mol = {}

        nH2O = float(inp.get(self.water_name))

        # send all non-water to gas
        for sp, n in inp.mol.items():
            n = float(n)
            if n <= EPS:
                continue
            if sp == self.water_name:
                continue
            gas.set(sp, n)

        # send all water to liquid
        if nH2O > EPS:
            liq.set(self.water_name, nH2O)

        # phases (useful even before stamping)
        gas.phase = "G"
        liq.phase = "L"

        if self.print_diagnostics:
            print(f"[{self.name}] knocked out H2O={nH2O:.6g} mol/s")
