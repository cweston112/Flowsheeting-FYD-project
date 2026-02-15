from __future__ import annotations

from ..flowsheet_tools import UnitOp, EPS


class UO3ToU3O8Reactor(UnitOp):
    """
    Convert all UO3(s) -> U3O8(s) with O2 evolution.

    Stoichiometry (balanced):
      3 UO3(s) -> U3O8(s) + 0.5 O2(g)

    Ports:
      in:      feed solid stream (expects UO3)
      solids:  product solid stream (U3O8)
      offgas:  gas stream (O2)
    """

    def __init__(
        self,
        name: str,
        *,
        uo3_name: str = "UO3",
        u3o8_name: str = "U3O8",
        o2_name: str = "O2",
        print_diagnostics: bool = False,
    ):
        super().__init__(name)
        self.uo3_name = str(uo3_name)
        self.u3o8_name = str(u3o8_name)
        self.o2_name = str(o2_name)
        self.print_diagnostics = bool(print_diagnostics)

    def apply(self) -> None:
        Sin = self.inlets["in"]
        Ssol = self.outlets["solids"]
        Sgas = self.outlets["offgas"]

        n_uo3 = Sin.get(self.uo3_name)

        # Clear outlets (explicitly write only produced species)
        Ssol.mol = {}
        Sgas.mol = {}

        if n_uo3 <= EPS:
            # nothing to do
            if self.print_diagnostics:
                print(f"[{self.name}] n(UO3)~0 -> no conversion")
            return

        # 3 UO3 -> U3O8 + 0.5 O2
        n_u3o8 = n_uo3 / 3.0
        n_o2 = (0.5 / 3.0) * n_uo3  # = n_uo3/6

        Ssol.set(self.u3o8_name, n_u3o8)
        Sgas.set(self.o2_name, n_o2)

        # Phases (keep explicit)
        Ssol.phase = "S"
        Sgas.phase = "G"

        # Carry T/p from inlet (will be overwritten by stamp_outlets_to_unit_conditions anyway)
        Ssol.T = Sin.T; Ssol.p = Sin.p
        Sgas.T = Sin.T; Sgas.p = Sin.p

        if self.print_diagnostics:
            print(
                f"[{self.name}] n(UO3)={n_uo3:.6g} -> "
                f"U3O8={n_u3o8:.6g}, O2={n_o2:.6g}"
            )
