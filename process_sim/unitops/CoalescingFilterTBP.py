from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

from ..flowsheet_tools import UnitOp, Stream, EPS


class CoalescingFilterTBP(UnitOp):
    """
    Coalescing filter to remove entrained TBP/diluent from an aqueous stream.
    Notes on parameterization:
      - TBP solubility in water ~ 0.28–0.40 g/L at ~25°C (≈ 1.05e-3 to 1.50e-3 mol/L). (PubChem/Wikipedia)
      - n-dodecane solubility is extremely low (~3.7e-3 mg/L reported by PubChem; order 1e-8 mol/L).
      - Coalescer "removal" here is applied only to the ENTRAINED portion (above solubility).
    """

    def __init__(
        self,
        name: str,
        *,
        tbp_name: str = "TBP",
        diluent_name: str = "Dodecane",
        water_name: str = "H2O",
        hno3_name: str = "HNO3",
        water_molarity_pure: float = 55.5,  # mol/L
        removal_eff: float = 0.99,          # fraction of ENTRAINED removed
        # dissolved-equilibrium solubilities (mol/L) at ~ambient T
        tbp_solubility_weak_mol_L: float = 1.5e-3,
        tbp_solubility_strong_mol_L: float = 5.0e-4,
        dil_solubility_weak_mol_L: float = 2.0e-8,
        dil_solubility_strong_mol_L: float = 1.0e-8,
        # "strong acid" threshold used to interpolate solubility based on HNO3 molarity
        hno3_strong_M: float = 3.0,
        assume_pure_water_for_weak_stream: bool = False,
        print_diagnostics: bool = False,
    ):
        super().__init__(name)
        self.tbp_name = tbp_name
        self.diluent_name = diluent_name
        self.water_name = water_name
        self.hno3_name = hno3_name
        self.water_molarity_pure = float(water_molarity_pure)

        self.removal_eff = float(removal_eff)
        if not (0.0 <= self.removal_eff <= 1.0):
            raise ValueError("removal_eff must be in [0,1]")

        self.tbp_solubility_weak_mol_L = float(tbp_solubility_weak_mol_L)
        self.tbp_solubility_strong_mol_L = float(tbp_solubility_strong_mol_L)
        self.dil_solubility_weak_mol_L = float(dil_solubility_weak_mol_L)
        self.dil_solubility_strong_mol_L = float(dil_solubility_strong_mol_L)

        self.hno3_strong_M = float(hno3_strong_M)
        self.assume_pure_water_for_weak_stream = bool(assume_pure_water_for_weak_stream)
        self.print_diagnostics = bool(print_diagnostics)

        # For reporting
        self.last_Vdot_L_s: float = 0.0
        self.last_HNO3_M: float = 0.0
        self.last_removed: Dict[str, float] = {}

    def _estimate_Vdot_L_s(self, s: Stream) -> float:
        """
        Estimate volumetric flow (L/s) from water molar flow using 55.5 mol/L basis.
        """
        nH2O = s.get(self.water_name)
        if self.water_molarity_pure <= EPS:
            return 0.0
        return max(nH2O / self.water_molarity_pure, 0.0)

    def _interp_solubility(self, sol_weak: float, sol_strong: float, HNO3_M: float) -> float:
        """
        Linearly interpolate solubility between weak and strong as acid increases 0 -> hno3_strong_M.
        Clamp outside range.
        """
        if self.hno3_strong_M <= EPS:
            return sol_strong
        f = max(0.0, min(1.0, HNO3_M / self.hno3_strong_M))
        # At HNO3=0 -> weak; at HNO3>=strong -> strong
        return (1.0 - f) * sol_weak + f * sol_strong

    def apply(self) -> None:
        sin = self.inlets["in"]
        saq = self.outlets["aq_out"]
        sorg = self.outlets["org_recovered"]

        # Start with copies (component-wise)
        saq.mol = dict(sin.mol)
        sorg.mol = {}

        # Volumetric flow estimate based on water only (consistent with your molarity helpers)
        Vdot = self._estimate_Vdot_L_s(sin)  # L/s

        # Acid molarity estimate (optionally overridden for weak stream case)
        if self.assume_pure_water_for_weak_stream:
            HNO3_M = 0.0
        else:
            nHNO3 = sin.get(self.hno3_name)
            HNO3_M = (nHNO3 / Vdot) if (Vdot > EPS) else 0.0

        # Determine dissolved (equilibrium) limits
        tbp_sol = self._interp_solubility(
            self.tbp_solubility_weak_mol_L,
            self.tbp_solubility_strong_mol_L,
            HNO3_M,
        )
        dil_sol = self._interp_solubility(
            self.dil_solubility_weak_mol_L,
            self.dil_solubility_strong_mol_L,
            HNO3_M,
        )

        tbp_allowed = tbp_sol * Vdot
        dil_allowed = dil_sol * Vdot

        # Remove entrained portion above allowed dissolved
        removed: Dict[str, float] = {}

        def _split(sp: str, allowed: float) -> None:
            n_in = sin.get(sp)
            if n_in <= EPS:
                saq.set(sp, 0.0)
                sorg.set(sp, 0.0)
                removed[sp] = 0.0
                return

            # dissolved portion (kept)
            n_keep = min(n_in, max(allowed, 0.0))

            # entrained (non-eq) portion
            n_entr = max(n_in - n_keep, 0.0)

            # removed to organic recovered
            n_to_org = self.removal_eff * n_entr

            # remainder stays in aqueous
            n_out_aq = n_keep + (n_entr - n_to_org)

            saq.set(sp, n_out_aq)
            sorg.set(sp, n_to_org)
            removed[sp] = n_to_org

        _split(self.tbp_name, tbp_allowed)
        _split(self.diluent_name, dil_allowed)

        # Stamp T/p/phase (keep aqueous phase on aq_out; organic recovered as L)
        saq.phase = "L"
        sorg.phase = "L"
        saq.T = sin.T
        saq.p = sin.p
        sorg.T = sin.T
        sorg.p = sin.p

        # Report
        self.last_Vdot_L_s = float(Vdot)
        self.last_HNO3_M = float(HNO3_M)
        self.last_removed = removed

        if self.print_diagnostics:
            print(f"[{self.name}] Vdot={Vdot:.4g} L/s, HNO3~{HNO3_M:.3g} M")
            for sp, nrem in removed.items():
                if nrem > 0:
                    print(f"  removed {sp}: {nrem:.4g} mol/s (entrained fraction only)")
