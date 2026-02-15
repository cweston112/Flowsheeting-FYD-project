from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

from ..flowsheet_tools import EPS, Stream, UnitOp


class IdealTSAColumnEMM17(UnitOp):
    """
    Ideal temperature swing adsorption (TSA) column model for iodine on EMM-17.

    Two operating modes:
      - mode="adsorb":
          * removes iodine from process gas (perfect capture until capacity is exceeded)
          * outputs cleaned gas (F50X)
      - mode="regen":
          * desorbs iodine (time-averaged) into a regeneration gas (air here)
          * outputs iodine-rich gas (F51X)

    Key idealisations:
      - Adsorption capacity uses a volumetric working capacity:
            cap_work_mol_per_m3 = (uptake_g_per_g * 1000 / MW_I2) * packing_density_kg_m3 * mtz_util
        where mtz_util = 0.70 implements your “70% of bed volume utilised” MTZ factor.
      - Regeneration is time-averaged: iodine captured during adsorption half-cycle
        is released uniformly during regeneration half-cycle.
      - “Minimise regen flow”: your sizing routine should pick the smallest air flow
        that keeps y_I2_out <= regen_yI2_max (to avoid impossible/supersaturated gas).
    """

    def __init__(
        self,
        name: str,
        *,
        mode: str,  # "adsorb" or "regen"
        iodine_name: str = "I_g",
        bed_volume_m3: float = 1.0,
        cap_work_mol_per_m3: float = 0.0,
        t_ads_s: float = 8.0 * 3600.0,
        t_regen_s: float = 4.0 * 3600.0,
        regen_yI2_max: float = 0.02,
        print_diagnostics: bool = False,
    ):
        super().__init__(name)
        self.mode = str(mode).lower().strip()
        if self.mode not in ("adsorb", "regen"):
            raise ValueError(f"{name}: mode must be 'adsorb' or 'regen'")

        self.iodine_name = iodine_name
        self.bed_volume_m3 = float(bed_volume_m3)
        self.cap_work_mol_per_m3 = float(cap_work_mol_per_m3)
        self.t_ads_s = float(t_ads_s)
        self.t_regen_s = float(t_regen_s)
        self.regen_yI2_max = float(regen_yI2_max)
        self.print_diagnostics = bool(print_diagnostics)

        # calculated / reported
        self.I_in_mol_s: float = 0.0
        self.I_removed_mol_s: float = 0.0
        self.I_slip_mol_s: float = 0.0
        self.I_desorb_mol_s: float = 0.0
        self.capacity_cycle_mol: float = 0.0
        self.loading_cycle_mol: float = 0.0

    def apply(self) -> None:
        if self.mode == "adsorb":
            self._apply_adsorb()
        else:
            self._apply_regen()

    def _apply_adsorb(self) -> None:
        gas_in = self.inlets["gas_in"]
        gas_out = self.outlets["gas_out"]

        # pass-through by default
        gas_out.mol = {}
        for sp, n in gas_in.mol.items():
            gas_out.set(sp, n)

        I_in = float(gas_in.get(self.iodine_name))
        self.I_in_mol_s = I_in

        # Available working capacity (per adsorption step)
        cap_cycle = max(self.cap_work_mol_per_m3, 0.0) * max(self.bed_volume_m3, 0.0)
        load_cycle = I_in * max(self.t_ads_s, 0.0)

        self.capacity_cycle_mol = float(cap_cycle)
        self.loading_cycle_mol = float(load_cycle)

        if cap_cycle <= EPS or self.t_ads_s <= EPS:
            # no capacity -> no removal
            I_removed = 0.0
        else:
            # perfect capture until capacity is exceeded over the adsorption half-cycle
            # if load_cycle > cap_cycle -> some breakthrough (slip) to satisfy capacity
            frac_captured = min(cap_cycle / max(load_cycle, EPS), 1.0)
            I_removed = frac_captured * I_in

        I_slip = max(I_in - I_removed, 0.0)

        # write outlet iodine
        gas_out.set(self.iodine_name, I_slip)

        self.I_removed_mol_s = float(I_removed)
        self.I_slip_mol_s = float(I_slip)

        if self.print_diagnostics:
            print(f"[{self.name}] ADSORB: I_in={I_in:.6g} mol/s, "
                  f"I_removed={I_removed:.6g}, I_slip={I_slip:.6g}, "
                  f"cap_cycle={cap_cycle:.6g} mol, load_cycle={load_cycle:.6g} mol")

    def _apply_regen(self) -> None:
        regen_in = self.inlets["regen_in"]
        regen_out = self.outlets["regen_out"]

        regen_out.mol = {}
        for sp, n in regen_in.mol.items():
            regen_out.set(sp, n)

        # The regeneration iodine release rate must be supplied by your sizing logic
        # via putting iodine into regen_in or by setting an attribute prior to apply.
        # Here we interpret ANY iodine already on regen_in as the release rate.
        I_rel = float(regen_in.get(self.iodine_name))
        self.I_desorb_mol_s = I_rel

        # Put it in outlet (already copied, but ensure)
        regen_out.set(self.iodine_name, I_rel)

        if self.print_diagnostics:
            n_tot = regen_out.total_molar_flow()
            y = I_rel / max(n_tot, EPS)
            print(f"[{self.name}] REGEN: I_out={I_rel:.6g} mol/s, total={n_tot:.6g} mol/s, y_I={y:.3e}")
