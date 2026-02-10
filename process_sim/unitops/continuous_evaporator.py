from __future__ import annotations

from typing import Dict, Optional
import numpy as np

from ..flowsheet_tools import UnitOp, EPS


class ContinuousEvaporatorAHA(UnitOp):
    """
    Continuous evaporator with:
      - Coupled evaporation controlled by H2O
      - Floors to prevent crystallisation
      - Full AHA destruction
      - Built-in diagnostics printing
    """

    def __init__(
        self,
        name: str,
        *,
        T_set_K: float,
        frac_vap: Dict[str, float],
        min_liq_mol_s: Optional[Dict[str, float]] = None,
        latent_kJ_per_mol: Optional[Dict[str, float]] = None,
        cp_liq_J_per_molK: float = 75.0,
        aha_name: str = "AHA",
        aha_products: Optional[Dict[str, float]] = None,
        aha_hno3_consumption: float = 0.0,
        coupled_by_water: bool = True,
        water_name: str = "H2O",
        # diagnostics
        print_diagnostics: bool = False,
        print_every: int = 1,
    ):
        super().__init__(name)

        self.T_set_K = float(T_set_K)
        self.frac_vap = {k: float(v) for k, v in frac_vap.items()}
        self.min_liq_mol_s = {k: float(v) for k, v in (min_liq_mol_s or {}).items()}

        self.latent_kJ_per_mol = {"H2O": 40.7, "HNO3": 39.0}
        if latent_kJ_per_mol:
            self.latent_kJ_per_mol.update(latent_kJ_per_mol)

        self.cp_liq_J_per_molK = float(cp_liq_J_per_molK)

        self.aha_name = aha_name
        self.aha_products = aha_products or {"AcOH": 1.0, "N2O": 0.5}
        self.aha_hno3_consumption = float(aha_hno3_consumption)

        self.coupled_by_water = coupled_by_water
        self.water_name = water_name

        # diagnostics
        self.print_diagnostics = print_diagnostics
        self.print_every = max(int(print_every), 1)
        self._call_counter = 0

        # reported fields
        self.theta = 0.0
        self.Q_sensible_kW = 0.0
        self.Q_latent_kW = 0.0
        self.Q_total_kW = 0.0
        self.vap_rate_mol_s: Dict[str, float] = {}
        self.liq_rate_mol_s: Dict[str, float] = {}
        self.aha_destroyed_mol_s = 0.0

    # ------------------------------------------------------------------
    def _maybe_print(self, msg: str) -> None:
        if self.print_diagnostics and (self._call_counter % self.print_every == 0):
            print(msg)

    # ------------------------------------------------------------------
    def apply(self) -> None:
        self._call_counter += 1

        s_in = self.inlets["in"]
        s_liq = self.outlets["liq"]
        s_vap = self.outlets["vap"]

        reg = s_in.reg
        idx = reg.index()
        x_in = s_in.to_dense().copy()

        # ==============================================================
        # AHA destruction
        # ==============================================================
        n_aha = float(s_in.get(self.aha_name))
        self.aha_destroyed_mol_s = n_aha

        if n_aha > EPS:
            if self.aha_name in idx:
                x_in[idx[self.aha_name]] = 0.0

            if self.aha_hno3_consumption > 0 and "HNO3" in idx:
                x_in[idx["HNO3"]] = max(
                    x_in[idx["HNO3"]] - self.aha_hno3_consumption * n_aha, 0.0
                )

            for sp, nu in self.aha_products.items():
                if sp not in reg.mw:
                    reg.add_species(sp, mw_g_mol=1.0)
                    idx = reg.index()
                x_in[idx[sp]] += nu * n_aha

        # ==============================================================
        # Coupled evaporation (H2O controls θ)
        # ==============================================================
        x_vap = np.zeros(reg.n())
        x_liq = x_in.copy()
        self.vap_rate_mol_s = {}
        self.liq_rate_mol_s = {}

        w = self.water_name
        fw = self.frac_vap.get(w, 0.0)
        nw = x_in[idx[w]] if w in idx else 0.0
        floor_w = self.min_liq_mol_s.get(w, 0.0)

        vap_nom_w = fw * nw
        vap_max_w = max(nw - floor_w, 0.0)
        vap_w = min(vap_nom_w, vap_max_w)

        self.theta = vap_w / vap_nom_w if vap_nom_w > EPS else 0.0

        if w in idx:
            x_vap[idx[w]] = vap_w
            x_liq[idx[w]] = nw - vap_w
            self.vap_rate_mol_s[w] = vap_w
            self.liq_rate_mol_s[w] = nw - vap_w

        for sp, f in self.frac_vap.items():
            if sp == w or sp not in idx:
                continue
            n = x_in[idx[sp]]
            if n <= EPS:
                continue
            floor = self.min_liq_mol_s.get(sp, 0.0)
            vap = min(self.theta * f * n, max(n - floor, 0.0))
            x_vap[idx[sp]] = vap
            x_liq[idx[sp]] = n - vap
            self.vap_rate_mol_s[sp] = vap
            self.liq_rate_mol_s[sp] = n - vap

        # ==============================================================
        # Energy
        # ==============================================================
        dT = max(self.T_set_K - s_in.T, 0.0)
        n_heat = sum(x_in[idx[sp]] for sp in ("H2O", "HNO3") if sp in idx)
        self.Q_sensible_kW = n_heat * self.cp_liq_J_per_molK * dT / 1000.0

        self.Q_latent_kW = sum(
            n * self.latent_kJ_per_mol.get(sp, 0.0)
            for sp, n in self.vap_rate_mol_s.items()
        )

        self.Q_total_kW = self.Q_sensible_kW + self.Q_latent_kW

        # ==============================================================
        # Write outlets
        # ==============================================================
        s_liq.from_dense(x_liq)
        s_vap.from_dense(x_vap)
        s_liq.phase = "L"
        s_vap.phase = "G"
        s_liq.T = self.T_set_K
        s_vap.T = self.T_set_K

        # ==============================================================
        # Diagnostics printout
        # ==============================================================
        self._maybe_print(
            f"""
--- {self.name} ---
T = {self.T_set_K:.1f} K | θ = {self.theta:.3f}
Inlet:  F = {s_in.total_molar_flow():.4e} mol/s
Water:  in={nw:.4e} vap={vap_w:.4e} liq={nw - vap_w:.4e} floor={floor_w:.2e}
HNO3:   in={s_in.get('HNO3'):.4e} vap={self.vap_rate_mol_s.get('HNO3',0):.4e}
AHA destroyed: {self.aha_destroyed_mol_s:.3e} mol/s
Q_sens = {self.Q_sensible_kW:.3f} kW
Q_lat  = {self.Q_latent_kW:.3f} kW
Q_tot  = {self.Q_total_kW:.3f} kW
--------------------------------
"""
        )
