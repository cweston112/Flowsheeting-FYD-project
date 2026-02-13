from __future__ import annotations

from typing import Dict, Optional
import numpy as np

from ..flowsheet_tools import UnitOp, EPS


class ContinuousEvaporatorAHA(UnitOp):
    """
    Continuous evaporator with:
      - Optional VLE "flash-like" evaporation at fixed T and P (representative of vacuum evaporator)
      - Floors to prevent crystallisation/dryout
      - Full AHA destruction
      - Built-in diagnostics printing

    Notes:
      * For a 50 °C, 70 torr vacuum evaporator, it is often more representative to
        compute evaporation from simple K-values (Psat/P) than to use fixed frac_vap.
      * This is still a simplified model (ideal K-values, no strong non-ideality).
    """

    def __init__(
        self,
        name: str,
        *,
        # Representative vacuum evaporator conditions
        T_set_K: float = 323.15,      # 50 °C
        P_set_torr: float = 70.0,     # 70 torr

        # If use_equilibrium=True, frac_vap is ignored for volatiles included in K model
        use_equilibrium: bool = True,

        # Legacy / fallback: fixed fractions evaporated, coupled by water
        frac_vap: Optional[Dict[str, float]] = None,
        coupled_by_water: bool = True,
        water_name: str = "H2O",

        # Floors
        min_liq_mol_s: Optional[Dict[str, float]] = None,

        # Thermo
        latent_kJ_per_mol: Optional[Dict[str, float]] = None,
        cp_liq_J_per_molK: float = 75.0,

        # Chemistry
        aha_name: str = "AHA",
        aha_products: Optional[Dict[str, float]] = None,
        aha_hno3_consumption: float = 0.0,

        # Diagnostics
        print_diagnostics: bool = False,
        print_every: int = 1,
    ):
        super().__init__(name)

        self.T_set_K = float(T_set_K)
        self.P_set_torr = float(P_set_torr)
        self.use_equilibrium = bool(use_equilibrium)

        # Fixed-fraction evaporation inputs (legacy mode)
        self.frac_vap = {k: float(v) for k, v in (frac_vap or {}).items()}
        self.coupled_by_water = coupled_by_water
        self.water_name = water_name

        self.min_liq_mol_s = {k: float(v) for k, v in (min_liq_mol_s or {}).items()}

        # Latent heats (kJ/mol) — defaults chosen for ~50 °C water; HNO3 left as prior default unless overridden
        self.latent_kJ_per_mol = {
            "H2O": 43.99,   # ~ΔHvap at 50 °C
            "HNO3": 39.0,   # keep prior default (user can override)
        }
        if latent_kJ_per_mol:
            self.latent_kJ_per_mol.update({k: float(v) for k, v in latent_kJ_per_mol.items()})

        self.cp_liq_J_per_molK = float(cp_liq_J_per_molK)

        self.aha_name = aha_name
        self.aha_products = aha_products or {"AcOH": 1.0, "N2O": 0.5}
        self.aha_hno3_consumption = float(aha_hno3_consumption)

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
    @staticmethod
    def _psat_torr_water(T_K: float) -> float:
        """
        Water saturation pressure using Antoine (valid ~1–100 °C):
            log10(P_mmHg) = A - B/(C + T_C)
        """
        T_C = T_K - 273.15
        A, B, C = 8.07131, 1730.63, 233.426
        return 10.0 ** (A - B / (C + T_C))  # mmHg ~ torr

    @staticmethod
    def _psat_torr_hno3(T_K: float) -> float:
        """
        Pure nitric acid vapor pressure correlation used in some VLE fits:
            log10(P_torr) = 8.25225 - 1918.3 / T_K
        """
        return 10.0 ** (8.25225 - 1918.3 / T_K)

    def _K_values(self, reg, idx, T_K: float, P_torr: float) -> np.ndarray:
        """
        Very simple K = Psat/P for selected volatiles.
        Everything else treated as non-volatile (K=0).
        """
        K = np.zeros(reg.n(), dtype=float)

        if "H2O" in idx:
            K[idx["H2O"]] = self._psat_torr_water(T_K) / max(P_torr, EPS)
        if "HNO3" in idx:
            K[idx["HNO3"]] = self._psat_torr_hno3(T_K) / max(P_torr, EPS)

        # If you later want NOx/HONO etc, add them here as additional volatiles.
        return K

    @staticmethod
    def _solve_beta(z: np.ndarray, K: np.ndarray) -> float:
        """
        Solve Rachford-Rice for vapor fraction beta in [0, 1].
        f(beta) = sum z_i*(K_i-1)/(1 + beta*(K_i-1)) = 0
        """
        ztot = float(np.sum(z))
        if ztot <= EPS:
            return 0.0

        # Normalize z to avoid scale issues
        z = z / ztot

        d = K - 1.0

        # If all K <= 1 -> no vapor
        if np.all(K <= 1.0 + 1e-12):
            return 0.0
        # If all K >= 1 and there are no nonvolatiles, beta -> 1;
        # but with nonvolatiles present, RR will give beta < 1.
        # We'll still solve on [0, 1-eps].
        lo, hi = 0.0, 1.0 - 1e-10

        def f(beta: float) -> float:
            denom = 1.0 + beta * d
            # avoid division by zero
            denom = np.maximum(denom, 1e-14)
            return float(np.sum(z * d / denom))

        flo = f(lo)
        fhi = f(hi)

        # If no sign change, clamp:
        if flo < 0.0 and fhi < 0.0:
            return 0.0
        if flo > 0.0 and fhi > 0.0:
            return float(hi)

        # Bisection
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            fmid = f(mid)
            if abs(fmid) < 1e-12:
                return mid
            if flo * fmid > 0:
                lo = mid
                flo = fmid
            else:
                hi = mid
        return 0.5 * (lo + hi)

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
        # Evaporation
        # ==============================================================
        x_vap = np.zeros(reg.n())
        x_liq = x_in.copy()
        self.vap_rate_mol_s = {}
        self.liq_rate_mol_s = {}

        if self.use_equilibrium:
            # --- Flash-like split based on K-values at (T_set, P_set) ---
            K = self._K_values(reg, idx, self.T_set_K, self.P_set_torr)
            beta = self._solve_beta(x_in, K)

            # Phase split
            d = K - 1.0
            denom = 1.0 + beta * d
            denom = np.maximum(denom, 1e-14)

            # liquid composition flow (not normalized): nL_i = z_i / denom
            nL = x_in / denom
            nV = x_in - nL  # by balance

            # Apply floors (prevent crystallisation/dryout)
            for sp, floor in self.min_liq_mol_s.items():
                if sp not in idx:
                    continue
                i = idx[sp]
                if nL[i] < floor:
                    # pull back from vapor into liquid to satisfy floor
                    deficit = floor - nL[i]
                    take = min(deficit, nV[i])
                    nL[i] += take
                    nV[i] -= take

            x_liq = nL
            x_vap = nV

            # For reporting, define a "theta-like" number based on water removal fraction
            w = self.water_name
            if w in idx:
                nw = x_in[idx[w]]
                vap_w = x_vap[idx[w]]
                self.theta = vap_w / max(nw, EPS)
            else:
                self.theta = 0.0

        else:
            # --- Legacy behaviour: coupled-by-water fixed fractions ---
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

        # Store rates dicts
        for sp, i in idx.items():
            v = float(x_vap[i])
            l = float(x_liq[i])
            if v > EPS:
                self.vap_rate_mol_s[sp] = v
            if l > EPS:
                self.liq_rate_mol_s[sp] = l

        # ==============================================================
        # Energy (simple sensible + latent)
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
        # -----------------------------
        # Guardrail: enforce non-negative component flows
        # while conserving total per-species: x_in = x_liq + x_vap
        # -----------------------------
        neg_v = x_vap < 0.0
        if np.any(neg_v):
            # move negative vapour amounts back into liquid
            x_liq[neg_v] = x_liq[neg_v] + x_vap[neg_v]
            x_vap[neg_v] = 0.0

        neg_l = x_liq < 0.0
        if np.any(neg_l):
            # move negative liquid amounts back into vapour (rare, but symmetric)
            x_vap[neg_l] = x_vap[neg_l] + x_liq[neg_l]
            x_liq[neg_l] = 0.0

        # final clip against tiny numerical noise
        x_liq[x_liq < 0.0] = 0.0
        x_vap[x_vap < 0.0] = 0.0

        s_liq.from_dense(x_liq)
        s_vap.from_dense(x_vap)

        s_liq.phase = "L"
        s_vap.phase = "G"
        s_liq.T = self.T_set_K
        s_vap.T = self.T_set_K

        # ==============================================================
        # Diagnostics
        # ==============================================================
        w = self.water_name
        nw_in = float(x_in[idx[w]]) if w in idx else 0.0
        nw_v = float(x_vap[idx[w]]) if w in idx else 0.0
        self._maybe_print(
            f"""
--- {self.name} ---
T = {self.T_set_K:.2f} K ({self.T_set_K-273.15:.1f} °C) | P = {self.P_set_torr:.1f} torr | eq={self.use_equilibrium}
θ (water removal frac proxy) = {self.theta:.3f}
Inlet:  F = {s_in.total_molar_flow():.4e} mol/s
Water:  in={nw_in:.4e} vap={nw_v:.4e} liq={nw_in - nw_v:.4e} floor={self.min_liq_mol_s.get(w,0.0):.2e}
HNO3:   in={float(x_in[idx['HNO3']]) if 'HNO3' in idx else 0.0:.4e} vap={self.vap_rate_mol_s.get('HNO3',0):.4e}
AHA destroyed: {self.aha_destroyed_mol_s:.3e} mol/s
Q_sens = {self.Q_sensible_kW:.3f} kW
Q_lat  = {self.Q_latent_kW:.3f} kW
Q_tot  = {self.Q_total_kW:.3f} kW
--------------------------------
"""
        )
