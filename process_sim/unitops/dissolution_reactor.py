from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Iterable, Tuple, List
import numpy as np

from ..flowsheet_tools import UnitOp, Stream


@dataclass(frozen=True)
class StoichReaction:
    name: str
    nu: Dict[str, float]


class DissolutionReactor(UnitOp):
    """
      - DOES NOT overwrite the acid inlet.
      - Checks sufficient HNO3 for full conversion of specified reactants.
      - UO2 split: a fixed fraction via R1, remainder via R2.
      - All other specified reactants convert to completion (extent = fuel reactant molar flow).
      - Noble gases to offgas; iodine 98% to offgas; NOx to offgas with produced NOx split
        as 41.5% NO2 / 58.5% NO.
    """

    def __init__(
        self,
        name: str,
        *,
        uo2_r1_frac: float = 0.18,
        iodine_total_species: str = "I",
        iodine_gas_species: str = "I_g",
        iodine_aq_species: str = "I_aq",
        iodine_offgas_frac: float = 0.98,
        noble_gases: Tuple[str, ...] = ("Xe", "Kr"),
        zircaloy_species: str = "Zircaloy",
        zircaloy_basket_mol_s: float = 0.0,
        zircaloy_escape_fines_mol_s: float = 0.0,
        fines_settle_frac: float = 0.20,
        extra_inlet_ports: Optional[Iterable[str]] = None,
        # Cleaning water inlet port (optional)
        cleaning_water_port: str = "cleaning_water",
        # NOx split (of produced NOx)
        nox_no2_frac: float = 0.415,
    ):
        super().__init__(name)

        if not (0.0 <= uo2_r1_frac <= 1.0):
            raise ValueError("uo2_r1_frac must be in [0,1]")
        self.uo2_r1_frac = float(uo2_r1_frac)

        self.iodine_total_species = iodine_total_species
        self.iodine_gas_species = iodine_gas_species
        self.iodine_aq_species = iodine_aq_species
        self.iodine_offgas_frac = float(iodine_offgas_frac)

        self.noble_gases = noble_gases

        self.zircaloy_species = zircaloy_species
        self.zirc_basket = float(zircaloy_basket_mol_s)
        self.zirc_escape = float(zircaloy_escape_fines_mol_s)
        if not (0.0 <= fines_settle_frac <= 1.0):
            raise ValueError("fines_settle_frac must be in [0,1]")
        self.fines_settle_frac = float(fines_settle_frac)

        self.extra_inlet_ports = tuple(extra_inlet_ports or ())

        self.cleaning_water_port = str(cleaning_water_port)

        if not (0.0 <= nox_no2_frac <= 1.0):
            raise ValueError("nox_no2_frac must be in [0,1]")
        self.nox_no2_frac = float(nox_no2_frac)

        self.reactions: List[StoichReaction] = [
            StoichReaction("R1", {"UO2": -1.0, "HNO3": -4.0, "UO2(NO3)2": 1.0, "NO2": 2.0, "H2O": 2.0}),
            StoichReaction("R2", {"UO2": -1.0, "HNO3": -(8.0 / 3.0), "UO2(NO3)2": 1.0, "NO2": (2.0 / 3.0), "H2O": (4.0 / 3.0)}),
            StoichReaction("R3", {"PuO2": -1.0, "HNO3": -4.0, "Pu(NO3)4": 1.0, "H2O": 2.0}),
            StoichReaction("R4", {"NpO2": -1.0, "HNO3": -4.0, "Np(NO3)4": 1.0, "H2O": 2.0}),
            StoichReaction("R5", {"AmO2": -1.0, "HNO3": -4.0, "Am(NO3)4": 1.0, "H2O": 2.0}),
            StoichReaction("R6", {"Cs2O": -1.0, "HNO3": -2.0, "CsNO3": 2.0, "H2O": 1.0}),
            StoichReaction("R7", {"SrO": -1.0, "HNO3": -2.0, "Sr(NO3)2": 1.0, "H2O": 1.0}),
            StoichReaction("R8", {"Nd2O3": -1.0, "HNO3": -6.0, "Nd(NO3)3": 2.0, "H2O": 3.0}),
            StoichReaction("R9", {"Sm2O3": -1.0, "HNO3": -6.0, "Sm(NO3)3": 2.0, "H2O": 3.0}),
            StoichReaction("R10", {"Eu2O3": -1.0, "HNO3": -6.0, "Eu(NO3)3": 2.0, "H2O": 3.0}),
            StoichReaction("R11", {"Gd2O3": -1.0, "HNO3": -6.0, "Gd(NO3)3": 2.0, "H2O": 3.0}),
            StoichReaction("R12", {"Tc": -1.0, "HNO3": -7.0, "HTcO4": 1.0, "NO2": 7.0, "H2O": 3.0}),
        ]

    def _build_Nu(self, reg) -> np.ndarray:
        idx = reg.index()
        nsp = reg.n()
        Nu = np.zeros((nsp, len(self.reactions)), dtype=float)
        for r, rxn in enumerate(self.reactions):
            for sp, nu in rxn.nu.items():
                if sp not in idx:
                    raise KeyError(f"{self.name}: species '{sp}' not present in registry")
                Nu[idx[sp], r] = float(nu)
        return Nu

    def _extents_from_fuel(self, fuel: Stream) -> np.ndarray:
        reg = fuel.reg
        idx = reg.index()
        F_fuel = fuel.to_dense()

        xi = np.zeros(len(self.reactions), dtype=float)

        # UO2 split across R1/R2
        if "UO2" in idx:
            n_uo2 = float(F_fuel[idx["UO2"]])
            xi[0] = self.uo2_r1_frac * n_uo2          # R1
            xi[1] = (1.0 - self.uo2_r1_frac) * n_uo2  # R2

        # Full conversion for the rest
        rmap = {
            2: "PuO2",
            3: "NpO2",
            4: "AmO2",
            5: "Cs2O",
            6: "SrO",
            7: "Nd2O3",
            8: "Sm2O3",
            9: "Eu2O3",
            10: "Gd2O3",
            11: "Tc",
        }
        for r, reactant in rmap.items():
            if reactant in idx:
                xi[r] = float(F_fuel[idx[reactant]])

        return xi

    def required_hno3_stoich(self, fuel: Stream) -> float:
        xi = self._extents_from_fuel(fuel)
        required = 0.0
        for r, rxn in enumerate(self.reactions):
            required += (-rxn.nu.get("HNO3", 0.0)) * xi[r]
        return float(required)

    def apply(self) -> None:
        fuel = self.inlets["fuel"]
        acid = self.inlets["acid"]

        extras = [self.inlets[p] for p in self.extra_inlet_ports if p in self.inlets]
        cleaning = self.inlets.get(self.cleaning_water_port, None)  # NEW

        offgas = self.outlets["offgas"]
        solids = self.outlets["solids"]
        aq = self.outlets["aq"]

        reg = fuel.reg
        idx = reg.index()
        nsp = reg.n()

        F_fuel = fuel.to_dense()
        F_acid = acid.to_dense()

        F_extra = np.zeros(nsp, dtype=float)
        for s in extras:
            F_extra += s.to_dense()
        if cleaning is not None:
            F_extra += cleaning.to_dense()  # NEW

        Nu = self._build_Nu(reg)
        xi = self._extents_from_fuel(fuel)

        required_hno3 = 0.0
        for r, rxn in enumerate(self.reactions):
            required_hno3 += (-rxn.nu.get("HNO3", 0.0)) * xi[r]

        # Count HNO3 in acid + extras + cleaning (safe if cleaning is pure water)
        n_hno3_in = acid.get("HNO3") + sum(s.get("HNO3") for s in extras)
        if cleaning is not None:
            n_hno3_in += cleaning.get("HNO3")

        if n_hno3_in + 1e-18 < required_hno3:
            raise ValueError(
                f"{self.name}: insufficient HNO3. Need >= {required_hno3:.6e} mol/s, "
                f"but inlet(s) have {n_hno3_in:.6e} mol/s."
            )

        delta = Nu @ xi
        total_after = F_fuel + F_acid + F_extra + delta

        # Produced NO2 from stoich (used as "produced NOx" basis for split)
        n_no2_prod = max(float(delta[idx["NO2"]]), 0.0) if "NO2" in idx else 0.0

        target_aq = total_after.copy()
        target_off = np.zeros(nsp, dtype=float)
        target_sol = np.zeros(nsp, dtype=float)

        # Noble gases to offgas
        for sp in self.noble_gases:
            if sp in idx:
                j = idx[sp]
                target_off[j] += target_aq[j]
                target_aq[j] = 0.0

        # Iodine bookkeeping
        nI_total = 0.0
        if self.iodine_total_species in idx:
            nI_total += target_aq[idx[self.iodine_total_species]]
            target_aq[idx[self.iodine_total_species]] = 0.0
        if self.iodine_aq_species in idx:
            nI_total += target_aq[idx[self.iodine_aq_species]]
            target_aq[idx[self.iodine_aq_species]] = 0.0
        if self.iodine_gas_species in idx:
            nI_total += target_aq[idx[self.iodine_gas_species]]
            target_aq[idx[self.iodine_gas_species]] = 0.0

        if nI_total > 0.0:
            if self.iodine_gas_species not in idx or self.iodine_aq_species not in idx:
                raise KeyError(f"{self.name}: iodine gas/aq species not present in registry")
            target_off[idx[self.iodine_gas_species]] += self.iodine_offgas_frac * nI_total
            target_aq[idx[self.iodine_aq_species]] += (1.0 - self.iodine_offgas_frac) * nI_total

        # Move any NO/NO2 present (inlets or produced) to offgas
        if "NO2" in idx:
            j = idx["NO2"]
            target_off[j] += target_aq[j]
            target_aq[j] = 0.0
        if "NO" in idx:
            j = idx["NO"]
            target_off[j] += target_aq[j]
            target_aq[j] = 0.0

        # Apply NOx split to NOx PRODUCED by reactions:
        # basis = produced NO2 moles from stoich (treat as total produced NOx)
        if n_no2_prod > 0.0:
            if "NO2" not in idx or "NO" not in idx:
                raise KeyError(f"{self.name}: need 'NO2' and 'NO' in registry for NOx split")

            j_no2 = idx["NO2"]
            j_no = idx["NO"]

            # remove "as-produced" NO2 from offgas, then re-add as NO2/NO split
            target_off[j_no2] = max(target_off[j_no2] - n_no2_prod, 0.0)

            n_no2 = self.nox_no2_frac * n_no2_prod
            n_no = (1.0 - self.nox_no2_frac) * n_no2_prod
            target_off[j_no2] += n_no2
            target_off[j_no] += n_no

        # Zircaloy handling (same as before)
        if self.zircaloy_species in idx:
            j = idx[self.zircaloy_species]
            target_aq[j] = 0.0
            target_off[j] = 0.0
            target_sol[j] = self.zirc_basket + self.fines_settle_frac * self.zirc_escape
            target_aq[j] = (1.0 - self.fines_settle_frac) * self.zirc_escape

        offgas.from_dense(target_off)
        solids.from_dense(target_sol)
        aq.from_dense(target_aq)
