from __future__ import annotations

"""
Configuration for the simplified PUREX-like process model.

Design intent
-------------
- Keep the simplified package easy to edit.
- Keep recycle-dependent *sizing logic* aligned with the old_complex package.
- Use fixed split fractions / conversions for the simplified units, with values
  taken from the old_complex base-case where appropriate.
- Where the PFD contains more detailed equipment than the old_complex package,
  represent that section with pass-through / dummy units that can be refined
  later without changing the recycle structure.
"""

# ---------------------------------------------------------------------------
# CHEMISTRY / MW SUPPORT
# ---------------------------------------------------------------------------
ATOMIC_WEIGHTS: dict[str, float] = {
    "U": 238.029, "Pu": 239.052, "Np": 237.048, "Am": 241.057,
    "Sr": 87.620, "Tc": 98.906, "Cs": 132.905, "Nd": 144.242,
    "Sm": 150.360, "Eu": 151.964, "Gd": 157.250,
    "Xe": 131.293, "Kr": 83.798, "I": 126.904,
    "O": 15.999, "H": 1.008, "N": 14.007, "Zr": 91.224,
    "C": 12.011, "P": 30.974,
}

SPECIES_FORMULAS: dict[str, str] = {
    # fuel / solids
    "UO2": "UO2",
    "PuO2": "PuO2",
    "NpO2": "NpO2",
    "AmO2": "AmO2",
    "Cs2O": "Cs2O",
    "SrO": "SrO",
    "Nd2O3": "Nd2O3",
    "Sm2O3": "Sm2O3",
    "Eu2O3": "Eu2O3",
    "Gd2O3": "Gd2O3",
    "Tc": "Tc",
    "Xe": "Xe",
    "Kr": "Kr",
    "I": "I",
    "I_g": "I2",
    "I_aq": "I",
    "Zircaloy": "Zr",
    # common process species
    "H2O": "H2O",
    "HNO3": "HNO3",
    "NO2": "NO2",
    "NO": "NO",
    "N2": "N2",
    "O2": "O2",
    "NH3": "NH3",
    "N2O": "N2O",
    "TBP": "C12H27O4P",
    "Dodecane": "C12H26",
    "AHA": "C2H5NO2",
    "AcOH": "C2H4O2",
    # dissolved nitrate salts
    "UO2(NO3)2": "UO2(NO3)2",
    "Pu(NO3)4": "Pu(NO3)4",
    "Np(NO3)4": "Np(NO3)4",
    "Am(NO3)4": "Am(NO3)4",
    "CsNO3": "CsNO3",
    "Sr(NO3)2": "SrNO32",
    "Nd(NO3)3": "Nd(NO3)3",
    "Sm(NO3)3": "Sm(NO3)3",
    "Eu(NO3)3": "Eu(NO3)3",
    "Gd(NO3)3": "Gd(NO3)3",
    "HTcO4": "HTcO4",
    # products
    "UO3": "UO3",
    "U3O8": "U3O8",
}

SPECIES_MW_OVERRIDE: dict[str, float] = {
    "Zircaloy": 91.224,
    "Sr(NO3)2": 211.630,
}

# ---------------------------------------------------------------------------
# EXOGENOUS FEEDS
# ---------------------------------------------------------------------------
EXOGENOUS_STREAMS: dict[str, dict] = {
    # Fuel feed to the effective continuous dissolver representation
    "F1": {
        "T": 298.15, "p": 1e5, "phase": "S",
        "mol": {
            "UO2": 0.810614973,
            "PuO2": 0.00765862,
            "NpO2": 0.000403571,
            "AmO2": 0.000321834,
            "Cs2O": 0.004272278,
            "SrO": 0.002094222,
            "Nd2O3": 0.00627727,
            "Sm2O3": 0.001264957,
            "Eu2O3": 0.000183208,
            "Gd2O3": 0.000106239,
            "Tc": 0.003139914,
            "Xe": 0.013807779,
            "Kr": 0.002438652,
            "I": 0.000552969,
            "Zircaloy": 0.131325068,
        },
    },
    # Wash / cleaning water to dissolver train
    "F2B": {"T": 298.15, "p": 1e5, "phase": "L", "mol": {"H2O": 1.034165}},
    # Recycle-sized feeds
    "F2M": {"T": 298.15, "p": 1e5, "phase": "L", "mol": {"HNO3": 0.5, "H2O": 0.2}},
    "F2W": {"T": 298.15, "p": 1e5, "phase": "L", "mol": {"H2O": 0.1}},
    "F16": {"T": 298.15, "p": 1e5, "phase": "L", "mol": {"TBP": 0.1, "Dodecane": 0.5}},
    "F9": {"T": 298.15, "p": 1e5, "phase": "L", "mol": {"AHA": 0.25, "HNO3": 0.5, "H2O": 1.0}},
    "F11": {"T": 298.15, "p": 1e5, "phase": "L", "mol": {"HNO3": 0.01, "H2O": 1.0}},
    "F40": {"T": 298.15, "p": 1e5, "phase": "G", "mol": {"O2": 0.21 * 2.135616222, "N2": 0.79 * 2.135616222}},
    "F43": {"T": 298.15, "p": 1e5, "phase": "G", "mol": {"NH3": 0.05}},
    "F47": {"T": 298.15, "p": 1e5, "phase": "L", "mol": {"H2O": 1.0}},
    "F52A": {"T": 298.15, "p": 1e5, "phase": "G", "mol": {"O2": 0.021, "N2": 0.079}},
    "F52B": {"T": 298.15, "p": 1e5, "phase": "G", "mol": {"O2": 0.021, "N2": 0.079}},
    # Updated final reactor-section air feed from the revised reactor stream table
    "F73": {"T": 414.0, "p": 1.0e5, "phase": "G", "mol": {"O2": 1.196823, "N2": 2.982937}},
    # Utility/service-side helper feeds used to reproduce the revised EV-103 train stream table
    "U67_src": {"T": 393.15, "p": 1.0e5, "phase": "L", "mol": {"H2O": 16.222300}},
}

# ---------------------------------------------------------------------------
# RECYCLE-DEPENDENT SIZING TARGETS (kept aligned with old_complex)
# ---------------------------------------------------------------------------
ACID_TARGETS = {"hno3_mol_s": 2.936403, "h2o_mol_s": 6.571934}

ACID_STOCK = {
    "hno3_vol_frac": 0.70,
    "h2o_vol_frac": 0.30,
    "mw_hno3": 63.012,
    "rho_hno3_g_ml": 1.51,
    "water_molarity": 55.5,
}
ACID_STOCK["hno3_molarity_pure"] = ACID_STOCK["rho_hno3_g_ml"] * 1000.0 / ACID_STOCK["mw_hno3"]
ACID_STOCK["water_per_hno3_mol"] = (
    (ACID_STOCK["h2o_vol_frac"] * ACID_STOCK["water_molarity"]) /
    (ACID_STOCK["hno3_vol_frac"] * ACID_STOCK["hno3_molarity_pure"])
)

SOLVENT = {"TBP_name": "TBP", "diluent_name": "Dodecane", "TBP_wt_frac": 0.30}

X101_SIZING = {"OA_target": 0.777772895903665, "OA_min": 0.5}
X102_SIZING = {
    "AO_target": 0.1,
    "water_molarity": 55.5,
    "molarities": {"AHA": 0.25, "HNO3": 0.5},
}
X103_SIZING = {
    "AO_target": 2.5,
    "water_molarity": 55.5,
    "molarities": {"HNO3": 0.01},
}

# ---------------------------------------------------------------------------
# DISSOLVER TRAIN / EFFECTIVE CONTINUOUS REPRESENTATION
# ---------------------------------------------------------------------------
DISSOLVER = {
    "uo2_r1_frac": 0.18,
    "iodine_total_species": "I",
    "iodine_gas_species": "I_g",
    "iodine_aq_species": "I_aq",
    "iodine_offgas_frac": 0.98,
    "noble_gases": ("Xe", "Kr"),
    "zircaloy_species": "Zircaloy",
    "zircaloy_basket_mol_s": 0.131325068 - 6.234e-3,
    "zircaloy_escape_fines_mol_s": 6.234e-3,
    "fines_settle_frac": 0.0,
    "nox_no2_frac": 0.415,
}

# ---------------------------------------------------------------------------
# FIXED SPLITS / CONVERSIONS FOR SIMPLIFIED UNITS
# Fractions are chosen to reproduce the old_complex base case for the current
# nominal feed / recycle state. They are easy to edit later.
# ---------------------------------------------------------------------------
X101 = {
    "frac_to_organic": {
        "TBP": 0.99871427976687,
        "Dodecane": 1.0,
        "UO2(NO3)2": 0.9989999973493199,
        "Pu(NO3)4": 0.9274835914098941,
        "Np(NO3)4": 0.9101737719162379,
        "HTcO4": 0.350657897312383,
    },
    "default_frac_to_organic": 0.0,
}

CF101 = {
    "frac_removed": {"TBP": 0.9595917602494292, "Dodecane": 0.0},
    "default_frac_removed": 0.0,
}

E101 = {
    "frac_to_vapour": {"H2O": 0.02206515731274, "HNO3": 0.5621469949705132},
    "default_frac_to_vapour": 0.0,
}

X102 = {
    "frac_to_aqueous": {
        "AHA": 1.0,
        "H2O": 1.0,
        "HNO3": 1.0,
        "Pu(NO3)4": 1.0,
        "Np(NO3)4": 1.0,
        "HTcO4": 0.19999795199475706,
        "UO2(NO3)2": 0.0024999999999999363,
        "TBP": 1.0e-4,
        "Dodecane": 0.0,
    },
    "default_frac_to_aqueous": 0.0,
}

CF102 = {"frac_removed": {"TBP": 0.86238997960775}, "default_frac_removed": 0.0}

E102 = {
    "frac_to_vapour": {"H2O": 0.010498937973213883, "HNO3": 0.3764594792366077},
    "default_frac_to_vapour": 0.0,
}

X103 = {
    "frac_to_aqueous": {
        "H2O": 1.0,
        "HNO3": 1.0,
        "UO2(NO3)2": 0.9999999099999991,
        "HTcO4": 0.999380673704373,
        "TBP": 0.0024999999999999294,
        "Dodecane": 0.0,
    },
    "default_frac_to_aqueous": 0.0,
}

CF103 = {"frac_removed": {"TBP": 0.9952853475981924}, "default_frac_removed": 0.0}

# Explicit EV-103A/B/C break-out.
#
# Intended simplified concentration logic
# --------------------------------------
# EV-103A
#   - main flash/concentration step
#   - removes most volatile H2O/HNO3 to F65 (sent to the hot-vapour header)
#   - keeps essentially all U in the liquid train
#   - bleeds a small heel of entrained organic / trace impurities to F69
# EV-103B
#   - hydraulic / internal stage split of the concentrated liquor into two feeds
#     to EV-103C (F66 and F68), with a small heel/waste bleed F72
#   - the upper feed F66 is intentionally lighter (more H2O/HNO3, less U)
#     than the lower feed F68
# EV-103C
#   - final concentration / product draw stage
#   - sends essentially all uranium nitrate to product F24
#   - sends a small mother-liquor / impurity bleed to F70
#
# The combined effect is kept close to the old_complex single-E103 behaviour:
# most H2O/HNO3 is removed from the U-product stream, but the staged PFD
# topology is now represented explicitly and the per-species stage routing is
# easy to edit here.

EV103A = {
    # Finalised performance data from the updated reactor / evaporator stream table.
    # Visible split:
    #   F57T -> F65 + F64_raw
    # Water is flashed here, while trace non-volatiles remain in the visible
    # liquid path so TBP handling is compositionally consistent through EV-103A/B/C.
    "frac_to_A": {
        "H2O": 6.340412 / 22.007472,
    },
    "default_frac_to_A": 0.0,
}

EV103B = {
    # Finalised visible split:
    #   F64 -> F66 + F68
    # F66 is the light side stream; F68 is the heavy U-bearing stream.
    # TBP is split visibly here so it remains traceable through EV-103B,
    # while HNO3 / HTcO4 continue to follow the fitted auxiliary route.
    "frac_to_A": {
        "H2O": 6.660651 / 15.667060,
        "TBP": 6.660651 / 15.667060,
    },
    "default_frac_to_A": 0.0,
}

EV103C_STAGE = {
    # Source-sensitive final evaporator stage.
    # Light feed F66 is mostly visible waste F70, but split TBP is retained to
    # product for compositionally consistent bookkeeping.
    # Heavy feed F68 is split to product F24 and vapour F25.
    "heavy_to_product": {
        "H2O": 2.206324,
        "UO2(NO3)2": 0.807780,
        "TBP": 0.000011643384884813555,
    },
    "heavy_to_vapour": {
        "H2O": 6.803606,
    },
    "light_to_product": {
        "TBP": 0.000008610815115186445,
    },
    "aux_to_product": {
        "HNO3": 0.000396602,
        "HTcO4": 0.000880431,
    },
    # Small fitted water generation term needed to reproduce the finalised
    # EV-103C outlet table exactly.
    "extra_h2o_to_vapour": 0.003520,
}

# Auxiliary revised reactor-section water/service streams
REACTOR_AUX = {
    "F67_H2O_mol_s": 16.222300,
    "F69_H2O_mol_s": 16.222300,
    "F72_H2O_mol_s": 6.340412,
}

R201 = {
    # Finalised reactor-section performance fit.
    "no2_per_u": 2.0,
    "o2_per_u": 0.0,
    "h2o_to_offgas_frac": 2.200747 / 2.206324,
    "route_hno3_to_offgas": 0.0,
    "route_tbp_to_offgas": 0.0,
    "route_htco4_to_offgas": 0.0,
}
R202 = {
    "reactions": [
        {
            "reactant": "UO3",
            "conversion": 1.0,
            "products": {"U3O8": 1.0 / 3.0, "O2": 1.0 / 6.0},
            "outlet_solid": ["U3O8"],
            "outlet_gas": ["O2"],
        }
    ],
    "default_to_solid": True,
}

K301 = {"frac_to_liquid": {"H2O": 1.0, "HNO3": 1.0}, "default_frac_to_liquid": 0.0}

D101 = {"NO2_capture_frac": 0.99, "NO_capture_frac": 0.95, "LG_water_per_gas": 6.0 / 2.68178}
E104 = {"target_hno3_mass_frac": 0.74}

EVAP_AHA = {"products": {"AcOH": 1.0, "N2O": 0.5}, "hno3_consumption": 0.0}
R103 = {"NOx_removal_frac": 0.9376, "NH3_excess_factor": 1.15}
KO101 = {"water_removal_frac": 1.0}
D102 = {"NH3_capture_frac": 0.99, "LG_water_per_gas": 6.77}
V104 = {"purge_frac": 0.05}
V105 = {"purge_frac": 0.05}
V133 = {"purge_frac": 0.05}

TSA = {
    "iodine_species": "I_g",
    "ads_capture_frac": 1.0,
    "regen_max_y_I2": 0.1,
    "regen_O2_frac": 0.21,
    "regen_N2_frac": 0.79,
    "active_column": "A",
    "t_ads_s": 8.0 * 3600.0,
    "t_regen_s": 4.0 * 3600.0,
}

# ---------------------------------------------------------------------------
# UNIT CONDITIONS
# ---------------------------------------------------------------------------
_EVAP_P_PA = 70.0 / 760.0 * 101325.0
UNIT_CONDITIONS = {
    # front end / dissolver train
    "M102_AcidMakeupMixer": {"T": 298.15, "p": 1e5},
    "V101_AcidTank": {"T": 298.15, "p": 1e5},
    "P107_AcidTransfer": {"T": 298.15, "p": 1.1e5},
    "E101_AcidHeater": {"T": 373.15, "p": 1e5},
    "DS101_104_DissolverTrain": {"T": 373.15, "p": 1e5},
    "V102_DissolverHold": {"T": 298.15, "p": 1e5},
    "P108_ToX101": {"T": 298.15, "p": 1.1e5},
    "V103_OffgasSurge": {"T": 298.15, "p": 1e5},
    "P104_OffgasBlower": {"T": 298.15, "p": 1.1e5},
    # solvent extraction section
    "M114_SolventMakeupMixer": {"T": 298.15, "p": 1e5},
    "X101_TBP_Extraction": {"T": 298.15, "p": 1e5},
    "CF101_Coalescer": {"T": 298.15, "p": 1e5},
    "E102_RaffinateEvaporator": {"T": 318.15, "p": _EVAP_P_PA},
    "X102_AHA_Strip": {"T": 298.15, "p": 1e5},
    "CF102_Coalescer": {"T": 298.15, "p": 1e5},
    "E103_RaffinateHeater": {"T": 318.15, "p": _EVAP_P_PA},
    "E104_PuNpFeedHeater": {"T": 318.15, "p": _EVAP_P_PA},
    "EV101_PuNpEvaporator": {"T": 318.15, "p": _EVAP_P_PA},
    "B101_OverheadBlower": {"T": 323.15, "p": 1e5},
    "E106_OrgStripPreheater": {"T": 343.15, "p": 1e5},
    "E105_AqStripPreheater": {"T": 343.15, "p": 1e5},
    "X103_Final_Strip": {"T": 343.15, "p": 1e5},
    "P109_UFeedPump": {"T": 343.15, "p": 1.1e5},
    "CF103_Coalescer": {"T": 379.20, "p": 0.8729e5},
    # uranium concentration / product section
    "E107_UFeedHeater": {"T": 379.20, "p": 0.8729e5},
    "EV103A_UConcentrator1": {"T": 383.15, "p": 1.0e5},
    "EV103B_UConcentrator2": {"T": 379.20, "p": 0.8729e5},
    "EV103C_UConcentrator3": {"T": 375.00, "p": 0.7526e5},
    "M270_EV103_OverheadMixer": {"T": 375.00, "p": 0.7526e5},
    "M271_EV103C_FeedMixer": {"T": 379.20, "p": 0.8729e5},
    "M272_EV103_WasteMixer": {"T": 375.00, "p": 0.7526e5},
    "P105_UProductPump": {"T": 375.00, "p": 1.0e5},
    "E109_UProductHeater": {"T": 475.00, "p": 1.0e5},
    "R101_ThermalReactor": {"T": 475.00, "p": 1.0e5},
    "P106_ReactorOffgasBlower": {"T": 475.00, "p": 1.1e5},
    "E111_ReactorOffgasCooler": {"T": 475.00, "p": 1.0e5},
    "E110_ProductCooler": {"T": 475.00, "p": 1.0e5},
    "K101_UProductSeparator": {"T": 650.00, "p": 1.0e5},
    "M301_HotVapourMixer": {"T": 373.15, "p": 1e5},
    "P101_HotVapourBlower": {"T": 373.15, "p": 1.1e5},
    "E112_HotVapourCooler": {"T": 298.15, "p": 1e5},
    "K301_Condenser": {"T": 298.15, "p": 1e5},
    # offgas treatment
    "M302_OffgasMixer": {"T": 298.15, "p": 1e5},
    "E113_OffgasCooler": {"T": 313.15, "p": 1e5},
    "D101_NO2_Absorber": {"T": 313.15, "p": 1e5},
    "M303_AcidMixer": {"T": 313.15, "p": 1e5},
    "E118_AcidTrimHeater": {"T": 313.15, "p": 1e5},
    "E104_HNO3_Concentrator": {"T": 373.15, "p": 1e5},
    "V104_AcidPurge": {"T": 298.15, "p": 1e5},
    "R103_NOxReduction": {"T": 673.15, "p": 1e5},
    "KO101_KnockoutPot": {"T": 323.15, "p": 1e5},
    "V105_KOWaterPurge": {"T": 323.15, "p": 1e5},
    "M105_WaterTrimMixer": {"T": 298.15, "p": 1e5},
    "E116_D102WaterHeater": {"T": 323.15, "p": 1e5},
    "D102_NH3_Absorber": {"T": 323.15, "p": 1e5},
    "V201_TSA_Split": {"T": 298.15, "p": 1e5},
    "TSA101A": {"T": 298.15, "p": 1e5},
    "TSA101B": {"T": 298.15, "p": 1e5},
    # solvent recycle
    "P110_SolventRecyclePump": {"T": 343.15, "p": 1.1e5},
    "M132_SolventRecoveryMixer": {"T": 298.15, "p": 1e5},
    "M133_SolventMixer": {"T": 298.15, "p": 1e5},
    "HX133_SolventCooler": {"T": 298.15, "p": 1e5},
    "V133_SolventPurge": {"T": 298.15, "p": 1e5},
}

# ---------------------------------------------------------------------------
# CONDITIONING FLAGS
# ---------------------------------------------------------------------------
CONDITION_FEEDS = {
    ("F2S", "DS101_104_DissolverTrain"): True,
    ("F53", "E102_RaffinateEvaporator"): True,
    ("F55", "EV101_PuNpEvaporator"): True,
    ("F57", "E107_UFeedHeater"): True,
    ("F30", "K301_Condenser"): True,
    ("F41", "E113_OffgasCooler"): True,
    ("F42", "R103_NOxReduction"): True,
    ("F43", "R103_NOxReduction"): True,
}

# ---------------------------------------------------------------------------
# SOLVER SETTINGS
# ---------------------------------------------------------------------------
SOLVER = {
    "tear_streams": ["F15", "F2R", "F45R"],
    "max_iter": 300,
    "tol": 1e-7,
    "relax": 0.45,
    "use_anderson": True,
    "anderson_m": 6,
    "verbose": True,
    "print_every": 10,
}
