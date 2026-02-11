from __future__ import annotations


class InputParameters:
    # -------------------------------------------------------------------------
    # GLOBAL DEFAULTS
    # -------------------------------------------------------------------------
    DEFAULT_T_K = 298.15
    DEFAULT_P_PA = 100000

    # -------------------------------------------------------------------------
    # FEED COMPOSITION (mol/s) - FIXED
    # -------------------------------------------------------------------------
    FEED_COMPOSITION = {
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
    }

    # -------------------------------------------------------------------------
    # ATOMIC WEIGHTS (g/mol) - for formula-based MW if used
    # -------------------------------------------------------------------------
    ATOMIC_WEIGHTS = {
        "U": 237.969, "Pu": 239.757, "Np": 237.000, "Am": 241.000,
        "Sr": 90.000, "Tc": 99.000, "Cs": 135.010, "Nd": 144.997,
        "Sm": 149.187, "Eu": 153.000, "Gd": 156.000,
        "Xe": 131.293, "Kr": 83.798, "I": 129.000,
        "O": 15.999, "H": 1.008, "N": 14.007, "Zr": 91.224,
        "C": 12.011, "P": 30.974,
    }

    # -------------------------------------------------------------------------
    # EXOGENOUS STREAM SPECS
    # -------------------------------------------------------------------------
    EXOGENOUS_STREAM_SPECS = {
        "F1": {"T": 298.15, "p": 1e5, "phase": "S"},  # Fuel feed
        "F2_M": {"T": 298.15, "p": 1e5, "phase": "L"},  # Nitric acid feed
        "F2_B": {"T": 298.15, "p": 1e5, "phase": "L"},  # Cleaning water feed
        "F16": {"T": 298.15, "p": 1e5, "phase": "L"},  # TBP make-up
        "F9": {"T": 298.15, "p": 1e5, "phase": "L"},  # AHA feed
        "F11": {"T": 298.15, "p": 1e5, "phase": "L"},  # Dilute nitric acid
    }

    # -------------------------------------------------------------------------
    # UNIT / VESSEL OPERATING CONDITIONS
    # -------------------------------------------------------------------------
    UNIT_CONDITIONS = {
        "R101_Dissolver": {"T": 373.15, "p": 1e5},
        "V101_Liquid_Buffer": {"T": 298.15, "p": 1e5},
        "X101_TBP_Extraction": {"T": 298.15, "p": 1e5},
        "X102_AHA_Strip": {"T": 298.15, "p": 1e5},
        "X103_Final_Strip": {"T": 323.15, "p": 1e5},
        "V133_SolventPurge": {"T": 298.15, "p": 1e5},
        "T101_Evaporator": {"T": 388.15, "p": 1e5},
        "T102_Evaporator": {"T": 388.15, "p": 1e5},
        "M201_EvapOverheadMixer": {"T": 388.15, "p": 1e5},
    }

    # -------------------------------------------------------------------------
    # DISSOLVER ACID FEED TARGET (F2)
    # -------------------------------------------------------------------------
    R1_F2_TARGET_HNO3_MOL_S = 2.936403
    R1_F2_TARGET_H2O_MOL_S = 6.571934

    # -------------------------------------------------------------------------
    # DISSOLVER PARAMETERS
    # -------------------------------------------------------------------------
    # Iodine & noble-gas behaviour
    R1_IODINE_OFFGAS_FRACTION = 0.98
    R1_NOBLE_GASES = ("Xe", "Kr")

    # UO2 pathway split: fraction via R1, remainder via R2
    R1_UO2_R1_FRACTION = 0.18

    # NOx speciation of NOx produced: 41.5% as NO2, remainder as NO
    R1_NOX_NO2_FRACTION = 0.415

    # Zircaloy handling
    R1_ZIRCALOY_BASKET_MOL_S = 0.131325068 - 6.234e-3
    R1_ZIRCALOY_ESCAPE_FINES_MOL_S = 6.234e-3
    R1_FINES_SETTLE_FRACTION = 0.0 #All zircaloy that escapes the basket goes into the downstream process

    # Species naming for iodine bookkeeping (match your reactor defaults unless you changed them)
    R1_IODINE_TOTAL_SPECIES = "I"
    R1_IODINE_GAS_SPECIES = "I_g"
    R1_IODINE_AQ_SPECIES = "I_aq"

    # Cleaning water time averaged flow rate
    F2B_flow_rate = 1.034165

    # -------------------------------------------------------------------------
    #  X1 TBP EXTRACTION PARAMETERS
    # -------------------------------------------------------------------------
    X1_DISTRIBUTION_COEFFICIENTS = {
        "UO2(NO3)2": 40.0,
        "Pu(NO3)4": 4.0,
        "Np(NO3)4": 3.5,
        "HTcO4": 0.5,
    }

    # minimum recovery-to-organic targets used by design_N_and_OA()
    X1_REQUIRED_RECOVERY = {
        "UO2(NO3)2": 0.999,
    }

    X1_N_MAX = 50
    X1_OA_MAX = 1000.0
    X1_N_BALANCE_CAP = 20
    X1_OA_MIN = 0.5  # So the optimiser does not use incredibly low flow rates

    X1_TBP_WT_FRACTION = 0.30
    X1_TBP_NAME = "TBP"
    X1_DILUENT_NAME = "Dodecane"

    # -------------------------------------------------------------------------
    # X2 AHA STRIP FEED (fixed ratios via molarities)
    # -------------------------------------------------------------------------
    X2_STRIP_VDOT_L_S = 1.0  # L/s (scales the total strip flow)
    X2_STRIP_MOLARITIES = {"AHA": 0.25, "HNO3": 0.5}  # mol/L
    X2_STRIP_WATER_MOLARITY = 55.5  # mol/L (pure water basis)
    X2_TARGET_AO = 0.5
    X2_AO_MIN = 0.5 #So the optimiser does not use incredibly low flow rates

    # -------------------------------------------------------------------------
    # X2 STRIP DISTRIBUTION COEFFICIENTS (org/aq)
    # -------------------------------------------------------------------------
    X2_DISTRIBUTION_COEFFICIENTS = {
        "UO2(NO3)2": 40.0,
        "Pu(NO3)4": 1e-2,
        "Np(NO3)4": 1e-3,  # assume Np(V)
        "HTcO4": 0.5,
    }

    # -------------------------------------------------------------------------
    # X2 REQUIRED RECOVERY TO AQUEOUS (percent or fraction)
    # -------------------------------------------------------------------------
    X2_REQUIRED_RECOVERY_TO_AQ = {
        "Pu(NO3)4": 99.9,
        "Np(NO3)4": 99.9,
    }

    X2_N_MAX = 50
    X2_AO_MAX = 1000.0
    X2_N_BALANCE_CAP = 20

    X2_NONTRANSFER_KEEP_IN_AQ = ["H2O", "HNO3", "AHA"]
    X2_NONTRANSFER_KEEP_IN_ORG = ["TBP", "Dodecane"]

    # -------------------------------------------------------------------------
    # X3 FINAL STRIP (dilute HNO3) – fixed composition via molarity ratios
    # -------------------------------------------------------------------------
    X3_STRIP_VDOT_L_S = 1.0
    X3_STRIP_MOLARITIES = {"HNO3": 0.01}  # mol/L
    X3_STRIP_WATER_MOLARITY = 55.5  # mol/L

    # D values are org/aq
    X3_DISTRIBUTION_COEFFICIENTS = { # Pessimistic values - find in more detail later
        "HTcO4": 0.93,
        "UO2(NO3)2": 0.25,
    }

    # Design targets for sizing A/O (pick something; adjust as desired)
    # NOTE: must be 0–1, not percent.
    X3_REQUIRED_RECOVERY_TO_AQ = {
        "UO2(NO3)2": 0.999,
    }

    X3_N_MAX = 50
    X3_AO_MAX = 1e4
    X3_AO_MIN = 0.5  # So the optimiser does not use incredibly low flow rates

    # -------------------------------------------------------------------------
    # Solvent purge fraction after X3 organic outlet
    # -------------------------------------------------------------------------
    PURGE_FRACTION_SOLVENT_LOOP = 0.05

    # -------------------------------------------------------------------------
    # EVAPORATORS T101 / T102 (waste concentration + AHA destruction)
    # -------------------------------------------------------------------------
    EVAP_T_K = 388.15  # 115C

    # Nominal vapour split fractions (H2O controls the coupled extent theta)
    T101_FRAC_VAP = {"H2O": 0.85, "HNO3": 0.60}
    T102_FRAC_VAP = {"H2O": 0.85, "HNO3": 0.60}

    # Floors in liquid product (anti-crystallisation guardrails)
    T101_MIN_LIQ_MOL_S = {"H2O": 0.25, "HNO3": 0.01}
    T102_MIN_LIQ_MOL_S = {"H2O": 0.05, "HNO3": 1e-4}

    # Simple energetics
    EVAP_LATENT_KJ_PER_MOL = {"H2O": 40.7, "HNO3": 39.0}
    EVAP_CP_LIQ_J_PER_MOLK = 75.0

    # AHA destruction (lumped). Full conversion each pass.
    EVAP_AHA_PRODUCTS = {"AcOH": 1.0, "N2O": 0.5}
    EVAP_AHA_HNO3_CONSUMPTION = 0.0

    # Coupling behaviour
    EVAP_COUPLED_BY_WATER = True
    EVAP_WATER_NAME = "H2O"
    # Diagnostics
    EVAP_PRINT_DIAGNOSTICS = True  # master on/off
    EVAP_PRINT_EVERY = 25  # print every N iterations
    EVAP_PRINT_STREAMS = True  # print inlet/outlet summaries





