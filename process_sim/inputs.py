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
        "TBP": 1e3, # Worst case scenario - low concentration of nitric acid and minimal metal salts
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
    X2_TARGET_AO = 0.1
    X2_AO_MIN = 0.1 #So the optimiser does not use incredibly low flow rates

    # -------------------------------------------------------------------------
    # X2 STRIP DISTRIBUTION COEFFICIENTS (org/aq)
    # -------------------------------------------------------------------------
    X2_DISTRIBUTION_COEFFICIENTS = {
        "UO2(NO3)2": 40.0,
        "Pu(NO3)4": 1e-4,
        "Np(NO3)4": 1e-6,  # assume Np(V)
        "HTcO4": 0.5,
        "TBP": 1e3,  # Worst case scenario - low concentration of nitric acid and minimal metal salts
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
    X2_NONTRANSFER_KEEP_IN_ORG = ["Dodecane"]

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
        "TBP": 1e3, # Worst case scenario - low concentration of nitric acid and minimal metal salts
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
    # EVAPORATORS E101 / E102 (vacuum: 50C, 70 torr)
    # -------------------------------------------------------------------------

    EVAP_T_K = 318.15  # 45C
    EVAP_USE_EQUILIBRIUM = True  # <-- NEW: flash-like split using K=Psat/P
    # EVAP_P_TORR already defined above (70.0)

    # frac_vap retained only for legacy mode if EVAP_USE_EQUILIBRIUM=False
    E101_FRAC_VAP = {"H2O": 0.85, "HNO3": 0.60}
    E102_FRAC_VAP = {"H2O": 0.85, "HNO3": 0.60}

    # Floors in liquid product (anti-crystallisation guardrails)
    E101_MIN_LIQ_MOL_S = {"H2O": 0.25, "HNO3": 0.01}
    E102_MIN_LIQ_MOL_S = {"H2O": 0.05, "HNO3": 1e-4}

    # Energetics
    EVAP_LATENT_KJ_PER_MOL = {"H2O": 43.99, "HNO3": 39.0}
    EVAP_CP_LIQ_J_PER_MOLK = 75.0

    # AHA destruction (lumped). Full conversion each pass.
    EVAP_AHA_PRODUCTS = {"AcOH": 1.0, "N2O": 0.5}
    EVAP_AHA_HNO3_CONSUMPTION = 0.0

    # Coupling behaviour (legacy mode only)
    EVAP_COUPLED_BY_WATER = True
    EVAP_WATER_NAME = "H2O"

    # Diagnostics
    EVAP_PRINT_DIAGNOSTICS = True
    EVAP_PRINT_EVERY = 25
    EVAP_PRINT_STREAMS = True

    # -------------------------------------------------------------------------
    # COALESCERS BETWEEN X101/X102 AND E101/E102 (TBP/DILUENT CARRYOVER CONTROL)
    # -------------------------------------------------------------------------
    COAL_TBP_NAME = "TBP"
    COAL_DILUENT_NAME = "Dodecane"
    COAL_WATER_NAME = "H2O"
    COAL_HNO3_NAME = "HNO3"
    COAL_WATER_MOLARITY_PURE = 55.5  # mol/L

    # Entrained-droplet removal efficiency (applied only to amount above solubility)
    COAL_REMOVAL_EFF = 0.999

    # Dissolved equilibrium solubilities (mol/L), ~ambient
    # TBP: 0.28–0.40 g/L => ~1.05e-3 to 1.50e-3 mol/L
    COAL_TBP_SOLUBILITY_WEAK_MOL_L = 1.5e-3
    COAL_TBP_SOLUBILITY_STRONG_MOL_L = 5.0e-4  # conservative lower dissolved TBP allowed in strong acid

    # Dodecane: extremely low water solubility (order 1e-8 mol/L)
    COAL_DIL_SOLUBILITY_WEAK_MOL_L = 2.0e-8
    COAL_DIL_SOLUBILITY_STRONG_MOL_L = 1.0e-8

    # "Strong acid" threshold for interpolation
    COAL_HNO3_STRONG_M = 3.0

    COAL_ASSUME_PURE_WATER_FOR_X102_AQ = True

    COAL_PRINT_DIAGNOSTICS = False

    COAL_X103_TO_E103_TBP_SOLUBILITY_MOL_L = 5.0e-5  # ~13 mg/L
    COAL_X103_TO_E103_REMOVAL_EFF = 0.9999

    # -------------------------------------------------------------------------
    # Fresh air feed for downstream offgas mixing
    # -------------------------------------------------------------------------
    F40_AIR_TOTAL_MOL_S = 2.135616222
    F40_AIR_O2_FRAC = 0.21
    F40_AIR_N2_FRAC = 0.79

    # -------------------------------------------------------------------------
    # E103 continuous evaporator @ 100C, 1 bar (species split fractions to vapour)
    # -------------------------------------------------------------------------
    E103_T_K = 373.15
    E103_P_PA = 1.0e5
    E103_FRAC_TO_VAP = {"H2O": 0.90, "HNO3": 0.90}

    # -------------------------------------------------------------------------
    # U nitrate -> UO3 reactor @ 383.15 K, 1 bar
    # -------------------------------------------------------------------------
    R201_T_K = 383.15
    R201_P_PA = 1.0e5
    R201_HYDRATE_N = 0.0

    # -------------------------------------------------------------------------
    # Hot-vapour mixing + condenser at 1 bar
    # -------------------------------------------------------------------------
    HOT_VAP_MIX_T_K = 373.15
    HOT_VAP_MIX_P_PA = 1.0e5

    COND_T_K = 298.15
    COND_P_PA = 1.0e5

    # -------------------------------------------------------------------------
    # Nitric acid recovery system
    # -------------------------------------------------------------------------

    D101_ABS_T_K = 313.15  # 40C
    D101_ABS_P_PA = 100000

    D101_NO2_CAPTURE_FRAC = 0.99
    D101_NO_CAPTURE_FRAC = 0.95
    D101_LG_WATER_PER_GAS = 6.0 / 2.68178  # mol water per mol gas

    E104_T_K = 373.15  # 100C
    E104_P_PA = 100000
    E104_TARGET_HNO3_MASS_FRAC = 0.74 # 45% mol fraction (taken from Kevin's table)

    # -------------------------------------------------------------------------
    # MAKEUP ACID SUPPLY: treat as 70 vol% HNO3 + 30 vol% water
    # with water molarity fixed at 55.5 mol/L
    # -------------------------------------------------------------------------
    ACID_STOCK_HNO3_VOL_FRAC = 0.70
    ACID_STOCK_H2O_VOL_FRAC = 0.30

    ACID_STOCK_MW_HNO3 = 63.012  # g/mol

    # Density of (near) anhydrous HNO3 @ 20 C ~ 1.51 g/cm3 = 1.51 g/mL
    ACID_STOCK_RHO_HNO3_G_ML = 1.51

    # Derived "pure liquid HNO3" molarity (mol/L)
    ACID_STOCK_HNO3_MOLARITY_PURE = (
            ACID_STOCK_RHO_HNO3_G_ML * 1000.0 / ACID_STOCK_MW_HNO3
    )

    # Water molarity in the aqueous fraction (mol/L)
    ACID_STOCK_WATER_MOLARITY = 55.5

    # Fixed mol water per mol HNO3 in the stock (used by sizing)
    ACID_STOCK_WATER_PER_HNO3_MOL = (
            (ACID_STOCK_H2O_VOL_FRAC * ACID_STOCK_WATER_MOLARITY) /
            (ACID_STOCK_HNO3_VOL_FRAC * ACID_STOCK_HNO3_MOLARITY_PURE)
    )

    # -------------------------------------------------------------------------
    # R103 NH3 NOx reduction column (SCR-like) on absorber offgas
    # -------------------------------------------------------------------------

    R103_T_K = 673.15  # 400C
    R103_P_PA = 1.0e5  # atmospheric
    R103_NOX_REMOVAL_FRAC = 0.9376  # applies to NO and NO2
    R103_NH3_EXCESS_FACTOR = 1.15  # Excess ammonia used (above the stoichiometric amount)

    # -------------------------------------------------------------------------
    # KO101 knock-out pot after R103 (remove all water)
    # -------------------------------------------------------------------------
    KO101_T_K = 323.15  # 50 C
    KO101_P_PA = 1.0e5  # 1 bar
    KO_WATER_PURGE_FRAC = 0.05  # 5% purge, 95% recycle

    # -------------------------------------------------------------------------
    # D102: NH3 absorption column (water wash) on F46
    # -------------------------------------------------------------------------
    D102_T_K = 323.15
    D102_P_PA = 1.0e5
    D102_NH3_CAPTURE_FRAC = 0.99

    D102_LG_WATER_PER_GAS = 6.77  # mol H2O per mol gas

    # -------------------------------------------------------------------------
    # TSA (EMM-17) IODINE POLISHING ON F49
    # -------------------------------------------------------------------------
    TSA_IODINE_NAME = "I_g"

    # Thesis-reported (order-of-magnitude) performance for EMM-17:
    TSA_EMM17_UPTAKE_G_PER_G = 0.42  # g I2 / g adsorbent (dynamic, ~150 ppmv tests)
    TSA_EMM17_PACKING_DENSITY_G_CM3 = 0.49  # g/cm^3
    TSA_MTZ_UTILIZATION = 0.70  # only 70% of bed utilised (finite MTZ)

    # Cycle assumptions (time-averaged TSA)
    TSA_T_ADS_H = 8.0
    TSA_T_REGEN_H = 4.0

    # Regeneration conditions
    TSA_REGEN_T_K = 423.15  # 150 C
    TSA_REGEN_P_PA = 1.0e5

    # "Minimise regen flow" via a max iodine mole fraction in regen outlet gas.
    # Choose a conservative bound to avoid unrealistic supersaturation/condensation.
    TSA_REGEN_MAX_Y_I2 = 0.1  # 10 mol% iodine max in regen exhaust

    # Which physical column is ADSORBING right now ("A" or "B")
    TSA_ACTIVE_COLUMN = "A"

    # Regen gas composition (AIR)
    TSA_REGEN_O2_FRAC = 0.21
    TSA_REGEN_N2_FRAC = 0.79

    # Safety factor on bed volume
    TSA_BED_VOL_SF = 1.10

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
        "F40": {"T": 298.15, "p": 1e5, "phase": "G"},  # Fresh air feed
        "F2_W": {"T": 298.15, "p": 1e5, "phase": "L"}, # Makeup water feed
        "F43": {"T": 298.15, "p": 1e5, "phase": "G"},  # pure NH3 feed
        "F47": {"T": 298.15, "p": 1e5, "phase": "L"},  # water wash into D102
        "F52A": {"T": 298.15, "p": 1e5, "phase": "L"},  # Regeneration feed in TSA system
        "F52B": {"T": 298.15, "p": 1e5, "phase": "L"},  # Regeneration feed in TSA system
    }

    # -------------------------------------------------------------------------
    # UNIT / VESSEL OPERATING CONDITIONS
    # -------------------------------------------------------------------------
    EVAP_P_TORR = 70.0
    EVAP_P_PA = 70.0 / 760.0 * 101325.0  # ≈ 9330 Pa

    UNIT_CONDITIONS = {
        # ------------------------------------------------------------------
        # Front-end dissolution & extraction
        # ------------------------------------------------------------------
        "R101_Dissolver": {"T": 373.15, "p": 1.0e5},
        "V101_Liquid_Buffer": {"T": 298.15, "p": 1.0e5},

        "X101_TBP_Extraction": {"T": 298.15, "p": 1.0e5},
        "CF101_Coalescer_X101_to_E101": {"T": 298.15, "p": 1.0e5},

        "E101_Evaporator": {"T": 323.15, "p": EVAP_P_PA},

        # ------------------------------------------------------------------
        # AHA strip + evaporation
        # ------------------------------------------------------------------
        "X102_AHA_Strip": {"T": 298.15, "p": 1.0e5},
        "CF102_Coalescer_X102_to_E102": {"T": 298.15, "p": 1.0e5},

        "E102_Evaporator": {"T": 323.15, "p": EVAP_P_PA},

        "M201_EvapOverheadMixer": {"T": 323.15, "p": EVAP_P_PA},

        # ------------------------------------------------------------------
        # Final strip + TBP polish
        # ------------------------------------------------------------------
        "X103_Final_Strip": {"T": 323.15, "p": 1.0e5},
        "CF103_Coalescer_X103_to_E103": {"T": 323.15, "p": 1.0e5},

        # ------------------------------------------------------------------
        # High-temperature evaporation & calcination
        # ------------------------------------------------------------------
        "E103_Evaporator100C": {"T": 373.15, "p": 1.0e5},

        "R201_UO3_Calciner": {"T": R201_T_K, "p": R201_P_PA},
        "R202_U3O8_Converter": {"T": 623.15, "p": 1.0e5},

        # ------------------------------------------------------------------
        # Hot vapour handling & condensation
        # ------------------------------------------------------------------
        "M301_HotVapMixer": {"T": HOT_VAP_MIX_T_K, "p": HOT_VAP_MIX_P_PA},

        "K301_Condenser": {"T": COND_T_K, "p": COND_P_PA},
        "M302_OffgasMixer": {"T": COND_T_K, "p": COND_P_PA},

        # ------------------------------------------------------------------
        # NOx handling & acid recovery
        # ------------------------------------------------------------------
        "V103_OffgasSurge": {"T": 298.15, "p": 1e5},
        "D101_NO2_Absorber": {"T": 313.15, "p": 1.0e5},  # 40 °C
        "E104_HNO3_Concentrator": {"T": 373.15, "p": 1.0e5},

        # ------------------------------------------------------------------
        # Solvent loop
        # ------------------------------------------------------------------
        "M133_OrgRecoveryMixer": {"T": 298.15, "p": 1.0e5},
        "V133_SolventPurge": {"T": 298.15, "p": 1.0e5},
        # ------------------------------------------------------------------
        # NOx reduction
        # ------------------------------------------------------------------
        "R103_NOxReduction": {"T": R103_T_K, "p": R103_P_PA},
        "KO101_KnockOutPot": {"T": KO101_T_K, "p": KO101_P_PA},
        "V105_KOWaterPurge": {"T": 323.15, "p": 1.0e5},
        "M105_WaterTrimMixer": {"T": 298.15, "p": 1.0e5},
        "V102_AcidSurge": {"T": 298.15, "p": 1.0e5},
        "D102_NH3_Absorber": {"T": D102_T_K, "p": D102_P_PA},
        # ------------------------------------------------------------------
        # Iodine polishing
        # ------------------------------------------------------------------
        "TSA101A_ColA": {"T": 298.15, "p": 1.0e5},
        "TSA101B_ColB": {"T": 298.15, "p": 1.0e5},



    }









