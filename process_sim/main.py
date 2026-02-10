from __future__ import annotations
from .build_flowsheet import build_flowsheet
from .flowsheet_tools import export_to_excel, export_units_to_excel
from .inputs import InputParameters



def main() -> None:
    fs = build_flowsheet()
    params = InputParameters()
    export_to_excel(fs, "stream_table.xlsx")
    export_units_to_excel(fs, params, "flowsheet_units.xlsx")
    print("Converged. Wrote flowsheet_results.xlsx")

    if hasattr(fs, "design"):
        print("\n--- Design summary ---")
        for k, v in fs.design.items():
            if k == "V3_recommended_cycle":
                print(f"{k}:")
                for kk, vv in v.items():
                    print(f"  {kk}: {vv}")
            else:
                print(f"{k}: {v}")

if __name__ == "__main__":
    main()
