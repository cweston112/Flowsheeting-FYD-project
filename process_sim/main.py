from __future__ import annotations

import os
import sys
import time

_here = os.path.dirname(os.path.abspath(__file__))

from .build_flowsheet import build_flowsheet
from .framework import export_to_excel, display_stream_name

KEY_STREAMS = [
    ("F2", "Acid makeup to dissolver train"),
    ("F3", "Dissolver offgas"),
    ("F5", "Dissolver liquor"),
    ("F7", "Solvent to X101"),
    ("F21", "Pu/Np product liquor"),
    ("F24PT", "U nitrate feed to calciner"),
    ("F27", "U3O8 product"),
    ("F2R", "Acid recycle"),
    ("F45R", "KO water recycle"),
    ("F15", "Solvent recycle"),
    ("F50A", "Final stack gas (active TSA outlet)"),
    ("F50B", "Final stack gas (standby TSA outlet)"),
]


def print_summary(fs) -> None:
    print("\n" + "=" * 80)
    print("SIMPLIFIED PROCESS MODEL SUMMARY")
    print("=" * 80)
    for sname, desc in KEY_STREAMS:
        s = fs.streams.get(sname)
        if s is None:
            continue
        dname = display_stream_name(sname)
        print(f"{dname:10s}  {s.phase:1s}  T={s.T:7.2f} K  p={s.p:9.1f} Pa  F={s.total_molar_flow():10.6f} mol/s  {desc}")
    print("-" * 80)
    print(f"Converged={getattr(fs, 'converged', '?')}  Iterations={getattr(fs, 'recycle_iter', '?')}  Error={getattr(fs, 'recycle_err', float('nan')):.3e}")


def main():
    t0 = time.time()
    fs = build_flowsheet()
    print_summary(fs)
    out_dir = os.path.join(_here, "results")
    os.makedirs(out_dir, exist_ok=True)
    xlsx_path = os.path.join(out_dir, "stream_results.xlsx")
    export_to_excel(fs, xlsx_path)
    print(f"\nWrote {xlsx_path} in {time.time() - t0:.2f} s")
    return fs


if __name__ == "__main__":
    main()
