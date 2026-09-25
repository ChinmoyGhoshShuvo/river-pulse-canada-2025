"""
Daily discharge to the sea from Canada's rivers, 2025 (data/discharge_to_sea_2025.csv,
written by code/render_river_pulse.py). Run:  python make_figures.py  -> ../images/
"""
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).parent
OUT = HERE.parent / "images"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 200,
    "savefig.bbox": "tight", "figure.facecolor": "white",
})

d = pd.read_csv(HERE / "data" / "discharge_to_sea_2025.csv", parse_dates=["date"]).set_index("date")
q = d["discharge_to_sea_m3s"] / 1000
fig, ax = plt.subplots(figsize=(7.2, 2.8))
ax.fill_between(q.index, q, color="#56B4E9", alpha=0.25, lw=0)
ax.plot(q.index, q, color="#0072B2", lw=1.2, label="Daily total")
ax.axhline(q.mean(), color="#777", lw=0.8, ls="--")
ax.text(q.index[5], q.mean() + 4, f"2025 mean {q.mean():,.1f}k", fontsize=7.5, color="#555")
for t, lab, va in [(q.idxmin(), "Low", "top"), (q.idxmax(), "Peak", "bottom")]:
    ax.scatter(t, q[t], color="#D55E00", zorder=3, s=18)
    ax.annotate(f"{lab}: {q[t]:,.1f}k m³/s\n{t:%d %b}", (t, q[t]), xytext=(8, -18 if va == "top" else 6),
                textcoords="offset points", fontsize=7.5)
ax.set_ylabel("Discharge to the sea (thousand m³/s)")
ax.set_ylim(0, 250)
ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
ax.set_xlim(q.index[0], q.index[-1])
ax.set_title(f"Canada's rivers to the sea, 2025: {q.max() / q.min():.1f}x swing from winter low to June peak", loc="left")
fig.text(0, -0.06, "Sum of modelled daily discharge at 3,833 river mouths (GloFAS v4 pattern scaled to HydroRIVERS means).",
         fontsize=7, color="#555")
fig.savefig(OUT / "discharge-to-sea-2025.png")
print("written")
