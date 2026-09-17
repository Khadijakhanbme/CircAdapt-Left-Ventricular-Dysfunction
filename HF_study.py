"""
In Silico Assessment of Left Ventricular Dysfunction Mechanisms
Using the CircAdapt Cardiovascular Model

This script simulates three independent pathological mechanisms and
compares their hemodynamic effects via PV loop analysis:
  1. Reduced contractility    → models systolic HF (HFrEF)
  2. Increased passive stiffness → models diastolic dysfunction (HFpEF)
  3. Increased afterload      → models systemic hypertension

Reference model: VanOsta2023 (CircAdapt v26.02)
Methodology informed by: van Loon et al., Eur Heart J Digital Health, 2020
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from circadapt.model.vanosta2023 import VanOsta2023

# ─── Output directory ─────────────────────────────────────────────
OUT_DIR = "results"
os.makedirs(OUT_DIR, exist_ok=True)

# ─── Column indices (from model cavity object order) ──────────────
# [0] SyArt  [1] SyVen  [2] PuArt  [3] PuVen  [4] Peri
# [5] La     [6] Ra     [7] cLv    [8] cRv
LV_IDX = 7  # Model.Peri.TriSeg.cLv
RV_IDX = 8  # Model.Peri.TriSeg.cRv

# Patch object indices:
# [0] La wall  [1] Ra wall  [2] LV wall  [3] Septum  [4] RV wall
LV_PATCHES = [2, 3]  # LV free wall + septum (both contribute to LV function)


# ─── Helper functions ─────────────────────────────────────────────

def extract_lv(model):
    """Extract LV pressure (mmHg) and volume (mL) from a solved model."""
    cav = model["Cavity"]
    V = np.array(cav["V"])[:, LV_IDX] * 1e6      # m³ → mL
    p = np.array(cav["p"])[:, LV_IDX] / 133.322   # Pa → mmHg
    return V, p


def compute_metrics(V, p):
    """Compute hemodynamic indices from one LV PV loop."""
    esv = V.min()
    edv = V.max()
    sv  = edv - esv
    ef  = 100.0 * sv / edv
    esp = p.max()                          # peak systolic pressure
    edp = p[np.argmax(V)]                  # pressure at max volume (end-diastole)
    # Stroke work = area inside PV loop (closed contour integral)
    if hasattr(np, "trapezoid"):
        sw = np.abs(np.trapezoid(p, V))
    else:
        sw = np.abs(np.trapz(p, V))
    sw_j = sw * 133.322e-6 * 1e-6          # convert to Joules (Pa·m³)
    return {
        "ESV (mL)":  round(esv, 1),
        "EDV (mL)":  round(edv, 1),
        "SV (mL)":   round(sv, 1),
        "EF (%)":    round(ef, 1),
        "ESP (mmHg)": round(esp, 1),
        "EDP (mmHg)": round(edp, 1),
        "SW (mmHg·mL)": round(sw, 0),
        "SW (J)":    round(sw_j, 4),
    }


def extract_time(model):
    """Extract time vector (seconds) for one cardiac cycle."""
    solver = model["Solver"]
    # Try direct time signal, fallback to dt-based reconstruction
    for sig in ["t", "time"]:
        if sig in solver.signals:
            t = np.array(solver[sig]).flatten()
            if len(t) > 1:
                return t - t[0]  # start from 0
    # Fallback
    n = np.array(model["Cavity"]["V"]).shape[0]
    t_cycle = float(np.array(model["General"]["t_cycle"]))
    return np.linspace(0, t_cycle, n)


def extract_flows_and_pressures(model):
    """Extract valve flows and aortic pressure for time-domain plots."""
    extras = {}

    # Aortic pressure = systemic arterial cavity pressure (index 0)
    p_all = np.array(model["Cavity"]["p"])
    extras["p_ao"] = p_all[:, 0] / 133.322  # SyArt → mmHg

    # Valve flows (q signal) — if available
    try:
        valve = model["Valve"]
        q = np.array(valve["q"])  # flow through each valve
        # Valve object order typically: mitral, aortic, tricuspid, pulmonary
        print("  Valve objects:", valve.objects)
        extras["q_valves"] = q * 1e6 * 60  # m³/s → mL/s → L/min (×60/1000)
        # Actually: m³/s → mL/s is ×1e6, then we keep mL/s for clarity
        extras["q_valves_mls"] = q * 1e6  # m³/s → mL/s
    except Exception as e:
        print(f"  Valve flow extraction: {e}")

    return extras


def run_condition(name, modify_fn=None):
    """
    Create a fresh model, optionally modify parameters, run to steady
    state, and return LV volume, pressure, computed metrics, and extras.
    """
    model = VanOsta2023()

    if modify_fn is not None:
        modify_fn(model)

    # Run solver — default call iterates until hemodynamic steady state
    model.run()

    V, p = extract_lv(model)
    metrics = compute_metrics(V, p)
    t = extract_time(model)
    extras = extract_flows_and_pressures(model)

    print(f"\n{'─'*50}")
    print(f"  {name}")
    print(f"{'─'*50}")
    for k, v in metrics.items():
        print(f"  {k:>16s}: {v}")

    return V, p, metrics, t, extras


# ─── Define the three pathological modifications ──────────────────

def reduce_contractility(model):
    """
    Reduce LV active stress (Sf_act) to 60% of baseline.
    
    Sf_act controls peak active fiber stress during contraction.
    Baseline = 120,000 Pa for LV patches. Reducing it simulates
    weakened systolic contraction as seen in HFrEF.
    
    Expected effect: ESV ↑, SV ↓, EF ↓, loop shifts right + narrows.
    """
    patch = model["Patch"]
    sf_act = np.array(patch["Sf_act"])
    for i in LV_PATCHES:
        sf_act[i] *= 0.60
    patch["Sf_act"] = sf_act


def increase_stiffness(model):
    """
    Increase LV passive stress (Sf_pas) to 300% of baseline.
    
    Sf_pas governs the passive elastic restoring force of myocardium.
    Baseline ≈ 420 Pa for LV patches. Tripling it simulates a stiff,
    fibrotic ventricle that resists filling — the hallmark of HFpEF.
    
    Expected effect: EDV ↓, EDP ↑, loop shifts left and upward at
    the bottom. EF may remain near-normal (preserved EF).
    
    Reference: van Loon et al. (2020) used end-diastolic elastance
    (Eed) ranging from 0.2 to 2.6 mmHg/mL to model HFpEF severity.
    """
    patch = model["Patch"]
    sf_pas = np.array(patch["Sf_pas"])
    for i in LV_PATCHES:
        sf_pas[i] *= 3.0
    patch["Sf_pas"] = sf_pas


def increase_afterload(model):
    """
    Increase systemic vascular reference pressure (p0) by 50%.
    
    ArtVen p0 is the reference pressure drop across the peripheral
    vasculature. Increasing it raises systemic vascular resistance,
    simulating chronic hypertension — the LV must generate higher
    pressure to eject against a stiffer arterial system.
    
    Expected effect: ESP ↑, ESV ↑, loop shifts upward + rightward.
    """
    artven = model["ArtVen"]
    p0 = np.array(artven["p0"])
    p0[0] *= 1.50  # index 0 = systemic circulation (Model.CiSy)
    artven["p0"] = p0


# ─── Run all conditions ──────────────────────────────────────────

conditions = [
    ("Baseline (healthy)",        None),
    ("Reduced contractility (HFrEF)", reduce_contractility),
    ("Increased stiffness (HFpEF)",   increase_stiffness),
    ("Increased afterload (HTN)",     increase_afterload),
]

colors = ["#2563eb", "#dc2626", "#16a34a", "#9333ea"]  # blue, red, green, purple
styles = ["-", "--", "-.", ":"]

results = {}
for (name, fn), color, ls in zip(conditions, colors, styles):
    V, p, metrics, t, extras = run_condition(name, fn)
    results[name] = {"V": V, "p": p, "metrics": metrics, "t": t,
                     "extras": extras, "color": color, "ls": ls}


# ═══════════════════════════════════════════════════════════════════
#  VISUALIZATIONS
# ═══════════════════════════════════════════════════════════════════

short_names  = ["Baseline", "HFrEF", "HFpEF", "HTN"]
result_list  = list(results.values())
name_list    = list(results.keys())


# ─── Figure 1: Overlaid PV loops (clean version) ─────────────────

fig1, ax1 = plt.subplots(figsize=(8, 6))
for name, data in results.items():
    ax1.plot(data["V"], data["p"], color=data["color"], linestyle=data["ls"],
             linewidth=2.2, label=name)
ax1.set_xlabel("LV Volume (mL)", fontsize=12)
ax1.set_ylabel("LV Pressure (mmHg)", fontsize=12)
ax1.set_title("LV Pressure–Volume Loops Under Pathological Conditions", fontsize=13)
ax1.legend(fontsize=9, loc="upper right")
ax1.set_xlim(left=0)
ax1.set_ylim(bottom=0)
ax1.grid(True, alpha=0.3)
plt.tight_layout()
fig1.savefig(os.path.join(OUT_DIR, "fig1_pv_overlay.png"), dpi=200)
plt.show()


# ─── Figure 2: 2×2 annotated panels (each condition vs baseline) ─

fig2, axes2 = plt.subplots(2, 2, figsize=(12, 10))
axes2 = axes2.flatten()
panel_labels = ["(a) Baseline (healthy)", "(b) Reduced contractility (HFrEF)",
                "(c) Increased stiffness (HFpEF)", "(d) Increased afterload (HTN)"]

for idx, (name, data) in enumerate(results.items()):
    ax = axes2[idx]
    V, p = data["V"], data["p"]
    m = data["metrics"]

    # Ghost baseline behind disease panels
    if idx > 0:
        bl = result_list[0]
        ax.plot(bl["V"], bl["p"], color="#cccccc", linewidth=1.5, label="Baseline")
        ax.fill(bl["V"], bl["p"], color="#f0f0f0", alpha=0.4)

    # Main loop
    ax.plot(V, p, color=data["color"], linewidth=2.2)
    ax.fill(V, p, color=data["color"], alpha=0.08)

    # Mark key points
    esv_idx = np.argmin(V)
    edv_idx = np.argmax(V)
    esp_idx = np.argmax(p)

    ax.plot(V[esv_idx], p[esv_idx], "o", color=data["color"], ms=7, zorder=5)
    ax.plot(V[edv_idx], p[edv_idx], "s", color=data["color"], ms=7, zorder=5)
    ax.plot(V[esp_idx], p[esp_idx], "^", color=data["color"], ms=7, zorder=5)

    # Annotation box with metrics
    info = (f"ESV = {m['ESV (mL)']:.0f} mL\n"
            f"EDV = {m['EDV (mL)']:.0f} mL\n"
            f"SV  = {m['SV (mL)']:.0f} mL\n"
            f"EF  = {m['EF (%)']:.1f}%\n"
            f"ESP = {m['ESP (mmHg)']:.0f} mmHg\n"
            f"EDP = {m['EDP (mmHg)']:.1f} mmHg")
    ax.text(0.03, 0.97, info, transform=ax.transAxes, fontsize=8.5,
            verticalalignment="top", fontfamily="monospace",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      edgecolor=data["color"], alpha=0.9))

    ax.set_title(panel_labels[idx], fontsize=11, fontweight="bold", loc="left")
    ax.set_xlabel("Volume (mL)", fontsize=10)
    ax.set_ylabel("Pressure (mmHg)", fontsize=10)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.2)

fig2.suptitle("LV Pressure–Volume Analysis: Baseline vs Pathological Conditions",
              fontsize=14, fontweight="bold", y=0.98)
fig2.subplots_adjust(top=0.88)
plt.tight_layout(rect=[0, 0, 1, 0.96])
fig2.savefig(os.path.join(OUT_DIR, "fig2_annotated_panels.png"), dpi=200,
             bbox_inches="tight")
plt.show()


# ─── Figure 3: Percent change from baseline (horizontal bars) ────

change_metrics = ["ESV (mL)", "EDV (mL)", "SV (mL)", "EF (%)", "ESP (mmHg)",
                  "EDP (mmHg)", "SW (mmHg·mL)"]
baseline_m = result_list[0]["metrics"]

fig3, ax3 = plt.subplots(figsize=(10, 6))
y_positions = np.arange(len(change_metrics))
bar_h = 0.25

for i in range(1, 4):  # skip baseline
    pct_changes = []
    for k in change_metrics:
        bv = baseline_m[k]
        cv = result_list[i]["metrics"][k]
        pct = 100.0 * (cv - bv) / abs(bv) if bv != 0 else 0
        pct_changes.append(pct)

    offset = (i - 2) * bar_h
    bars = ax3.barh(y_positions + offset, pct_changes, bar_h,
                    color=result_list[i]["color"], alpha=0.85,
                    label=short_names[i])
    # Value labels
    for bar, val in zip(bars, pct_changes):
        x_pos = bar.get_width()
        ha = "left" if x_pos >= 0 else "right"
        ax3.text(x_pos + (2 if x_pos >= 0 else -2), bar.get_y() + bar.get_height()/2,
                 f"{val:+.0f}%", va="center", ha=ha, fontsize=7.5, fontweight="bold")

ax3.axvline(0, color="black", linewidth=0.8)
ax3.set_yticks(y_positions)
ax3.set_yticklabels(change_metrics, fontsize=10)
ax3.set_xlabel("Change from Baseline (%)", fontsize=11)
ax3.set_title("Hemodynamic Deviations from Healthy Baseline", fontsize=13)
ax3.legend(fontsize=9, loc="lower right")
ax3.grid(True, axis="x", alpha=0.3)
plt.tight_layout()
fig3.savefig(os.path.join(OUT_DIR, "fig3_pct_change.png"), dpi=200)
plt.show()


# ─── Figure 4: Time-domain waveforms (pressure + flow) ───────────
#
# This shows WHEN things happen during the cardiac cycle, which the
# PV loop compresses away. Top row: LV and aortic pressure traces.
# Bottom row: valve flows (mitral inflow / aortic outflow).

fig4, axes4 = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

# ── Top: LV pressure + aortic pressure over time ──
ax_p = axes4[0]
for name, data in results.items():
    t_ms = data["t"] * 1000  # seconds → ms
    ax_p.plot(t_ms, data["p"], color=data["color"], linestyle=data["ls"],
              linewidth=1.8, label=f"LV – {name}")
    # Aortic pressure (dashed, thinner, same color)
    if "p_ao" in data["extras"]:
        ax_p.plot(t_ms, data["extras"]["p_ao"], color=data["color"],
                  linestyle=data["ls"], linewidth=0.9, alpha=0.5)

ax_p.set_ylabel("Pressure (mmHg)", fontsize=11)
ax_p.set_title("Cardiac Cycle: Pressure and Flow Waveforms", fontsize=13)
ax_p.legend(fontsize=7.5, ncol=2, loc="upper right")
ax_p.set_ylim(bottom=0)
ax_p.grid(True, alpha=0.2)
# Add a text label for aortic traces
ax_p.text(0.01, 0.92, "Bold = LV pressure, faint = aortic pressure",
          transform=ax_p.transAxes, fontsize=7.5, color="gray")

# ── Bottom: Valve flows (if available) ──
ax_q = axes4[1]
has_flows = False
for name, data in results.items():
    t_ms = data["t"] * 1000
    if "q_valves_mls" in data["extras"]:
        has_flows = True
        q = data["extras"]["q_valves_mls"]
        # Typically: columns map to valve objects in order
        # We plot total flow through mitral (index 0) and aortic (index 1)
        # Mitral flow is positive during filling, aortic during ejection
        if q.ndim == 2 and q.shape[1] >= 2:
            ax_q.plot(t_ms, q[:, 0], color=data["color"], linestyle=data["ls"],
                      linewidth=1.5, alpha=0.8)

if has_flows:
    ax_q.set_ylabel("Mitral Flow (mL/s)", fontsize=11)
    ax_q.axhline(0, color="black", linewidth=0.5)
    ax_q.grid(True, alpha=0.2)
else:
    # Fallback: plot LV volume over time (always available)
    for name, data in results.items():
        t_ms = data["t"] * 1000
        ax_q.plot(t_ms, data["V"], color=data["color"], linestyle=data["ls"],
                  linewidth=1.8, label=name)
    ax_q.set_ylabel("LV Volume (mL)", fontsize=11)
    ax_q.legend(fontsize=8, ncol=2)
    ax_q.grid(True, alpha=0.2)

ax_q.set_xlabel("Time (ms)", fontsize=11)
plt.tight_layout()
fig4.savefig(os.path.join(OUT_DIR, "fig4_waveforms.png"), dpi=200)
plt.show()


# ─── CSV export ───────────────────────────────────────────────────

header = ["Condition"] + list(baseline_m.keys())
csv_path = os.path.join(OUT_DIR, "hemodynamic_table.csv")
with open(csv_path, "w") as f:
    f.write(",".join(header) + "\n")
    for name, data in results.items():
        row = [name] + [str(v) for v in data["metrics"].values()]
        f.write(",".join(row) + "\n")


# ─── Console summary ─────────────────────────────────────────────

print(f"\n{'='*80}")
print("  RESULTS SUMMARY")
print(f"{'='*80}")
print(f"\n{'Condition':<35s}", end="")
for k in ["EF (%)", "SV (mL)", "ESP (mmHg)", "EDP (mmHg)"]:
    print(f"{k:>15s}", end="")
print()
print("─" * 95)
for name, data in results.items():
    print(f"{name:<35s}", end="")
    for k in ["EF (%)", "SV (mL)", "ESP (mmHg)", "EDP (mmHg)"]:
        print(f"{data['metrics'][k]:>15}", end="")
    print()

print(f"\n→ All figures saved to: {OUT_DIR}/")
print("  fig1_pv_overlay.png       — Overlaid PV loops")
print("  fig2_annotated_panels.png — 2×2 annotated comparison")
print("  fig3_pct_change.png       — Percent deviation bar chart")
print("  fig4_waveforms.png        — Time-domain pressure & flow")
print("  hemodynamic_table.csv     — Raw data export")