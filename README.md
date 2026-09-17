## Simulation of Left Ventricular Dysfunction Mechanisms via Pressure–Volume and Transmitral Flow Analysis Using the CircAdapt Lumped-Parameter Model

> **A computational cardiovascular modelling study isolating three independent mechanisms of left ventricular dysfunction and quantifying their hemodynamic effects through pressure–volume loop and transmitral flow analysis.**

---

## Motivation

Heart failure involves simultaneous changes in contractility, passive myocardial stiffness, and arterial afterload but their individual contributions to hemodynamic deterioration are difficult to isolate clinically because they co-exist in patients. This project uses the CircAdapt lumped-parameter cardiovascular model to simulate each mechanism in isolation and compare the resulting shifts in LV pressure–volume relationships and diastolic filling dynamics.

## Model

**Platform:** [CircAdapt](https://www.circadapt.org/) (v26.02) - a closed-loop, multi-scale computational model of the heart and circulation that solves coupled ODEs representing cardiac chambers, valves, pericardium, and the systemic/pulmonary vasculature.

**Reference parameterization:** VanOsta2023, representing a healthy adult cardiovascular system at rest (HR = 70 bpm, CO ≈ 5 L/min).

## Simulated Conditions

| Condition | Parameter Modified | Change | Clinical Analogue |
|---|---|---|---|
| **Baseline** | — | — | Healthy adult |
| **Reduced contractility** | `Sf_act` (active fiber stress, LV wall + septum) | ×0.60 | HFrEF (systolic heart failure) |
| **Increased stiffness** | `Sf_pas` (passive fiber stress, LV wall + septum) | ×6.0 | HFpEF (diastolic dysfunction) |
| **Increased afterload** | `p0` (systemic vascular reference pressure) | ×1.50 | Systemic hypertension |

## Results

### 1. PV Loop Overlay

Each mechanism produces a distinct, physiologically consistent shift in the LV pressure–volume relationship:

<img width="1600" height="1200" alt="fig1_pv_overlay" src="https://github.com/user-attachments/assets/02c87bb3-f49c-4bbd-8e30-9324e44bb38d" />


- **HFrEF** (red): loop shifts rightward - the weakened muscle cannot empty the ventricle, raising ESV.
- **HFpEF** (green): loop shifts leftward and narrows - the stiff wall resists filling, reducing EDV and raising EDP.
- **HTN** (purple): loop shifts upward - the ventricle generates higher pressure to eject against increased resistance.

---

### 2. Annotated PV Panels

Individual PV loops with hemodynamic metrics annotated. Baseline (grey) is ghosted behind each disease condition for direct comparison:

<img width="2379" height="1969" alt="fig2_annotated_panels" src="https://github.com/user-attachments/assets/23139bed-099f-4b1c-8a2a-5fbaa20d2315" />


---

### 3. Hemodynamic Deviations

Percent change from the healthy baseline across all measured indices:

<img width="2000" height="1200" alt="fig3_pct_change" src="https://github.com/user-attachments/assets/1409dd3a-931c-4252-bc5d-403bc3272691" />


| Metric | HFrEF | HFpEF | HTN |
|---|---|---|---|
| EF | −23% | −11% | −4% |
| SV | −20% | −23% | −4% |
| ESP | −13% | −10% | +4% |
| EDP | −15% | **+114%** | +0% |
| SW | −31% | −32% | +1% |

**Key findings:**
- **HFrEF** shows the largest EF reduction (66.5% → 51.3%) and ESV increase (+50%), consistent with systolic pump failure.
- **HFpEF** preserves near-normal ESV (+6%) but dramatically elevates end-diastolic pressure (+114%) - the hallmark of diastolic dysfunction with preserved ejection fraction.
- **HTN** shows compensated hemodynamics with elevated ESP (+4%) but minimal volume or EF changes, representing early-stage pressure overload before structural remodeling.

---

### 4. Cardiac Cycle: Pressure and Transmitral Flow Waveforms

Time-domain analysis showing LV/aortic pressure traces (top) and mitral valve flow with E/A wave identification (bottom):

<img width="2000" height="1400" alt="fig4_waveforms" src="https://github.com/user-attachments/assets/9d6cc846-6f77-4661-8fc5-4c9f0ad76842" />


The **E wave** (early passive filling) and **A wave** (atrial contraction) are marked on the transmitral flow. 
E/A ratio - a primary echocardiographic index for grading diastolic dysfunction; showed limited sensitivity to isolated parameter changes at moderate perturbation levels, remaining >1 across all conditions. 
HTN produced the largest relative E/A reduction (−8% vs baseline), consistent with afterload elevation being an early driver of diastolic impairment. This is consistent with clinical observations that E/A alone is insufficient for diagnosing early-stage diastolic dysfunction.

---

## Reference

van Loon T, Knackstedt C, Cornelussen R, et al. *Increased myocardial stiffness more than impaired relaxation function limits cardiac performance during exercise in heart failure with preserved ejection fraction: a virtual patient study.* European Heart Journal – Digital Health, 1(1):40–50, 2020. [DOI: 10.1093/ehjdh/ztaa009](https://doi.org/10.1093/ehjdh/ztaa009)

## Requirements

- Python 3.10+
- [CircAdapt](https://www.circadapt.org/) v26.02 (`circadapt-2602` wheel)
- NumPy
- Matplotlib
- SciPy

### Installation

```bash
pip install numpy matplotlib scipy
pip install circadapt-2602-py3-none-win_amd64.whl
```

### Run

```bash
python circadapt_hf_study.py
```

Figures and metrics are saved to `results/`.

## Repository Structure

| File | Description |
|---|---|
| `circadapt_hf_study.py` | Complete simulation, analysis, and visualization script |
| `README.md` | Project documentation |
| `results/fig1_pv_overlay.png` | Overlaid PV loops across all conditions |
| `results/fig2_annotated_panels.png` | 2×2 annotated PV comparison with hemodynamic metrics |
| `results/fig3_pct_change.png` | Percent deviation from healthy baseline |
| `results/fig4_waveforms.png` | Time-domain LV pressure and transmitral flow waveforms |
| `results/hemodynamic_table.csv` | Raw hemodynamic metrics export |


## Author 

Khadija Khan (khadijakhanbme@gmail.com)

MSc Biomedical Engineering
