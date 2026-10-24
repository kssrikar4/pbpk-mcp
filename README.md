# PBPK MCP

This repository contains a re-implementation of the core Physiologically Based Pharmacokinetic (PBPK) ordinary differential equation (ODE) execution engine extracted directly from the Open Systems Pharmacology (OSP) repositories (`PK-Sim`, `OSPSuite.Core`, `OSPSuite.SimModel`).

## System Requirements
- Operating System: Linux / macOS / Windows
- Python: Version 3.10, 3.11, 3.12, or 3.13
- C/C++ Compiler / LLVM

## Dependencies
- `numpy >= 1.26.0`
- `scipy >= 1.12.0`
- `numba >= 0.60.0`
- `matplotlib >= 3.8.0`
- `fastmcp >= 0.1.0`
- `mcp >= 1.0.0`
- `rdkit >= 2023.9.0`

## Installation

```bash
git clone https://github.com/kssrikar4/pbpk-mcp.git
cd pbpk-mcp
pip install -r requirements.txt
pip install -e .
```

## Quick Start

The test suite validates mass balance conservation ($dC/dt = 0$ at steady state and total mass conservation to within $10^{-6}$ relative tolerance), Rodgers & Rowland tissue partition coefficients against analytical benchmarks, parallel population simulations, and MCP tools:

```bash
python3 -m unittest discover tests
```

To run individual benchmark suites:

```bash
python3 -m unittest tests/test_mass_balance.py
python3 -m unittest tests/test_partition.py
python3 -m unittest tests/test_simulation.py
python3 -m unittest tests/test_server.py
python3 -m unittest tests/test_properties.py
python3 -m unittest tests/test_performance.py
```

## Running the FastMCP Server

The MCP server exposes flat primitive tools (`get_physiological_params`, `calculate_partition_coefficients`, `estimate_compound_properties`, `convert_clearance_units`, `run_pbpk_simulation`, `run_smiles_pbpk_simulation`, `compare_renal_impairment`, `generate_plots`, and `generate_comparison_plots`):

### Standard I/O (STDIO) Mode
```bash
python3 -m pbpk.server
```

Or via FastMCP CLI:
```bash
fastmcp run pbpk/server.py
```

### Inspecting MCP Server Schema
```bash
fastmcp inspect pbpk/server.py
```

## Usage Guide

### 1. Scaling Physiology
Extract exact anatomical organ volumes and blood flow rates for an adult human scaled to a custom body weight:

```python
from pbpk.physiology import scale_physiology

phys = scale_physiology("Human", weight_kg=75.0)

print("Arterial Blood (L):", phys["volumes"]["ArterialBlood"])
print("Venous Blood (L):", phys["volumes"]["VenousBlood"])
print("Liver (L):", phys["volumes"]["Liver"])
print("Cardiac Output (L/h):", phys["flows_l_h"]["CardiacOutput"])
```

### 2. Calculating Rodgers & Rowland Tissue Partition Coefficients ($K_p$)
Calculate tissue-to-plasma partition coefficients across all 15 human organs for neutral, basic, acidic, or zwitterionic compounds:

```python
from pbpk.partition import calculate_rodgers_rowland_kp

# Propranolol: Strong Monoprotic Base (pKa = 9.45, LogP = 2.6, fu = 0.13, B/P = 1.15)
kp_dict = calculate_rodgers_rowland_kp(
    compound_type="monoprotic_base",
    log_p=2.6,
    pka1=9.45,
    fup=0.13,
    b2p=1.15,
    species_data=phys
)

for organ, val in kp_dict.items():
    print(f"{organ:15s}: {val:.3f}")
```

### 3. Automated In Silico Property Estimation from SMILES
Automatically derive physicochemical parameters (`log_p`, `pka1`, `compound_type`, `fup`, `b2p`, `ka_h`, `clearance_hep_l_h`, `clearance_ren_l_h`) for novel drugs without experimental wet-lab measurements:

```python
from pbpk.properties import estimate_parameters_from_smiles

# Metformin: SMILES CN(C)C(=N)NC(=N)N
metformin_props = estimate_parameters_from_smiles("CN(C)C(=N)NC(=N)N")
print("Compound Type:", metformin_props["compound_type"])
print("pKa1:", metformin_props["pka1"])
print("LogP:", metformin_props["log_p"])
print("Hepatic Clearance (L/h):", metformin_props["clearance_hep_l_h"])
print("Renal Clearance (L/h):", metformin_props["clearance_ren_l_h"])
```

### 4. Running Single-Subject Simulations

#### A. IV Bolus Dose
```python
from pbpk.model import WholeBodyPBPKModel
from pbpk.simulation import PBPKSimulation

model = WholeBodyPBPKModel(phys, kp_dict)

sim = PBPKSimulation(
    model=model,
    dose_mg=100.0,
    dose_route="iv_bolus",
    simulation_time_h=24.0,
    clearance_hep_l_h=18.0,
    clearance_ren_l_h=4.5
)

results = sim.run()

print("Cmax (mg/L):", results["summary"]["c_max_mg_l"])
print("Tmax (h):", results["summary"]["t_max_h"])
print("AUC (mg*h/L):", results["summary"]["auc_plasma_mg_h_l"])
print("Max Mass Discrepancy (mg):", results["max_mass_balance_error_mg"])
```

#### B. Intravenous Infusion
```python
sim_inf = PBPKSimulation(
    model=model,
    dose_mg=200.0,
    dose_route="iv_infusion",
    infusion_duration_h=2.0,
    simulation_time_h=24.0,
    clearance_hep_l_h=18.0,
    clearance_ren_l_h=4.5
)
results_inf = sim_inf.run()
```

#### C. Oral Administration
```python
sim_oral = PBPKSimulation(
    model=model,
    dose_mg=100.0,
    dose_route="oral",
    simulation_time_h=24.0,
    clearance_hep_l_h=18.0,
    clearance_ren_l_h=4.5,
    ka_h=1.2
)
results_oral = sim_oral.run()
```

### 4. Running Multi-Subject Virtual Populations
Simulate inter-individual anatomical variability across virtual subjects using multi-threaded parallel execution:

```python
from pbpk.simulation import run_virtual_population

compound_params = {
    "compound_type": "monoprotic_base",
    "log_p": 2.6,
    "pka1": 9.45,
    "fup": 0.13,
    "b2p": 1.15
}

dose_params = {
    "dose_mg": 100.0,
    "dose_route": "iv_bolus",
    "simulation_time_h": 24.0,
    "clearance_hep_l_h": 18.0,
    "clearance_ren_l_h": 4.5,
    "num_time_points": 300
}

pop_results = run_virtual_population(
    compound_params=compound_params,
    dose_params=dose_params,
    num_subjects=50,
    mean_weight_kg=73.0,
    weight_cv=0.15,
    n_workers=4
)

print("Population Mean AUC:", pop_results["mean_auc"])
print("Population Mean Cmax:", pop_results["mean_cmax"])
```

### 5. Plotting
Generate plots comparing plasma and target organ concentrations over time:

```python
from pbpk.visualization import plot_pbpk_curves

plot_pbpk_curves(
    sim_result=results,
    output_path="propranolol_curves.png",
    target_organ="Liver"
)
```

## MCP Server Setup & AI Client Configuration

### 6. Installing the Package 
Install the package locally in editable mode:

```bash
cd pbpk-mcp
pip install -e .
```

This creates the registered command line utility:
```bash
pbpk-mcp --help
```

### 7. Starting the MCP Server Manually
To run the server in standard I/O (STDIO) mode:

```bash
pbpk-mcp
```

Or run directly with Python:
```bash
python3 -m pbpk.server
```

### 8. Configuring AI Clients to Use the MCP Server

#### Claude Desktop Configuration
Add the server entry to your `claude_desktop_config.json`:
- **Linux**: `~/.config/Claude/claude_desktop_config.json`
- **macOS**: `~/Library/Application Support/Claude/claude_desktop_config.json`
- **Windows**: `%APPDATA%\Claude\claude_desktop_config.json`

```json
{
  "mcpServers": {
    "pbpk-mcp": {
      "command": "~/.local/bin/pbpk-mcp"
    }
  }
}
```

*Or use the direct python module invocation:*
```json
{
  "mcpServers": {
    "pbpk-mcp": {
      "command": "python3",
      "args": ["-m", "pbpk.server"],
      "cwd": "~/pbpk-mcp"
    }
  }
}
```

#### Cursor / Antigravity Configuration
In Cursor or Antigravity MCP settings, add:
- **Name**: `pbpk-mcp`
- **Type**: `command` (stdio)
- **Command**: `pbpk-mcp` (or `~/.local/bin/pbpk-mcp`)

---

## 9. MCP JSON Schemas

Small LLM agents invoke these deterministic tools using standard primitives:

#### Tool 1: `get_physiological_params`
```json
{
  "species": "Human",
  "weight_kg": 73.0
}
```

#### Tool 2: `calculate_partition_coefficients`
```json
{
  "compound_type": "monoprotic_base",
  "log_p": 2.6,
  "pka1": 9.45,
  "pka2": 0.0,
  "fup": 0.13,
  "b2p": 1.15
}
```

#### Tool 3: `estimate_compound_properties`
Predicts $pK_a$, $f_{up}$, $\log P$, $k_a$, hepatic and renal clearance from SMILES:
```json
{
  "smiles": "CN(C)C(=N)NC(=N)N"
}
```

#### Tool 4: `run_pbpk_simulation`
```json
{
  "dose_mg": 100.0,
  "dose_route": "iv_bolus",
  "infusion_duration_h": 0.0,
  "simulation_time_h": 24.0,
  "clearance_hep_l_h": 18.0,
  "clearance_ren_l_h": 4.5,
  "log_p": 2.6,
  "pka1": 9.45,
  "compound_type": "monoprotic_base",
  "fup": 0.13,
  "b2p": 1.15,
  "weight_kg": 73.0,
  "num_subjects": 1
}
```

#### Tool 5: `run_smiles_pbpk_simulation`
Runs full PBPK ODE simulation directly from SMILES without manual parameters:
```json
{
  "smiles": "CN(C)C(=N)NC(=N)N",
  "dose_mg": 500.0,
  "dose_route": "oral",
  "simulation_time_h": 48.0
}
```

#### Tool 6: `compare_renal_impairment`
Simulates Healthy, Mild, Moderate, and Severe renal impairment, outputting multi-curve plots and metrics:
```json
{
  "smiles": "CN(C)C(=N)NC(=N)N",
  "dose_mg": 500.0,
  "simulation_time_h": 48.0,
  "output_dir": "."
}
```

#### Tool 7: `generate_plots`
Pass either the `simulation_id` or `simulation_file` returned by simulation tools:
```json
{
  "simulation_json": "sim_8e6bafab",
  "output_path": "simulation_plot.png",
  "target_organ": "Liver"
}
```

#### Tool 8: `generate_comparison_plots`
Compares multiple simulations across groups:
```json
{
  "simulation_ids": "sim_01,sim_02,sim_03,sim_04",
  "labels": "Healthy,Mild,Moderate,Severe",
  "output_path": "comparison_plot.png",
  "plot_type": "multicurve"
}
```

#### Tool 9: `convert_clearance_units`
Converts clearance metrics between `L/h` and `mL/min`:
```json
{
  "value": 28.0,
  "from_unit": "L/h",
  "to_unit": "mL/min"
}
```

## 10. Verified End-to-End Test Prompt & Benchmark Output

Users can copy and paste the following prompt into any AI assistant connected to the `pbpk-mcp` MCP server to verify that all tools, simulations, and visualizations execute properly.

### Test Prompt to Copy & Paste

```text
Develop a metformin PBPK model using its SMILES string "CN(C)C(=N)NC(=N)N".
Simulate plasma concentration–time profiles after oral administration.
Compare healthy individuals with mild, moderate and severe renal impairment.
Determine changes in peak concentration, total exposure and elimination.
Produce a final multi-curve graph demonstrating how renal function affects drug exposure.

2. Simulation design
Parameter: Proposed setup
Drug: Metformin
SMILES: CN(C)C(=N)NC(=N)N
Administration: Oral
Dosage: 500 mg single dose
Simulation duration: 48 hours
Reference population: Healthy adults
Comparison populations: Mild, moderate and severe renal impairment
Primary output: Plasma concentration versus time
Secondary outputs: C_max, T_max, AUC and half-life

Use a consistent virtual adult population and change renal physiology appropriately for each scenario. The severity categories and model parameters should be defined using published clinical and physiological data or estimated from the compound's structure rather than arbitrary settings.

Your final figure should contain four simulated curves:
- Healthy renal function
- Mild renal impairment
- Moderate renal impairment
- Severe renal impairment

4. Additional plots for your report
Figure 2 — AUC comparison:
A bar chart showing total drug exposure across the four renal-function groups.

Figure 3 — Renal function versus drug exposure:
Plot renal function against AUC or apparent elimination half-life to illustrate the relationship.

5. Workflow
Create a metformin compound using appropriate physicochemical and pharmacokinetic information derived from its chemical structure (SMILES: "CN(C)C(=N)NC(=N)N").
Build or select a suitable adult reference individual.
Configure oral administration and the single-dose protocol.
Parameterize renal impairment using evidence-based physiological changes, particularly renal function and renal drug clearance.
Run each scenario using the same dosing and observation schedule.
Export plasma concentration–time results and overlay all four curves.
Calculate and compare pharmacokinetic parameters, then assess the model against available clinical data.
```

### Verified Benchmark Output (Using qwen3-4b at 32k context)

```text
## **Metformin PBPK Model Development Summary**

### **1. Compound Properties (from SMILES "CN(C)C(=N)NC(=O)N")**
- Molecular Weight: ~130 g/mol
- logP: -0.849 (highly water-soluble, polar compound)
- pKa1: 12.4, pKa2: 0.0 (monoprotic base at physiological pH)
- Compound Type: monoprotic_base
- Fraction unbound in plasma (fup): 0.878
- b2p ratio: 1.0

### **2. Physiological Reference Population**
- Species: Human, Weight: 73 kg
- Kidney blood flow: ~79.5 L/h
- Plasma pH: 7.4 (physiologic)

### **3. Simulated Pharmacokinetic Parameters by Renal Function Group**

| Renal Function | eGFR (mL/min) | CrCl (L/h) | c_max (mg/L) | T_max (h) | AUC (mg·h/L) | t₁/₂ (h) | Relative AUC Fold |
|----------------|----------------|------------|--------------|-----------|--------------|----------|-------------------|
| Healthy | 100 | 27.4 | 2.02 | 1.73 | 18.16 | 4.31 | 1.0 (baseline) |
| Mild Impairment | 75 | 20.5 | 2.22 | 2.02 | 24.14 | 5.36 | 1.33 |
| Moderate Impairment | 45 | 12.3 | 2.55 | 2.69 | 39.47 | 8.21 | 2.17 |
| Severe Impairment | 20 | 5.5 | 2.99 | 3.94 | 76.07 | 17.32 | 4.19 |

### **4. Key Findings**

**Peak Concentration (c_max):**
- Minimal increase in mild impairment (+10%)
- Moderate: +27% increase
- Severe: +48% increase

**Total Exposure (AUC):**
- Mild impairment: 33% increase (clinically insignificant)
- Moderate impairment: 118% increase (dose consideration needed)
- Severe impairment: 319% increase (significant accumulation risk)

**Elimination Half-life:**
- Healthy to severe impairment: 4.3 h → 17.3 h (4× prolongation)

### **5. Generated Plots**

**Figure 1 — Multi-curve Plasma Concentration-Time Profiles:**
- Shows all four renal function groups overlaid
- Demonstrates progressive increase in exposure with worsening renal impairment
- Plot location: `/tmp/metformin_pk_comparison_alt.svg`

**Figure 2 — AUC Comparison Bar Chart:**
- Quantitative comparison of total drug exposure across groups
- Clear visualization showing exponential increase from severe to healthy impairment
- Plot location: `/tmp/renal_impairment_analysis/figure2_auc_comparison.png`

**Figure 3 — Renal Function vs Drug Exposure:**
- Correlates eGFR/CrCl with AUC and half-life changes
- Illustrates non-linear relationship between renal function and drug exposure
- Plot location: `/tmp/renal_impairment_analysis/figure3_renal_vs_exposure.png`

**Individual Concentration-Time Plots:**
- All four scenarios visualized separately for detailed analysis
- Locations:
  - Healthy: `/tmp/metformin_healthy.svg`
  - Mild: `/tmp/metformin_mild_impairment.svg`
  - Moderate: `/tmp/metformin_moderate_impairment.svg`
  - Severe: `/tmp/metformin_severe_impairment.svg`

### **6. Clinical Relevance**

The model predictions align with FDA recommendations for metformin:
- eGFR ≥60 mL/min (mild impairment): No dose adjustment needed (model shows <35% AUC increase)
- eGFR 30-59 mL/min (moderate impairment): Dose consideration (model shows ~2× AUC increase)
- eGFR <30 mL/min (severe impairment): Contraindicated (model shows >4× AUC increase with significant accumulation risk)
```

## Mathematical Formulations

### 1. Rodgers & Rowland Tissue Partition Coefficients ($K_p$)

#### Ionization Terms ($X$, $Y$, $Z$)
Evaluated at tissue intracellular pH ($pH_{iw} = 7.0$), plasma pH ($pH_{pls} = 7.4$), and red blood cell pH ($pH_{rbc} = 7.22$):

- **Neutrals**:
  $$X = 0, \quad Y = 0, \quad Z = 1$$
- **Monoprotic Acids**:
  $$X = 10^{pH_{iw} - pK_a}, \quad Y = 10^{pH_{pls} - pK_a}, \quad Z = 1$$
- **Diprotic Acids**:
  $$X = 10^{pH_{iw} - pK_{a1}} + 10^{2\cdot pH_{iw} - pK_{a1} - pK_{a2}}$$
  $$Y = 10^{pH_{pls} - pK_{a1}} + 10^{2\cdot pH_{pls} - pK_{a1} - pK_{a2}}, \quad Z = 1$$
- **Monoprotic Bases**:
  $$X = 10^{pK_a - pH_{iw}}, \quad Y = 10^{pK_a - pH_{pls}}, \quad Z = 10^{pK_a - pH_{rbc}}$$
- **Diprotic Bases**:
  $$X = 10^{pK_{a2} - pH_{iw}} + 10^{pK_{a1} + pK_{a2} - 2\cdot pH_{iw}}$$
  $$Y = 10^{pK_{a2} - pH_{pls}} + 10^{pK_{a1} + pK_{a2} - 2\cdot pH_{pls}}$$
  $$Z = 10^{pK_{a2} - pH_{rbc}} + 10^{pK_{a1} + pK_{a2} - 2\cdot pH_{rbc}}$$
- **Zwitterions / Ampholytes**:
  $$X = 10^{pK_{a,base} - pH_{iw}} + 10^{pH_{iw} - pK_{a,acid}}$$
  $$Y = 10^{pK_{a,base} - pH_{pls}} + 10^{pH_{pls} - pK_{a,acid}}$$
  $$Z = 10^{pK_{a,base} - pH_{rbc}} + 10^{pH_{rbc} - pK_{a,acid}}$$

#### Blood Cell Partitioning & Acidic Phospholipid Affinity ($Ka_{AP}$)
For moderate-to-strong bases ($pK_a \ge 7$):
$$Kpu_{bc} = \frac{HCT - 1 + B2P}{f_u \cdot HCT}$$
$$\Delta_{AP} = Kpu_{bc} - \left(\frac{1+Z}{1+Y}\right) f_{iw,bc} - \frac{10^{\log P} f_{nl,bc} + (0.3 \cdot 10^{\log P} + 0.7) f_{np,bc}}{1+Y}$$
$$Ka_{AP} = \max(0, \Delta_{AP}) \cdot \frac{1+Y}{AP_{bc} \cdot Z}$$

Where:
- $HCT = 0.47$ (Hematocrit)
- $f_{iw,bc} = 0.63$ (RBC intracellular water)
- $f_{nl,bc} = 0.0012$ (RBC neutral lipids)
- $f_{np,bc} = 0.0033$ (RBC neutral phospholipids)
- $AP_{bc} = 0.57\text{ mg/g}$ (RBC acidic phospholipids)

#### Protein Affinity Constant ($Ka_{PR}$)
For neutrals, acids, and weak bases:
$$Ka_{PR} = \max\left(0, \frac{1}{f_u} - 1 - \frac{10^{\log P} f_{nl,p} + (0.3 \cdot 10^{\log P} + 0.7) f_{np,p}}{1+Y}\right)$$

Where plasma lipid fractions are:
- $f_{nl,p} = 0.0032$ (Plasma neutral lipids)
- $f_{np,p} = 0.0021$ (Plasma neutral phospholipids)

#### Tissue-to-Plasma Partition Coefficient ($K_{p,T}$)
For Adipose (Fat):
$$\log P_{ow} = 1.115 \cdot \log P - 1.35, \quad P_{eff} = 10^{\log P_{ow}}$$
For all other tissues:
$$P_{eff} = 10^{\log P}$$

Binding term:
$$\text{Binding} = \begin{cases} 
\frac{Ka_{AP} \cdot AP_T \cdot X}{1+Y} & \text{if moderate/strong base} \\
Ka_{PR} \cdot PR_{T/P,Lipoprotein} & \text{if neutral} \\
Ka_{PR} \cdot PR_{T/P,Albumin} & \text{if acid or weak base}
\end{cases}$$

Final Partition Coefficient:
$$K_{p,T} = \left[ \left(\frac{1+X}{1+Y}\right) f_{iw,T} + f_{ew,T} + \frac{P_{eff} \cdot f_{nl,T} + (0.3 P_{eff} + 0.7) f_{np,T}}{1+Y} + \text{Binding} \right] \cdot f_u$$

### 2. Whole-Body PBPK Ordinary Differential Equation System

Let $A_k$ denote the drug amount (mg) in compartment $k$, $V_k$ the physiological compartment volume (L), and $Q_k$ the blood flow rate (L/h).

#### Non-Eliminating Parallel Organs
For organ $T \in \{\text{Bone}, \text{Brain}, \text{Fat}, \text{Gonads}, \text{Heart}, \text{Muscle}, \text{Skin}\}$:
$$C_{v,T} = \frac{A_T}{V_T \cdot K_{p,T}}$$
$$\frac{dA_T}{dt} = Q_T \left( C_{art} - C_{v,T} \right)$$

#### Splanchnic Organs & Portal Vein
For $j \in \{\text{Pancreas}, \text{Spleen}, \text{Stomach}, \text{LargeIntestine}\}$:
$$\frac{dA_j}{dt} = Q_j \left( C_{art} - C_{v,j} \right)$$
For Small Intestine (with oral absorption rate $k_a$ and gut lumen amount $A_{lumen}$):
$$\frac{dA_{lumen}}{dt} = -k_a \cdot A_{lumen}$$
$$\frac{dA_{small\_int}}{dt} = Q_{small\_int} \left( C_{art} - C_{v,small\_int} \right) + k_a \cdot A_{lumen}$$

Portal vein blood concentration entering the liver:
$$Q_{pv} = Q_{pancreas} + Q_{spleen} + Q_{stomach} + Q_{small\_int} + Q_{large\_int}$$
$$C_{pv} = \frac{\sum_{j} Q_j \cdot C_{v,j}}{Q_{pv}}$$

#### Liver Compartment
Receives blood from the hepatic artery ($Q_{ha}$) and the portal vein ($Q_{pv}$):
$$Q_{liv} = Q_{ha} + Q_{pv}$$
$$C_{v,liv} = \frac{A_{liv}}{V_{liv} \cdot K_{p,liv}}$$
$$R_{elim,hep} = CL_{hep} \cdot C_{v,liv}$$
$$\frac{dA_{liv}}{dt} = Q_{ha} \cdot C_{art} + Q_{pv} \cdot C_{pv} - Q_{liv} \cdot C_{v,liv} - R_{elim,hep}$$
$$\frac{dA_{elim,hep}}{dt} = R_{elim,hep}$$

#### Kidney Compartment
$$C_{v,kid} = \frac{A_{kid}}{V_{kid} \cdot K_{p,kid}}$$
$$R_{elim,ren} = CL_{ren} \cdot C_{art}$$
$$\frac{dA_{kid}}{dt} = Q_{kid} \left( C_{art} - C_{v,kid} \right) - R_{elim,ren}$$
$$\frac{dA_{elim,ren}}{dt} = R_{elim,ren}$$

#### Venous Pool
$$Q_{co} = \sum_{i} Q_i + Q_{liv} + Q_{kid}$$
$$\frac{dA_{ven}}{dt} = \sum_{i} Q_i C_{v,i} + Q_{liv} C_{v,liv} + Q_{kid} C_{v,kid} - Q_{co} C_{ven} + R_{infusion}(t)$$

#### Lung
Receives the complete cardiac output ($Q_{co}$) from the venous pool:
$$\frac{dA_{lung}}{dt} = Q_{co} \left( C_{ven} - C_{v,lung} \right)$$

#### Arterial Pool
Receives blood from the lungs and supplies all systemic organs:
$$\frac{dA_{art}}{dt} = Q_{co} \left( C_{v,lung} - C_{art} \right)$$

### 3. Mass Conservation Proof
Summing all differential equations:
$$\frac{d}{dt} \left[ A_{ven} + A_{lung} + A_{art} + \sum_{T} A_T + A_{elim,hep} + A_{elim,ren} + A_{lumen} \right] = R_{infusion}(t)$$
Integrating with respect to time yields:
$$\text{TotalMass}(t) = \text{Dose}_{\text{administered}}(t)$$
conserved identically to within machine precision.

## References

1. **Rodgers T, Leahy D, Rowland M.** (2005). Physiologically based pharmacokinetic modeling 1: predicting the tissue distribution of moderate-to-strong bases. *Journal of Pharmaceutical Sciences*, 94(6), 1259–1276. https://doi.org/10.1002/jps.20322
2. **Rodgers T, Rowland M.** (2006). Physiologically based pharmacokinetic modelling 2: predicting the tissue distribution of acids, very weak bases, neutrals and zwitterions. *Journal of Pharmaceutical Sciences*, 95(6), 1238–1257. https://doi.org/10.1002/jps.20502
3. **Rodgers T, Rowland M.** (2007). Mechanistic approaches to volume of distribution predictions: understanding the processes. *Pharmaceutical Research*, 24(5), 918–933. https://doi.org/10.1007/s11095-006-9210-3
4. **International Commission on Radiological Protection (ICRP).** (2002). Basic Anatomical and Physiological Data for Use in Radiological Protection: Reference Values. *ICRP Publication 89*, Annals of the ICRP, 32(3-4).
5. **Willmann S, Lippert J, Sevestre M, Solodenko J, Fois F, Schmitt W.** (2003). PK-Sim: A novel tool for physiologically based pharmacokinetic modeling and simulation. *Biosilico*, 1(4), 121–124.
6. **Willmann S, Lippert J, Schmitt W.** (2005). From physicochemistry to absorption and distribution: Predictive PBPK modeling and its role in drug development. *Expert Opinion on Drug Metabolism & Toxicology*, 1(1), 159–168.
7. **Jones H, Rowland-Yeo K.** (2013). Basic concepts in physiologically based pharmacokinetic modeling in drug discovery and development. *CPT: Pharmacometrics & Systems Pharmacology*, 2(8), e63. https://doi.org/10.1038/psp.2013.41

## Repository Sources

Formulas and parameters were extracted directly from the Open Systems Pharmacology suite:
- **`Open-Systems-Pharmacology/PK-Sim`**:
  - `src/Db/Dump/PKSimDB_dump.txt` and `src/Db/PKSimDB.sqlite`
  - Table: `tab_calculation_method_rates`
  - Table: `tab_container_parameter_curves` (European ICRP 2002 population curves)
  - Table: `tab_container_parameter_values` (tissue lipid and water volume fractions)
- **`Open-Systems-Pharmacology/OSPSuite.SimModel`**:
  - `src/OSPSuite.SimModelNative/src/CppODEExporter.cpp`
  - `src/OSPSuite.SimModelNative/src/MatlabODEExporter.cpp`

## Acknowledgments

We gratefully acknowledge the open-source community behind the [Open Systems Pharmacology (OSP) Suite](https://open-systems-pharmacology.org/) for their foundational work, transparent database structures, and mathematical implementations that make rigorous, reproducible PBPK modeling possible.

## License

This project is licensed under the **GNU General Public License v3.0 (GNU GPL v3.0)**.

You may use, modify, and distribute this software under the terms of the GPL v3.0. See the [LICENSE](LICENSE).
