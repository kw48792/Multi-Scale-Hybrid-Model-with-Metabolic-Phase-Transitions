🧬 Multi-Scale Hybrid Model with Metabolic Phase Transitions

Multi-scale hybrid modeling framework to predict CHO cell culture dynamics with metabolic phase transitions

This repository contains Python scripts and data required to reproduce the modeling, prediction, and visualization results described in the manuscript “Multi-Scale Hybrid Modeling to Predict Cell Culture Process with Metabolic Phase Transitions.”


📂 Repository Overview
DataViz.py

Generates figures related to cell culture dynamic profiles and model predictions.

OTR.py

Applies the Savitzky–Golay (SG) filter to smooth the oxygen uptake rate (OUR) trajectory.

pH.py

Applies the Savitzky–Golay (SG) filter to smooth the pH trajectory.

State_Shift_Pred_A.py

Uses params_A1.pkl to predict dynamic trajectories for Case A.

State_Shift_Pred_B.py

Uses params_B.pkl to predict dynamic trajectories for Case B.

State_Shift_Pred_C.py

Uses params_C.pkl to predict dynamic trajectories for Case C.

State_Shift_Pred_PI.py

Generates prediction errors (WAPE) and 95% prediction intervals (PI) across all cases.


⚙️ Requirements

Python ≥ 3.8

numpy, scipy, pandas, matplotlib, lmfit
