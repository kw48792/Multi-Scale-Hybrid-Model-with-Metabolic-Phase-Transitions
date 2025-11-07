# Multi-Scale Hybrid Model with Metabolic Phase Transitions

This repository contains the code and data used in the journal paper: **"Multi-Scale Hybrid Modeling to Predict Cell Culture Process with Metabolic Phase Transitions"**

## Overview

This project implements a multi-scale hybrid modeling approach for predicting cell culture processes with metabolic phase transitions. The model combines mechanistic and data-driven approaches to capture complex cellular behaviors during bioprocessing.


## Installation

1. Clone this repository:
```bash
git clone https://github.com/kw48792/Multi-Scale-Hybrid-Model-with-Metabolic-Phase-Transitions.git
cd Multi-Scale-Hybrid-Model-with-Metabolic-Phase-Transitions
```

2. Install the required dependencies:
```bash
pip install -r requirements.txt
```

## Usage

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

## Data

The datasets generated and analyzed during this study are available from the corresponding author upon reasonable request for academic research purposes.

## Citation

If you use this code in your research, please cite:

```bibtex
@article{wang2024multi,
  title={Multi-scale kinetics modeling for cell culture process with metabolic state transition},
  author={Wang, Keqi and Harcum, Sarah W and Xie, Wei},
  journal={arXiv preprint arXiv:2412.03883},
  year={2024}
}
```

## Contact

For questions or issues, please open an issue on GitHub or contact the authors.

## Acknowledgments

This work was developed to support the journal publication and is made available in response to reviewer requirements for reproducibility.
