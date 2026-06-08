# Hect2026_voicelikeness

**Acoustic regularities define perceptual and cortical representations of voice-likeness**

Jasmine L. Hect, Kyle M. Rupp, Avniel Singh Ghuman, Lori L. Holt, Taylor J. Abel

*Current Biology* (2026) — CURRENT-BIOLOGY-D-26-00179

---

## Overview

This repository contains all analysis code used to produce the figures and tables in Hect et al. (2026). The study examines whether voice perception is organized along a continuous, category-defining acoustic dimension and whether this structure is hierarchically encoded across the human auditory ventral stream, using intracranial EEG (iEEG/sEEG) recordings from 39 participants and perceptual ratings from 253 online listeners.

**Key findings:**
- A linear discriminant axis derived from acoustic features of natural sounds captures a continuous, perceptually valid dimension of voice-likeness that generalizes to novel synthetic stimuli without true category membership
- Neural population responses along the auditory ventral stream mirror this acoustic-perceptual organization, with graded voice-likeness structure emerging hierarchically from primary auditory cortex through belt, parabelt, and superior temporal cortex
- Acoustic-neural alignment is significant across all six auditory ventral stream regions, increasing systematically in magnitude from primary to association cortex
- Continuous perceptual ratings better explain neural responses than binary category labels in mid- and higher-level auditory cortex

---

## Repository structure

```
Hect2026_voicelikeness/
├── README.md
├── main_analysis.m                         — full analysis pipeline (run this to reproduce all figures and tables)
├── savePublicDataArchive.m                 — function to save results archive for OSF/Zenodo deposit
├── functions/
│   ├── loadSettings.m                      — returns gcs struct with local file paths
│   ├── importGorilla.m                     — loads and parses Gorilla behavioral data
│   ├── getDemographics.m                   — returns participant demographic table from IDs
│   ├── getStimAcFts.m                      — extracts TCD/FCD + GeMAPS features for all stimuli
│   ├── defineRois.m                        — assigns HCPex parcels to auditory ventral stream ROIs
│   ├── plotChansSurf.m                     — plots channel values on MNI inflated surface
│   ├── sig2star.m                          — converts FDR q-value to significance star string
│   ├── plot_acoustic_feature_hist_grid.m   — grid of overlaid histograms for acoustic features
│   └── run_yamnet_specs.m                  — runs YAMNet classifier on stimulus directory
├── preprocessing/
│   ├── process_sEEG.m                      — epochs, baseline-normalizes, and saves ld_vl and ld_gsp
│   └── extract_acoustic_features.py        — openSMILE Python wrapper for GeMAPS feature extraction
├── external/                               — vendored copies of File Exchange dependencies
│   ├── cbrewer2/                           — colorbrewer2 colormaps (https://github.com/scottclowe/cbrewer2)
│   └── cmocean/                            — perceptually uniform ocean colormaps (https://github.com/chadagreene/cmocean)
└── data/
    └── README_data.md                      — instructions for downloading data from OSF/Zenodo
```

---

## Data

Preprocessed neural data, behavioral ratings, acoustic features, and all analysis results are deposited on OSF/Zenodo:

> **DOI: [to be assigned upon acceptance]**
> **URL: [to be assigned upon acceptance]**

Download the archive and place the `.mat` files in a `data/` directory at the repository root before running the analysis pipeline. The archive contains the following files:

| File | Contents |
|---|---|
| `Hect2026_behavioral.mat` | Perceptual ratings, stimulus sort indices, category labels, YAMNet predictions |
| `Hect2026_acoustic.mat` | Acoustic features (GeMAPS + TCD + FCD), LDA model, projections, permutation results |
| `Hect2026_neural_lfp.mat` | Preprocessed LFP data matrices + channel metadata (MNI ROI labels) |
| `Hect2026_results_roi.mat` | ROI-level LDA accuracy, Spearman rho, bootstrap CIs, FDR q-values, model comparison |
| `Hect2026_results_singlechan.mat` | Single-channel rho, LDA accuracy, model comparison, surface coordinates |
| `Hect2026_results_timedomain.mat` | Sliding-window LDA accuracy, time-resolved rho, time-varying delta AIC |
| `Hect2026_lme.mat` | LME long-format tables and fitted coefficients |
| `Hect2026_metadata.mat` | Provenance: date, MATLAB version, variable manifest |

**Note on raw iEEG recordings:** Raw neural recordings are deposited separately on OSF/Zenodo with IRB-compliant anonymization (UPMC IRB protocol STUDY20030060). The preprocessed `data_matrix` fields in `Hect2026_neural_lfp.mat` are baseline-normalized, repeat-averaged, and sorted by stimulus rating rather than by participant, and contain no protected health information.

**Natural sound stimuli:** The Belin voice localizer corpus (144 natural sounds) is available from the original authors (Belin et al., 2000; https://doi.org/10.1038/35002078). These stimuli are not redistributed here.

**Synthetic stimuli:** The 250 synthetic sound textures generated for this study are deposited on OSF/Zenodo alongside the neural data.

---

## Requirements

### MATLAB (primary analysis environment)

- MATLAB R2024b
- Statistics and Machine Learning Toolbox (for `fitcdiscr`, `fitlme`, `cvpartition`)
- Bioinformatics Toolbox (for `mafdr`)
- Audio Toolbox (for `yamnet`, `yamnetPreprocess`)
- Signal Processing Toolbox

All toolboxes were tested with MATLAB R2024b. Earlier versions may work but are not guaranteed.

### Python (acoustic feature extraction only)

- Python 3.10
- openSMILE >= 3.0 (https://github.com/audeering/opensmile)

Python is only required to re-extract GeMAPS acoustic features from raw audio. If you are working from the deposited `.mat` archive, Python is not needed.

### External MATLAB dependencies (vendored in `external/`)

| Package | Version | Source |
|---|---|---|
| cbrewer2 | 1.0.0 | https://github.com/scottclowe/cbrewer2 |
| cmocean | 2.02 | https://github.com/chadagreene/cmocean |

These are included in the `external/` directory and do not require separate installation. Add them to your MATLAB path with `addpath(genpath('external'))`.

---

## Setup

**1. Clone the repository**

```bash
git clone https://github.com/pbe-lab/Codeshare_Hect2026.git
cd Hect2026_voicelikeness
```

**2. Download the data archive**

Download all `.mat` files from the OSF/Zenodo repository (DOI above) and place them in `data/`:

```
data/
├── Hect2026_behavioral.mat
├── Hect2026_acoustic.mat
├── Hect2026_neural_lfp.mat
├── Hect2026_results_roi.mat
├── Hect2026_results_singlechan.mat
├── Hect2026_results_timedomain.mat
├── Hect2026_lme.mat
└── Hect2026_metadata.mat
```

**3. Configure local paths**

Edit `functions/loadSettings.m` to point to your local data directory:

```matlab
function gcs = loadSettings()
    gcs.dataDir    = '/path/to/data/';           % directory containing .mat archive files
    gcs.fOneDrive  = '/path/to/stimuli/';        % directory containing stimulus .wav files (if re-running preprocessing)
    gcs.outputDir  = '/path/to/output/figures/'; % directory for saving figures
end
```

**4. Add paths and run**

```matlab
addpath(genpath('functions'))
addpath(genpath('external'))
main_analysis
```

The full pipeline takes approximately 4-6 hours to run on a modern workstation due to permutation testing (10,000 iterations per ROI for the neural LDA analyses). To reproduce figures only without re-running permutations, load the pre-computed results from the archive and run the visualization sections of `main_analysis.m` directly.

---

## Reproducing specific figures and tables

The `main_analysis.m` script is organized into labeled sections using MATLAB's `%%` cell structure. Each section corresponds to a specific figure or table. Navigate to the relevant section using MATLAB's cell mode (Ctrl+Enter to run a single cell).

| Figure / Table | Script section |
|---|---|
| Figure 1 (experimental design, ratings) | `%% Behavioral data and ratings` |
| Figure 2 (acoustic LDA) | `%% Acoustic Linear discriminant analysis` |
| Figure 3 (neural LDA, brain maps) | `%% Neural Linear discriminant analysis` and `%% RANK ORDER ACOUSTIC-NEURAL CORRELATIONS` |
| Figure 4 (time-resolved LDA) | `%% TIME VARYING Neural LDA` |
| Figure 5 (rank-based rho bar plots) | `%% Rank-based trend strength across ROIs` |
| Figure 6 (model comparison) | `%% Model comparison: binary vs continuous` |
| Supplementary Fig. 1 (acoustic feature distributions) | `%% plot big grid of histograms` |
| Supplementary Fig. 2 (LDA weights) | Commented block in acoustic LDA section |
| Supplementary Fig. 3 (channel coverage) | `%% Brain maps` |
| Supplementary Fig. 4 (representative LFP electrodes) | `%% example channels` (natural sounds block) |
| Supplementary Fig. 5 (ROI-averaged LFP) | `%% make LME tables` — ROI mean waveform plots |
| Supplementary Fig. 6 (BHA analyses) | BHA analysis script (separate file, see below) |
| Supplemental Table 5 (LDA accuracy) | Printed to command window by neural LDA section |
| Supplemental Table 6 (LME coefficients) | `%% fit LMEs and plot coefficients` |
| Supplemental Table 7 (rho, within-category) | Printed to command window by bootstrap section |
| Supplemental Table 8 (model comparison AIC) | Printed to command window by model comparison section |
| Supplemental Table 9 (BHA accuracy) | BHA analysis script |
| Supplemental Table 10 (BHA model comparison) | BHA analysis script |

**BHA analyses** are run from a separate script (`main_analysis_bha.m`) that is structured identically to `main_analysis.m` but loads the BHA data files (`dataIn_daniel_hga_*.mat`) instead of the LFP files. Switch between LFP and BHA by changing the load call at the top of the script — the commented-out load lines are already in place.

---

## Preprocessing pipeline

If you wish to reproduce the analysis from raw iEEG recordings rather than from the preprocessed `.mat` archive, run the preprocessing pipeline first:

```matlab
% From MATLAB, after downloading raw data from OSF/Zenodo:
process_sEEG('voiceLocalizer')   % produces ld_vl struct
process_sEEG('gsp')              % produces ld_gsp struct
```

Acoustic feature extraction from raw `.wav` files requires Python and openSMILE:

```bash
# From the repository root:
python preprocessing/extract_acoustic_features.py \
    --stim_dir /path/to/stimuli/ \
    --output_dir /path/to/data/
```

This produces the GeMAPS feature matrices that are then loaded and combined with TCD/FCD coefficients by `getStimAcFts.m`.

---

## Citation

If you use this code or data, please cite:

> Hect, J.L., Rupp, K.M., Ghuman, A.S., Holt, L.L., & Abel, T.J. (2026). Acoustic regularities define perceptual and cortical representations of voice-likeness. *Current Biology*. https://doi.org/[to be assigned]

---

## License

Code in this repository is released under the MIT License. See `LICENSE` for details.

External dependencies (`cbrewer2`, `cmocean`) are subject to their own licenses, which are included in their respective subdirectories under `external/`.

---

## Contact

**Lead contact:** Jasmine L. Hect — jasminehect@gmail.com

Department of Neurological Surgery, University of Pittsburgh

For questions about the iEEG data collection or clinical protocols, contact the corresponding author: Taylor J. Abel, MD — abeltj@upmc.edu

---

## Acknowledgements

This work was funded by 1F30DC021342-01 (JLH), R21DC019217-01A1 (TJA and LLH), R01DC013315-07 (TJA), R01MH132225 (ASG), and T32GM008208. The content is solely the responsibility of the authors and does not represent the official views of the National Institutes of Health. We thank Emily Harford, Sreekrishna Ramakrishnapillai, Arish Alreja, and Mary Kate Richey for their support of this project and assistance with data collection.
