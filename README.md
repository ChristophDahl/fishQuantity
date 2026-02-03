# FishQuantity: shoal-size choice and stimulus-control analysis (MATLAB)

This repository contains the MATLAB pipeline used to analyse zebrafish shoaling-choice
behaviour and to quantify image-derived stimulus metrics (coverage, density, spacing, etc.)
used as control covariates in the paper.

The code is organised as a small set of **stage scripts** (reproducible entry points) and
a set of **core functions**.

## What this repo contains

### Pipeline runner
- `runFishQuantityPipeline.m`

### Stage scripts (called by the runner)
- `stage1_smallFish_profiles_AB.m`
- `stage1_smallFish_laterality_export.m`
- `stage2_smallFish_laterality_LME.m`
- `stage1_largeFish_profiles_AB.m`
- `stage2_largeFish_laterality_LME.m`
- `stage3_singleContrast_smallFish.m`
- `stage3_singleContrast_largeFish.m`

### Core functions
- `fitOTSvsANSModels.m`
- `summarize_small_large_runs.m`
- `fit_laterality_condition_models.m`
- `interactionFishQuantity.m`
- `analyze_fish_pages_range.m`
- `make_stimulus_pair_figure.m`

## Requirements
- MATLAB R2021b+ (recommended)
- Toolboxes typically required:
  - Statistics and Machine Learning Toolbox (mixed models / regression)
  - Image Processing Toolbox (stimulus page metrics)

## Data and folder structure

This repo is **code only**. Data are not stored in GitHub.

By default the pipeline expects a local working directory structured like:

