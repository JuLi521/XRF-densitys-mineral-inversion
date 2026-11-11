# XRF-densitys-mineral-inversion
# ICDP-DIVE 5071_1_B: XRF–density mineral inversion (10-sample example)

This repository provides a compact demonstration of the original MATLAB workflow for Bayesian mineral-fraction inversion of kinzigite using XRF oxide concentrations and bulk core density.

## Files

- `XRF_Kinzigite_Publication_Main.m`: main inversion and visualization script.
- `mcmces_xrf1.m`: **unchanged** from the author-supplied original. 
- `DIVE_5071_1_B_KINZ_10samples.csv`: 10 real representative samples from 100–200 m in borehole 5071_1_B. Depth alignment and density interpolation have been performed in advance;

## Execution

Requires MATLAB with Statistics and Machine Learning Toolbox (`ksdensity`). Place all three files in the current directory and run:

```matlab
XRF_Kinzigite_Publication_Main
```

The script outputs four figure windows (posterior PDF depth tracks; XRF measured versus reconstructed; density measured versus reconstructed; grouped mineral volume fractions), and `Inversion_results_10samples.csv`.

## Inputs and model

XRF input columns (wt.%): SiO2, K2O, Al2O3, CaO, MgO, MnO, Na2O, FeO, Fe2O3. The first eight oxides enter the likelihood, exactly as in the original working script; Fe2O3 is reconstructed but not independently fitted. Core density is supplied in g/cm3.

The original seven end-members are phlogopite, annite, quartz, albite, anorthite, garnet and K-feldspar. Their oxide reference compositions and densities are preserved. The inverse parameters are mineral **weight fractions**, constrained to the simplex by normalization of the initial walkers and affine proposals. Volume fractions are computed afterwards by dividing each weight fraction by the corresponding mineral density and normalizing.

The original likelihood convention, 100 walkers, 2000 steps, 60% burn-in, StepSize=5 and ThinChain=1 are maintained. As in the original main script, XRF residual squares are divided by the *sample standard deviations* of the example dataset, not their squared values. This is a legacy model choice, not a validated Gaussian measurement-error model. Thus, quantitative uncertainty interpretation needs independent validation before publication.

## Sampler attribution and methodological caveat

Methodological background:
Cheng, L., Jin, G., Michelena, R., & Tura, A. (2022). 
Practical Bayesian Inversions for Rock Composition and Petrophysical Endpoints in Multimineral Analysis. 
SPE Reservoir Evaluation & Engineering, 25(04), 849–865. https://doi.org/10.2118/210576-PA


## Figures

1. **Posterior distributions:** KDE of each mineral's retained weight-fraction samples evaluated on 100 bins, mapped as a depth track with the original white-to-color palette and black posterior-mean stars. Color intensity is the KDE discretely normalized per sample, not a continuous probability density with physical units.
2. **Oxide reconstruction:** black stars show predicted concentrations, red circles the 10 measured samples.
3. **Bulk density reconstruction:** black stars show predicted mixture density, red circles the observed density.
4. **Mineral volume fractions:** the original grouping Bt (Phl+Ann), Qtz, Pl (Ab+An), Grt and Kfs, with the original black filled-marker depth-track style.

No comparison is made to logging inversion or thin-section abundance files, since these data were not included in the 10-sample release. No unavailable measurements are fabricated.

## Reproducibility / limitations

The random number generator is seeded with `rng(2026,'twister')`. This script has been statically checked for data dimensions, column names and nonzero oxide standard deviations, but **has not been executed in MATLAB** in this environment. A local MATLAB test and MCMC convergence assessment remain necessary before publication. The 10 points support a methods demonstration, not full-paper figure reproduction.

