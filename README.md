# tensor_fUS



This repository contains the MATLAB code for tensor decomposition applied to functional ultrasound (fUS) imaging data for feature extraction. The code is organized into three main scripts that demonstrate the Canonical Polyadic Decomposition (CPD), Block Term Decomposition (BTD), and an improved BTD approach.

## Prerequisites

Before you run the scripts, ensure you have MATLAB installed with the following specific tools and functions:

- [TensorLab](https://www.tensorlab.net/) library (for `assignment_4_fUS_BTD_improved` script) (already included in our upload)

## Usage

To run the code, follow these steps:

1. Replace the `given/display_brain_img.m` script with the `given/display_brain_img_sub.m`. This step is crucial as it modifies the plotting function to improve visualization with transparent background and resolve the issue with negative correlations appearing as black on the background.

2. Start MATLAB and navigate to the directory containing the scripts.

3. Run the scripts in the following order:
   - `assignment_4_fUS_CPD.m` for CPD decomposition.
   - `assignment_4_fUS_BTD.m` for BTD decomposition.
   - `assignment_4_fUS_BTD_improved.m` for the improved BTD approach.

Note: Some `.mat` files are provided to bypass repetitive experiments. These files contain precomputed results that can be loaded directly into the workspace.

## Files

- `assignment_4_fUS_CPD.m`: Script for CPD decomposition.
- `assignment_4_fUS_BTD.m`: Script for BTD decomposition.
- `assignment_4_fUS_BTD_improved.m`: Script for the improved BTD approach.
- `given/`
  - `display_brain_img_sub.m`: Modified plotting function.
- `*.mat`: Precomputed results for quick loading.

## Reference

1. Hunyadi, B., Camps, D., Sorber, L., et al. (2014). Block term decomposition for modelling epileptic seizures. *EURASIP Journal on Advances in Signal Processing, 2014*(1), 1-19.

2. De Lathauwer, L. (2008). Decompositions of a higher-order tensor in block terms—Part II: Definitions and uniqueness. *SIAM Journal on Matrix Analysis and Applications, 30*(3), 1033-1066.

3. Chatzichristos, C., Kofidis, E., Morante, M., et al. (2019). Blind fMRI source unmixing via higher-order tensor decompositions. *Journal of Neuroscience Methods, 315*, 17-47.

4. Vervliet, N., Debals, O., & De Lathauwer, L. (2016). Tensorlab 3.0—Numerical optimization strategies for large-scale constrained and coupled matrix/tensor factorization. In *2016 50th Asilomar Conference on Signals, Systems and Computers* (pp. 1733-1738).