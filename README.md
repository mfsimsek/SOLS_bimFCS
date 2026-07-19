# SOLS_bimFCS
Binned-imaging fluorescence correlation spectroscopy measurements on single objective light sheet microscope

REFERENCES:
1- Simsek, Thesis, 2015. Proquest ID: 10141671.
2- Weixian, Simsek, Pralle, Methods, 2018. DOI: 10.1016/j.ymeth.2018.02.019.
3- Huang, Simsek, et al. PLoS One, 2015. DOI: 10.1371/journal.pone.0121777.

Main analysis file is SOLS_bimFCS_Analysis_SimsekLab.
Outputs are:
1- png plot for autocorrelation functions for each super-pixel overlaid in a single plot as well as two-component diffusive behaviour (D1 fast, D2 slow) fits.
2- png plot for linear fitted FCS difussion law plot for all analyzed super-pixels
3- png plot for linear fitted FCS diffusion law plot for super-pixels after smallest 5.
4- spreadsheet file with all the linear fit information and fitted results.

This file relies on following supplementary files to run:
1- lintolog.m : this file resamples acquired time data in logarithmic scale.
2- bin2d_overlap.m : this file calculates convoluted intensities for each super-pixel by scanning the 2D pixel array.

Next two files are used to measure w2 values and plug those in to the main analysis file.
3- beamwaister_FCS.m : tis file calculates w2 effective beam waist values for each super-pixel using measured PSF value and camera pixel size.
4- beprofi.m : this is the square detection profile PSF convolution function used in beamwaister_FCS.m.
