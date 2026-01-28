# Visualization tool for interferometric data

# Table of Contents
* [About](#About)
* [Dependencies](#dependencies)
* [Using the tool](#Using-the-tool)
* [Acknowledgements](#acknowledgements)


## About
Visualization tool for interferometric data do standard calibration and data visualization of interfereometric data. 

## Dependencies

- python 3
  - casacore 3.7.1
  - numpy 7.2.0
  - matplotlib 3.10.8
  - astropy 7.2.0
  - tqdm 4.67.1
  - baseband 4.3.0

## Using the tool

| **Scripts**  | **Description**                                                                     |
|--------------|-------------------------------------------------------------------------------------|
| aips.py      | Applay baseband, flux and elevation calibration to MS files. Show dynamical spectra |
| plot_vdif.py | Perform FFT on VDIF file and show results.                                          |

# Acknowledgements
This software was written by Jānis Šteinbergs. If you make use of this software to get results that appear in a publication or presentation please include this acknowledgement: &quot;We have made use of Visualization tool for interferometric data, a tool developed by Jānis Šteinbergs. &quot;

This work was supported by Latvian Council of Science Project "Multi-wavelength Study of Quasi-Periodic Pulsations in Solar and Stellar Flares" Nr.: Izp-2022/1-0017.

