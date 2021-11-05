# MultiModes
Author: David Pamos Ortega - University of Granada (UGR) -

# What it is
MultiModes is a python code to extract the most significant frequencies of a sample of variable stars

# Input
- Directory with light curves in format .fits
- ini.txt with the initial parameters: 
  sim_fit_n: Number of simultaneous peaks to be fit before extracting to the original light curve for obtaining the residual: 20 by default
  max_freq: Maximum value of the analysed frequencies domain: 100 muHz by default
  os_ratio: oversampling factor: 5 by default
  stop: Stop criterion, FAP or SNR: SNR by default
  min_snr: Minimum signal to noise ratio: 4 by default (Breger 1993)
  max_fap: Maximum value of the False Alarm Probability: 0.01 by default (Balona et al. 2014)
  tail_per: Minimum  frequency of the tail of the periodogram: 80 muHz by default
  
# Output
For every light curve, a directory with:
- file best_modes.dat, containing the values of the most significant frequencies, amplitudes, phases, corresponding errors and FAPs/SNRs
- file lc.dat, the light curve in format .dat for using with other codes, such as SigSpec (Reegen 2007)
- file pg.dat, the periodogram of the original light curve
- LC.png, the plot of the light curve
- LS.png, the periodogram
- LS_n.png, the periodogram obtained every number n of extracted peaks
- res.dat, with the final residual after extracting the most significant frequencies

# Pre-Installed Packages Requirements
- python 3.8.5
- numpy 1.19.2
- matplotlib 3.3.2
- pandas 1.1.2
- astropy 4.0.2
- lmfit 1.0.2
- scipy 1.5.2

