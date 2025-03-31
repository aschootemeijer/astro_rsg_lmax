##### in this project, we compare observational data of red supergiants (RSGS) in different galaxies to stellar models. Specifically, we perform RSG model population synthesis to compare the predicted nr of bright RSGs to the observed nr.

###### OBSERVATIONAL DATA
We have the following sample of six galaxies with individually resolved stars, including RSGs. 
NGC300 from Dalcanton+2009; NGC5352 and NGC4395 from Sabbi+2019; SMC and LMC from Davies+2018; IZw18 from Hirschauer+2024
In NGC300, NGC5352, and NGC4395 we still need to identify the RSGs. We use a neural network to do this. Also, as a test we compare the values we get for the luminosity of RSGs in the Small Magellanic Cloud (SMC) to those in earlier work by Davies+2018. We do this with the script compare_methods.py and compare_methods.R. The R script does the same thing as the python script.
- neural_network.py: assign probability that a star is an RSG (trained on SMC data).
- compare_methods.py: calculate luminosity using visible wavelength and IR data, and compare to earlier results of Davies+2018.
- compare_methods.R: same as the python script, but is written from scratch in R as a test. It obtains identical results.

###### STELLAR MODELS:
We obtain the interpolated BOOST stellar evolution tracks from http://galaxy.asu.cas.cz/pages/boost
- convert_to_usable_data.py: converts BOOST data to a csv that is more convenient to use for us
- condense_data.py: data are condensed because we do not need all columns and initial masses separated as densely as by 0.001 dex.

######  OBSSERVATIONS vs MODEL comparison
We calculate number of observed stars in the luminosity range 5.0<logL<5.4 and make a theoretical population of stars with the same number of stars in the range range 5.0<logL<5.4.
- model_vs_obs_hists.py compares the luminosity distributions
- z_vs_brrsg.py compares the luminosity of the brighest star in the theory population to the observed population, and investigates how this depends on metallicity (Z).
