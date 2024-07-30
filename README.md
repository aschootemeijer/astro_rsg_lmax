### in this project, we compare OBSERVATIONAL DATA of red supergiants (RSGS) in different galaxies to STELLAR MODELS.
### specifically, we perform RSG model population synthesis to compare the predicted nr of bright RSGs to the observed nr

## OBSERVATIONAL DATA
We have the following sample of six galaxies with individually resolved stars, including RSGs. 
NGC300 from Dalcanton+2009; NGC5352 and NGC4395 from Sabbi+2019; SMC and LMC from Davies+2018; IZw18 from Hirschauer+2024
In NGC300, NGC5352, and NGC4395 we still need to identify the RSGs. We use a neural network to do this. 

## STELLAR MODELS:
We obtain the interpolated BOOST stellar evolution tracks from http://galaxy.asu.cas.cz/pages/boost
- convert_to_usable_data.py: converts BOOST data to a csv that is more convenient to use for us
- condense_data.py/: data are condensed because we do not need all columns and initial masses separated as densely as by 0.001 dex.

