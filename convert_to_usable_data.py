import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from matplotlib import rcParams
rcParams['font.family'] = 'STIXGeneral'
c1,c2,c3,c4,c5,c6,c7,c8,c9,c10 = '#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf'
fig,ax=plt.subplots(1,1,figsize=(6,3.6))
import time
start = time.time()

gal = 'MW'    # Data of which galaxy?   (MW, LMC, SMC, IZw18)

# THIS WE USE TO CREATE A USABLE DATA FILE OUT OF THE ASCII FILE
d = pd.read_csv('interpolatedtracks-%s.dat'%gal,sep=r'\s{2,}', engine='python')
print( d, '\n', d.columns )
d = d[['mZAMS', 'time(02)', 'mass(03)','Lbol(07)', 'Rstar(08)', 'Teff(09)']]  # we don't need all columns
d.to_csv( 'boost_%s_int.csv'%gal,index=False )  # save as csv
print(d)
print( 'THAT took %s seconds.'%int( time.time()-start + 0.5) )
