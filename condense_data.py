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

gal = 'MW'    #  MW, LMC, SMC, IZw18 

# with this part we convert to normal units
d = pd.read_csv('boost_%s_int.csv'%gal)
d = d.apply( pd.to_numeric,errors='coerce' )
d = d.dropna( subset=['mZAMS'] )
Msol = 1.9884777777777777e33  # this returns exactly 9 Msol for the supposedly 9Msol star 
d.mZAMS, d.mass = d.mZAMS/Msol, d.mass/Msol
d.age   = d.age/(365.25*86400)*1e-6  # from sec to Myr
d.Rstar = d.Rstar/69634000000.       # from cm to Rsol
d.logL  = np.log10(d.logL/3.846e33)  # from erg to logLsol
print(d,'\n',d.columns)

# with this part we make steps of 0.02 dex instead of 0.001 dex
u = d.mZAMS.unique()
print(u)
la = list(range(0,1856, 20)) # steps of 20
lala = u[la]                 # ZAMS masses in steps of 20
d = d[ d.mZAMS.isin(lala) ]  # we take only data of those ZAMS masses
print(d)
print( d.mZAMS.unique() )
d.to_csv( 'boost_%s_int_usable.csv'%gal, index=False )
