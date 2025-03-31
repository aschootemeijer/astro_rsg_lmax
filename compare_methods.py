import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rcParams
rcParams['font.family'] = 'STIXGeneral'
import time
start = time.time()

save_fig  = False 
save_as   = '../figs/meth_test_loglonly.png'
visorirs = ['Vis','IR']                       # we use data in the visible part of the spectrum as well as IR data
gal      = 'SMC'                              # Galaxy: Small Magallanic Cloud

## FUNCTIONS TO LOAD DATA AND TO GET LOGL
path = '/home/aschoot/Desktop/etteren_projects'
def get_gal_props( gal,visorir ):
    if gal == 'SMC':
        #d = pd.read_csv('%s/tidy_lrsg_vs_z/lrsg_meth/data/d18smc_x_skymapper1p1.csv'%path) # both for IR and vis
        d = pd.read_csv('%s/tidy_lrsg_vs_z/lrsg_meth/data/yang19_x_d18smc.csv'%path) # both for IR and vis
        print(d.columns)
        DM, AV = 18.977, 0.35   # Graczyk+2020, Schootemeijer+2021
        if visorir == 'IR':
            d = d[ d.Jmag-d.Kmag > 0.6 ]
            cb_poly = [5.00513314,-18.39883674,25.83141869,-18.2448306,8.62224833,-0.19544753]    # 20250317 jwst/boco    #for 2mass
            mags, colors = d['Kmag'], d['Jmag']-d['Kmag']
        if visorir == 'Vis':
            cb_poly = [0.12977808, -1.07136445, 3.16558804, -3.92916783, 1.26515746, 0.1510785 ]
            cb_fit  = np.poly1d(cb_poly)  
            mags, colors = d['imagSkyM'], d['rmagSkyM']-d['imagSkyM']
    cb_fit  = np.poly1d(cb_poly)
    return mags, colors, d, cb_fit, DM, AV

def get_logls( mags, colors, DM, AV,cb_fit): # function to calculate luminosities    
    Ablue, Ared = Ablue_div_AV * AV, Ared_div_AV * AV
    mabs0       = mags - Ared- DM
    colors0     = colors - (Ablue-Ared)
    logls= -0.4*( mabs0-4.74 + cb_fit(colors0) )
    return logls 

def get_mag0_color0( mags, colors, AV ):     # to get intrinsic magnitude, colors (i.e., correct for extinction)
    Ablue, Ared = Ablue_div_AV * AV, Ared_div_AV * AV
    mags0       = mags - Ared
    colors0     = colors - (Ablue-Ared)
    return mags0, colors0

fig,(ax1b,ax1a)=plt.subplots(1,2,figsize=(7,3.5) )# , dpi=250)
axes = [ax1a,ax1b]
xfit = np.linspace(0, 3, 101)
lss, lws = ['solid','--','-.',':'], [2.5,2,1.5,1]

# OK HERE WE REALLY START
for i,visorir in enumerate(visorirs):
    # OBTAIN MAGNITUDES AND COLORS
    ax = axes[i]
    Ablue_div_AV, Ared_div_AV = 1, 0.567                             # Gordon+2003    V I band, SMC. in Ax/AV
    if visorir == 'IR':  Ablue_div_AV, Ared_div_AV = 0.243, 0.078    # Wang,Chen+2019 J,K band 
    if visorir == 'Vis': Ablue_div_AV, Ared_div_AV = 0.843, 0.628    # Wang,Chen+2019 r,i band PANSTARRS 
    mags, colors, d, cb_fit, DM, AV = get_gal_props( gal,visorir )
    print('GAL:', gal, visorir)
    rsg_mags0, rsg_colors0 = get_mag0_color0( mags,colors,AV )       # correct for extinction

    # CALCULATE THE RSG LUMINOSITIES
    logls = get_logls( rsg_mags0, rsg_colors0, DM, 0, cb_fit )       # use AV=0 since extinction is already taken into account
    loglsorted = sorted(logls, reverse=True)
    print('sorted logLs', loglsorted[:10])
    print('\n\n\n')
    
    # MAKE THE PLOTS
    ax.scatter( d.logL, logls, s=15, lw=0, c='#ff0000')#'#ddaa33' )
    ax.plot([4.5,5.6],[4.5,5.6], ls='--',c='#333333', lw=1.8)
    diff = d.logL - logls         # calculate the difference between our method and Davies+2018
    diff = diff[ abs(diff) < 1 ]  # to remove 104 dex outlier
    std = np.round( np.std( diff ),2)
    if visorir == 'Vis': visorir='VIS'  # for text in plot
    ax.text(0.05,0.95,'%s, %s\nstdev($\Delta$logL) = %s'%(gal,visorir,std),transform=ax.transAxes,va='top',ha='left',size=12)
    ax.set_xlim(4.5,5.7)
    ax.set_ylim(4.5,5.8)
    ax.set_xlabel('$\log( L / L_\odot$) [Davies+2018]', size=12)
    ax.set_ylabel('$\log( L / L_\odot$) [This work]',   size=12)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params( axis='both',direction='in')
        
plt.tight_layout()
if save_fig == True:
    print('SAVED.   %s'%save_as)
    plt.savefig(save_as)
plt.show()
