import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import poisson
from matplotlib import rcParams
rcParams['font.family'] = 'STIXGeneral'
c1,c2,c3,c4,c5,c6,c7,c8,c9,c10 = '#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf'
import time
start = time.time()

weight_factor = 200                                   # lower value --> longer simulation in while loop. 10 = ~100s
save_fig = False                                      # save figure? 
save_as  = 'figs/model_vs_obs_hists_zwzw.png'  

nbins,xmin,xmax = 15, 4, 7                            # values used for binning
fig,((ax1a,ax1b),(ax2a,ax2b),(ax3a,ax3b))=plt.subplots(3,2,figsize=(7,5))
axes = [ax1a,ax1b,ax2a,ax2b,ax3a,ax3b]

def pro_binning(var_list, nbins,xmin,xmax,weights):   # binning function
        hist, bins = np.histogram(var_list, bins=nbins, weights=weights, range=(xmin,xmax))
        width = 1. * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        return hist, center, width

gals = { 'n4395':{'Z':1,    'Zrep':'MW'},             # make a dictionary of GALAXIES to get relevant metallicities Z
         'LMC':  {'Z':0.5,  'Zrep':'LMC' },
         'n300': {'Z':0.4,  'Zrep':'LMC'  },
         'n5253':{'Z':0.33, 'Zrep':'SMC'},
         'SMC':  {'Z':0.2,  'Zrep':'SMC'},
         'izw18':{'Z':0.02, 'Zrep':'IZw18'} }

for i, gal in enumerate( gals ):                      # loop over the 6 galaxies we consider
    ax = axes[i]
    z, zrep = gals[gal]['Z'], gals[gal]['Zrep']       # get metallicity and closest model metallicity for each galaxy
    if gal in ['n5253','n4395','n300']:
        logl_rsg = pd.read_csv( 'data/pergal/%s_vis_logls.csv'%gal, header=None )   # calculated previously with lrsg_meth/logl_from_cmd.py
    if gal == 'izw18':
        logl_rsg = pd.read_csv( 'data/pergal/%s_ir_logls.csv'%gal, header=None )
    elif gal in ['SMC','LMC']:
        d = pd.read_csv('/home/aschoot/Desktop/etteren_projects/tidy_lrsg_vs_z/data/pergal/davies18_%s.csv'%gal)
        logl_rsg = d.logL

    # OBSERVATIONS
    ohist,ocenter,owidth = pro_binning(logl_rsg,nbins,xmin,xmax,None)
    ax.step( ocenter,np.log10(ohist+1e-20),where='mid',label='Obs.')  # show observed logL distribution
    nobs_anchor = sum(ohist[ (ocenter>4.99) & (ocenter<5.41) ])       # count observed N_star in range 5.0 < logL < 5.4
    print( '\n%s nobs_anchor: %s'%(gal,nobs_anchor))

    # MODELS
    b     = pd.read_csv( '/home/aschoot/Desktop/etteren_projects/tidy_lrsg_vs_z/data/boost_%s_int_usable.csv'%zrep )
    if i==0: print(b.columns)
    b     = b[ (b.Teff < 7e3) & (b.age < 99) & (b.Teff>1) ]   # roughly kill Hburn, and also BHs
    minis = b.mZAMS.unique()
    minis = minis[ minis<200 ]                                # used below to select models by their mass at ZAMS
    props = minis**-1.35                                      # weight by imf   -1.35 instead of -2.35 because of equal spacing in LOG space
    props = props / sum(props)                                # normalize
    nmodel_anchor, model_logLs, model_weights = 0,[],[]
    while nmodel_anchor < nobs_anchor:                        # randomly draw models 
        randmini = np.random.choice( minis,p=props )
        bi = b[ b.mZAMS == randmini ]  
        bi = bi[ bi.age > 0.92*max(bi.age) ]
        randj=np.random.randint( 1,len(bi) )                  # Draw integer for random model. Start at 1 to avoid trouble when we draw 0
        randlogL, randdt = bi.logL.iloc[randj], bi.age.iloc[randj]-bi.age.iloc[randj-1]
        if randdt < 0: continue                               # because of a mistake in the interpolated BOOST models (Szecsi+2022)
        #if bi.index[randj]-bi.index[randj-1] > 1.5: randdt = 0   # if the previous model is many models ago, dt is troubled
        if bi.index[randj]-bi.index[randj-1] > 1.5: continue  # test: works it?
        model_logLs.append(   randlogL )
        model_weights.append( randdt*weight_factor)           # give models a weight based on time step and a constant regulating sim. time
        if 5 < randlogL < 5.4:
            nmodel_anchor += randdt*weight_factor             # use N_star in range 5.0<logL<5.4 as an anchor to match to N_observed
    mhist,mcenter,mwidth = pro_binning(model_logLs,nbins,xmin,xmax,model_weights)
    ax.step( mcenter,np.log10(mhist+1e-20),where='mid',label='Models', ls='--', lw=3 )
    model_logLs, model_weights = np.array(model_logLs), np.array(model_weights)
    q = np.array(model_logLs) > 5.6                           # investigate number of bright theory models
    occ_model = np.round( sum(model_weights[q]),2 )
    occ_obs = (logl_rsg.values > 5.6).sum()
    print( '%s integrated model N_RSG brighter than logL = 5.6: %s. Observed: %s'%(gal,occ_model, occ_obs) )
    cum_pois_prob = poisson.cdf( occ_obs, occ_model )         # calculate cumulative Poisson probability
    formatted_cum = "{:.1e}".format( cum_pois_prob )          # transform that into the format in which we want to write it
    print( 'cumulative poisson possibility to have %s or less occurrences: %s'%(occ_obs,formatted_cum) )

    # COSMETICS
    ax.set_ylim(-0.4,3.4)
    ax.set_xlim(4.5,7.25)
    ax.set_ylabel('$\log N_\mathrm{RSG}$')
    ax.text(0.05,0.925,'%s, %s$\,$Z$_\odot$'%(gal,z),va='top', ha='left', transform=ax.transAxes,size=11)
    suffix = ''
    if gal == 'IZw18':                                        # then we write which models we show
        suffix = '%s models'%zrep
    ax.text(0.925,0.925,'at logL>5.6:\nObs: %s\nModel: %s\n$p\\thinspace_\mathrm{Pois.}$=%s\n%s'%(occ_obs,occ_model,formatted_cum,suffix),
            va='top', ha='right', transform=ax.transAxes,size=9)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params( axis='both',direction='in')
    if i>3.5:
        ax.set_xlabel('$\log ( L / L_\odot )$')
ax1a.legend( loc='upper center', framealpha=1, prop={'size':9} )
plt.tight_layout()
print( 'THAT took %s seconds.'%int( time.time()-start + 0.5) )
if save_fig == True:
    print('SAVED.   %s'%save_as)
    plt.savefig(save_as)
plt.show()
