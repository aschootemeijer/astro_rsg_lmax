import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.optimize import curve_fit
from matplotlib import rcParams
rcParams['font.family'] = 'STIXGeneral'
c1,c2,c3,c4,c5,c6,c7,c8,c9,c10 = '#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf'
fig,ax=plt.subplots(1,1,figsize=(6,4.75), dpi=100)
import time
start = time.time()

simspergal = 10                            # it will take almost 10 seconds per simulation nr entered here
nbrightest = 3                             # do we consider the brighest star, or e.g. 3rd brightest?
save_fig = False 
save_as  = 'figs/z_vs_brrsg_violin_rsgnr%s.png'%nbrightest

# make galaxy dictionary to use later on
gals = { 'n4395':{'Z':1,     'write':'NGC$\,$4395',  'Zrep':'MW' },   
         'n300': {'Z':0.5,   'write':'NGC$\,$300',   'Zrep':'LMC'},
         'LMC':  {'Z':0.4,   'write':'LMC',          'Zrep':'LMC'},      
         'n5253':{'Z':0.33,  'write':'NGC$\,$5253',  'Zrep':'SMC'},
         'SMC':  {'Z':0.2,   'write':'SMC',          'Zrep':'SMC'},
         'izw18':{'Z':0.025, 'write':'I$\,$Zw$\,$18','Zrep':'IZw18'}, 
         'm31':  {'Z':1.5,   'write':'M$\,$31',      'Zrep':'MW'} }

ss, lw = 200, 2.5                          # for scatter plot

# LOOP OVER 6 GALAXIES
for i,gal in enumerate(gals):           
    print('\n%s'%gal)
    zrep, z, write_name, mark= gals[gal]['Zrep'], gals[gal]['Z'], gals[gal]['write'], 'x'
    ax.text(np.log10(z),7.1, '- %s'%write_name, rotation=90, ha='center', va='top', size=10)
    
    # OBSERVATIONS
    if gal in ['n4395','n300','n5253', 'izw18']:            # diff galaxies have different data files that read
        logl_obs = pd.read_csv( 'data/pergal/%s_vis_logls.csv'%gal, header=None)
        ax.scatter( np.log10( gals[gal]['Z']), logl_obs.iloc[ nbrightest-1 ], c=c2, s=ss, marker='3', lw=lw )
        if gal == 'izw18':                                  # so for izw18 we overwrite logls from vis by ir logls
            logl_obs = pd.read_csv( 'data/pergal/%s_ir_logls.csv'%gal, header=None)
            ax.scatter( np.log10( gals[gal]['Z']), logl_obs.iloc[ nbrightest-1 ], c=c1, s=ss, marker='4', lw=lw )
    if gal in ['LMC','SMC']:
        logl_obs = pd.read_csv( 'data/pergal/davies18_%s.csv'%gal, usecols=['logL'] )
        ax.scatter( np.log10( gals[gal]['Z']), logl_obs.iloc[ nbrightest-1 ], c=c3, s=ss, marker='2', lw=lw ) 
    if gal == 'm31':
        logl_obs = [5.75, 5.53, 5.49, 5.46, 5.44]           # from McDonald+2022
        ax.scatter( np.log10( gals[gal]['Z']), logl_obs[ nbrightest-1 ], c=c4, s=ss, marker='1', lw=lw )
        continue                                            # we have no theory models of M31 hence this move
    logl_obs = np.array( logl_obs )
    print( logl_obs[:5] )
    nobs_anchor = len( logl_obs[ (logl_obs>5.0)&(logl_obs<5.4) ] )
    print( 'nr %s brighest:'%nbrightest,logl_obs[ nbrightest-1 ] )
    print( 'nobs_anchor',nobs_anchor )

    # THEORY POPULATIONS
    b     = pd.read_csv( '/home/aschoot/Desktop/etteren_projects/tidy_lrsg_vs_z/data/boost_%s_int_usable.csv'%zrep )
    if i==5: print(b.columns)
    b     = b[ (b.age < 99) & (b.Teff>1) ]                  # remove dead stars (black holes, neutron stars)

    minis = b.mZAMS.unique()                                # get initial masses and their statistical weights
    minis = minis[ minis<200 ]                              # adopt upper mass limit of 200M
    props = minis**-1.35                                    # weight by IMF by -1.35 instead of -2.35 because of equal spacing in LOG space
    props = props / sum(props)                              # normalize
    drawn_max_logls = []
    #continue                                               # if we want to show just observations. If also theory, uncomment
    for j in range( simspergal ):                     
        nmodel_anchor, drawn_logls, ndrawn = 0, [], 0

        while nmodel_anchor < nobs_anchor:                  # start drawing models with random age and initial mass
            randmini = np.random.choice( minis,p=props )
            randage  = np.random.random()*22                # random age in Myr to take lifetime effect into account
            if randage < 11*(randmini/12)**-1 or randage > 20*(randmini/12)**-0.7: continue # 20241213 to speed by factor 3. checked in tests/
            bi = b[ b.mZAMS == randmini ]
            bi = bi[ bi.age > 0.91*max(bi.age) ]            # select only models with the drawn initial mass      
            if randage < bi.age.iloc[0] or randage > bi.age.iloc[-1]: continue # do not consider ages outside of He-burning age of model
            randint  = (bi.age - randage).abs().idxmin()    # find closest age in model sequence
            drawn_teff, drawn_logl =  bi.Teff.loc[randint], bi.logL.loc[randint] # take Teff and logL of that model

            ndrawn += 1
            if drawn_teff > 7e3: continue
            drawn_logls.append(drawn_logl)
            if 5 < drawn_logl < 5.4:
                nmodel_anchor += 1                          # count how many theory model in range 5.0 < logL < 5.4

        # now we are done drawing and making theoretical population: evaluate theoretical population.
        drawn_logls  = np.array( drawn_logls )
        bright_logls = drawn_logls[ drawn_logls >5.6 ]
        print (  j, max(drawn_logls), ndrawn, 'drawn. Of these %s have logL > 5.6'%len(bright_logls)  )
        drawn_logls = sorted( drawn_logls, reverse=True )
        logl_of_interest = drawn_logls[ nbrightest ]
        drawn_max_logls.append( logl_of_interest )
        ax.scatter( np.log10(z), logl_of_interest, c='k', marker='_', alpha=0.25, zorder=12 )

    # show violin plots to show density of results from simulations
    parts = plt.violinplot( drawn_max_logls, positions=[np.log10(z)], widths=0.25, showextrema=False )
    for pc in parts['bodies']:
        pc.set_facecolor('k')
        pc.set_alpha(0.15)

# COSMETICS
# These 6 lines are to have extra control over the legend
ax.scatter( -100, -100, c='k', marker='_', alpha=0.35, zorder=12, label='Simulation' )
ax.scatter( -100, -100, c=c2, s=ss, marker='3', lw=lw, label='HST' )
ax.scatter( -100, -100, c=c1, s=ss, marker='4', lw=lw, label='JWST' )
ax.scatter( -100, -100, c=c3, s=ss, marker='2', lw=lw, label='Davies+2018' ) 
ax.scatter( -100, -100, c=c4, s=ss, marker='1', lw=lw, label='McDonald+2022' )
ax.legend(  loc='center right', framealpha=0.8, prop={'size':10.5}, bbox_to_anchor=(0.98,0.65))
labs = 13
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.tick_params( axis='both',direction='in')
ax.set_xlabel('$\log( \, Z / Z_\odot \, )$', size=labs)
if nbrightest == 1: ax.set_ylabel('$\log( \, L / L_\odot )$ of brightest RSG', size=labs)
if nbrightest == 3: ax.set_ylabel('$\log( \, L / L_\odot )$ of 3$\mathrm{rd}$ brightest RSG', size=labs)
ax.set_xlim(-1.8,0.6)
ax.set_ylim(5,7.15)
plt.tight_layout()
print( 'THAT took %s seconds.'%int( time.time()-start + 0.5) )
if save_fig == True:
    print('SAVED.   %s'%save_as)
    plt.savefig(save_as)
plt.show()
