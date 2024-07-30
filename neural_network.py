import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from tensorflow import keras
from tensorflow.keras import layers
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

# READ SMC TRAINING DATA
d = pd.read_csv( '/kaggle/input/vi-data/smc_training_data.csv' )
rsg_flag_ir_bool = d.rsg_flag_ir.astype( bool )

# PREPROCESS DATA
d = d.dropna( subset=['color','mabs'] )
qout1 = (d.rsg_flag_ir == 1) & (d.color<0)   # find outlier RSGs
qout2 = (d.rsg_flag_ir == 1) & (d.color>1.1) # find more outlier RSGs
qout  = qout1 | qout2
d = d[~qout]   # remove the outliers, 5 in this case
scaler = StandardScaler()
# add artificial non-RSG points to help training where data is sparse
xgrid, ygrid = np.linspace(-1,2,51), np.linspace(-9.5,0.5,51)
Xgrid, Ygrid = np.meshgrid( xgrid, ygrid )
art_data = np.vstack( [Xgrid.flatten(),Ygrid.flatten()] ).T
e = pd.DataFrame( art_data, columns=['color','mabs'])
e['rsg_flag_ir'] = 0
# Remove artificial points close to RSG branch
slope, cept = -13, 2
e = e[ abs( e.mabs - slope*e.color - cept) > 6 ]
de = pd.concat( [d,e], axis=0 )  # stack real stars and training stars vertically 

# DO A TRAIN - TEST SPLIT
X = de[['color','mabs']]
X = scaler.fit_transform(X)
y = de.rsg_flag_ir
X_train,X_val,y_train,y_val = train_test_split( X,y,train_size=0.7,random_state=1 )
lala = y_train > 0.5
X_train_rsg  = X_train[lala]
X_train_rest = X_train[~lala]

# TRAIN NEURAL NETWORK WITH TRAINING DATA
model = keras.Sequential( 
    [layers.Dense(units=16,input_shape=[2],activation='relu'),
    layers.Dense(units=16,activation='sigmoid'), # relu sigmoid same performance
    layers.Dense(units=1, activation='sigmoid') ])   # with sigmoid it works
model.compile( loss='binary_crossentropy',optimizer='adam')
history=model.fit( X_train,y_train,batch_size=256, epochs=120 )
print(model)
pred = model.predict( X_val )

# TEST THE PREDICTIONS
lala = pred > 0.5   # here the neural network predicts a high RSG probability
if lala.ndim != 1:
    lala = lala.reshape(-1)
X_val_rsg  = X_val[lala]
X_val_rest = X_val[~lala]
print( 'shapes rsg and rest:', X_val_rsg.shape[0],X_val_rest.shape[0] )
# plot to validate method
#plt.scatter( X_val_rsg[:,0], X_val_rsg[:,1], alpha=0.2, s=5, c='r' )
#plt.scatter( X_val_rest[:,0], X_val_rest[:,1], alpha=0.2, s=5, c='k' )

# FINAL STEP: USE TRAINED NEURAL NETWORK TO FIND RSGs IN OTHER GALAXIES
gal = 'n300'
di = pd.read_csv( '/kaggle/input/vi-data/%s.csv'%gal )
di = di.dropna( subset=['color','mag'] )
di = di[ di.mag < 24 ] # remove dim sources (high magnitude --> dim)
# standard scaling like we did for the training data
Xi = scaler.fit_transform(di)
predi = model.predict(Xi)
di = di.assign( p_rsg_nn=predi )
lili = predi > 0.5
if lili.ndim != 1:
    lili = lili.reshape(-1)
print( lili[0] )
di.to_csv( '/kaggle/working/%s_nnprob.csv'%gal )  # SAVE DATA WITH PREDICTED RSG PROBABILITY
di_rsg = di[ di.p_rsg_nn > 0.5 ]
di_rest= di[ di.p_rsg_nn < 0.5 ]
plt.scatter( di_rsg.color,  di_rsg.mag, alpha=0.2, s=5, c='r' )
plt.scatter( di_rest.color, di_rest.mag, alpha=0.2, s=5, c='k' )
plt.ylim(25,15)
plt.show()

