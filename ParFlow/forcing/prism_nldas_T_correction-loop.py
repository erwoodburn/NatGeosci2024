
import numpy as np
import pandas as pd
import glob
import os
import sys
import pickle as pk
import pdb
import itertools
from scipy import stats


from parflow.tools.io import write_pfb, read_pfb
import parflow as pf
import parflow.tools as pftools

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, MultipleLocator, AutoMinorLocator
plt.rcParams['font.size']=14
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import lines



# --- Read in DEM ---

# --- Read in Parflow static fields ---
slopex   = pf.read_pfb('./ER_slope_x.pfb')[0,:,:] # shape [170,150]
slopey   = pf.read_pfb('./ER_slope_y.pfb')[0,:,:] # shape [170,150]
mannings = pf.read_pfb('./ER_mannings_100m.pfb')[0,:,:] # shape [170,150]
dx,dy    = 100.0, 100.0
nz       = 5

#dem = np.loadtxt('ER_dem.sa', skiprows=1)
#dem = dem.reshape(170,150)
dem = pf.read_pfb('ER_dem.pfb')[0,:,:]


# Setup mask for only East River domain
dom   = np.loadtxt('./subbasin_indicator.txt', skiprows=1)
mask  = np.loadtxt('./subbasin_indicator.txt', skiprows=1)  # 170x150
mask_ = np.zeros([nz, mask.shape[0], mask.shape[1]]) # needs to be 3D array
mask_[-1] = mask




# --- Utility Functions ---

def read_nldas(dirname, wy, time, var):
    f  = './baseline/{}/wy{}/NLDAS.{}.{:06}.pfb'.format(dirname, wy, var, time)
    dd = read_pfb(f)[0,:,:] # 2d surface array
    return dd

def read_prism(dirname, wy, time, var):
    f  = '{}/PRISM.{}.{:06}.pfb'.format(dirname, var, time)
    dd = read_pfb(f)[0,:,:] # 2d surface array
    return dd

def find_times(dirname, wy):
    fnames = glob.glob('./baseline/{}/wy{}/*Temp*'.format(dirname, wy))
    times  = []
    for f in fnames:
        times.append(int(f.split('.')[-2]))
    times = np.array(times)
    return np.sort(times)




def build_nldas_T_wy(dirname, wy, times=None):
    # Note, dictionary is about 6 GB per WY
    # times needs to be list or array of .pfb files to look for
    if times is None:
        times = find_times(dirname, wy)
    else:
        times = times
        
    dates_hourly = pd.date_range(start='{}-10-01'.format(wy-1), freq='1H', periods=len(times))
    
    # arrays to hold raw-hourly output
    T_  = np.zeros([len(times), dem.shape[0],dem.shape[1]])
    for i,t in enumerate(times):
        T_[i]    = read_nldas(dirname=dirname, wy=wy, time=t, var='Temp')-273.15   # NLDAS Temp. (C)
    
    return T_, dates_hourly



def build_prism_T_wy(dirname, wy, times, var):
    # Note, dictionary is about 6 GB per WY
    # times needs to be list or array of .pfb files to look for
    if times is None:
        times = find_times(dirname, wy)
    else:
        times = times
        
    dates_hourly = pd.date_range(start='{}-10-01'.format(wy-1), freq='1H', periods=len(times))

    # arrays to hold raw-hourly output
    T_  = np.zeros([len(times), dem.shape[0],dem.shape[1]])
    for i,t in enumerate(times):
        T_[i]    = read_prism(dirname=dirname, wy=None, time=t, var=var)   # NLDAS Temp. (C)
    
    return T_, dates_hourly




# -----------------------------------------------------------------------------------
# --- PRISM Mean Daily T Analysis ---
# --- Loop Through and Perform Correction Factor on all data X,Y points and Times ---
# -----------------------------------------------------------------------------------

#wy = 2015
wy = int(sys.argv[1])
print ('working on {}'.format(wy))


# --- NLDAS ---
dirname = 'nldas'
Tn, datesn = build_nldas_T_wy(dirname, wy)


    
# --- PRISM Mean Temperatures ---
dirname = '/pscratch/sd/p/pjdf/Projects/EastRiver/runsParflow/prism/wy2015-2021_bilinear'

prism_dates_all = pd.date_range('2014-10-01 01:00:00', '2021-10-01 00:00:00', freq='H')
wy_msk = np.where((prism_dates_all>'{}-10-01 00:00:00'.format(wy-1))&(prism_dates_all<'{}-10-01 01:00:00'.format(wy)), True, False)
prism_dates_wy = prism_dates_all[wy_msk]
prism_times    = np.arange(len(prism_dates_all))[wy_msk]
    
Tp, datesp = build_prism_T_wy(dirname, wy, times=prism_times, var='Temp')


# In[6]:


def discrete_sum(dates_arr, T_arr, y, x):
    if len(T_arr.shape) == 3:  # 3d array of T at [time,x,y]
        _df1   = pd.DataFrame(index=dates_arr, data=T_arr[:,y,x], columns=['T1'])
    else: # 1d array of T at [time]
        _df1   = pd.DataFrame(index=dates_arr, data=T_arr, columns=['T1'])
        
    _t1    = _df1.groupby(by=pd.Grouper(freq='1D')).mean()
    _df    = _t1.resample('1H').mean()
    _df.rename(columns={'T1':'T2'}, inplace=True)
    _df_avg = _df1.join(_df)['T2'].ffill()

    return _df_avg.to_numpy().ravel()

def ma(x, w):
    '''moving average on 1D timeseries x with window size w'''
    return np.convolve(x, np.ones(w), 'same') / w


# In[7]:


# update temperatures now
# mean PRISM only here
Tn_cor = np.ones_like((Tn))*-9999

for y in range(Tn_cor.shape[1]):
    for x in range(Tn_cor.shape[2]):
        _Tn_davg = discrete_sum(datesn, Tn, y, x)
        
        cf = Tp[:,y,x] - _Tn_davg
        cf = ma(cf,12)

        Tn_cor[:,y,x] = Tn[:,y,x]+cf


# In[8]:


# break up into water years and write pfb files
_dir = './NLDAS_PRISM_meantemp_pfb/wy{}'.format(wy)

if not os.path.exists(_dir):
    os.makedirs(_dir)
    
for j,l in enumerate(datesn):
    #print ('{}/{}'.format(j,len(datesn[_datesmsk])))
    pf.write_pfb('{}/NLDAS.Temp.{:06d}.pfb'.format(_dir,j), Tn_cor[j]+273.15, dist=True, p=8 , q=8, r=1)
    pass


# In[ ]:





# In[16]:


# --------------------------------------------------------------------
# --- PRISM Tmax and Tmin Daily T Analysis ---
# --------------------------------------------------------------------

# --- NLDAS ---
dirname = 'nldas'
Tn, datesn = build_nldas_T_wy(dirname, wy)

    

# --- PRISM Min/Max Temperatures ---
dirname = '/pscratch/sd/p/pjdf/Projects/EastRiver/runsParflow/prism/wy2015-2021_bilinear'


prism_dates_all = pd.date_range('2014-10-01 01:00:00', '2021-09-30 23:00:00', freq='H')

wy_msk = np.where((prism_dates_all>'{}-10-01 00:00:00'.format(wy-1))&(prism_dates_all<'{}-10-01 01:00:00'.format(wy)), True, False)
prism_dates_wy = prism_dates_all[wy_msk]
prism_times    = np.arange(len(prism_dates_all))[wy_msk]
    
Tpmax, datespmax = build_prism_T_wy(dirname, wy, times=prism_times, var='Tmax')
Tpmin, datespmin = build_prism_T_wy(dirname, wy, times=prism_times, var='Tmin')


# In[17]:


def min_max_norm(y, x, Tn, datesn, Tpmin, Tpmax, datespmin):
    # NLDAS Min and Max Temperatures before corrections
    _T = Tn[:,y,x].reshape(int(len(Tn[:,y,x])/24), 24)
    _d = datesn.to_numpy().reshape(int(len(datesn)/24), 24)

    _dmin_ = _d.ravel()[_T.argmin(axis=1) + np.arange(len(_T))*24] # NLDAS dates for tmin
    _dmax_ = _d.ravel()[_T.argmax(axis=1) + np.arange(len(_T))*24] # NLDAS dates for tmax
    #_Tmin_ = _T.ravel()[_T.argmin(axis=1) + np.arange(len(_T))*24] # NLDAS tmin
    #_Tmax_ = _T.ravel()[_T.argmax(axis=1) + np.arange(len(_T))*24] # NLDAS tmax

    _Tn = Tn[:,y,x]
    _Tn = _Tn.reshape(int(len(_Tn)/24), 24)

    _Tmin = Tpmin[np.intersect1d(datespmin, _dmin_, return_indices=True)[1], y, x]
    _Tmax = Tpmax[np.intersect1d(datespmax, _dmax_, return_indices=True)[1], y, x]

    T_std = (_Tn - _Tn.min(axis=1)[:,np.newaxis]) / (_Tn.max(axis=1)-_Tn.min(axis=1))[:,np.newaxis]
    T_scl = T_std * (_Tmax[:,np.newaxis]-_Tmin[:,np.newaxis]) + _Tmin[:,np.newaxis]
    T_n_pminmax = T_scl.ravel()
    #dates_n_pminmax = datesn.copy()
    
    return T_n_pminmax #, dates_n_pminmax


# In[20]:


# update temperatures now
# min/max PRISM analysis
Tn_cor = np.ones_like((Tn))*-9999

for y in range(Tn_cor.shape[1]):
    for x in range(Tn_cor.shape[2]):
        Tn_cor[:,y,x] = min_max_norm(y, x, Tn, datesn, Tpmin, Tpmax, datespmin)


# In[23]:


_dir = './NLDAS_PRISM_minmaxtemp_pfb/wy{}'.format(wy)

if not os.path.exists(_dir):
    os.makedirs(_dir)
    
for j,l in enumerate(datesn):
    #print ('{}/{}'.format(j,len(datesn)))
    pf.write_pfb('{}/NLDAS.Temp.{:06d}.pfb'.format(_dir,j), Tn_cor[j]+273.15, dist=True, p=8 , q=8, r=1)
    #pass


# In[ ]:




