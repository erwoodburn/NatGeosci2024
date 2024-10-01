#bl_wy2016.out.clm_output.04464.sa_PRISM_hetero_regridll.nc

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
import glob
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import sys 



target = sys.argv[1]

pf_2016   = xr.open_dataset("/pscratch/sd/w/woodburn/er_deltat_fall22/aso_comparison/bl_wy2016.out.clm_output.04464.sa_%s_regridll.nc"%target)
pf_2018a  = xr.open_dataset("/pscratch/sd/w/woodburn/er_deltat_fall22/aso_comparison/bl_wy2018.out.clm_output.04344.sa_%s_regridll.nc"%target)
pf_2018b  = xr.open_dataset("/pscratch/sd/w/woodburn/er_deltat_fall22/aso_comparison/bl_wy2018.out.clm_output.05640.sa_%s_regridll.nc"%target)
pf_2019a  = xr.open_dataset("/pscratch/sd/w/woodburn/er_deltat_fall22/aso_comparison/bl_wy2019.out.clm_output.04512.sa_%s_regridll.nc"%target)
pf_2019b  = xr.open_dataset("/pscratch/sd/w/woodburn/er_deltat_fall22/aso_comparison/bl_wy2019.out.clm_output.06048.sa_%s_regridll.nc"%target)

aso_2016  = xr.open_dataset("/pscratch/sd/w/woodburn/er_deltat_fall22/aso_comparison/aso_files/ASO_50M_SWE_USCOCB_20160404_llregrid_reg.nc")
aso_2018a = xr.open_dataset("/pscratch/sd/w/woodburn/er_deltat_fall22/aso_comparison/aso_files/ASO_50M_SWE_USCOCB_20180330_llregrid_reg.nc")
aso_2018b = xr.open_dataset("/pscratch/sd/w/woodburn/er_deltat_fall22/aso_comparison/aso_files/ASO_50M_SWE_USCOGE_20180524_llregrid_reg.nc")
aso_2019a = xr.open_dataset("/pscratch/sd/w/woodburn/er_deltat_fall22/aso_comparison/aso_files/ASO_50M_SWE_USCOGE_20190407_llregrid_reg.nc")
aso_2019b = xr.open_dataset("/pscratch/sd/w/woodburn/er_deltat_fall22/aso_comparison/aso_files/ASO_50M_SWE_USCOGE_20190610_llregrid_reg.nc")


dem = xr.open_dataset("/pscratch/sd/w/woodburn/er_deltat_fall22/aso_comparison/aso_files/upper_erw_dem_ll_0005regrid_reg.nc")



cmap_swe = plt.cm.viridis
cmap_swe.set_over(color='orange')
cmap_swe.set_under(color='lightblue')

cmap_swe_diff = plt.cm.coolwarm
cmap_swe_diff.set_over(color='red')
cmap_swe_diff.set_under(color='purple')



fig, ax = plt.subplots(5,3, figsize=(12,18))

vmax = 1400
vmin = 50 

vmaxdiff = 900
vmindiff = -900

def make_a_plot(aso_dat, pfdat, dem, axg, row):
    aso_dat = (aso_dat.Band1 * 1000).where(dem.Band1)
    pf_dat  =  pfdat.swe.where(dem.Band1)
    diff    =  pf_dat - aso_dat

    cb0=    aso_dat.plot(vmax=vmax, vmin=vmin, ax=axg[row,0], cmap=cmap_swe, add_colorbar=False)
    cb1=    pf_dat.plot(vmax=vmax, vmin=vmin,  ax=axg[row,1], cmap=cmap_swe, add_colorbar=False)
    cb2=    diff.plot(vmax=vmaxdiff, vmin=vmindiff,    ax=axg[row,2], cmap=cmap_swe_diff, add_colorbar=False)

    cb0=fig.colorbar(cb0, extend='both')
    cb1=fig.colorbar(cb1, extend='both')
    cb2=fig.colorbar(cb2, extend='both')
    cb2.set_ticks([-900, -600, -300, 0, 300, 600, 900])

    # set the labels
    cb0.set_label('SWE (mm)')
    cb1.set_label('SWE (mm)')
    cb2.set_label('SWE (mm)')


make_a_plot(aso_2016,  pf_2016,  dem, ax, 0)
make_a_plot(aso_2018a, pf_2018a, dem, ax, 1)
make_a_plot(aso_2018b, pf_2018b, dem, ax, 2)
make_a_plot(aso_2019a, pf_2019a, dem, ax, 3)
make_a_plot(aso_2019b, pf_2019b, dem, ax, 4)



ax[0,0].set_title("ASO")  
ax[0, 1].set_title("PF-CLM")
ax[0, 2].set_title("PF-CLM  - ASO")

for i in range(4):
    ax[i,0].set_xlabel("")
    ax[i,1].set_xlabel("")
    ax[i,2].set_xlabel("")

    ax[i,0].set_xticklabels([])
    ax[i,1].set_xticklabels([])
    ax[i,2].set_xticklabels([])

ax[1,0].set_xticklabels(ax[1,0].get_xticklabels(), rotation=45)
ax[1,1].set_xticklabels(ax[1,1].get_xticklabels(), rotation=45)


for z,lab in zip(range(5), ["2016", "2018a", "2018b", "2019a", "2019b"]):
    ax[z,0].set_ylabel('lat')
    ax[z,0].text(-.5, .5, lab, va='center', rotation='vertical', transform=ax[z,0].transAxes)
    ax[z,1].set_ylabel('')
    ax[z,2].set_ylabel('')

for yx in range(5):
    ax[yx,1].set_yticklabels([])
    ax[yx,2].set_yticklabels([])

plt.savefig("%s_map_plots.png"%target, dpi=300)




### now make the boxplots...abs
pflist = [pf_2016,
          pf_2018a,
          pf_2018b, 
          pf_2019a, 
          pf_2019b] 


asolist = [aso_2016,
           aso_2018a,
           aso_2018b,
           aso_2019a,
           aso_2019b]

eblist = [2000, 2500, 3000, 3500]


total_swe_aso   = []
total_swe_pf    = []
total_bias_list = []
total_pbias_list = []
total_rmse_list = []

k = 0 
fig,ax = plt.subplots()
for pf, aso in zip(pflist, asolist):
    pfmnx    = pf.swe.where(dem.Band1 > 0)
    pfmn    = pfmnx.values.flatten()
    pfmn    = pfmn[~np.isnan(pfmn)]
    
    asomnx    = (aso.Band1*1000).where(dem.Band1 > 0)
    asomn    = asomnx.values.flatten()
    asomn    = asomn[~np.isnan(asomn)]


    bp1=ax.boxplot(pfmn,  positions=[k-.1], showfliers=False, patch_artist=True)
    for patch in bp1['boxes']:
        patch.set_facecolor('red')


    bp2=ax.boxplot(asomn, positions=[k+.1], showfliers=False, patch_artist=True)
    for patch in bp2['boxes']:
        patch.set_facecolor('green')

    k+=1

    # calculate total swe bias 
    bias = np.mean(pfmn) - np.mean(asomn)
    total_bias_list.append(bias)
    
    # calculate total swe bias %
    pbias = bias/np.mean(asomn) * 100
    total_pbias_list.append(pbias)
    
    # calculate RMSE
    rmse = np.sqrt(np.mean((pfmnx - asomnx)**2))
    total_rmse_list.append(rmse.values)


    total_swe_aso.append(np.mean(asomn))
    total_swe_pf.append(np.mean(pfmn))

ax.legend([bp1["boxes"][0], bp2["boxes"][0]], ['PF-CLM', 'ASO'], loc='upper left')


ax.set_xticks(range(5))
ax.set_xticklabels(["2016", "2018a", "2018b", "2019a", "2019b"])

ax.set_title("SWE (mm) in East River Watershed -- Upper ERW Only")
ax.set_ylabel("SWE (mm)")

plt.savefig("%s_boxplots.png"%target, dpi=300)



## now make the table
df = pd.DataFrame(data=np.array([total_swe_pf, total_swe_aso, total_bias_list, total_pbias_list, total_rmse_list])).T
df.columns = ["PF-CLM-mn-swe", "ASO-mn-swe", "Bias", "Biasp", "RMSE"]
df.index = ["2016", "2018a", "2018b", "2019a", "2019b"]
df.to_csv("%s_swe_compare_table.csv"%target)

