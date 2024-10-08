import xarray as xr
import numpy as np
import xesmf as xe
import pandas as pd 
import pathlib as pl
import sys

file_path = pl.Path(sys.argv[1])
parent_folder = str(file_path.parent).split("/")[-1]
print(parent_folder)


#file_path="/pscratch/sd/w/woodburn/er_deltat_fall22/aso_comparison/NLDAS/bl_wy2016.out.clm_output.04464.sa"

# Assuming the file is named 'numbers.txt' and is in the current directory
#file_path = 'numbers.txt'
data = np.loadtxt(file_path, skiprows=1)


#array_list = [] 
#
#for z in range(1,18):
#    array = np.zeros(150 * 170)
#    for y in range(0,170):
#        for x in range(0,150):
#            xxx = z*y*x
#            array[y*x] = data[xxx]  
#    array_list.append(array)


s = 170*150
i= 0
swe = data[i*s:(i+1)*s].reshape(170,150)


latlons = pd.read_csv("/pscratch/sd/w/woodburn/er_deltat_fall22/aso_comparison/modelCellCenters_latLon_csv.csv")
the_lats = np.array(latlons["Lat"])
the_lons = np.array(latlons["Lon"])

latgrid = the_lats[::-1].reshape(170,150)
longrid = the_lons.reshape(170,150)


ds = xr.Dataset(data_vars={"snowfall": (("x","y"), swe),
                           "lat": (("x","y"), latgrid),
                           "lon": (("x","y"), longrid)}) 


ds = ds.assign_coords({"lat": (("x","y"), latgrid),
                  "lon": (("x","y"), longrid)})
# now let's regrid the data to a regular grid...


lats = np.arange(ds.lat.min().values, ds.lat.max(), .001)
lons = np.arange(ds.lon.min().values, ds.lon.max(), .001)


target = xr.DataArray(np.zeros((len(lats), len(lons))), coords={'lat': lats, 'lon': lons})

regridder = xe.Regridder(ds, target, 'bilinear', ignore_degenerate=True)

ds_regrid = regridder(ds)

ds_regrid.to_netcdf(file_path.name+"_"+parent_folder+"_regridll.nc")
