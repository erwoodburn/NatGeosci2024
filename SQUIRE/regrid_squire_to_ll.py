import xarray as xr
import numpy as np
import xesmf as xe
import pathlib as pl 
import sys 


aso_file = pl.Path(sys.argv[1])

aso = xr.open_dataset(aso_file)

# open up the target file 
ds = xr.open_dataset("SQUIRE_PERIOD.NLDAS_snowfall_sum_map.h.txt_._regridll.nc")

# regrid the to the regular lat lon grid in the parflow file 
regridder = xe.Regridder(aso, ds, 'bilinear', ignore_degenerate=True)

# the output name 
the_name = aso_file.name.split(".nc")[0] + "regrid_reg.nc"

aso_regrid = regridder(aso)
aso_regrid.to_netcdf(the_name)
