import pandas as pd
import pathlib as pl
import shutil
import act
import pandas as pd
import xarray as xr
import numpy as np
from datetime import timedelta

token='7915e9fe3c987b12'

start_date = '2022-03-04'
end_date = '2022-04-01'


base_dest = pl.Path("/global/homes/r/rudisill/gshare/sail_data_will/data_store_sail_period/squire/")
final_dest = base_dest.joinpath("daily_acc_files")

# output_file = act.discovery.download_data(
#     "willrudisill", token, "gucxprecipradarsquire", sd, ed, output=dest)


date_range = pd.date_range(start=start_date, end=end_date, freq='D')

# loop through dates and download data
for date in date_range:
    sd = date.strftime("%Y%m%d")
    ed = (date + timedelta(days=1)).strftime("%Y%m%d")

    # temp location for hte files
    month_dest = final_dest.joinpath(f"{date.year}{date.month:02d}")
    month_dest.mkdir(parents=True, exist_ok=True)

    # download the files
    output_file = act.discovery.download_data(
        "willrudisill", token, "gucxprecipradarsquire", sd, ed, output=month_dest)

    # now process them ...
    dsxx=xr.open_mfdataset(month_dest.glob("*nc"))

    # output file name
    ofilename = final_dest.joinpath("squire_%s_snow_rate_ws88diw.nc"%sd)

    # make the mean
    dsxx.snow_rate_ws88diw.mean(dim='time').to_netcdf(ofilename)

    # now delete the rest of the files that we dont need
    shutil.rmtree(month_dest)