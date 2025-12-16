import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import cdflib
import pandas as pd
import datetime as dt


import xarray as xr
import fsspec


from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)


a = 2 + 4
print(a)