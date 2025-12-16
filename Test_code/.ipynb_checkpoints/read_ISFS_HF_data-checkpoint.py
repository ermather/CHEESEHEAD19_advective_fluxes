# from setup import * #Import setup module
#Import packages
import fsspec
import xarray as xr
import pandas as pd
import numpy as np
import datetime as dt

#Filepath for saving data
filepath = '.\hh_means.csv'

#Set up start and end dates for reading in data
#Format: df.datetime(year, month, date, hour)
#Start and end dates for the full CHEESEHEAD dataset 
# startdate = dt.datetime(2019, 6, 20, 0)
# enddate = dt.datetime(2019, 10, 14, 23, 30)

#Start and end timestamps for a short test period 
startdate = dt.datetime(2019, 8, 1, 0) 
enddate = dt.datetime(2019, 8, 1, 3)



hours_f = np.arange(0,24) #Array of hours as floats
hours_s = [str(x).zfill(2) for x in hours_f] #Make list of hours as strings to use in calling files

#Set up dataframe to hold desired variables
#Array of vertical wind variable names to read in
#WLEF data is in a different directory
w_vars = ['w_30m_nw1', 
          'w_10m_nw2', 
          'w_2m_nw3', 
          'w_30m_nw4',
          'w_30m_ne1', 
          'w_30m_ne2', 
          'w_30m_ne3', 
          'w_30m_ne4', 
          'w_30m_sw1',
          'w_25m_sw2', 
          'w_30m_sw3', 
          'w_30m_sw4', 
          'w_30m_se2', 
          'w_30m_se3',
          'w_3m_se4', 
          'w_10m_se5', 
          'w_30m_se6']

dtindex =  pd.date_range(start=startdate, end=enddate, freq='30min')
df_w = pd.DataFrame(index = dtindex, columns = w_vars)


#Series of string timestamps to use 
ts_s = pd.Series(data = df_w.index.strftime('%m%d_%H'), index = dtindex)

#This creates an array to hold the high-frequency data. Since I just used 30-minute averages, I never actually used this, the code to add data to the array is commented below (at the bottom of the for loop). Would need to be modified to handle multiple variables at once.
flat_dat = np.empty(0)

#Open files and pull out desired data
for hour in dtindex:
    print(hour)

    url = 'http://co2.aos.wisc.edu/data/CHEESEHEAD-incoming/ISFS/high-rate/isfs_geo_hr_2019' +ts_s.loc[hour] + '.nc'
    with fsspec.open(url) as fobj:
        ds = xr.open_dataset(fobj) #This is an xarray dataset with 20 hz data. The dimensions for each variable is 3600 s/h x 20 samples/s
        
    #Takes 30-minute averages of each variable in w_vars    
    if hour.minute ==0:    
        for var in w_vars[1:]: #Skip WLEF - not in deataset
            df_w.loc[hour, var] = ds[var].mean(axis = 1).values[0:1800].mean(); 
            #first mean is average (over 20 measurements/sec), second mean is 30 min avg
            
    elif hour.minute ==30:
        for var in w_vars[1:]:
            df_w.loc[hour, var] = ds[var].mean(axis = 1).values[1800:3600].mean();
            
    df_w.to_csv(filepath)
            
    #This line can create a flat array of high-frequency data for a single variable ((n+1)th variable in w_vars array)
    # n = 1
    # flat_dat = np.append(flat_dat, ds[w_vars.iloc[n]].values.ravel()) 
