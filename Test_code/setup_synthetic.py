#This file has all the basic setup needed for various notebooks and py files used in Cheesehead advection project
#Includes: Imports for commonly used packages, Constants, setup of tc (var containing basic tower info), and data

#Imports and settings
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime as dt
import cdflib
from scipy.interpolate import RBFInterpolator
import pickle



#Use this to ignore warning in inerpolation script, maybe get rid of later to check for other issues
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

pd.set_option('display.max_columns', None) #Change settings in Pandas so it shows all columns
 
idx = pd.IndexSlice #Used for slicing dataframes
    


    
#Constants
#Met
R = 8.314 #kgm^2s^-2K^-1mol^-1
h = 30 #Measurement height in m
L = 2.5*10**3 #Latent heat of vapoization of water at 0 deg C [J/g], should try to calc based on dependency on C
cp = 1005 #J/kg*K
mm = 0.02896 #kg/mol molar mass dry air
mm_w = 1801 # molar mass of water [g/mol]
CtoK = 273.15 #Conversion from Celcius to Kelvin
k = 0.4 #von Karmann constant
   
    
#Tower info   
tc = pd.read_csv('../UnitTest/Inputs/tower_coords.csv', index_col = 'Tower') #datafrane of tower coordinates


#Read in run cases
cases = pd.read_csv('../Inputs/Cases.csv', index_col = 'case') #dataframe of tower coordinates


#Output base folder filepaths
input_filepath = '../UnitTest/Inputs/'
intermed_filepath = '../UnitTest/Intermed_data/'
output_filepath = '../UnitTest/Outputs/'

'''---------------------------------------------------------------------------------------'''



#Function to read in universal data
#Reads in all data unless start and end date are specified
def readdata(var, startdate = dt.datetime(2019, 6, 1), enddate = dt.datetime(2019, 11, 1)):
    unidata_paths = {
        'TA':'../UnitTest/Inputs/TA_orig_all.pickle',
        'H2O':'../UnitTest/Inputs/H2O_orig_all.pickle',
        'hrz_wind': '../UnitTest/Inputs/wind_dat.pickle',
        'PA_top': '../UnitTest/Inputs/PA_top.pickle',
        'vert_wind_sonic': '../UnitTest/Inputs/w_all.pickle', 
        'dtindex': '../Inputs/dtindex.pickle'} 
    
    if var == 'PA_top' or var =='vert_wind_sonic':
        #No levels for PA, dtindex, select start to end date
        data = pd.read_pickle(unidata_paths[var]).loc[startdate:enddate] 
    elif var == 'dtindex':
        full_dtindex = pd.read_pickle(unidata_paths[var]) #read in full dtindex
        data = full_dtindex[(full_dtindex>=startdate) & (full_dtindex <= enddate)] #cut to desired times
    
    else:
        #Select all levels, start to end date, all tows
        data = pd.read_pickle(unidata_paths[var]).loc[idx[:, startdate:enddate], :] 

    return data


def make_multi_df(col_lists, level_names, index):
    # dtindex = readdata('dtindex', startdate, enddate)
    col_ind = pd.MultiIndex.from_product(col_lists, names=level_names)
    df = pd.DataFrame(index = index, columns = col_ind)
    return df




NSEWdirlist = ['N', 'S', 'E', 'W'] #list of directions
c2_dirlist = ['NS','EW'] #Cardinal 2 dir (N/S, E/W list)
u2_dirlist = ['up', 'down'] #Up/downwind 2 dir list
levlist = [2, 10, 30]
varlist = ['TA', 'H2O']
dimlist = ['x', 'y']

towlist = tc.index #List of towers

#Multiindices for dataframes, names give nesting order 
varlevtow4dir_cols = pd.MultiIndex.from_product([varlist, levlist, towlist, NSEWdirlist], names=["var", "lev", "tow", "direct"])
varlevtow_cols = pd.MultiIndex.from_product([varlist, levlist, towlist], names=["var", "lev", "tow"])
varlev_cols = pd.MultiIndex.from_product([varlist, levlist], names=["var", "lev"])
tow4dir_cols = pd.MultiIndex.from_product([towlist, NSEWdirlist], names=["tower", "direct"]) #for columns of interp data df
varlevtowc2dir_cols = pd.MultiIndex.from_product([varlist, levlist, towlist, c2_dirlist], names=["var", "lev", "tow", "direct"])
towc2dir_cols = pd.MultiIndex.from_product([towlist, c2_dirlist], names=["tower", "direct"]) 
towu2dirdim_cols = pd.MultiIndex.from_product([towlist, u2_dirlist, dimlist], names=["tow", "dir", "dim"]) #for updown corrds, into dict of directions


