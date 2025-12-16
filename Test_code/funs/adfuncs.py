#Read in all data

#Read in tower coords and set up dataframe
tc = pd.read_csv('../Inputs/tower_coords.csv') #datafrane of tower coordinates
tc.set_index('Tower', inplace = True)

#Convert lat and lon to coordinates in km, origin currently set to PFA coordinates, should change this since PFk is out of range
tc['x'] = (tc.Lon - tc.Lon['PFA'])*111*np.cos(tc.Lat['PFA']*np.pi/180)
tc['y'] = (tc.Lat - tc.Lat['PFA'])*111 


#Read in TA, H2O, and wind data
TA = pd.read_pickle('../Inputs/TA_orig_all.pickle')
H2O = pd.read_pickle('../Inputs/H2O_orig_all.pickle')
TA_grad = pd.read_pickle('../Inputs/TA_int_all.pickle')
H2O_grad = pd.read_pickle('../Inputs/H2O_int_all.pickle')
hrz_wind = pd.read_pickle('../Inputs/wind_dat.pickle')
vert_wind = pd.read_pickle('../Inputs/w_all.pickle')
PA_top = pd.read_pickle('../Inputs/PA_top.pickle')

#Rename vert_wind columns
tow_names =  dict(zip(tc.wind_var[1:], tc.index[1:]))
vert_wind = vert_wind.rename(columns = tow_names)






def WS_top_calc(tower, d_rat = 0.65):
    global tc #Use tc defined in main notebook
    
    k = 0.4 #
    
    #Measurement height (z1 for all towers except PFA)

    h_tow = tc.z1.loc[tower]
        
    h_veg = tc.veg_h.loc[tower] #Vegetation height
    WS_calc = (hrz_wind.loc['Ustar'][tower]/k)*np.log((h_tow - d_rat*h_veg)/(0.1*h_veg))
    return WS_calc



def wind_prof(tower, step = 0.1):
    
    global tc #Use tc defined in main notebook
    
    towdat = tc.loc[tower] #Tower data
    
    #Step is step size to be used for integration
    k = 0.4 #Von Karmann constant
    h_veg = tc.veg_h.loc[tower] #vegetation height
    h_tow = tc.z1.loc[tower] #tower height
    z_range = np.arange(0, h_tow + step, step)
    z_m = 0.1*h_veg #momentum friction parameter
    z_m_soil = 0.006 #From biophysic notes, check book Table_, ch.5
    
    d = 0.65*h_veg #zero plane displacement
    
    # d_mean = h_tow - z_m*np.exp(k*WS_top[tower]/Ustar[tower])#zero plane displacement, mean of val calculated for all timestamps
    
    WS_prof = pd.DataFrame() #windspeed profile
    
    
    WS_calc = WS_top_calc(tower)
    WS_calc_rat = hrz_wind.loc['WS_top'][tower]/WS_calc
    u_of_h = hrz_wind.loc['Ustar'][tower]/k*np.log((h_veg-d)/z_m)*WS_calc_rat #wind speed at top of canopy
    u_of_soil = u_of_h*np.exp(0.1 - 1) #WS at 0.1h, using exponential 
    ustar_soil = u_of_soil*k/np.log(0.1*h/z_m_soil) #Soil friction velocity, calculated from WS at 01.h
    
    
    for z in z_range: 
        if z >= h_veg:
            WS_prof[z] = hrz_wind.loc['Ustar'][tower]/k*np.log((z-d)/z_m)*WS_calc_rat
        elif z>=0.1*h_veg:
            WS_prof[z] = u_of_h*np.exp(z/h - 1)
        elif z>0:
            WS_prof[z] = ustar_soil/k*np.log((z)/z_m_soil)
        else:
            WS_prof[0] = 0
            
    return WS_prof



    
def prof_calcs(tower, step = 0.1):
    
    global tc #Use tc defined in main notebook
    
    towdat = tc.loc[tower] #Tower data
    WS_prof = wind_prof(tower)
    
    h_tow = tc.z1.loc[tower] #tower height
    z_range = np.arange(0, h_tow + step, step)
    
    T_grad_NS, T_grad_EW = pd.DataFrame(), pd.DataFrame()
    H2O_grad_NS, H2O_grad_EW = pd.DataFrame(), pd.DataFrame()
    
   
    for z in z_range:     
        if towdat.z3 == 0:
            pass
        else:
            #should maybe make this more flexible for different tower heights
            T_grad_slope_NS_up = (TA_grad.loc['NS'].loc[30][tower] - TA_grad.loc['NS'].loc[10][tower])/20 #Slope of temp gradient from 10 to 30 m, NS
            T_grad_slope_EW_up = (TA_grad.loc['EW'].loc[30][tower] - TA_grad.loc['EW'].loc[10][tower])/20 #Slope of temp gradient from 10 to 30 m, EW
            T_grad_slope_NS_low = (TA_grad.loc['NS'].loc[10][tower] - TA_grad.loc['NS'].loc[2][tower])/8 #Slope of temp gradient from 2 to 10 m, NS
            T_grad_slope_EW_low = (TA_grad.loc['EW'].loc[10][tower] - TA_grad.loc['EW'].loc[2][tower])/8 #Slope of temp gradient from 2 to 10 m, EW
            
            H2O_grad_slope_NS_up = (H2O_grad.loc['NS'].loc[30][tower] - H2O_grad.loc['NS'].loc[10][tower])/20 #Slope of temp gradient from 10 to 30 m, NS
            H2O_grad_slope_EW_up = (H2O_grad.loc['EW'].loc[30][tower] - H2O_grad.loc['EW'].loc[10][tower])/20 #Slope of temp gradient from 10 to 30 m, EW
            H2O_grad_slope_NS_low = (H2O_grad.loc['NS'].loc[10][tower] - H2O_grad.loc['NS'].loc[2][tower])/8 #Slope of temp gradient from 2 to 10 m, NS
            H2O_grad_slope_EW_low = (H2O_grad.loc['EW'].loc[10][tower] - H2O_grad.loc['EW'].loc[2][tower])/8 #Slope of temp gradient from 2 to 10 m, EW
            
            if z >= towdat.z2:
                T_grad_NS[z] = TA_grad.loc['NS'].loc[10][tower] + T_grad_slope_NS_up*(z-towdat.z2)
                T_grad_EW[z] = TA_grad.loc['EW'].loc[10][tower] + T_grad_slope_EW_up*(z-towdat.z2)
                
                H2O_grad_NS[z] = H2O_grad.loc['NS'].loc[10][tower] + H2O_grad_slope_NS_up*(z-towdat.z2)
                H2O_grad_EW[z] = H2O_grad.loc['EW'].loc[10][tower] + H2O_grad_slope_EW_up*(z-towdat.z2)
            else:
                T_grad_NS[z] = TA_grad.loc['NS'].loc[10][tower] - T_grad_slope_NS_low*(towdat.z2 - z)
                T_grad_EW[z] = TA_grad.loc['EW'].loc[10][tower] - T_grad_slope_EW_low*(towdat.z2 - z)
                
                H2O_grad_NS[z] = H2O_grad.loc['NS'].loc[10][tower] - H2O_grad_slope_NS_low*(towdat.z2 - z)
                H2O_grad_EW[z] = H2O_grad.loc['EW'].loc[10][tower] - H2O_grad_slope_EW_low*(towdat.z2 - z)

                

            
    return WS_prof, T_grad_NS, T_grad_EW, H2O_grad_NS, H2O_grad_EW





#Hrz H and LE advection calcs

def hrz_ad_calc(tower, step = 0.1):
    

    WS_prof, T_grad_NS, T_grad_EW, H2O_grad_NS, H2O_grad_EW = prof_calcs(tower, step)
    
    TA_K = TA.loc[30][tower] + CtoK
    
    PA_Pa = PA_top[tower]*1000 #kPa to Pa **this is at top of tower(~400m), need PA val for 30m and below**
    rho = mm*PA_Pa/TA_K/R #Dry air density at each tower [kg/m^3] **should probably use virtual temperature here** 
    z = step #np.asarray(tc.loc[tower].z1) #Height of each tower

    SN_wind = WS_prof.mul(np.sin(hrz_wind.loc['WD'][tower]*np.pi/180), axis = 0)
    WE_wind = WS_prof.mul(np.cos(hrz_wind.loc['WD'][tower]*np.pi/180), axis = 0)
    
    #H advection
    H_SN = T_grad_NS*SN_wind.mul(-rho*z*cp/2000, axis = 0) #Advection of sensible heat by southerly wind [W/m^2]
    H_WE = T_grad_EW*WE_wind.mul(-rho*z*cp/2000, axis = 0) #Advection of sensible heat by westerly wind [W/m^2]
    H_ad = H_SN + H_WE #Total advection

    #LE advection
    LE_SN = H2O_grad_NS*SN_wind.mul(-z*L/2000) #Advection of latend heat by southerly wind [W/m^2]
    LE_WE = H2O_grad_EW*WE_wind.mul(-z*L/2000 )#Advection of latent hear by westerly wind [W/m^2]
    LE_ad = (LE_SN + LE_WE)#*rho/rho #LE advection (the *rho/rho just makes it a df with the same tower and timestamp as indices)


    # hrz_ad = (H_ad + LE_ad)
    H_ad_tot = H_ad.sum(axis = 1)
    LE_ad_tot = LE_ad.sum(axis = 1)
    
    return H_ad_tot, LE_ad_tot
    #Deal with water vapor later
    # md_w = e/(R*TA_K) #Molar density of water vapor [mol/m^3], used in SLE calc