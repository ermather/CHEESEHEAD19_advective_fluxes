import os
from setup import * #Import setup module



''' Wind Profile Functions'''


#Stability correction
def stab_cor(stab_meth, z, z0, MO_tow, dtindex):
    #Calculates stability correction using Benoit 1977
    #Inputs: 
    #stab_meth = stability method
    #z = height
    #z0 = roughness length
    #MO_tow = MO length for given tower
    #dtindex - datetime index for dataset
    
    if pd.isnull(stab_meth): 
        Psi = np.zeros(shape = len(dtindex)) #Psi is zeros if not using stability correction
    
    elif stab_meth == 'Benoit': 
        #Only calculate correction factor for 
        x0 = (1-(15*z0/MO_tow.where(MO_tow < 0)))**(1/4)
        x = (1-(15*z/MO_tow.where(MO_tow < 0)))**(1/4)
        Psi = np.log((x0**2 + 1)*(x0 + 1)**2/((x**2 + 1)*(x + 1)**2)) + 2*(np.arctan(x) - np.arctan(x0))
        Psi = Psi.map(lambda n: 0 if pd.isnull(n) == True else n) #Set equal to zero where Psi is undefined
        #(where MO_tow >= 0 and where MO_tow is nan), note that n is an element in Psi
        # Psi = Psi.where(pd.isnull(MO_tow)==False)
        Psi = Psi.where(np.isfinite(MO_tow))
    
    return Psi



def wind_prof(tow, stab_meth, WS_use, MO_tow, dtindex, step = 0.1, z_Umes = 'use_mes'):
    
    #Step is step size to be used for integration
    
    #Tower constants
    h_veg = tc.veg_h.loc[tow] #vegetation height
    
    #If these aren't specified, use vals in tc for wind and sonic measurement heights
    if z_Umes == 'use_mes':
        z_Umes= tc.z1.loc[tow] #Wind measurement height
    z_max = round(tc.z_son.loc[tow], 1) #sonic height, rounded to nearest 0.1
    z_range = np.arange(0, z_max + step, step) #need wind profile up to level of sonic
    z0 = 0.1*h_veg #friction parameter
    z0_soil = 0.006 #From biophysic notes, check book Table_, ch.5
    d = tc.d.loc[tow] #zero plane displacement
    a = pd.read_pickle(intermed_filepath + 'a_from_lai.pickle')[dtindex[0]:dtindex[-1]][tow]

    
    WS_prof = pd.DataFrame(index = dtindex, columns = z_range) #windspeed profile
    
    Psi_z1 = stab_cor(stab_meth, z_Umes, z0, MO_tow, dtindex)
    ustar_calc = WS_use*k/np.log((z_Umes - d)/z0 + Psi_z1) #ustar to use in wind profile
    u_of_h = ustar_calc/k*np.log((h_veg-d)/z0 + stab_cor(stab_meth, h_veg, z0, MO_tow, dtindex)) #Wind speed at top of canopy
    u_of_soil = u_of_h*np.exp(a*(z0/h_veg - 1)) #WS at 0.1h, using exponential sub-canopy model
    ustar_soil = u_of_soil*k/np.log(z0/z0_soil) #Soil friction velocity, calculated from WS at 01.h
    

    for z in z_range: 
        WS_prof[z] = WS_of_z(stab_meth, z, z0, d, h_veg, a, z0_soil, u_of_h, ustar_calc, ustar_soil, MO_tow, dtindex)
            
    return WS_prof
    # return a
    

    WS_of_z(stab_meth, z, z0, d, h_veg, z0_soil, u_of_h, ustar_calc, ustar_soil, MO[tow], dtindex)
    
def WS_of_z(stab_meth, z, z0, d, h_veg, a, z0_soil, u_of_h, ustar_calc, ustar_soil, MO_tow, dtindex):
    
    if z >= h_veg:
        Psi = stab_cor(stab_meth, z, z0, MO_tow, dtindex)
        WS_z = ustar_calc/k*np.log(((z-d)/z0) + Psi)
    elif z>=z0:
        WS_z = u_of_h*np.exp(a*(z/h_veg - 1)) 
    elif z>0:
        WS_z = ustar_soil/k*np.log((z)/z0_soil)
    else:
        WS_z = 0
        
    return WS_z


'''get measurement heights to use for profiles'''
def get_zs(dat_type, tow):
    if dat_type == 'interp':
        zl, zm, zt = mlh, mmh, mth #For interpolated data, assume heights are mean heights
    elif dat_type == 'mes': #For measured data, use measurement heights
        towdat = tc.loc[tow]
        if tow == 'PFd' or tow == 'PFr': #Short towers
            zl, zm, zt = towdat.z1, np.nan, np.nan
        elif tow == 'PFc'or tow =='PFs': #Mid-height towers
            zl, zm, zt = towdat.z2, towdat.z1, np.nan
        elif tow =='PFA':
            zl, zm, zt = np.nan, np.nan, towdat.z1 
        elif tow == 'PFl': #25m tower, top mes were shifted up to mth (33m)
            zl, zm, zt = towdat.z3, towdat.z2, mth
        elif tow == 'PFe': #mid level is near 16m, interpolated to get at mmh (12m)
            zl, zm, zt = towdat.z3, mmh, towdat.z1
        else: #Regular towers
            zl, zm, zt = towdat.z3, towdat.z2, towdat.z1
       
    return zl, zm, zt #Low mid and hight levels
        
    
    
    
'''Calculate coordinates of upwind/downwind interpolation locatins'''
def updowncoords(dist, WD, tc):
    
    coords_dist = pd.DataFrame(index = WD.index, columns = towu2dirdim_cols)
    WD = WD.astype(float)
    #Locations dist km up and down wind from 
    coords_dist.loc[:, (slice(None), 'up', 'x')] = tc.x.values + dist*np.sin(WD*np.pi/180).values
    coords_dist.loc[:, (slice(None), 'up', 'y')] = tc.y.values + dist*np.cos(WD*np.pi/180).values
    coords_dist.loc[:, (slice(None), 'down', 'x')] = tc.x.values - dist*np.sin(WD*np.pi/180).values
    coords_dist.loc[:, (slice(None), 'down', 'y')] = tc.y.values - dist*np.cos(WD*np.pi/180).values
    
    #Create dictionary 
    
    return coords_dist



'''Calculate profiles of variables'''

def prof_calcs(tow, data, dat_type, dtindex, step = 0.1):
    towdat = tc.loc[tow] #Tower data
    
    zl, zm, zt = get_zs(dat_type, tow)
    
    h_tow = towdat.z_son #height of sonic measurements
    z_range = np.arange(0, h_tow + step, step)
    
    #Create dataframe to hold hrz gradient profile data
    data_prof = pd.DataFrame(columns = z_range, index = dtindex)
    
    # data_tow = data.loc[dtindex[0]:dtindex[-1], (slice(None), tow)]
    data_tow = data.loc[:, (slice(None), tow)]
    
    # should maybe make this more flexible for different tower heights #fixthis
    # Vertical gradients of horizontal gradients
    grad_up = (data_tow[30].values - data_tow[10].values)/(zt-zm) #(33-12)
    grad_low = (data_tow[10].values - data_tow[2].values)/(zm-zl) #(12-2)
    
    for z in z_range:
        if z >= zm:

            data_prof[z] = data_tow[10] + grad_up*(z-zm)

        else:
            data_prof[z] = data_tow[10] - grad_low*(zm-z)
                
            
    return data_prof