from setup import *

'''------------------------------------------Coordinates-------------------------------------------------'''


#Function to calculate coordinates at different dirs NSEW from each tower
def NSEWcoords(dist, tc = tc): 
    tc_int = tc.copy() #Make copy of tc to put interpolation coords in

    coords_dist = {} #Dictionary of NSEW coords at the specified distance
    
    #N, S, E, W coords (in km)
    tc_int['y_N'] = tc_int['y'] + dist
    tc_int['y_S'] = tc_int['y'] - dist
    tc_int['x_E'] = tc_int['x'] + dist
    tc_int['x_W'] = tc_int['x'] - dist

    #Dictionary holding 2D array of coords for each tower in each direction 
    #Dims are tower, time, and direction (x/y)
    coords_dist['N'] = np.stack((tc_int.x.values, tc_int.y_N.values)).T
    coords_dist['S'] = np.stack((tc_int.x.values, tc_int.y_S.values)).T
    coords_dist['E'] = np.stack((tc_int.x_E.values, tc_int.y.values)).T
    coords_dist['W'] = np.stack((tc_int.x_W.values, tc_int.y.values)).T
    
    return coords_dist


'''---------------------------------------------Interpolation--------------------------------------------------'''

def interp_NSEW(var, lev, coords_dist, kernel = 'thin_plate_spline'):
     
    #Select needed data
    data = unidata[var].loc[lev] #level is part of index (maybe better as column but probably not worth changing)
    
    #Create interp_dat dataframe with nested columns (tower>direction)
    interp_dat = make_multi_df([towlist, NSEWdirlist], ["tower", "direct"])
    
    #Loop through timestep and interpolate
    for i in dtindex[0:1]:
        dat_use = data.loc[data.index == i].values[0, :] #For each loop interpolate with just one row of data
        
        dat_nn = dat_use[nanmask[var][lev].loc[i]] #dat_use without nan values
        ar_nn = ar[nanmask[var][lev].loc[i]] #coordinates of towers that have finite values
        
        if nancount[var][lev].loc[i] >= 16:
            pass
            
        else:

            #Run interpolator
            interpfunction = RBFInterpolator(ar_nn, dat_nn, kernel = kernel)
            
            #Apply interpolation function to coords of points in each direction
            for direct in NSEWdirlist:
                interp_dat.loc[i, idx[:, direct]] = interpfunction(coords_dist[direct])
                
            
    return interp_dat


'''---------------------------------------------Gradients--------------------------------------------------'''

#Gradient calcs

def NSEWgradcalc(interpdist = 1, kernel = 'thin_plate_spline'):
    
    #Calc coords of interp locations
    coords_dist = NSEWcoords(interpdist)
    
    #Create output grad df
    grad_updown = pd.DataFrame(index = dtindex, columns = varlevtow_cols) 
    
    for var in ['TA', 'H2O']: #Loop through both variable
        for lev in [2, 10, 30]: #Loop through all three levels
            pointdat = interp_NSEW(var, lev, coords_dist, kernel) #gradients for particular var, lev, all tow, di
            grad_NSEW = NSEW_point_to_grad(pointdat, interpdist)
            grad_updown[var, lev] = grad_NSEW_to_updown(grad_NSEW)
            
    
    return grad_updown


#Converts point data to gradients
def NSEW_point_to_grad(pointdat, dist):
    grad = pd.DataFrame(index = dtindex, columns = towc2dir_cols)
    grad.loc[:, (slice(None), 'NS')] = \
    (pointdat.loc[:, (slice(None), 'N')].values - pointdat.loc[:, (slice(None), 'S')].values)/dist
    grad.loc[:, (slice(None), 'EW')] = \
    (pointdat.loc[:, (slice(None), 'E')].values - pointdat.loc[:, (slice(None), 'W')].values)/dist

    return grad

def grad_NSEW_to_updown(grad_NSEW):
    WD = unidata['hrz_wind'].loc['WD']*np.pi/180 #Convert to radians
    NS_grad = grad_NSEW.xs("NS", level="direct", axis=1, drop_level=True) #grad_NSEW.loc[:, (slice(None), 'NS')]
    EW_grad = grad_NSEW.xs("NS", level="direct", axis=1, drop_level=True)
    #Note this is the same as grad_NSEW.loc[:, (slice(None), 'EW')].droplevel(leve = 'direct', axis = 1)
    grad_updown = NS_grad*np.cos(WD) + EW_grad*np.sin(WD)
    return grad_updown