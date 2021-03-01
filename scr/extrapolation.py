# -*- coding: utf-8 -*-
"""Interpolation and extrapolation.

This module contains all the functions required to fit a power law to Wind 
Ninja's data and to extrapolate them to the required altitudes.

@author: Q.Abeille (2020)
"""
import numpy as np
import scipy.optimize as opt
import scipy.interpolate as intp


def get_Zlist_pos(x_pos,y_pos,points):
    """get_Zlist_pos.
    
    A function retrieving the altitudes of the different data points at a given
    position.

    Parameters
    ----------
    x_pos : float
        Coordinate of the point along the x-axis.
    y_pos : float
        Coordinate of the point along the y-axis.
    points : Numpy array
        Coordinates in (x,y,z) convention of the data points.

    Returns
    -------
    Zarg_pos : Numpy array of int
        Arguments of the points with (x_pos, y_pos) as (x, y) coordinates in 
        the data array.
    Zxy : Numpy array of float
        Z coordinates of points with (x_pos, y_pos) as (x, y) coordinates.

    """
    cond_X = points[:,0]==x_pos  # Elts situés en x_pos
    cond_Y = points[:,1]==y_pos  # Elts situés en y_pos
    cond = np.multiply(cond_X, cond_Y)  # Elts situés en x_pos, y_pos
    
    Zarg_pos = np.argwhere(cond==True)  # Récupération de la position des elts
    Zxy = []
    
    # Récupération des positions en z du point (x_pos, y_pos)
    for k in range (0,len(Zarg_pos)):
        cur_arg = Zarg_pos[k]
        Zxy = np.append(Zxy,points[cur_arg,2])
        
    return Zarg_pos, Zxy

def alt_to_elev(points_surf,points_wind):
    """Switching from altitude to elevation.
    
    This function gets the wind field in the "elevation from the ground"
    convention.

    Parameters
    ----------
    points_surf : Numpy array of float
        (X,Y,Z) coordinates of the surface at each (X,Y) measurement point. The
        Z coordinate must be the altitude above the sea.
    points_wind : Numpy array of float
        (X,Y,Z) coordinates of each measurement point of the wind field. The Z
        coordinate must be the altitude above the sea.

    Returns
    -------
    points_wind_elev : Numpy array of float
        Coordinates of each measurement point of the wind field with the Z 
        coordinate being the elevation from the ground.

    """
    # Getting the coordinates of all the ground points.
    X = np.array(points_surf[:,0])
    Y = np.array(points_surf[:,1])
    
    Zelev_wind = points_wind[:,2]
    points_wind_elev = points_wind  # wind field in elevation coonvention
    
    for k in range(len(X)):
        x = X[k]
        y = Y[k]
        
        # Retrieving the altitude of the surface and the measurement points of
        # current ground position.
        Zwind_argpos, Zwind_pos = get_Zlist_pos(x,y,points_wind)
        Zsurf_pos = np.unique(get_Zlist_pos(x,y,points_surf)[1])
        
        # Calculating the elevation from the surface of each measurement point.
        for i in range(len(Zwind_argpos)):
            cur_pos = Zwind_argpos[i]
            Zelev_wind[cur_pos] += -Zsurf_pos
            
    points_wind_elev[:,2] = Zelev_wind
            
    return points_wind_elev

def get_vert_profile(x_pos,y_pos,field):
    """Get vertical profile.
    
    This function gets the vertical profile of a scalar field at a given
    ground position.
    
    Parameters
    ----------
    x_pos : float
        X coordinate of the ground position.
    y_pos : float
        Y coordinate of the ground position.
    field : Numpy array of float
        The scalar field from wich the vertical profile must be extracted. Must
        be in the format (X,Y,Z in elevation convention, Field value).

    Returns
    -------
    Zxy : array of float
        Z coordinates of the vertical profile.
    Fieldxy : array of float
        Field values along the z-axis at the ground point.
    """
    cond_X = field[:,0]==x_pos  # Elts situés en x_pos.
    cond_Y = field[:,1]==y_pos  # Elts situés en y_pos.
    cond = np.multiply(cond_X, cond_Y)  # Elts situés en x_pos, y_pos.
    
    ZUarg_pos = np.argwhere(cond==True)  # Récupération de la position des elts.
    
    # Tableaux de sortie
    Zxy = np.array([])
    Fieldxy = np.array([])
    
    # Remplissage des tableaux
    for k in range (0,len(ZUarg_pos)):
        cur_arg = ZUarg_pos[k]
        Zxy = np.append(Zxy,field[cur_arg,2])
        Fieldxy = np.append(Fieldxy,field[cur_arg,3])
        
    # Elimination des points sous la surface
    check = Zxy<0
    Zneg_argpos = np.argwhere(check == True)
    if len(Zneg_argpos)>0:
        arg_surf = np.max(Zneg_argpos)
        Zxy = Zxy[arg_surf+1:]
        Fieldxy = Fieldxy[arg_surf+1:]
        
    return(Zxy,Fieldxy)

def get_extrap_law(Z,wind,law):
    """Interpolation along z-axis.
    
    Fitting of a given law to the data along z-axis. The data points for the
    fitting are automatically chosen by the function.

    Parameters
    ----------
    Z : Array of float
        Z positions of the data points.
    wind : Array of float
        Wind value at the data points.
    law : function
        The function to fit the data to.

    Returns
    -------
    Zextrap_points : Array of float
        Z positions of the data points on which the fitting was done.
    Vextrap_points : Array of float
        Wind value of the data points on which the fitting was done.
    best_val : Array of float
        Law parameters that makes it fit the selected data points.

    """
    # Initialization
    Vextrap_points = wind
    Zextrap_points = Z
    success = False
    
    # Error indicators.
    r2_new = np.inf  
    r2_old = np.inf
    while len(Vextrap_points)>1 and success == False:  # Minimum of 2 points to do the fitting.
        try:
            # Fitting of the law to the selected data points.
            best_val, covar = opt.curve_fit(law, Zextrap_points, Vextrap_points, p0=[2.5,0.128])
        except RuntimeError:  # Fitting did not succeed.
            # Selection of the highest points
            Vextrap_points = Vextrap_points[1:]
            Zextrap_points = Zextrap_points[1:]
        else:
            # Calculation of the error to the data points.
            r2_old = r2_new
            res = law(Zextrap_points,best_val[0],best_val[1]) - Vextrap_points
            r2_new = np.sum(np.multiply(res,res))
            
            if r2_new>r2_old:  # Minimal error has been reached.
                success = True
            else:
                # Selection of the highest points.
                Vextrap_points = Vextrap_points[1:]
                Zextrap_points = Zextrap_points[1:]
    return Zextrap_points,Vextrap_points, best_val

def main(points,wind,points_surf,elev_max, step):
    """Compute the wind cube.
    
    This function uses all the above functions to create an extrapolated wind
    cube from the initial data. It is the main function of the file.
    
    Parameters
    ----------
    points : Narray of floats (x,y,z)
        Contains the x,y,z coordinates of each data point of the initial wind 
        field.
    wind : Narray of floats (u,v,w)
        Contains the values of the wind field components at each of the initial
        data points.
    points_surf : Narray of floats (x,y,z)
        Contains the location of the ground surface in the fiels.
    elev_max : float
        Maximal elevation from the ground that is wanted for the extrapolated 
        field.
    step : float
        Size of the steps (in m) required for the extrapolation.

    Returns
    -------
    extrap_field : Narray of floats
        Contains the extrapolated wind field in format (x,y,z,U,V,W).
    Z_surf : Narray of floats
        Array containing the altitude of the ground surface for each of the 
        (x,y,z) point of extrap_field.

    """    
    # Decomposition along the axes
    X = np.array(points[:,0])
    Y = np.array(points[:,1])
    
    elev_wind = alt_to_elev(points_surf,points)    
    
    # Retrieving the wind components
    U = np.array(wind[:,0])
    V = np.array(wind[:,1])
    W = np.array(wind[:,2])
    
    # Mesh
    X_tick = np.unique(X)
    Y_tick = np.unique(Y)
    Z_tick = get_Zlist_pos(X_tick[0], Y_tick[0], points)[1]
    Z_tick = np.extract(Z_tick>=0, Z_tick)
    Z_tick = np.round(Z_tick)
    
    # Wind field for each component
    U_field = np.array([X,Y,elev_wind[:,2],U])
    U_field = np.transpose(U_field)
    
    V_field = np.array([X,Y,elev_wind[:,2],V])
    V_field = np.transpose(V_field)
    
    W_field = np.array([X,Y,elev_wind[:,2],W])
    W_field = np.transpose(W_field)
    
    # Extrapolation array
    X_extrap = []
    Y_extrap = []
    Z_extrap = []
    U_extrap = []
    V_extrap = []
    W_extrap = []
    Z_surf = []
    
    # For each point of the mesh
    for i in range(len(X_tick)):
        for j in range(len(Y_tick)):
            # Current position
            x_pos = X_tick[i]
            y_pos = Y_tick[j]
            
            # Current surface altitude
            z_surf = get_Zlist_pos(x_pos,y_pos,points_surf)[1]
            
            # Retrieving the vertical profile of each component
            z_pos, U_pos = get_vert_profile(x_pos, y_pos, U_field)
            V_pos = get_vert_profile(x_pos,y_pos,V_field)[1]
            W_pos = get_vert_profile(x_pos,y_pos,W_field)[1]
            
            # Interpolation over Z_tick to have a regular mesh over z-axis
            Ucalc = intp.interp1d(z_pos,U_pos, bounds_error = False, 
                                  fill_value='extrapolate')
            Vcalc = intp.interp1d(z_pos,V_pos, bounds_error = False, 
                                  fill_value='extrapolate')
            Wcalc = intp.interp1d(z_pos,W_pos, bounds_error = False, 
                                  fill_value='extrapolate')
            
            U_pos = Ucalc(Z_tick)
            V_pos = Vcalc(Z_tick)
            W_pos = Wcalc(Z_tick)
            
            # Adding the coordinates of the current point
            X_extrap = np.append(X_extrap,np.ones(len(Z_tick))*x_pos)
            Y_extrap = np.append(Y_extrap,np.ones(len(Z_tick))*y_pos)
            Z_extrap = np.append(Z_extrap,Z_tick)  # In altitude conv
            
            Z_surf = np.append(Z_surf,np.ones(len(Z_tick))*z_surf)
            
            # Adding the data points of the initial fields
            U_extrap = np.append(U_extrap,U_pos)
            V_extrap = np.append(V_extrap,V_pos)
            W_extrap = np.append(W_extrap,W_pos)
            
            # Extrapolation of the different fields
            Z_pos_extrap = np.linspace(Z_tick[-1],elev_max,step)
            Z_pos_extrap = Z_pos_extrap[1:]
            param_pow = get_extrap_law(Z_tick,U_pos,power_law)[2]
            U_pos_extrap = power_law(Z_pos_extrap,param_pow[0], param_pow[1])
            
            param_pow = get_extrap_law(Z_tick,V_pos,power_law)[2]
            V_pos_extrap = power_law(Z_pos_extrap,param_pow[0], param_pow[1])
            
            param_pow = get_extrap_law(Z_tick,W_pos,power_law)[2]
            W_pos_extrap = power_law(Z_pos_extrap,param_pow[0], param_pow[1])
            
            # Adding the coordinates of the extrapolated points
            X_extrap = np.append(X_extrap,np.ones(len(Z_pos_extrap))*x_pos)
            Y_extrap = np.append(Y_extrap,np.ones(len(Z_pos_extrap))*y_pos)
            Z_extrap = np.append(Z_extrap, Z_pos_extrap)
            
            Z_surf = np.append(Z_surf,np.ones(len(Z_pos_extrap))*z_surf)
            
            # Adding the extrapolated points
            U_extrap = np.append(U_extrap,U_pos_extrap)
            V_extrap = np.append(V_extrap,V_pos_extrap)
            W_extrap = np.append(W_extrap,W_pos_extrap)     
              
    # Creation of the extrapolated field
    extrap_field = np.array([X_extrap, Y_extrap, Z_extrap, U_extrap, 
                                   V_extrap, W_extrap, Z_surf])
    
    # Transposition into a vertical array
    extrap_field = np.transpose(extrap_field)
    
    return extrap_field

def power_law(z,U10,alpha):
    """Power law.
    
    Empiric power law of the evolution of wind with the altitude.

    Parameters
    ----------
    z : float
        Height from the ground.
    U10 : float
        Velocity at 10m height.
    alpha : float
        Alpha parameter.

    Returns
    -------
    float
        value given by the power law at height z with parameters U10 and alpha.

    """
    return U10*np.power((z)/10., alpha)