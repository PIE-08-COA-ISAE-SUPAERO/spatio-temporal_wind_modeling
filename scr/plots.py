# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 17:10:12 2021

@author: abeil
"""

#%% Imports
import numpy as np
import matplotlib.pyplot as plt
#plt.style.use('presentation.mplstyle')
from mpl_toolkits.mplot3d import axes3d
from scipy.interpolate import LinearNDInterpolator
import extrapolation as interp

#%% Récupération fichier texte
data = np.loadtxt('extrap_field.txt')

#Récupération du maillage horizontal
X = data[:,0]
Y = data[:,1]
Z_surf = data[:,6]
X_tick = np.unique(X)
Y_tick = np.unique(Y)
Z_tick = interp.get_Zlist_pos(X_tick[0], Y_tick[0], data[:,0:3])[1]

#%% Cube de vent

def get_tick_argminmax(cond_min,cond_max):
    """
    
    Get the boundary arguments for the interpolation domain.

    Parameters
    ----------
    cond_min : Narray of bool
        Numpy array of booleans representing the condition value >= min_value.
    cond_max : Narray of bool
        Numpy array of booleans representing the condition value <= min_value.

    Returns
    -------
    arg_min : int
        argument of the lower boundary condition.
    arg_max : TYPE
        argument of the upper boundary condition.

    """
    #Initialization
    arg_max = len(cond_max)-1
    arg_min = 0
    
    # retrieving the lower boundary condition going from the bottom
    while cond_min[arg_min]==False:
        arg_min += 1
    
    #retrieving the upper boundary condition going from the top
    while cond_max[arg_max]==False:
        arg_max += -1
    
    if arg_min>arg_max:  # The domain is between two consecutive points
        # Switching the values
        temp = arg_min
        arg_min = arg_max
        arg_max = temp
        
    return arg_min, arg_max

def get_interv(xlim,ylim,zlim):
    """
    Get the interpolation ranges.

    Parameters
    ----------
    xlim : Narray of floats
        Numpy array of size 2 containing the limits of the x-axis range.
    ylim : Narray of floats
        Numpy array of size 2 containing the limits of the y-axis range.
    zlim : Narray of floats
        Numpy array of size 2 containing the limits of the z-axis range.

    Returns
    -------
    X_range : Narray of floats
        Numpy array containing the x-coord of the points inside the 
        interpolation domain.
    Y_range : Narray of floats
        Numpy array containing the y-coord of the points inside the 
        interpolation domain.
    Z_range : Narray of floats
        Numpy array containing the z-coord of the points inside the 
        interpolation domain.

    """
    
    # Retrieving the min and the max values along each axis
    x_min = xlim[0]
    x_max = xlim[1]
    y_min = ylim[0]
    y_max = ylim[1]
    z_min = zlim[0]
    z_max = zlim[1]
    
    # Getting the arguments for the limits along x-axis
    cond_X1 = x_min <= X_tick
    cond_X2 = X_tick <= x_max
    arg_min, arg_max = get_tick_argminmax(cond_X1, cond_X2)
    # Array with all the arguments
    Xarg_pos = np.arange(arg_min, arg_max+1, step =1)
    
    # Getting the arguments for the limits along y-axis
    cond_Y1 = y_min <= Y_tick 
    cond_Y2 = Y_tick<= y_max
    arg_min, arg_max = get_tick_argminmax(cond_Y1, cond_Y2)
    # Array with all the arguments
    Yarg_pos = np.arange(arg_min, arg_max+1, step =1)
    
    # Getting the arguments for the limits along z-axis
    cond_Z1 = z_min <= Z_tick 
    cond_Z2 = Z_tick<= z_max
    arg_min, arg_max = get_tick_argminmax(cond_Z1, cond_Z2)
    # Array with all the arguments
    Zarg_pos = np.arange(arg_min, arg_max+1, step =1)
    
    # Expending the interpolation range to avoid interpolation the bounds of the
    # selected data field to be strictly inside the interpolation domain
    # (interpolation wouldn't work in tha case)
    xarg_min = int(max(Xarg_pos[0]-1, 0))
    xarg_max = int(min(Xarg_pos[-1]+1, len(X_tick)-1))
    yarg_min = int(max(Yarg_pos[0]-1, 0))
    yarg_max = int(min(Yarg_pos[-1]+1, len(Y_tick)-1))
    zarg_min = int(max(Zarg_pos[0]-1, 0))
    zarg_max = int(min(Zarg_pos[-1]+1, len(Z_tick)-1))
    
    # Cutting the interpolation ranges from the axis
    X_range = X_tick[xarg_min : xarg_max+1]
    Y_range = Y_tick[yarg_min : yarg_max+1]
    Z_range = Z_tick[zarg_min : zarg_max+1]

    return X_range, Y_range, Z_range


def get_interp_data(xlim, ylim, zlim):
    """
    
    Calculates an interpolation from the wind_cube data to a defined domain.

    Parameters
    ----------
    xlim : Narray of floats
        Numpy array of size 2 containing the limits of the x-axis range.
    ylim : Narray of floats
        Numpy array of size 2 containing the limits of the y-axis range.
    zlim : Narray of floats
        Numpy array of size 2 containing the limits of the z-axis range.

    Returns
    -------
    X_mesh : Narray of floats
        3D-mesh for the interpolated x coordinates.
    Y_mesh : Narray of floats
        3D-mesh for the interpolated x coordinates.
    Z_mesh : Narray of floats
        3D-mesh for the interpolated x coordinates.
    Uinterp : Narray of floats
        3D-mesh for the interpolated wind speed component along x-axis.
    Vinterp : Narray of floats
        3D-mesh for the interpolated wind speed component along y-axis.
    Winterp : Narray of floats
        3D-mesh for the interpolated wind speed component along z-axis.
    Sinterp : Narray of floats
        3D-mesh for the interpolated surface altitude.

    """
    
    # Initialization
    if zlim[0]==zlim[1]:   # Interpolation at a constant altitude -> faster to 
        data_interp = data # use the whole domain
    else:
        # Calculating the ranges along each axis
        X_range, Y_range, Z_range = get_interv(xlim,ylim,zlim)
        
        if len(X_range)>len(X_tick)/2: # More than half of the cube required
            data_interp = data         # faster to use the whole domain
        else:
            data_interp = [[]]
            
            # Retrieving the lines of data that will be used in the interpoaltion
            for x_pos in X_range:
                for y_pos in Y_range:
                    for z_pos in Z_range:
                        cond_X = data[:,0]==x_pos  # Elts located at x_pos.
                        cond_Y = data[:,1]==y_pos  # Elts located at y_pos.
                        cond_Z = data[:,2]==z_pos  # Elts located at z_pos.
                        cond = np.multiply(cond_X, cond_Y) # Elts at (x_pos,y_pos)
                        cond = np.multiply(cond, cond_Z)  # Elts at (x_pos,y_pos,z_pos)
                    
                        # Arguments of elements at the point
                        arg_pos = np.argwhere(cond==True)
                        
                        # Appending the lines to data_interp
                        for k in arg_pos:
                            if np.size(data_interp)==0:
                                data_interp = np.append(data_interp,data[k,:], 1)
                            else:
                                data_interp = np.append(data_interp,data[k,:], 0)
    
    # Defined horizontal ranges for the interpolation
    X_interp = np.linspace(xlim[0],xlim[1],10)
    Y_interp = np.linspace(ylim[0],ylim[1],10)
    
    # Vertical range for the interpolation
    if zlim[0]==zlim[1]:  # Constant altitude
        Z_interp = zlim[1]
    else:
        cond_Z1 = zlim[0] <= Z_tick
        cond_Z2 = Z_tick <= zlim[1]
        arg_min, arg_max = get_tick_argminmax(cond_Z1, cond_Z2)
        Z_interp = Z_tick[arg_min : arg_max+1]
    
    # Creating interpolators
    Ucalc = LinearNDInterpolator(data_interp[:,0:3],data_interp[:,3])
    Vcalc = LinearNDInterpolator(data_interp[:,0:3],data_interp[:,4])
    Wcalc = LinearNDInterpolator(data_interp[:,0:3],data_interp[:,5])
    Scalc = LinearNDInterpolator(data_interp[:,0:3],data_interp[:,6])
    
    # Creating mesh
    X_mesh, Y_mesh, Z_mesh = np.meshgrid(X_interp,Y_interp,Z_interp)
    
    # Interpolation
    Uinterp = Ucalc(X_mesh, Y_mesh, Z_mesh)
    Vinterp = Vcalc(X_mesh, Y_mesh, Z_mesh)
    Winterp = Wcalc(X_mesh, Y_mesh, Z_mesh)
    Sinterp = Scalc(X_mesh, Y_mesh, Z_mesh)

    return X_mesh, Y_mesh, Z_mesh, Uinterp, Vinterp, Winterp, Sinterp

def plot_wind_cube(xlim, ylim, zlim, plot):
    """
    
    Calculates a wind_cube inside the wind_cube and plots it if needed.

    Parameters
    ----------
    xlim : Narray of floats
        Numpy array of size 2 containing the limits of the x-axis range.
    ylim : Narray of floats
        Numpy array of size 2 containing the limits of the y-axis range.
    zlim : Narray of floats
        Numpy array of size 2 containing the limits of the z-axis range.
    plot : Bool
        Activates the plot or not.
        
    Returns
    -------
    X_mesh : Narray of floats
        3D-mesh for the interpolated x coordinates.
    Y_mesh : Narray of floats
        3D-mesh for the interpolated x coordinates.
    Z_mesh : Narray of floats
        3D-mesh for the interpolated x coordinates.
    Uinterp : Narray of floats
        3D-mesh for the interpolated wind speed component along x-axis.
    Vinterp : Narray of floats
        3D-mesh for the interpolated wind speed component along y-axis.
    Winterp : Narray of floats
        3D-mesh for the interpolated wind speed component along z-axis.
    Sinterp : Narray of floats
        3D-mesh for the interpolated surface altitude.

    """
    
    # Interpolation
    X_mesh, Y_mesh, Z_mesh, Uinterp, Vinterp, Winterp, Sinterp = get_interp_data(xlim,ylim,zlim)
    
    # Plotting the wind-cube if required
    if plot == True:
        fig = plt.figure(figsize=(16,12)) 
        ax = fig.gca(projection='3d') 
        ax.quiver(X_mesh, Y_mesh, Z_mesh+Sinterp, Uinterp, Vinterp, Winterp, length=max(X_mesh[0,:,0])/400) 
        ax.plot_surface(X_mesh[:,:,0],Y_mesh[:,:,0], Sinterp[:,:,0],color="#FBEEE6")
        plt.ion()
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')
        
    return X_mesh, Y_mesh, Z_mesh, Uinterp, Vinterp, Winterp, Sinterp

def plot_wind_surface(axis, coord, alt, plot):
    """
    
    Calculates a wind profile or a wind surface from the wind_cube.

    Parameters
    ----------
    axis : string
        Type of interpoltion to do. "z" for a wind surface, "x" or "y" for a
        wind profile.
    coord : Narry of float
        Coordinates of the point for the wind profile.
    alt : float
        Altitude for the wind surface.
    plot : Bool
        Activates the plot option.

    Returns
    -------
    X_mesh : Narray of floats
        3D-mesh for the interpolated x coordinates.
    Y_mesh : Narray of floats
        3D-mesh for the interpolated x coordinates.
    Z_mesh : Narray of floats
        3D-mesh for the interpolated x coordinates.
    Uinterp : Narray of floats
        3D-mesh for the interpolated wind speed component along x-axis.
    Vinterp : Narray of floats
        3D-mesh for the interpolated wind speed component along y-axis.
    Winterp : Narray of floats
        3D-mesh for the interpolated wind speed component along z-axis.
    Sinterp : Narray of floats
        3D-mesh for the interpolated surface altitude.

    """
    global Z_tick 
    global data
    # Wind profile
    if axis == "x" or axis == "y":
        #Decomposing the coordinates
        x = coord[0]
        y = coord[1]
        
        # Creating the limit arrays
        xlim = np.array([x,x])
        ylim = np.array([y,y])
        zlim = np.array([Z_tick[0], alt])
        
        # Interpolation
        X_mesh, Y_mesh, Z_mesh, Uinterp, Vinterp, Winterp, Sinterp = get_interp_data(xlim,ylim,zlim)
        
        # Resizing the result meshes
        Z_mesh =Z_mesh[0,0,:]
        Uinterp = Uinterp[0,0,:]
        Vinterp = Vinterp[0,0,:]
        Winterp = Winterp[0,0,:]
        
        # Calculating the norm
        Norm = np.sqrt(Uinterp*Uinterp + Vinterp*Vinterp + Winterp*Winterp)
        
        # Plot if required
        if plot == True:
            # U-component
            fig = plt.figure(figsize=(16,12))
            plt.plot(Uinterp,Z_mesh,'b-o',label='U')
            plt.legend(fontsize=16)
            plt.xlabel('Vitesse algébrique (m/s)')
            plt.ylabel('altitude (m)')
            plt.grid(which='both')
            plt.title('Wind profile at x=%6.2fm y=%6.2fm' %(x,y))
            
            # V component
            fig = plt.figure(figsize=(16,12))
            plt.plot(Vinterp,Z_mesh,'g-o',label='V')
            plt.legend(fontsize=16)
            plt.xlabel('Vitesse algébrique (m/s)')
            plt.ylabel('altitude (m)')
            plt.grid(which='both')
            plt.title('Wind profile at x=%6.2fm y=%6.2fm' %(x,y))
            
            # W component
            fig = plt.figure(figsize=(16,12))
            plt.plot(Winterp,Z_mesh,'r-o',label='W')
            plt.legend(fontsize=16)
            plt.xlabel('Vitesse algébrique (m/s)')
            plt.ylabel('altitude (m)')
            plt.grid(which='both')
            plt.title('Wind profile at x=%6.2fm y=%6.2fm' %(x,y))
            
            # Norm
            plt.figure(figsize=(16,12))
            plt.plot(Norm,Z_mesh,'m-o',label='W')
            plt.legend(fontsize=16)
            plt.xlabel('Norme (m/s)')
            plt.ylabel('altitude (m)')
            plt.grid(which='both')
            plt.title('Wind profile at x=%6.2fm y=%6.2fm' %(x,y))
            
    else:
        # Wind surface
        if axis == "z":
            # Creating the limit arrays
            xlim = np.array([X_tick[0], X_tick[-1]])
            ylim = np.array([Y_tick[0], Y_tick[-1]])
            
            Z_elev = np.array(data[:,2])
            data[:,2] += data[:,6]
            
            Z_tick = np.unique(data[:,2])
            #zlim = np.array([Z_tick[0], Z_tick[-1]])
            zlim = np.array([alt,alt])
            
            # Interpolation
            X_mesh, Y_mesh, Z_mesh, Uinterp, Vinterp, Winterp, Sinterp = get_interp_data(xlim,ylim,zlim)
            
            # Resizing the meshes
            X_mesh = X_mesh[:,:,0]
            Y_mesh = Y_mesh[:,:,0]
            Uinterp = Uinterp[:,:,0]
            Vinterp = Vinterp[:,:,0]
            Sinterp = Sinterp[:,:,0]
            
            # Calculating the norm
            M = np.hypot(Uinterp,Vinterp)
            
            # Plot if required
            if plot == True:
                
                plt.figure(figsize=(14,12))
                plt.quiver(X_mesh, Y_mesh, Uinterp, Vinterp, M, pivot='mid', units='xy')
                plt.colorbar()
                plt.scatter(X_mesh, Y_mesh, color='r', s=10)
                plt.title('Wind plot at elevation %im a.g.l (m/s)' %alt)
                plt.xlabel('x (m)')
                plt.ylabel('y (m)')
                plt.grid()
            
            data[:,2] = Z_elev
            Z_tick = interp.get_Zlist_pos(X_tick[0], Y_tick[0], data[:,0:3])[1]
            
    return X_mesh, Y_mesh, Z_mesh, Uinterp, Vinterp, Winterp

#%% Test
xlim = np.array([X_tick[0], X_tick[-1]])
ylim = np.array([Y_tick[0], Y_tick[-1]])
zlim = np.array([Z_tick[0], Z_tick[-1]])

X_mesh, Y_mesh, Z_mesh, Uinterp, Vinterp, Winterp, Sinterp = plot_wind_cube(xlim, ylim, zlim, True)

#X_mesh, Y_mesh, Z_mesh, Uinterp, Vinterp, Winterp = plot_wind_surface("z", [X_tick[0], Y_tick[0]], 1000, True)