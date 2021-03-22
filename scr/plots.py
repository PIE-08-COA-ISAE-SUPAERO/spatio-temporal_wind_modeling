# -*- coding: utf-8 -*-
"""Plots.

This module contains all the functions required to plot the wind cube data

@author: Q.Abeille (2020)
"""

#%% Imports
import numpy as np
import matplotlib.pyplot as plt
#plt.style.use('presentation.mplstyle')
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
from scipy.interpolate import LinearNDInterpolator
import extrapolation as interp

#%% Wind cube
def get_surf(X_range, Y_range):
    """
    Create the arrays for the surface plot.
    
    Parameters
    ----------
    X_range : Narray of float
        Range along x_axis for the surface.
    Y_range : Narray of float
        Range along z axis for the surface.

    Returns
    -------
    X_surf : Narray of float
        Meshgrid for x-axis.
    Y_surf : Narray of float
        Meshgrid for y-axis.
    Surf : Narray of float
        Surface altitude grid.

    """
    Surf = np.zeros((len(Y_range),len(X_range)))
    X_surf, Y_surf = np.meshgrid(X_range,Y_range)
    
    for i in range(len(X_range)):
        for j in range(len(Y_range)):
            zpos = np.min(interp.get_Zlist_pos(X_range[i], Y_range[j], np.concatenate((data[:,0:2],np.reshape(Zsurf,(len(Zsurf),1))),axis=1))[1])
            Surf[j,i] = zpos;
    return X_surf, Y_surf, Surf

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
    
    # Checking if limits inside the wind_cube boundaries
    if x_min < X_tick[0]:
        raise ValueError("Lower x limit out of bounds")
    
    if x_max > X_tick[-1]:
        raise ValueError("Upper x limit out of bounds")
    
    if y_min < Y_tick[0]:
        raise ValueError("Lower y limit out of bounds")
    
    if y_max > Y_tick[-1]:
        raise ValueError("Upper y limit out of bounds")
    
    if z_min < Z_tick[0]:
        raise ValueError("Lower z limit out of bounds")
    
    if z_max > Z_tick[-1]:
        raise ValueError("Upper z limit out of bounds")
    
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


def get_interp_data(xlim, ylim, zlim, nb_points, surf):
    """
    
    Calculate an interpolation from the wind_cube data to a defined domain.

    Parameters
    ----------
    xlim : Narray of floats
        Numpy array of size 2 containing the limits of the x-axis range.
    ylim : Narray of floats
        Numpy array of size 2 containing the limits of the y-axis range.
    zlim : Narray of floats
        Numpy array of size 2 containing the limits of the z-axis range.
    nb_points : int
        Number of output points along x and y axes
    surf : Boolean
        True if calculation of a wind_surface

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
    # Case of a wind surface
    if surf == True:
        z_elev_range = np.ones(len(data[:,6]))*zlim[0] - data[:,6]  # from alt to elevation range
        zmin = max(0,min(z_elev_range))
        zmax = min(Z_tick[-1],max(z_elev_range))
        zlim = np.array([zmin,zmax])
        
    # Calculating the ranges along each axis
    X_range, Y_range, Z_range = get_interv(xlim,ylim,zlim)
    
    nb_points_required = len(X_range) * len(Y_range) * len(Z_range)
    nb_points_cube = len(X_tick) * len(Y_tick) * len(Z_tick)
    
    if nb_points_required > nb_points_cube/2: 
        # More than half of the cube required, faster to use the whole domain
        data_interp = data         
    else:
        data_interp = [[]]
        
        # Retrieving the lines of data that will be used in the interpolation
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
    X_interp = np.linspace(xlim[0],xlim[1],nb_points)
    Y_interp = np.linspace(ylim[0],ylim[1],nb_points)
    
    # Vertical range for the interpolation
    if surf==True:
        Z_interp = np.unique(z_elev_range)  # To have the wind at the right altitude
    else:
        cond_Z1 = zlim[0] <= Z_tick
        cond_Z2 = Z_tick <= zlim[1]
        arg_min, arg_max = get_tick_argminmax(cond_Z1, cond_Z2)
        Z_interp = Z_tick[arg_min : arg_max+2]
    
    # Creating interpolators
    Ucalc = LinearNDInterpolator(data_interp[:,0:3],data_interp[:,3])
    Vcalc = LinearNDInterpolator(data_interp[:,0:3],data_interp[:,4])
    Wcalc = LinearNDInterpolator(data_interp[:,0:3],data_interp[:,5])
    Scalc = LinearNDInterpolator(data_interp[:,0:2],data_interp[:,6])
    
    # Creating mesh
    X_mesh, Y_mesh, Z_mesh = np.meshgrid(X_interp,Y_interp,Z_interp)
    
    # Interpolation
    Uinterp = Ucalc(X_mesh, Y_mesh, Z_mesh)
    Vinterp = Vcalc(X_mesh, Y_mesh, Z_mesh)
    Winterp = Wcalc(X_mesh, Y_mesh, Z_mesh)
    Sinterp = Scalc(X_mesh, Y_mesh)

    return X_mesh, Y_mesh, Z_mesh, Uinterp, Vinterp, Winterp, Sinterp

def plot_wind_cube(wind_cube, xlim, ylim, zlim, nb_points, plot):
    """
    
    Calculate a wind_cube inside the wind_cube and plots it if needed.

    Parameters
    ----------
    wind_cube : wind_cube
        The wind cube on which the interpolation will be done.
    xlim : Narray of floats
        Numpy array of size 2 containing the limits of the x-axis range.
    ylim : Narray of floats
        Numpy array of size 2 containing the limits of the y-axis range.
    zlim : Narray of floats
        Numpy array of size 2 containing the limits of the z-axis range. In elevation convention.
    nb_points : int
        Number of output points along x and y axes. Default = 10.
    plot : Bool
        Activates the plot option. Default = False.
        
    Returns
    -------
    X_mesh : Narray of floats
        3D-mesh for the interpolated x coordinates.
    Y_mesh : Narray of floats
        3D-mesh for the interpolated x coordinates.
    Z_mesh : Narray of floats
        3D-mesh for the interpolated x coordinates. In elevation convention.
    Uinterp : Narray of floats
        3D-mesh for the interpolated wind speed component along x-axis.
    Vinterp : Narray of floats
        3D-mesh for the interpolated wind speed component along y-axis.
    Winterp : Narray of floats
        3D-mesh for the interpolated wind speed component along z-axis.
    Sinterp : Narray of floats
        3D-mesh for the interpolated surface altitude.

    """
    global data_points
    data_points = wind_cube["Position"]
    global data_wind
    data_wind = wind_cube["Wind_speed"]
    global Zsurf
    Zsurf = wind_cube["Surface_altitude"].reshape(-1,1)
    global data
    data = np.concatenate((data_points, data_wind, Zsurf),axis=1)

    # Retrieving the horizontal mesh
    global X
    X = data[:,0]
    global Y
    Y = data[:,1]
    global X_tick
    X_tick = np.unique(X)
    global Y_tick
    Y_tick = np.unique(Y)
    global Z_tick
    Z_tick = interp.get_Zlist_pos(X_tick[0], Y_tick[0], data[:,0:3])[1]
    
    # Interpolation
    X_mesh, Y_mesh, Z_mesh, Uinterp, Vinterp, Winterp, Sinterp = get_interp_data(xlim,ylim,zlim,nb_points,False)
    
    # Plotting the wind-cube if required
    if plot == True:

        X_range, Y_range, Z_range = get_interv(xlim, ylim, zlim)
        
        X_surf, Y_surf, Surf = get_surf(X_range, Y_range)
        
        fig = plt.figure(figsize=(16,10)) 
        ax = fig.gca(projection='3d')
        ax.plot_surface(X_surf,Y_surf, Surf, cmap=cm.terrain)
        ax.quiver(X_mesh, Y_mesh, Z_mesh+Sinterp, Uinterp, Vinterp, Winterp, colors='#7F8C8D', length=max(X_mesh[0,:,0])/400)
        ax.set_zlim([max(0,min(Sinterp[0,0,:])),max(2000+min((Z_mesh+Sinterp)[0,0,:]),max((Z_mesh+Sinterp)[0,0,:]))])
        plt.ion()
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')
        plt.show()
        
    return X_mesh, Y_mesh, Z_mesh, Uinterp, Vinterp, Winterp, Sinterp

def plot_wind_surface(wind_cube, axis, coord, alt, nb_points, plot):
    """
    
    Calculate a wind profile or a wind surface from the wind_cube.

    Parameters
    ----------
    wind_cube : wind_cube
        The wind cube on which the interpolation will be done.
    axis : string
        Type of interpoltion to do. "z" for a wind surface, "x" or "y" for a
        wind profile.
    coord : Narray of float
        Coordinates of the point for the wind profile. In GPS coordinates.
    alt : float
        Altitude for the wind surface. Must be the altitude above sea level.
        If wind profile required, elevation above ground max for the plot. 
    nb_points : int
        Number of output points along x and y axes. Default = 10.
    plot : Bool
        Activates the plot option. Default = False.

    Returns
    -------
    X_mesh : Narray of floats
        Interpolated x coordinates. 1D list.
    Y_mesh : Narray of floats
        Interpolated y coordinates. 1D list
    Z_mesh : Narray of floats
        Interpolated z coordinates. 1D list in elevation convention for 
        the wind_profile and altitude convention for the wind_surface.
    Uinterp : Narray of floats
        1D array for the interpolated wind speed component along x-axis.
    Vinterp : Narray of floats
        1D array for the interpolated wind speed component along y-axis.
    Winterp : Narray of floats
        1D array for the interpolated wind speed component along z-axis.
    Sinterp : Narray of floats
        1D array for the interpolated surface altitude.

    """
    global data_points
    data_points = wind_cube["Position"]
    global data_wind 
    data_wind = wind_cube["Wind_speed"]
    global Zsurf
    Zsurf = wind_cube["Surface_altitude"].reshape(-1,1)
    global data 
    data = np.concatenate((data_points, data_wind, Zsurf),axis=1)

    # Retrieving the horizontal mesh
    global X
    X = data[:,0]
    global Y
    Y = data[:,1]
    global X_tick
    X_tick = np.unique(X)
    global Y_tick
    Y_tick = np.unique(Y)
    global Z_tick
    Z_tick = interp.get_Zlist_pos(X_tick[0], Y_tick[0], data[:,0:3])[1]
    
    # Wind profile --> interpolation in elevation convention
    if axis == "x" or axis == "y":
        #Decomposing the coordinates
        x = coord[0]
        y = coord[1]
        
        # Creating the limit arrays
        xlim = np.array([x,x])
        ylim = np.array([y,y])
        zlim = np.array([Z_tick[0], alt])
        
        # Interpolation
        X_mesh, Y_mesh, Z_mesh, Uinterp, Vinterp, Winterp, Sinterp = get_interp_data(xlim, ylim, zlim, nb_points, False)
        
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
            plt.ylabel('elevation a.g.l (m)')
            plt.grid(which='both')
            plt.title('Wind profile at x=%6.2fm y=%6.2fm' %(x,y))
            
            # V component
            fig = plt.figure(figsize=(16,12))
            plt.plot(Vinterp,Z_mesh,'g-o',label='V')
            plt.legend(fontsize=16)
            plt.xlabel('Vitesse algébrique (m/s)')
            plt.ylabel('elevation a.g.l (m)')
            plt.grid(which='both')
            plt.title('Wind profile at x=%6.2fm y=%6.2fm' %(x,y))
            
            # W component
            fig = plt.figure(figsize=(16,12))
            plt.plot(Winterp,Z_mesh,'r-o',label='W')
            plt.legend(fontsize=16)
            plt.xlabel('Vitesse algébrique (m/s)')
            plt.ylabel('elevation a.g.l (m)')
            plt.grid(which='both')
            plt.title('Wind profile at x=%6.2fm y=%6.2fm' %(x,y))
            
            # Norm
            plt.figure(figsize=(16,12))
            plt.plot(Norm,Z_mesh,'m-o',label='W')
            plt.legend(fontsize=16)
            plt.xlabel('Norme (m/s)')
            plt.ylabel('elevation a.g.l (m)')
            plt.grid(which='both')
            plt.title('Wind profile at x=%6.2fm y=%6.2fm' %(x,y))

            plt.show()
            
    else:
        # Wind surface
        if axis == "z":
            # Checking altitude inside bounds
            if alt < min(Zsurf):
                raise ValueError("Altitude below ground surface")
            if alt > (max(Zsurf) + max(Z_tick)):
                raise ValueError("Altitude exceed the maximum altitude of the domain")
            if alt > (min(Zsurf) + max(Z_tick)):
                raise Warning("Altitude exceed the maximum calculated elevation over the lowest point of the domain. Some points might not be interpolated")
            # Creating the limit arrays
            xlim = np.array([X_tick[0], X_tick[-1]])
            ylim = np.array([Y_tick[0], Y_tick[-1]])
            zlim = np.array([alt,alt])
            
            # Interpolation
            X_mesh, Y_mesh, Z_mesh, Uinterp, Vinterp, Winterp, Sinterp = get_interp_data(xlim,ylim,zlim, nb_points, True)
            
            # Retrieving the position of the points at the required altitude
            select = np.argwhere(np.round(Sinterp+Z_mesh)==alt)
            
            X = np.zeros(len(select))
            Y = np.zeros(len(select))
            U = np.zeros(len(select))
            V = np.zeros(len(select))
            S = np.zeros(len(select))
            
            # Building the surface
            for i in range(len(select)-1):
                xarg = select[i,0]
                yarg = select[i,1]
                zarg = select[i,2]
                
                X[i] = X_mesh[xarg,yarg,zarg]
                Y[i] = Y_mesh[xarg,yarg,zarg]
                U[i] = Uinterp[xarg,yarg,zarg]
                V[i] = Vinterp[xarg,yarg,zarg]
                S[i] = Sinterp[xarg,yarg,zarg]
            
            # Calculating the norm
            M = np.hypot(U,V)
            
            # Plot if required
            if plot == True:
                
                X_surf, Y_surf, Surf = get_surf(X_tick, Y_tick)

                plt.figure(figsize=(14,12))
                # wind surface
                plt.quiver(X, Y, U, V, M, pivot='mid', units='xy')
                plt.colorbar()
                #iso-altitude contours
                plt.scatter(X, Y, color='r', s=5)
                CS = plt.contour(X_surf, Y_surf, Surf, colors='black')
                plt.clabel(CS, inline=1, fontsize=10, fmt='%1.f')
                plt.title('Wind plot at altitude %im a.s.l (m/s)' %alt)
                plt.xlabel('x (m)')
                plt.ylabel('y (m)')
                plt.grid()
                plt.show()
    return X_mesh, Y_mesh, Z_mesh, Uinterp, Vinterp, Winterp, Sinterp



def spectre_prin_p(frequency, speed, altitude):
    
    """
    Return the DSP of the principal wind speed
    Parameters
    ----------
    frequency : double
         The frequency of the component expected.
    speed : double
         The stationnary speed of the wind at the point wanted.
    altitude : double
         The altitude of the point wanted.  
    Returns
    -------
    S : double 
         The value of the DSP for these parameters
    """
    
    lu = altitude/((A_Z+B_Z*altitude)**(C_Z))
    S = 4*lu/speed * 1/((1 + D_PRIN*(frequency*lu/speed)**2)**(PENTE_PRIN))
    
    return(S)


def spectre_late_p(frequency, speed, altitude):
    
    """
    Return the DSP of the lateral wind speed
    Parameters
    ----------
    frequency : double
         The frequency of the component expected.
    speed : double
         The stationnary speed of the wind at the point wanted.
    altitude : double
         The altitude of the point wanted.  
    Returns
    -------
    S : double 
         The value of the DSP for these parameters
    """
    
    lv = altitude/((A_Z+B_Z*altitude)**(C_Z))
    S = 4*lv/speed * (1 + D_LATW*(2*frequency*lv/speed)**2)/((1 + D_PRIN*(2*frequency*lv/speed)**2)**(PENTE_LATW))
    
    return(S)


def spectre_w_p(frequency, speed, altitude):
    
    """
    Return the DSP of the vertical wind speed
    Parameters
    ----------
    frequency : double
         The frequency of the component expected.
    speed : double
         The stationnary speed of the wind at the point wanted.
    altitude : double
         The altitude of the point wanted.  
    Returns
    -------
    S : double 
         The value of the DSP for these parameters
    """
    
    lw = altitude
    S = 4*lw/speed * (1 + D_LATW*(2*frequency*lw/speed)**2)/((1 + D_PRIN*(2*frequency*lw/speed)**2)**(PENTE_LATW))
    
    return(S)


def plot_wind_cube_turbulent(wind_cube, xlim, ylim, zlim, T, dt, nb_points,  plot):
    """
    
    Calculate a wind_cube with turbulence inside the wind_cube and plots it if 
    needed (if the number of iteration is reasonnable).

    Parameters
    ----------
    wind_cube : wind_cube
        The wind cube on which the interpolation will be done.
    xlim : Narray of floats
        Numpy array of size 2 containing the limits of the x-axis range.
        It must be in meters from the bottom left corner.
    ylim : Narray of floats
        Numpy array of size 2 containing the limits of the y-axis range.
        It must be in meters from the bottom left corner.
    zlim : Narray of floats
        Numpy array of size 2 containing the limits of the z-axis range. In elevation convention.
    nb_points : int
        Number of output points along x and y axes. Default = 10.
    T : float
        duration of observation of the wind speed.
    dt : float
        timestep of the computation
    plot : Bool
        Activates the plot option. Default = False.
        
    Returns
    -------
    cube_list : list of 6 Narray containing the evolution of the wind speed in time :
        X_mesh : Narray of floats
            3D-mesh for the interpolated x coordinates.
        Y_mesh : Narray of floats
            3D-mesh for the interpolated x coordinates.
        Z_mesh : Narray of floats
            3D-mesh for the interpolated x coordinates. In elevation convention.
        Uinterp : Narray of floats
            3D-mesh for the interpolated wind speed component along x-axis.
        Vinterp : Narray of floats
            3D-mesh for the interpolated wind speed component along y-axis.
        Winterp : Narray of floats
            3D-mesh for the interpolated wind speed component along z-axis.

    """
    global data_points
    data_points = wind_cube["Position"]
    global data_wind
    data_wind = wind_cube["Wind_speed"]
    global Zsurf
    Zsurf = wind_cube["Surface_altitude"].reshape(-1,1)
    global data
    data = np.concatenate((data_points, data_wind, Zsurf),axis=1)

    #Récupération du maillage horizontal
    
    global X
    X = data[:,0]
    global Y
    Y = data[:,1]
    global X_tick
    X_tick = np.unique(X)
    global Y_tick
    Y_tick = np.unique(Y)
    global Z_tick
    Z_tick = interp.get_Zlist_pos(X_tick[0], Y_tick[0], data[:,0:3])[1]
    
    # Interpolation
    
    X_mesh, Y_mesh, Z_mesh, Uinterp, Vinterp, Winterp, Sinterp = get_interp_data(xlim,ylim,zlim,nb_points,False)
    
    # Preparation of the values for turbulence
    
    fs = 10
    N = 120
    Ntheta = 25
    df, dtheta = fs/2 / N, np.pi / Ntheta
    nb_iter = int(T/dt) + 1
    if nb_iter > 15 :
        print("Too much iteration ! No images will be generated")
        plot = False
    
    # Getting the mean and direction on each elevation
        
    cube_list = [0 for k in range(nb_iter)]
    list_Uz = []
    list_dz = []
    tableau_Sprin = []
    tableau_Slate = []
    tableau_Svert = []
    
    for z in range(np.shape(Z_mesh)[0]) :
        hauteur = Z_mesh[0,0,z]
        Umoy = 0
        dmoy = 0
        count = 0
        for x in range(np.shape(X_mesh)[0]) :
            for y in range(np.shape(Y_mesh)[0]) :
                u, v = Uinterp[x,y,z], Vinterp[x,y,z]
                norme = np.sqrt(u**2 + v**2)
                direction = np.arccos(abs(u)/norme)
                if v < 0 :
                    direction = - direction
                count = count + 1
                Umoy = Umoy + norme
                dmoy = dmoy + direction
        Umoy = Umoy / count
        dmoy = dmoy / count
        list_Uz.append(Umoy)
        list_dz.append(dmoy)
        
        altitude_correction = max((CORREC_ALTI - 1)/(1000 - 0)*hauteur + 1, CORREC_ALTI)
        s_base = Umoy * altitude_correction * CORREC_GLOBAL
        s_prin = RATIO_PRIN * s_base
        s_late = RATIO_LAT * s_base
        s_w = RATIO_VERT * s_base
        X_prin, X_late, X_vert = [], [], []
        
        for k in range(N):
            frequency = k/(2*N) * fs
            s_prink = np.sqrt(TIME_MAX/(2*np.pi) * s_prin**2 * spectre_prin_p(frequency,\
                                                        Umoy, hauteur))
            s_latek = np.sqrt(TIME_MAX/(2*np.pi) * s_late**2 * spectre_late_p(frequency,\
                                                        Umoy, hauteur))
            s_wk = np.sqrt(TIME_MAX/(2*np.pi) * s_w**2 * spectre_w_p(frequency, \
                                                        Umoy, hauteur))
            X_prin.append(s_prink)
            X_late.append(s_latek)
            X_vert.append(s_wk)
        
        tableau_Sprin.append(X_prin)
        tableau_Slate.append(X_late)
        tableau_Svert.append(X_vert)
    
    Psi = rd.uniform(0, 2*np.pi, N * Ntheta + 1)
    
    # Computing for each timestep the wind cube
    
    for t in range(1, nb_iter + 1):
        
        Ut = Uinterp.copy()
        Vt = Vinterp.copy()
        Wt = Winterp.copy()
        
        for x in range(np.shape(X_mesh)[0]) :
            for y in range(np.shape(Y_mesh)[0]) :
                for z in range(np.shape(Z_mesh)[0]) :
                    coordx = X_mesh[x,y,z]
                    coordy = Y_mesh[x,y,z]
                    coordz = Z_mesh[x,y,z]
                    u = Uinterp[x, y, z]
                    v = Vinterp[x, y, z]
                    w = Winterp[x, y, z]
                    
                    coordx_prin = np.cos(list_dz[z]) * coordx - np.sin(list_dz[z]) * coordy
                    coordx_late = np.sin(list_dz[z]) * coordx + np.cos(list_dz[z]) * coordy
                    uprin, ulate, uvert = 0, 0, 0
                    windmoy = list_Uz[z]
                    for k in range(N) :
                        for kk in range(Ntheta):
                            uprin = uprin + np.sqrt(4/np.pi * tableau_Sprin[z][k] * 2*np.pi* df * dtheta * (np.cos(-np.pi/2 + kk * dtheta))**2) * np.cos(2*np.pi*k*df/windmoy * coordx_prin * np.cos(-np.pi/2 + kk * dtheta) + 2*np.pi*k*df/windmoy * coordx_late * np.sin(-np.pi/2 + kk * dtheta) - 2*np.pi*k*df * t*dt + Psi[k*Ntheta + kk]) 
                            ulate = ulate + np.sqrt(4/np.pi * tableau_Slate[z][k] * 2*np.pi* df * dtheta * (np.cos(-np.pi/2 + kk * dtheta))**2) * np.cos(2*np.pi*k*df/windmoy * coordx_prin * np.cos(-np.pi/2 + kk * dtheta) + 2*np.pi*k*df/windmoy * coordx_late * np.sin(-np.pi/2 + kk * dtheta) - 2*np.pi*k*df * t*dt + Psi[k*Ntheta + kk]) 
                            uvert = uvert + np.sqrt(4/np.pi * tableau_Svert[z][k] * 2*np.pi* df * dtheta * (np.cos(-np.pi/2 + kk * dtheta))**2) * np.cos(2*np.pi*k*df/windmoy * coordx_prin * np.cos(-np.pi/2 + kk * dtheta) + 2*np.pi*k*df/windmoy * coordx_late * np.sin(-np.pi/2 + kk * dtheta) - 2*np.pi*k*df * t*dt + Psi[k*Ntheta + kk]) 
                    
                    unew = u + np.cos(list_dz[z]) * uprin - np.sin(list_dz[z]) * ulate
                    vnew = v + np.sin(list_dz[z]) * uprin + np.cos(list_dz[z]) * ulate
                    wnew = w + uvert
                    
                    Ut[x, y, z] = unew
                    Vt[x, y, z] = vnew
                    Wt[x, y, z] = wnew
        
        cube_list.append([X_mesh, Y_mesh, Z_mesh, Ut, Vt, Wt])
        
        if plot == True:
            X_range, Y_range, Z_range = get_interv(xlim, ylim, zlim)
            X_surf, Y_surf, Surf = get_surf(X_range, Y_range)
            fig = plt.figure(figsize=(16,10)) 
            ax = fig.gca(projection="3d")
            ax.plot_surface(X_surf,Y_surf, Surf, cmap=cm.terrain)
            ax.quiver(X_mesh, Y_mesh, Z_mesh+Sinterp, Uinterp, Vinterp, Winterp, color = "limegreen", length=max(X_mesh[0,:,0])/400, label = "Steady")
            ax.quiver(X_mesh, Y_mesh, Z_mesh+Sinterp, Ut, Vt, Wt, color = "red", length=max(X_mesh[0,:,0])/400, label = "Turbulent")
            ax.set_zlim([max(0,min(Sinterp[0,0,:])),max(2000+min((Z_mesh+Sinterp)[0,0,:]),max((Z_mesh+Sinterp)[0,0,:]))])
            ax.legend() 
            plt.ion()
            plt.xlabel("x (m)")
            plt.ylabel("y (m)")
            plt.show()
    
    return(cube_list)
