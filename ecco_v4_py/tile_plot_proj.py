#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 16:11:15 2017

@author: ifenty
"""
from __future__ import division,print_function
import numpy as np
import matplotlib.pylab as plt
import matplotlib.path as mpath
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import ecco_v4_py.resample_to_latlon as resample_to_latlon

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
def plot_proj_to_latlon_grid(lons, lats, data, 
                             projection_type = 'robin', 
                             plot_type = 'pcolormesh', 
                             user_lon_0 = -66,
                             lat_lim = 50, 
                             levels = 20, 
                             cmap='jet', 
                             dx=.25, 
                             dy=.25,
                             show_colorbar = False, 
                             show_grid_lines = True,
                             show_grid_labels = True,
                             **kwargs):
    """Generate a plot of llc data, resampled to lat/lon grid, on specified 
    projection.

    Parameters
    ----------
    lons    : 
    lats    :
    data    : 

    """
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # default projection type = robinson
    # default central longitude = 60W
    # default no colorbar, no grid labels, no grid lines.
    # default color limits take the min and max of the values
    # default plot_type is pcolormesh.
    # default lat/lon spacing in lat/lon grid is 0.25 degrees
    # default number of levels for contourf is 20 (levels)
    # default latitude limit for polar stereographic plots is 50N (lat_lim)
    # default colormap is 'jet'

    #%%    
    cmin = np.nanmin(data)
    cmax = np.nanmax(data)
    subplot_grid = None

    for key in kwargs:
        if key == "cmin":
            cmin = kwargs[key]
        elif key == "cmax":
            cmax =  kwargs[key]
        elif key == "subplot_grid":
            subplot_grid = kwargs[key]
        else:
            print("unrecognized argument ", key)

    #%% 
    # set lat/lon bounds
    left_lon    = -180
    right_lon   =  180
    top_lat     =  89.5 
    bottom_lat  = -89.5 


    #%%
    # Make projection axis
    (ax,show_grid_labels) = _create_projection_axis(
            projection_type,user_lon_0,lat_lim,subplot_grid)
    

    #%%
    # do interpolation
    f = plt.gcf()
    new_grid_lon, new_grid_lat, data_latlon_projection = \
        resample_to_latlon(lons, lats, data, 
                           bottom_lat, top_lat, dy, 
                           left_lon, right_lon, dx,
                           mapping_method='nearest_neighbor')
        
    #%% 
    # make the plot
    if isinstance(ax.projection, ccrs.NorthPolarStereo) or \
       isinstance(ax.projection, ccrs.SouthPolarStereo) :
        p, gl, cbar = \
            plot_pstereo(new_grid_lon,
                         new_grid_lat, 
                         data_latlon_projection,
                         4326, lat_lim, 
                         cmin, cmax, ax,
                         plot_type = plot_type,
                         show_colorbar=False, 
                         circle_boundary=True,
                         cmap=cmap, 
                         show_grid_lines=False)

    else: # not polar stereo
        p, gl, cbar = \
            plot_global(new_grid_lon,
                        new_grid_lat, 
                        data_latlon_projection,
                        4326, 
                        cmin, cmax, ax,
                        plot_type = plot_type,                                       
                        show_colorbar = False,
                        cmap=cmap, 
                        show_grid_lines = False,
                        show_grid_labels = False)
                
    #%% 
    # add plot features
    if show_grid_lines :
        ax.gridlines(crs=ccrs.PlateCarree(), 
                              linewidth=1, color='black', 
                              alpha=0.5, linestyle='--', 
                              draw_labels = show_grid_labels)
    
    if show_colorbar:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(cmin,cmax))
        sm._A = []
        cbar = plt.colorbar(sm,ax=ax)        
    
    #%%
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.COASTLINE)

    ax= plt.gca()

    #%%
    return f, ax, p, cbar
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    

def plot_pstereo(xx,yy, data, 
                 data_projection_code, \
                 lat_lim, 
                 cmin, cmax, ax, 
                 plot_type = 'pcolormesh', 
                 show_colorbar=False, 
                 circle_boundary = False, 
                 cmap='jet', 
                 show_grid_lines=False,
                 levels = 20):

                            
    if isinstance(ax.projection, ccrs.NorthPolarStereo):
        ax.set_extent([-180, 180, lat_lim, 90], ccrs.PlateCarree())
    elif isinstance(ax.projection, ccrs.SouthPolarStereo):
        ax.set_extent([-180, 180, -90, lat_lim], ccrs.PlateCarree())
    else:
        raise ValueError('ax must be either ccrs.NorthPolarStereo or ccrs.SouthPolarStereo')

    if circle_boundary:
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)
        ax.set_boundary(circle, transform=ax.transAxes)

    if show_grid_lines :
        gl = ax.gridlines(crs=ccrs.PlateCarree(), 
                          linewidth=1, color='black', 
                          alpha=0.5, linestyle='--')
    else:
        gl = []

    if data_projection_code == 4326: # lat lon does nneed to be projected
        data_crs =  ccrs.PlateCarree()
    else:
        # reproject the data if necessary
        data_crs=ccrs.epsg(data_projection_code)
    

    p=[]    
    if plot_type == 'pcolormesh':
        p = ax.pcolormesh(xx, yy, data, transform=data_crs, \
                          vmin=cmin, vmax=cmax, cmap=cmap)

    elif plot_type =='contourf':
        p = ax.contourf(xx, yy, data, levels, transform=data_crs,  \
                 vmin=cmin, vmax=cmax, cmap=cmap)

    else:
        raise ValueError('plot_type  must be either "pcolormesh" or "contourf"')

         
    ax.add_feature(cfeature.LAND)
    ax.coastlines('110m', linewidth=0.8)

    cbar = []
    if show_colorbar:
        sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(cmin,cmax))
        sm._A = []
        cbar = plt.colorbar(sm,ax=ax)
    
    return p, gl, cbar

#%%    

def plot_global(xx,yy, data, 
                data_projection_code,
                cmin, cmax, ax, 
                plot_type = 'pcolormesh', 
                show_colorbar=False, 
                cmap='jet', 
                show_grid_lines = True,
                show_grid_labels = True,
                levels=20):

    if show_grid_lines :
        gl = ax.gridlines(crs=ccrs.PlateCarree(), 
                          linewidth=1, color='black', 
                          draw_labels = show_grid_labels,
                          alpha=0.5, linestyle='--')
    else:
        gl = []
        
    if data_projection_code == 4326: # lat lon does nneed to be projected
        data_crs =  ccrs.PlateCarree()
    else:
        data_crs =ccrs.epsg(data_projection_code)
        
    if plot_type == 'pcolormesh':
        p = ax.pcolormesh(xx, yy, data, transform=data_crs, 
                          vmin=cmin, vmax=cmax, cmap=cmap)
    elif plot_type =='contourf':
        p = ax.contourf(xx, yy, data, levels, transform=data_crs,
                        vmin=cmin, vmax=cmax, cmap=cmap)
    else:
        raise ValueError('plot_type  must be either "pcolormesh" or "contourf"') 
                         
    ax.coastlines('110m', linewidth=0.8)
    ax.add_feature(cfeature.LAND)

    cbar = []
    if show_colorbar:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(cmin,cmax))
        sm._A = []
        cbar = plt.colorbar(sm,ax=ax)
    
    return p, gl, cbar

# -----------------------------------------------------------------------------

def _create_projection_axis(projection_type,user_lon_0,lat_lim,subplot_grid):
    """Set appropriate axis for projection type.
        Note: optional subplot makes this messy, but has to happen here because
              axis projection cannot be "reset" once created.

    Parameters
    ----------
    projection_type :   string
                        user specified projection

    user_lon_0      :   double
                        center longitude value

    lat_lim         :   double
                        limiting latitude value

    subplot_grid    :   (optional) dict or list
                        specifying placement on subplot as 
                        dict: 
                        {'nrows': rows_val, 'ncols': cols_val, 'index': index_val}

                        or list: 
                        [nrows_val, ncols_val, index_val]

                        equates to

                        matplotlib.pyplot.subplot(
                            row=nrows_val, col=ncols_val,index=index_val)


    Returns
    -------
    ax              :   matplotlib axis object

    show_grid_labels:   logical
                        True = show the grid labels, only currently
                        supported for PlateCarree and Mercator projections
    """

    # initialize (optional) subplot variables
    row = []
    col = []
    ind = []

    if subplot_grid is not None:

        if type(subplot_grid) is dict:
            row = subplot_grid['nrows']
            col = subplot_grid['ncols']
            ind = subplot_grid['index']

        elif type(subplot_grid) is list:
            row = subplot_grid[0]
            col = subplot_grid[1]
            ind = subplot_grid[2]

        else:
            raise TypeError('Unexpected subplot_grid type: ',type(subplot_grid))

    if projection_type == 'Mercator':
        if subplot_grid is not None:
            ax = plt.subplot(row, col, ind, 
                    projection=ccrs.Mercator(central_longitude=user_lon_0))
        else:
            ax = plt.axes(projection=ccrs.Mercator(central_longitude=user_lon_0))

        show_grid_labels = True

    # ---

    elif projection_type == 'PlateCaree':
        if subplot_grid is not None:
            ax = plt.subplot(row, col, ind, 
                    projection=ccrs.PlateCarree(central_longitude=user_lon_0))
        else:
            ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=user_lon_0))
        show_grid_labels = True

    # ---

    elif projection_type == 'cyl':
        if subplot_grid is not None:
            ax = plt.subplot(row, col, ind, 
                    projection=ccrs.LambertCylindrical(central_longitude=user_lon_0))
        else:
            ax = plt.axes(projection=ccrs.LambertCylindrical(central_longitude=user_lon_0))
        show_grid_labels = False

    # ---

    elif projection_type == 'robin':    
        if subplot_grid is not None:
            ax = plt.subplot(row, col, ind, 
                    projection=ccrs.Robinson(central_longitude=user_lon_0))
        else:
            ax = plt.axes(projection=ccrs.Robinson(central_longitude=user_lon_0))
        show_grid_labels = False

    # ---

    elif projection_type == 'ortho':
        if subplot_grid is not None:
            ax = plt.subplot(row, col, ind, 
                    projection=ccrs.Orthographic(central_longitude=user_lon_0))
        else:
            ax = plt.axes(projection=ccrs.Orthographic(central_longitude=user_lon_0))
        show_grid_labels = False

    # ---

    elif projection_type == 'stereo':    
        if lat_lim > 0:
            stereo_proj = ccrs.NorthPolarStereo()
        else:
            stereo_proj = ccrs.SouthPolarStereo()

        if subplot_grid is not None:
            ax = plt.subplot(row, col, ind, 
                    projection=stereo_proj)
        else:
            ax = plt.axes(projection=stereo_proj)

        show_grid_labels = False

    # ---

    elif projection_type == 'InterruptedGoodeHomolosine':
        if subplot_grid is not None:
            ax = plt.subplot(row, col, ind, 
                    projection=ccrs.InterruptedGoodeHomolosine(central_longitude=user_lon_0))
        else:
            ax = plt.axes(projection=ccrs.InterruptedGoodeHomolosine(central_longitude=user_lon_0))
        show_grid_labels = False

    # ---
        
    else:
        raise ValueError('projection type must be either "Mercator", "PlateCaree",  "cyl", "robin", "ortho", or "stereo"')

    #print ('projection type ', projection_type)
    return (ax,show_grid_labels)
