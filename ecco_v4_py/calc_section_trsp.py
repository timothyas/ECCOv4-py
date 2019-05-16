"""
Compute transport (freshwater, heat, salt) across a section, e.g. Drake Passage
"""

import numpy as np
import xarray as xr
# xarray compatibility
try:
    from xarray.core.pycompat import OrderedDict
except ImportError:
    from collections import OrderedDict

from .ecco_utils import get_llc_grid
from .get_section_masks import get_available_sections, \
        get_section_endpoints, get_section_line_masks

# -------------------------------------------------------------------------------
# Main functions for computing standard transport quantities
# -------------------------------------------------------------------------------

# Define constants
METERS_CUBED_TO_SVERDRUPS = 10**-6
WATTS_TO_PETAWATTS = 10**-15
RHO_CONST = 1000
HEAT_CAPACITY = 4000

def calc_section_vol_trsp(ds,
                          pt1=None, pt2=None,
                          section_name=None,
                          maskW=None, maskS=None, 
                          grid=None):
    """Compute volumetric transport across section in Sverdrups
    There are 3 ways to call this function:

    1. Provide pre-defined section_name, e.g.

        >> trsp = calc_section_vol_trsp(ds,'Drake Passage')

            * Computes volumetric trsp across predefined Drake Passage line
            * See get_available_sections for available definitions

    2. Provide lat/lon pairs to compute transport across, e.g. 
        
        >> pt1 = [lon1, lat1]
        >> pt2 = [lon2, lat2]
        >> trsp = calc_section_vol_trsp(ds,pt1,pt2)

            * Computes volumetric transport across a band between pt1 -> pt2
            * If section name is provided, it gets added to returned DataArray

    3. Provide maskW, maskS, e.g.

        >> _, maskW, maskS = get_section_line_masks(pt1, pt2, ds)
        >> trsp = calc_section_vol_trsp(ds,maskW,maskS)

            * Compute trsp across band defined by masks    
            * If section name is provided, it gets added to returned DataArray

    Parameters
    ----------
    ds : xarray Dataset
        must contain UVELMASS,VVELMASS, drF, dyG, dxG
    pt1, pt2 : list or tuple with two floats, optional
        end points for section line as [lon lat] or (lon, lat)
    maskW, maskS : xarray DataArray, optional
        masks denoting the section, created by get_section_line_masks
    section_name: string, optional
        name for the section. If predefined value, section mask is defined 
        via get_section_endpoints
        otherwise, adds name to returned DataArray
    grid : xgcm Grid object, optional
        denotes LLC90 operations for xgcm, see ecco_utils.get_llc_grid
        see also the [xgcm documentation](https://xgcm.readthedocs.io/en/latest/grid_topology.html)

    Returns
    -------
    vol_trsp_ds : xarray Dataset
        includes variables as xarray DataArrays
            vol_trsp
                freshwater transport across section in Sv
                with dimensions 'time' (if in given dataset) and 'lat' 
            maskW, maskS
                defining the section
        and the section_name as an attribute if it is provided
    """

    maskW, maskS = _parse_section_trsp_inputs(ds,pt1,pt2,maskW,maskS,section_name)

    # Define volumetric transport
    x_vol = ds['UVELMASS'] * ds['drF'] * ds['dyG'] 
    y_vol = ds['VVELMASS'] * ds['drF'] * ds['dxG'] 

    # Computes salt transport in m^3/s at each depth level
    vol_trsp = section_trsp_at_depth(x_vol,y_vol,maskW,maskS,
                                     cds=ds.coords.to_dataset(),
                                     grid=grid)
    # Sum over depth
    vol_trsp = vol_trsp.sum('k')

    # Convert to Sv
    vol_trsp = METERS_CUBED_TO_SVERDRUPS * vol_trsp
    vol_trsp.attrs['units'] = 'Sv'

    # Add this to a Dataset
    vol_trsp_ds = vol_trsp.to_dataset(name='vol_trsp')

    # Add section name and masks to Dataset
    vol_trsp_ds['maskW'] = maskW
    vol_trsp_ds['maskS'] = maskS
    if section_name is not None:
        vol_trsp_ds.attrs['name'] = section_name

    return vol_trsp_ds

def calc_section_heat_trsp(ds,
                           pt1=None, pt2=None,
                           section_name=None,
                           maskW=None, maskS=None, 
                           grid=None):
    """Compute heat transport across section in PW
    Inputs and usage are same as calc_section_vol_trsp. 
    The only differences are:

    Parameters
    ----------
    ds : xarray Dataset
        must contain ADVx_TH, ADVy_TH, DFxe_TH, DFyE_TH

    Returns
    -------
    heat_trsp_ds : xarray Dataset
        includes variables as xarray DataArrays
            heat_trsp
                heat transport across section in PW
                with dimensions 'time' (if in given dataset) and 'lat' 
            maskW, maskS
                defining the section
        and the section_name as an attribute if it is provided
    """

    maskW, maskS = _parse_section_trsp_inputs(ds,pt1,pt2,maskW,maskS,section_name)

    # Define heat transport
    x_heat = ds['ADVx_TH'] * ds['DFxE_TH']
    y_heat = ds['ADVy_TH'] * ds['DFyE_TH']

    # Computes salt transport in degC * m^3/s at each depth level
    heat_trsp = section_trsp_at_depth(x_heat,y_heat,maskW,maskS,
                                      cds=ds.coords.to_dataset(),
                                      grid=grid)
    # Sum over depth
    heat_trsp = heat_trsp.sum('k')

    # Convert to PW
    heat_trsp = WATTS_TO_PETAWATTS * RHO_CONST * HEAT_CAPACITY * heat_trsp
    heat_trsp.attrs['units'] = 'PW'

    # Add this to a Dataset
    heat_trsp_ds = heat_trsp.to_dataset(name='heat_trsp')

    # Add section name and masks to Dataset
    heat_trsp_ds['maskW'] = maskW
    heat_trsp_ds['maskS'] = maskS
    if section_name is not None:
        heat_trsp_ds.attrs['name'] = section_name

    return heat_trsp_ds

def calc_section_salt_trsp(ds,
                           pt1=None, pt2=None,
                           section_name=None,
                           maskW=None, maskS=None, 
                           grid=None):
    """Compute salt transport across section in psu*Sv
    Inputs and usage are same as calc_section_vol_trsp. 
    The only differences are:

    Parameters
    ----------
    ds : xarray Dataset
        must contain ADVx_SLT, ADVy_SLT, DFxe_SLT, DFyE_SLT

    Returns
    -------
    salt_trsp_ds : xarray Dataset
        includes variables as xarray DataArrays
            salt_trsp
                salt transport across section in psu*Sv
                with dimensions 'time' (if in given dataset) and 'lat' 
            maskW, maskS
                defining the section
        and the section_name as an attribute if it is provided
    """

    maskW, maskS = _parse_section_trsp_inputs(ds,pt1,pt2,maskW,maskS,section_name)

    # Define salt transport
    x_salt = ds['ADVx_SLT'] * ds['DFxE_SLT']
    y_salt = ds['ADVy_SLT'] * ds['DFyE_SLT']

    # Computes salt transport in psu * m^3/s at each depth level
    salt_trsp = section_trsp_at_depth(x_salt,y_salt,maskW,maskS,
                                      cds=ds.coords.to_dataset(),
                                      grid=grid)
    # Sum over depth
    salt_trsp = salt_trsp.sum('k')

    # Convert to PW
    salt_trsp = METERS_CUBED_TO_SVERDRUPS * salt_trsp
    salt_trsp.attrs['units'] = 'psu.Sv'

    # Add this to a Dataset
    salt_trsp_ds = salt_trsp.to_dataset(name='salt_trsp')

    # Add section name and masks to Dataset
    salt_trsp_ds['maskW'] = maskW
    salt_trsp_ds['maskS'] = maskS
    if section_name is not None:
        salt_trsp_ds.attrs['name'] = section_name

    return salt_trsp_ds

# -------------------------------------------------------------------------------
# Main function for computing standard transport quantities
# -------------------------------------------------------------------------------

def section_trsp_at_depth(xfld, yfld, maskW, maskS, cds, 
                          grid=None):
    """
    Compute transport of vector quantity at each depth level 
    across latitude(s) defined in lat_vals

    Parameters
    ----------
    xfld, yfld : xarray DataArray
        3D spatial (+ time, optional) field at west and south grid cell edge
    maskW, maskS : xarray DataArray
        defines the section to define transport across
    cds : xarray Dataset
        with all LLC90 coordinates, including: maskW/S, YC
    grid : xgcm Grid object, optional
        denotes LLC90 operations for xgcm, see utils.get_llc_grid

    Returns
    -------
    lat_trsp : xarray DataArray
        transport of vector quantity across denoted latitude band at
        each depth level with dimensions 'time' (if in given dataset),
        'k' (depth), and 'lat' 
    """

    if grid is None:
        grid = get_llc_grid(cds)

    # Initialize empty DataArray with coordinates and dims
    sec_trsp = _initialize_section_trsp_data_array(cds)

    # Apply section mask and sum horizontally
    sec_trsp_x = (xfld * maskW).sum(dim=['i_g','j','tile'])
    sec_trsp_y = (yfld * maskS).sum(dim=['i','j_g','tile'])

    return sec_trsp_x + sec_trsp_y


# -------------------------------------------------------------------------------
#
# All functions below are non-user facing
#
# -------------------------------------------------------------------------------
# Helper functions for the computing volume, heat, and salt transport 
# -------------------------------------------------------------------------------

def _parse_section_trsp_inputs(ds,pt1,pt2,maskW,maskS,section_name):
    """Handle inputs for computing volume, heat, or salt transport across
    a section

    Parameters
    ----------
    see calc_section_vol_trsp

    Returns
    -------
    maskW, maskS : xarray DataArray
        masks defining the section
    """

    use_predefined_section = False
    use_endpoints = False
    use_masks = False

    # Test if section name is in available basins
    if get_section_endpoints(section_name) is not None:
        use_predefined_section = True

    # Test if endpoints provided
    if (pt1 is not None and pt2 is not None):
        use_endpoints = True

    # Test if masks provided
    if (maskW is not None and maskS is not None):
        use_masks = True

    # Test to make sure section is defined by at least one method
    if not use_predefined_section and not use_endpoints and not use_masks:
        raise TypeError('Must provide one method for defining section')

    # First, try to use predefined section
    if use_predefined_section:
        if use_endpoints or use_masks:
            raise TypeError('Cannot provide more than one method for defining section')
        pt1, pt2 = get_section_endpoints(section_name)
        _, maskW, maskS = get_section_line_masks(pt1, pt2, ds)
    else:
        # Secondly, try to use endpoints or mask
        if use_endpoints and use_masks:
            raise TypeError('Cannot provide more than one method for defining section')
        elif use_endpoints:
            _, maskW, maskS = get_section_line_masks(pt1, pt2, ds)

    return maskW, maskS

def _initialize_section_trsp_data_array(cds):
    """Create an xarray DataArray with time, depth, and latitude dims

    Parameters
    ----------
    cds : xarray Dataset
        contains LLC coordinates 'k' and (optionally) 'time'

    Returns
    -------
    da : xarray DataArray
        zero-valued DataArray with time and depth dimensions
    """

    coords = OrderedDict()
    dims = ()

    if 'time' in cds:
        coords.update( {'time': cds['time'].values} )
        dims += ('time',)
        zeros = np.zeros((len(cds['time'].values),
                          len(cds['k'].values)))
    else:
        zeros = np.zeros((len(cds['k'].values)))

    coords.update( {'k': cds['k'].values} )

    dims += ('k',)

    return xr.DataArray(data=zeros, coords=coords, dims=dims)
