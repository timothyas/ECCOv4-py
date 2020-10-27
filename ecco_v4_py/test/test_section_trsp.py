"""
Test routines for computing meridional transport
"""
import warnings
import numpy as np
import xarray as xr
import pytest
import ecco_v4_py

from .test_common import llc_mds_datadirs, get_test_ds, get_test_vectors
from .test_vector_calc import get_fake_vectors

def test_trsp_ds(get_test_ds):
    """stupid simple"""
    exp = get_test_ds
    test = ecco_v4_py.calc_section_trsp._initialize_section_trsp_data_array(exp)
    assert np.all(test.time==exp.time)
    assert np.all(test.k == exp.k)

@pytest.mark.parametrize("name, pt1, pt2, maskW, maskS, expArr",
        [
            ("drakepassage",None,None,None,None,None),
            (None,[-173,65.5],[-164,65.5],None,None,None),
            (None,None,None,True,True,None),
            (None,None,None,None,None,TypeError),
            ("drakepassage",[-173,65.5],[-164,65.5],None,None,TypeError),
            ("drakepassage",None,None,True,True,TypeError),
            (None,[-173,65.5],[-164,65.5],True,True,TypeError),
            ("noname",None,None,None,None,TypeError)
        ])
def test_vol_trsp(get_test_vectors,name,pt1,pt2,maskW,maskS,expArr):
    """compute a volume transport"""

    ds = get_test_vectors
    grid = ecco_v4_py.get_llc_grid(ds)

    ds['U'],ds['V'] = get_fake_vectors(ds['U'],ds['V'])
    ds = ds.rename({'U':'UVELMASS','V':'VVELMASS'})

    if maskW is not None and maskS is not None:
        if maskW and maskS:
            maskW,maskS = ecco_v4_py.vector_calc.get_latitude_masks(30,ds['YC'],grid)

    if expArr is None:
        trsp = ecco_v4_py.calc_section_vol_trsp(ds,
                        pt1=pt1,pt2=pt2,
                        maskW=maskW,maskS=maskS,
                        section_name=name,
                        grid=grid)

        maskW,maskS = ecco_v4_py.calc_section_trsp._parse_section_trsp_inputs(ds,
                        pt1=pt1,pt2=pt2,maskW=maskW,maskS=maskS,
                        section_name=name)

        trspx = (ds['drF']*ds['dyG']*np.abs(maskW)).where(ds['maskW']).sum(dim=['i_g','j','tile'])
        trspy = (ds['drF']*ds['dxG']*np.abs(maskS)).where(ds['maskS']).sum(dim=['i','j_g','tile'])
        test = trsp.vol_trsp_z.reset_coords(drop=True)
        expected = (1e-6*(trspx+trspy)).reset_coords(drop=True)
        xr.testing.assert_allclose(test,expected)

    else:
        with pytest.raises(expArr):
            trsp = ecco_v4_py.calc_section_vol_trsp(ds,
                            pt1=pt1,pt2=pt2,
                            maskW=maskW,maskS=maskS,
                            section_name=name,
                            grid=grid)

            maskW,maskS = ecco_v4_py.calc_section_trsp._parse_section_trsp_inputs(ds,
                            pt1=pt1,pt2=pt2,maskW=maskW,maskS=maskS,
                            section_name=name)

@pytest.mark.parametrize("name, pt1, pt2, maskW, maskS, expArr",
        [
            ("drakepassage",None,None,None,None,None),
            (None,[-173,65.5],[-164,65.5],None,None,None),
            (None,None,None,True,True,None),
            (None,None,None,None,None,TypeError),
            ("drakepassage",[-173,65.5],[-164,65.5],None,None,TypeError),
            ("drakepassage",None,None,True,True,TypeError),
            (None,[-173,65.5],[-164,65.5],True,True,TypeError),
            ("noname",None,None,None,None,TypeError)
        ])
def test_heat_trsp(get_test_vectors,name,pt1,pt2,maskW,maskS,expArr):
    """compute heat transport"""

    ds = get_test_vectors
    grid = ecco_v4_py.get_llc_grid(ds)

    ds['U'],ds['V'] = get_fake_vectors(ds['U'],ds['V'])
    ds = ds.rename({'U':'ADVx_TH','V':'ADVy_TH'})
    ds['DFxE_TH'] = ds['ADVx_TH'].copy()
    ds['DFyE_TH'] = ds['ADVy_TH'].copy()

    if maskW is not None and maskS is not None:
        if maskW and maskS:
            maskW,maskS = ecco_v4_py.vector_calc.get_latitude_masks(30,ds['YC'],grid)

    if expArr is None:
        trsp = ecco_v4_py.calc_section_heat_trsp(ds,
                        pt1=pt1,pt2=pt2,
                        maskW=maskW,maskS=maskS,
                        section_name=name,
                        grid=grid)

        maskW,maskS = ecco_v4_py.calc_section_trsp._parse_section_trsp_inputs(ds,
                        pt1=pt1,pt2=pt2,maskW=maskW,maskS=maskS,
                        section_name=name)

        trspx = (2*np.abs(maskW)).where(ds['maskW']).sum(dim=['i_g','j','tile'])
        trspy = (2*np.abs(maskS)).where(ds['maskS']).sum(dim=['i','j_g','tile'])
        test = trsp.heat_trsp_z.reset_coords(drop=True)
        expected = (1e-15*1029*4000*(trspx+trspy)).reset_coords(drop=True)
        xr.testing.assert_allclose(test,expected)

    else:
        with pytest.raises(expArr):
            trsp = ecco_v4_py.calc_section_heat_trsp(ds,
                            pt1=pt1,pt2=pt2,
                            maskW=maskW,maskS=maskS,
                            section_name=name,
                            grid=grid)

            maskW,maskS = ecco_v4_py.calc_section_trsp._parse_section_trsp_inputs(ds,
                            pt1=pt1,pt2=pt2,maskW=maskW,maskS=maskS,
                            section_name=name)

@pytest.mark.parametrize("name, pt1, pt2, maskW, maskS, expArr",
        [
            ("drakepassage",None,None,None,None,None),
            (None,[-173,65.5],[-164,65.5],None,None,None),
            (None,None,None,True,True,None),
            (None,None,None,None,None,TypeError),
            ("drakepassage",[-173,65.5],[-164,65.5],None,None,TypeError),
            ("drakepassage",None,None,True,True,TypeError),
            (None,[-173,65.5],[-164,65.5],True,True,TypeError),
            ("noname",None,None,None,None,TypeError)
        ])
def test_salt_trsp(get_test_vectors,name,pt1,pt2,maskW,maskS,expArr):
    """compute salt transport"""

    ds = get_test_vectors
    grid = ecco_v4_py.get_llc_grid(ds)

    ds['U'],ds['V'] = get_fake_vectors(ds['U'],ds['V'])
    ds = ds.rename({'U':'ADVx_SLT','V':'ADVy_SLT'})
    ds['DFxE_SLT'] = ds['ADVx_SLT'].copy()
    ds['DFyE_SLT'] = ds['ADVy_SLT'].copy()

    if maskW is not None and maskS is not None:
        if maskW and maskS:
            maskW,maskS = ecco_v4_py.vector_calc.get_latitude_masks(30,ds['YC'],grid)

    if expArr is None:
        trsp = ecco_v4_py.calc_section_salt_trsp(ds,
                        pt1=pt1,pt2=pt2,
                        maskW=maskW,maskS=maskS,
                        section_name=name,
                        grid=grid)

        maskW,maskS = ecco_v4_py.calc_section_trsp._parse_section_trsp_inputs(ds,
                        pt1=pt1,pt2=pt2,maskW=maskW,maskS=maskS,
                        section_name=name)

        trspx = (2*np.abs(maskW)).where(ds['maskW']).sum(dim=['i_g','j','tile'])
        trspy = (2*np.abs(maskS)).where(ds['maskS']).sum(dim=['i','j_g','tile'])
        test = trsp.salt_trsp_z.reset_coords(drop=True)
        expected = (1e-6*(trspx+trspy)).reset_coords(drop=True)
        xr.testing.assert_allclose(test,expected)

    else:
        with pytest.raises(expArr):
            trsp = ecco_v4_py.calc_section_salt_trsp(ds,
                            pt1=pt1,pt2=pt2,
                            maskW=maskW,maskS=maskS,
                            section_name=name,
                            grid=grid)

            maskW,maskS = ecco_v4_py.calc_section_trsp._parse_section_trsp_inputs(ds,
                            pt1=pt1,pt2=pt2,maskW=maskW,maskS=maskS,
                            section_name=name)

@pytest.mark.parametrize("myfunc, fld",
        [   (ecco_v4_py.calc_section_vol_trsp,"vol_trsp"),
            (ecco_v4_py.calc_section_heat_trsp,"heat_trsp"),
            (ecco_v4_py.calc_section_salt_trsp,"salt_trsp")])
@pytest.mark.parametrize("section_name","drakepassage") # unnecessary to do more...
def test_separate_coords(get_test_vectors,myfunc,fld,section_name):
    ds = get_test_vectors
    grid = ecco_v4_py.get_llc_grid(ds)

    ds['U'],ds['V'] = get_fake_vectors(ds['U'],ds['V'])
    ds = ds.rename({'U':'UVELMASS','V':'VVELMASS'})
    for f in ['ADVx_TH','DFxE_TH','ADVx_SLT','DFxE_SLT']:
        ds[f] = ds['UVELMASS'].copy()
    for f in ['ADVy_TH','DFyE_TH','ADVy_SLT','DFyE_SLT']:
        ds[f] = ds['VVELMASS'].copy()

    expected = myfunc(ds,section_name=section_name,grid=grid)
    coords = ds.coords.to_dataset().reset_coords()
    ds = ds.reset_coords(drop=True)

    test = myfunc(ds,section_name,coords=coords,grid=grid)
    xr.test.assert_allclose(test[fld],expected[fld])

@pytest.mark.parametrize("myfunc, fld, xflds, yflds",
        [   (ecco_v4_py.calc_section_vol_trsp,"vol_trsp",
                ['UVELMASS'],['VVELMASS']),
            (ecco_v4_py.calc_section_heat_trsp,"heat_trsp",
                ['ADVx_TH','DFxE_TH'],['ADVy_TH','DFyE_TH']),
            (ecco_v4_py.calc_section_salt_trsp,"salt_trsp",
                ['ADVx_SLT','DFxE_SLT'],['ADVy_SLT','DFyE_SLT'])])
@pytest.mark.parametrize("section_name",["drakepassage"]) # more is unnecessary
def test_separate_coords(get_test_vectors,myfunc,fld,xflds,yflds,section_name):
    ds = get_test_vectors
    grid = ecco_v4_py.get_llc_grid(ds)

    ds['U'],ds['V'] = get_fake_vectors(ds['U'],ds['V'])
    for fx,fy in zip(xflds,yflds):
        ds[fx] = ds['U']
        ds[fy] = ds['V']

    expected = myfunc(ds,section_name=section_name,grid=grid)
    coords = ds.coords.to_dataset().reset_coords()
    ds = ds.reset_coords(drop=True)

    test = myfunc(ds,section_name=section_name,coords=coords,grid=grid)
    xr.testing.assert_allclose(test[fld].reset_coords(drop=True),
                               expected[fld].reset_coords(drop=True))
