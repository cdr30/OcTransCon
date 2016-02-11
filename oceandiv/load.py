"""
Module containing routines to load netcdf data into <iris.cube.Cube>

"""


import iris
import numpy as np


class ShapeError(Exception):
    pass


def test_shape(cube, shape, label):
    """ Raise ShapeError if cube.shape != shape """
    
    if cube.shape != shape:
        raise ShapeError('Shape of %s must be %s, not %s ' % 
                         (label, repr(shape), repr(cube.shape)))


def load_cube(f, ncvar):
    """ Return specified netcdf variable as an <iris.cube.Cube> """
    
    var_constraint = iris.Constraint(cube_func=lambda c:
                                     c.var_name == ncvar)
    cube = iris.load_cube(f, var_constraint)
    
    return cube


def mask_cube(cube, mdi):
    """ Return cube following application of mask where data == mdi """
    
    mask = (cube.data == mdi)
    cube.data = np.ma.MaskedArray(cube.data, mask=mask, fill_value=mdi)
    
    return cube


def load_geodata(config, geotype='areas'):
    """ Return geographic data as an <iris.cube.Cube> """
    
    datadir = config.get(geotype, 'dir')
    f = config.get(geotype, 'f')
    mdi = config.getfloat(geotype, 'mdi')
    ncvar = config.get(geotype, 'var')
    
    nx = config.getint('metadata', 'nx')
    ny = config.getint('metadata', 'ny')
    shape = (ny, nx)
    
    cube = load_cube(datadir + f, ncvar)
    test_shape(cube, shape, geotype)
    cube = mask_cube(cube, mdi)

    return cube


def load_data(config, dtype='ohc'):
    """ Return data as list of <iris.cube.Cube> objects"""
    
    cubes = []
    ncubes = config.getint('metadata', 'n%s' % (dtype))
    nx = config.getint('metadata', 'nx')
    ny = config.getint('metadata', 'ny')
    nt = config.getint('metadata', 'nt')
    shape = (nt, ny, nx)
    
    for ncube in range(ncubes):
        n = ncube + 1    
        datadir = config.get('%s%i' % (dtype, n), 'dir')
        f = config.get('%s%i' % (dtype, n), 'f')
        ncvar = config.get('%s%i' % (dtype, n), 'var')
        t0 = config.getint('%s%i' % (dtype, n), 't0')
        mdi = config.getfloat('%s%i' % (dtype, n), 'mdi')
        
        cube = load_cube(datadir + f, ncvar)
        cube = cube[t0:t0+nt]
        test_shape(cube, shape,'%s%i' % (dtype, n))
        cube = mask_cube(cube, mdi)
        cubes.append(cube)

    return cubes



