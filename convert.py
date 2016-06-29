# useful atmospheric conversions and stuff
import iris
import numpy as np
import iris.analysis.maths as im
import IrisMetConversions as IMC
import sys

con = IMC.constants()

#-------------------
# TEMPERATURE/DENSITY CONVERSIONS
#-------------------

# POTENTIAL TEMPERATURE FROM TEMPERATURE
def theta2temp(TH, P=None, P0 = 1000.):

    '''
    Convert potential temperature to temperature.

    temperature = theta2temp(TH, P=None, P0=1000.)

    TH = cube of potential temperature.
    P (optional) = cube of pressure (only required if pressure is not a coordinate in TH)
    P0 = reference pressure in hPa (default=1000.)
    '''

    # convert p0 to a coordinate with appropriate units
    pref = iris.coords.AuxCoord(P0,
                            long_name='reference_pressure',
                            units='hPa')

    pressure = _get_pressure(TH, P)
    pref.convert_units(pressure.units)

    # if p is specified, use iris arithmetic
    if (type(pressure) is iris.cube.Cube):
        temperature = TH * ( (pressure / pref) ** (con.Rd / con.cpd).data )
    # if pressure is a coordinate object, pressure/pref is not possible, 
    # so calculate using numpy arrays first (it's unitless anyway)
    else: 
        s = TH.shape
        pind = TH.coord_dims(pressure)
        pfrac = iris.util.broadcast_to_shape(pressure.points / pref.points,
                                            s, pind)

        temperature = TH * ( pfrac ** (con.Rd / con.cpd).data )

    temperature.rename('air_temperature')

    return temperature

# TEMPERATURE FROM POTENTIAL TEMPERATURE
def temp2theta(T, P=None, P0 = 1000.):

    '''
    Convert temperature to potential temperature.

    theta = theta2temp(T, p=None, p0=1000.)

    T = cube of temperature.
    P (optional) = cube of pressure (only required if pressure is not a coordinate in T)
    P0 = reference pressure in hPa (default=1000.)
    '''

    # convert p0 to a coordinate with appropriate units
    pref = iris.coords.AuxCoord(P0,
                            long_name='reference_pressure',
                            units='hPa')

    pressure = _get_pressure(T, P)
    pref.convert_units(pressure.units)
       
    # if p is specified, use iris arithmetic
    if (type(pressure) is iris.cube.Cube):
        theta = T * (pref * (pressure)**-1)** (con.Rd / con.cpd).data

    # if pressure is a coordinate object, pressure/pref is not possible, 
    # so calculate using numpy arrays first (it's unitless anyway)
    else: 
        s = T.shape
        pind = T.coord_dims(pressure)

        pfrac = iris.util.broadcast_to_shape(pref.points / pressure.points,
                                            s, pind)

        theta = T * (pfrac** (con.Rd / con.cpd).data)

    theta.rename('potential_temperature')

    return theta

def temp2tvirtual(T, MIXR=None, SH=None):
    '''
    Compute virtual temperature from temperature and humidity

    TV = temp2tvirtual(T, MR)

    T = cube of temperature
    One of the following humidity variables MUST be set
    MIXR = cube of water vapour mixing ratio
    SH = cube of specific humidity
    '''

    # Check humidity has been set, and convert SH to MIXR if necessary
    if MIXR:
        pass
    elif SH:
        MIXR = sh2mixr(SH)
    else:
        print('Need to set either MIXR or SH. No humidity cube found.')
        return

    # force loading of data first, or the function doesn't work
    T.data
    MIXR.data

    TV = T*((1 + (MIXR/(con.Rd/con.Rv)))/(1+MIXR))

    TV.rename('virtual_temperature')

    return TV

def theta2thvirtual(TH, MIXR=None, SH=None):
    '''
    Compute virtual potential temperature from potential temperature and humidity

    THV = theta2thvirtual(TH, MIXR=None, SH=None)

    TH = cube of potential temperature
    One of the following humidity variables MUST be set
    MIXR = cube of water vapour mixing ratio
    SH = cube of specific humidity
    '''

    # Check humidity has been set, and convert SH to MIXR if necessary
    if MIXR:
        pass
    elif SH:
        MIXR = sh2mixr(SH)
    else:
        print('Need to set either MIXR or SH. No humidity cube found.')
        return

    # force loading of data first, or the function doesn't work
    TH.data
    MIXR.data

    THV = TH*((1 + (MIXR/(con.Rd/con.Rv)))/(1+MIXR))

    THV.rename('virtual_potential_temperature')

    return THV

def rho(T, P=None, SH=0.):
    '''
    Compute air density of dry (default) or moist air.

    rho = rho(T, P=None, SH=0.)

    T = cube of temperatures
    P (optional) = cube of pressure (only required if pressure is not a coordinate in T)
    SH = cube of specific humidity 
    '''

    pressure = _get_pressure(T, P)

    # water vapour contribution
    R = ((-1.*SH + 1) * con.Rd)  + (SH*con.Rv)

    # have to invert the division because coord/cube is not possible (but cube/coord is ok)
    rho = ((T*R)/pressure)**-1

    rho.rename('air_density')
    rho.convert_units('kg m-3')

    return rho

    
#---------------------------------
# HUMIDITY CONVERSIONS
#---------------------------------

# conversion relative humidity to mixing ratio
def rh2mixr(RH, T, P=None):
    '''
    Compute mixing ratio from relative humidity, pressure and temperature.

    mixr = rh2mixr(RH, T, P=None)

    RH = cube of relative humidity
    T = cube of temperature
    P (optional) = cube of pressure (only required if pressure is not a coordinate in RH)
    '''

    pressure = _get_pressure(RH, P)

    es = esat(T)
    fact = RH * es
    fact.convert_units(pressure.units)

    mixr = (fact*(con.Mw/con.Md))/(-1.*fact+pressure)
    mixr.long_name = 'water vapour mixing ratio'
    mixr.convert_units('kg kg**-1')

    return mixr

def mixr2rh(MIXR, T, P=None):
    '''
    Compute relative humidity from mixing ratio, pressure and temperature.

    rh = mixr2rh(MIXR, T, P=None)

    MIXR = cube of water vapour mixing ratio
    T = cube of temperature
    P (optional) = cube of pressure (only required if pressure is not a coordinate in MIXR)
    '''

    pressure = _get_pressure(MIXR, P)

    # calculation must be done in kg/kg, so store units first and convert.
    mixrunits = MIXR.units
    MIXR.convert_units('kg kg**-1')

    es = esat(T)
    fact = MIXR * con.Md/con.Mw
    rh = (pressure*es**-1)*fact/(1+fact)
    rh.long_name = 'Relative Humidity'
    rh.convert_units('%')

    # convert back to original units
    MIXR.convert_units(mixrunits)


    return rh


# conversion mixing ratio to specific humidity 
def mixr2sh(MIXR):
    '''
    Convert mixing ratio to specific humidity keeping the same units.

    sh = mixr2sh(MIXR)

    MIXR = cube of mixing ratio
    '''

    # calculation must be done in kg/kg, so store units first and convert.
    mixrunits = MIXR.units
    MIXR.convert_units('kg kg-1')

    sh = MIXR/(1.+MIXR)
    sh.long_name = 'specific humidity'

    # and convert back to original units
    sh.convert_units(mixrunits)
    MIXR.convert_units(mixrunits)

    return sh

def sh2mixr(SH):
    '''
    Convert specific humidity to mixing ratio keeping the same units.

    mixr = mixr2sh(SH)

    SH = cube of specific humidity
    '''

    # calculation must be done in kg/kg, so store units first and convert.
    shunits = SH.units
    SH.convert_units('kg kg-1')

    mixr = SH/(-1.*SH+1)
    mixr.long_name = 'water vapour mixing ratio'

    # convert back to original units
    mixr.convert_units(shunits)
    SH.convert_units(shunits)

    return mixr

def rh2tdew(RH, T):
    '''
    Calculate dew point temperature from relative humidity and temperature 
    according (nearly) to WMO standard procedure.
    The parameters are slightly adapted however to fit the
    formula of Bolton (1980).
    (Code adapted from idl scripts written by Dominik Brunner)


    tdew = rh2tdew(T, RH)

    T = cube of temperature
    RH = cube of relative humidity

    Uses parameters from Bolton (Mon. Wea. Rev., 108, 1047-1053, 1980.)
    '''

    # values from Bolton (Mon. Wea. Rev., 108, 1047-1053, 1980.)
    b = 17.67
    c = 243.5 # c has units degC

    # T must be in degC - store original units so output is the same as input
    Tunits = T.units
    T.convert_units('degC')

    logrh = im.log(RH)
    logrh.convert_units('ln(re 1)')

    # The following calculations are done as numpy arrays before recreating the cube
    # because the operation with cubes is very slow/memory intensive (for some reason)
    fact = logrh.data+(b*T.data)/(c+T.data)
    tdewdata = c*(fact/(-1.*fact+b))

    # create tdew cube taking coordinates from T
    dimcoords = [(coord, i) for i, coord in enumerate(T.coords())]

    tdew = iris.cube.Cube(tdewdata,
                dim_coords_and_dims = dimcoords,
                units = 'degC', long_name = 'dew point temperature')
    
    tdew.convert_units(Tunits)
    T.convert_units(Tunits)

    return tdew

def sh2tdew(SH, T, P=None):
    '''
    Convert specific humidity to dew point temperature

    tdew = sh2tdew(SH, T, P=None)

    SH = cube of specific humidity
    T = cube of temperature
    P (optional) = cube of pressure (only required if pressure is not a coordinate in SH)
    '''

    pressure = _get_pressure(SH, P)

    mixr = IMC.convert.sh2mixr(SH)
    rh = IMC.convert.mixr2rh(mixr, T, P=pressure)

    tdew = rh2tdew(rh, T)

    return tdew

def mixr2tdew(MIXR, T, P=None):
    '''
    Convert water vapour mixing ratio to dew point temperature

    tdew = mixr2tdew(MIXR, T, P=None)

    MIXR = cube of specific humidity
    T = cube of temperature
    P (optional) = cube of pressure (only required if pressure is not a coordinate in MIXR)
    '''

    pressure = _get_pressure(MIXR, P)

    rh = mixr2rh(MIXR, T, P=pressure)
    tdew = rh2tdew(rh, T)

    return tdew


# saturation pressure
def esat(T):
    '''
    Get saturation pressure (units [hPa]) for a given air temperature
    
    es = esat(T)

    T = cube of temperatures
    '''

    # force the loading of T data, as otherwise get some errors
    T.data

    # units must be in K. Note the conversion may change the inputted cube
    if (T.units != 'K'): 
        print('converting temperature units to Kelvin')
        T.convert_units('K')

    TK = iris.coords.AuxCoord(273.15, units='K')
    e1 = iris.coords.AuxCoord(1013.25, units='hPa')

    logTTK = im.log10(T/TK)
    logTTK.units = '1'

# pwr10(cube) is equivalent to 10**cube
    pwr10 = im.IFunc(cube_power_10, cube_power_10_unitcheck)

    es =  e1*pwr10( 10.79586*(1+(-1.*TK*T**-1))-
                    5.02808*logTTK+ 
                    1.50474*1e-4*(-1.*pwr10(-8.29692*(T/TK-1))+1)+ 
                    0.42873*1e-3*(pwr10(4.76955*((-TK*T**-1)+1))-1)-
                    2.2195983 ) 
    es.long_name = 'saturation vapour pressure'
    es.convert_units('hPa')

    return es

#-------------------
# 'UTILITY' FUNCTIONS USED ABOVE
#-------------------

# determine whether to use the cube's pressure coordinate or optional input 'P'
def _get_pressure(cube, P):

    if not P: # optional P input takes priority
        try: 
            pressure = cube.coord(axis='Z')
        except iris.exceptions.CoordinateNotFoundError:
            raise ValueError(
                    'P not set and the input cube does not have a z coordinate')
    else:
        pressure = P

    # check pressure is really pressure
    if not pressure.units.is_convertible('hPa'):
        raise ValueError(
                    'Neither P nor the z coordinate in the input cube are a pressure')
        
    return pressure


# functions to raise 10 to the power of a cube
def cube_power_10(cube_data):
    return 10**(cube_data)

def cube_power_10_unitcheck(cube):
    if cube.units != '1':
        raise ValueError("Units must be '1'")

    return cube.units

