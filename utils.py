from math import sin, cos, tan, atan2, sqrt, pi, degrees, radians, fabs
from angle import Angle, Latitude, Longitude
from transforms import SphCoord
from interpolation import Inter3polate


"""
Compute the correction due to refraction to altitude
@param alt_ The *apparent* altitude as an Angle object
"""
def refraction(alti_):
    assert isinstance(alti_, Angle), 'alti_ should be an Angle'
    alti_degs = degrees(alti_.rads)
    denom = tan(radians(alti_degs + 7.31/(alti_degs+4.4)))
    R = 1/denom # units of arc-minutes
    #return Angle(radians((R + 0.001351521673756295)/60))
    R -= 0.06*sin(radians(14.7*R/60+13)) # units of arc-minutes
    return Angle(radians(R/60))


"""
Compute the angular separation between two spherical coordinates
Uses Thierry Pauwel's formula
@param sph1 SphCoord for object 1
@param sph2 SphCoord for object 2
"""
def angular_separation(sph1,sph2):
    assert isinstance(sph1, SphCoord), 'sph1 must be a SphCoord'
    assert isinstance(sph2, SphCoord), 'sph2 must be a SphCoord'

    alpha1, delta1 = sph1.a.rads, sph1.b.rads
    alpha2, delta2 = sph2.a.rads, sph2.b.rads

    x = cos(delta1)*sin(delta2) - sin(delta1)*cos(delta2)*cos(alpha2-alpha1)
    y = cos(delta2)*sin(alpha2-alpha1)
    z = sin(delta1)*sin(delta2) + cos(delta1)*cos(delta2)*cos(alpha2-alpha1)

    d = atan2(sqrt(x*x+y*y),z)
    return Angle(d)


"""
Computes the least angular separation between two celestial objects
using interpolation
@param coords1 List of 3 SphCoords for object 1 (equidistant times)
@param coords2 List of 3 SphCoords for object 2 (equidistant times)
@param precision The precision at which to terminate the iteration
"""
def least_angular_separation(coords1, coords2, precision):
    assert isinstance(coords1, list) and isinstance(coords2, list), \
           'coords must be lists'
    # CAUTION: Pythonism below!
    assert len(coords1) == 3 == len(coords2), \
           'coords must be 3 element lists'
    for sph1, sph2 in zip(coords1, coords2):
        assert isinstance(sph1, SphCoord) and isinstance(sph2, SphCoord), \
               'Elements of coords must be SphCoords'
    assert precision != 0, 'precision should not be 0'

    u_list = []
    v_list = []
    for idx in range(len(coords1)):
        sph1, sph2 = coords1[idx], coords2[idx]
        alpha1, delta1 = sph1.a.rads, sph1.b.rads
        alpha2, delta2 = sph2.a.rads, sph2.b.rads
        dalpha = alpha2 - alpha1
        ddelta = delta2 - delta1
        K = (degrees(1)*3600) / \
            (1 + (sin(delta1)**2)*tan(dalpha)*tan(dalpha/2))
        u = -K*(1-tan(delta1)*sin(ddelta))*cos(delta1)*tan(dalpha)
        v = K*(sin(ddelta)+sin(delta1)*cos(delta1)*tan(dalpha)*tan(dalpha/2))
        u_list.append(u)
        v_list.append(v)
    #print 'u_list:', u_list
    #print 'v_list:', v_list

    # For interpolation, treat the u_list as the values of y corresponding to
    # x = -1,0,1. Ditto for v_list
    u_ipol = Inter3polate( zip([-1,0,1], u_list) )
    v_ipol = Inter3polate( zip([-1,0,1], v_list) )
    n = 0.0  # Interpolating factor around the central value
    while (True):
        # Interpolate for u and v around the central value of x=0
        u = u_ipol.compute(0+n)
        v = v_ipol.compute(0+n)
        # Determine variations in u, v
        u_ = (u_list[2]-u_list[0])/2.0 + n*u_ipol.c
        v_ = (v_list[2]-v_list[0])/2.0 + n*v_ipol.c
        #print 'n:', n, 'u:', u, 'v:', v, 'u_:', u_, 'v_:', v_
        # Correction to n
        n_ = -(u*u_ + v*v_)/float(u_**2 + v_**2)
        #print 'n_:', n_
        n += n_
        if fabs(n_) <= fabs(precision):
            break
    # Compute the final values of u and v
    u = u_ipol.compute(0+n)
    v = v_ipol.compute(0+n)
    #print 'final n:', n, 'u:', u, 'v:', v
    #@todo return final 'n' as well
    return Angle(radians(sqrt(u**2+v**2)/3600))


"""
Relative Position Angle
@param ref The SphCoord of the reference body
@param obj The SphCoord of the object for which the RPA is to be computed
@return Angle object corresponding to the RPA
"""
def relative_position_angle(ref, obj):
    assert isinstance(ref, SphCoord) and isinstance(obj, SphCoord), \
           'Arguments should be a SphCoord'
    alpha1 = obj.a.rads
    delta1 = obj.b.rads
    alpha2 = ref.a.rads
    delta2 = ref.b.rads
    dalpha = alpha1-alpha2

    rpa = atan2(sin(dalpha), cos(delta2)*tan(delta1)-sin(delta2)*cos(dalpha))
    return Angle(rpa)


def proper_motion_correction():
    assert False, '@todo To Be Implemented'


if __name__ == "__main__":
    ## Refraction tests
    print refraction(Angle(0))
    print refraction(Angle(pi/4))
    print refraction(Angle(pi/2))

    ## Angular Separation Tests
    # List of (alpha1, delta1, alpha2, delta2) tuples
    data = [
        # sph1 a number
        (
            pi/2, SphCoord(Longitude(pi/6), Latitude(pi/3))
        ),
        # sph2 a tuple
        (
            SphCoord(Longitude(pi/6), Latitude(pi/3)), (0, pi/2)
        ),
        # (0,45), (30,60). Should be just over 23 degrees
        (
            SphCoord(Longitude(0), Latitude(pi/4)),
            SphCoord(Longitude(pi/6), Latitude(pi/3))
        ),
        # Example 17.a, Arcturus and Spica. (32.7930 degrees)
        (
            SphCoord(Longitude(radians(213.9154)), Latitude(radians(19.1825))),
            SphCoord(Longitude(radians(201.2983)), Latitude(radians(-11.1614)))
        ),
    ]
    for item in data:
        try:
            sph1, sph2 = item
            print angular_separation(sph1,sph2)
        except Exception as e:
            print 'Error:', e
    print

    ## Least Angular Separation
    precision = 1e-5          # Corresponds to 0.9sec for interval of 1 day
    data = [
        (0, []),              # coord1 not a list
        ([], 0),              # coord2 not a list
        ([1], [1,2,3]),       # coord1 not a 3 element list
        ([1,2,3], [1]),       # coord2 not a 3 element list
        ([1,2,3], [1,2,3,4]), # coord2 not a 3 element list
        ([1,2,3], [1,2,3])    # coord elements are not SphCoords
    ]
    for tup in data:
        try:
            least_angular_separation(tup[0], tup[1], precision)
        except AssertionError as e:
            print 'Error:', e

    # For 1978 0h TD (pg 110)
    # Mercury
    mercury_coords = [
        # Sep 13
        SphCoord (
            Longitude(radians((10+(29+44.27/60.0)/60.0)*15)),
            Latitude(radians(11+(02+5.9/60)/60.0))
        ),
        # Sep 14
        SphCoord (
            Longitude(radians((10+(36+19.63/60.0)/60.0)*15)),
            Latitude(radians(10+(29+51.7/60)/60.0))
        ),
        # Sep 15
        SphCoord (
            Longitude(radians((10+(43+1.75/60.0)/60.0)*15)),
            Latitude(radians(9+(55+16.7/60)/60.0))
        )
    ]
    saturn_coords = [
        # Sep 13
        SphCoord (
            Longitude(radians((10+(33+29.64/60.0)/60.0)*15)),
            Latitude(radians(10+(40+13.2/60)/60.0))
        ),
        # Sep 14
        SphCoord (
            Longitude(radians((10+(33+57.97/60.0)/60.0)*15)),
            Latitude(radians(10+(37+33.4/60)/60.0))
        ),
        # Sep 15
        SphCoord (
            Longitude(radians((10+(34+26.22/60.0)/60.0)*15)),
            Latitude(radians(10+(34+53.9/60)/60.0))
        )
    ]
    try:
        least_angular_separation(mercury_coords, saturn_coords, 0)
    except AssertionError as e:
        print 'Error:', e

    print least_angular_separation(mercury_coords, saturn_coords, precision) # 3'44"
    print least_angular_separation(saturn_coords, mercury_coords, precision) # 3'44"

    for m,s in zip(mercury_coords,saturn_coords):
        print relative_position_angle(m,s)

    try:
        proper_motion_correction()
    except AssertionError as e:
        print 'Error:', e
