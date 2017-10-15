from math import sin, cos, tan, asin, acos, atan2, sqrt, pi, degrees, radians, fabs, copysign
from types import IntType
from numbers import Number
from mymath import pfmod
from angle import Angle, Latitude, Longitude
from coordinates import SphCoord, Equatorial
from interpolation import Inter3polate
from calendar import Date, Time, JulianDayNumber
from constants import epoch_j2000, julian_year, julian_century
from precession import Precession

def refraction(alti_):
    """Compute the correction to altitude due to refraction.
    @param alt_ The *apparent* altitude as an Angle object.
    @return The correction to alt_ as an Angle object.
    """
    assert isinstance(alti_, Angle), 'alti_ should be an Angle'
    alti_degs = degrees(alti_.rads)
    denom = tan(radians(alti_degs + 7.31/(alti_degs+4.4)))
    R = 1/denom # units of arc-minutes
    #return Angle(radians((R + 0.001351521673756295)/60))
    R -= 0.06*sin(radians(14.7*R/60+13)) # units of arc-minutes
    return Angle(radians(R/60))

def angular_separation(sph1,sph2):
    """Compute the angular separation between two spherical coordinates using
    Thierry Pauwel's formula.
    @param sph1 SphCoord for object 1.
    @param sph2 SphCoord for object 2.
    @return The angular separation as an Angle object.
    """
    assert isinstance(sph1, SphCoord), 'sph1 must be a SphCoord'
    assert isinstance(sph2, SphCoord), 'sph2 must be a SphCoord'

    alpha1, delta1 = sph1.a.rads, sph1.b.rads
    alpha2, delta2 = sph2.a.rads, sph2.b.rads

    x = cos(delta1)*sin(delta2) - sin(delta1)*cos(delta2)*cos(alpha2-alpha1)
    y = cos(delta2)*sin(alpha2-alpha1)
    z = sin(delta1)*sin(delta2) + cos(delta1)*cos(delta2)*cos(alpha2-alpha1)

    d = atan2(sqrt(x*x+y*y),z)
    return Angle(d)

def least_angular_separation(coords1, coords2, precision):
    """Compute the least angular separation between two celestial objects
    using iteration.
    @param coords1 List of 3 SphCoords for object 1 (equidistant times).
    @param coords2 List of 3 SphCoords for object 2 (equidistant times).
    @param precision The precision at which to terminate the iteration.
    @return The least angular separation as an Angle object
    """
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
        # Correction to n
        n_ = -(u*u_ + v*v_)/float(u_**2 + v_**2)
        n += n_
        if fabs(n_) <= fabs(precision):
            break
    # Compute the final values of u and v
    u = u_ipol.compute(0+n)
    v = v_ipol.compute(0+n)
    #@todo return final 'n' as well
    return Angle(radians(sqrt(u**2+v**2)/3600))

def relative_position_angle(ref, obj):
    """Compute the Relative Position Angle.
    @param ref The SphCoord of the reference body.
    @param obj The SphCoord of the object for which the RPA is to be computed.
    @return The RPA as an Angle object.
    """
    assert isinstance(ref, SphCoord) and isinstance(obj, SphCoord), \
           'Arguments should be a SphCoord'
    alpha1 = obj.a.rads
    delta1 = obj.b.rads
    alpha2 = ref.a.rads
    delta2 = ref.b.rads
    dalpha = alpha1-alpha2

    rpa = atan2(sin(dalpha), cos(delta2)*tan(delta1)-sin(delta2)*cos(dalpha))
    return Angle(rpa)

def proper_motion(coord, r_parsecs, v_parsecs_per_year, annual_pm, epoch_yrs):
    """Compute the proper motion of a star over the given number of years.
    @param coord The coordinates of the star as a Equatorial object at the epoch.
    @param r_parsecs The radial distance to the star in parsecs.
    @param v_parsecs_per_year The radial velocity of the star in parsecs/year.
    @param annual_pm The annual proper motion of the star as a tuple of Angles.
    @param epoch_yrs The number of years from the starting epoch.
    @return The updated coordinates for the star as an Equatorial object.
    """
    assert isinstance(coord, Equatorial), 'coord should be a Equatorial'
    assert isinstance(annual_pm, tuple), 'annual_pm should be a tuple'
    assert len(annual_pm) == 2, 'annual_pm should be a 2-element tuple'
    assert isinstance(annual_pm[0], Angle) and isinstance(annual_pm[1],Angle), \
           'annual_pm elements should be Angles'
    alpha_0, delta_0 = coord.a.rads, coord.b.rads
    dalpha, ddelta = annual_pm[0].rads, annual_pm[1].rads

    x = r_parsecs*cos(delta_0)*cos(alpha_0) # parsecs
    y = r_parsecs*cos(delta_0)*sin(alpha_0) # parsecs
    z = r_parsecs*sin(delta_0)              # parsecs

    dx = (x/r_parsecs)*v_parsecs_per_year - z*ddelta*cos(alpha_0) - y*dalpha
    dy = (y/r_parsecs)*v_parsecs_per_year - z*ddelta*sin(alpha_0) + x*dalpha
    dz = (z/r_parsecs)*v_parsecs_per_year + r_parsecs*ddelta*cos(delta_0)
    # dx,dy,dz are in parsecs/year

    x_new = x + epoch_yrs*dx # parsecs
    y_new = y + epoch_yrs*dy # parsecs
    z_new = z + epoch_yrs*dz # parsecs

    alpha_new = atan2(y_new, x_new)
    delta_new = atan2(z_new, sqrt(x_new**2+y_new**2))

    return Equatorial(Longitude(alpha_new), Latitude(delta_new))

def proper_motion_classical(coord, annual_pm, epoch_yrs):
    """Compute the proper motion using the classical method of uniform changes
    in RA and declination.
    @param coord The coordinates of the star as a Equatorial object at the epoch.
    @param annual_pm The annual proper motion of the star as a tuple of Angles.
    @param epoch_yrs The number of years from the starting epoch.
    @return The updated coordinates for the star as an Equatorial object.
    """
    assert isinstance(coord, Equatorial), 'coord should be a Equatorial'
    assert isinstance(annual_pm, tuple), 'annual_pm should be a tuple'
    assert len(annual_pm) == 2, 'annual_pm should be a 2-element tuple'
    assert isinstance(annual_pm[0], Angle) and isinstance(annual_pm[1],Angle), \
           'annual_pm elements should be Angles'
    alpha_0, delta_0 = coord.a.rads, coord.b.rads
    dalpha, ddelta = annual_pm[0].rads, annual_pm[1].rads
    alpha_new = alpha_0 + dalpha*epoch_yrs
    delta_new = delta_0 + ddelta*epoch_yrs
    return Equatorial(Longitude(alpha_new), Latitude(delta_new))

def deprecated_precession(coord, epoch_start, epoch_end):
    """Calculate the precession in the RA and Declination at the ending epoch
    of a body given the mean coordinates at a starting epoch.
    @param coord The mean coordinates referred to epoch_start as a Equatorial object
    @param epoch_start The epoch to which the coord is referred as a
    Julian Day Number in TD
    @param epoch_end The epoch for which the new coord is to be computed as a
    Julian Day Number in TD
    @return The RA and Declination for epoch_end as a SphCoord object
    @caution We assume that the coord is already corrected for the proper
    motion of the body over the epoch interval in question
    """
    assert isinstance(coord, Equatorial), 'coord must be a Equatorial'
    assert isinstance(epoch_start, JulianDayNumber)    \
           and isinstance(epoch_end, JulianDayNumber), \
           'epochs much be JulianDayNumbers'

    T = (epoch_start.jdn - epoch_j2000.jdn) / julian_century
    t = (epoch_end.jdn - epoch_start.jdn) / julian_century

    K = (2306.2181 + (1.39656 - 0.000139*T)*T)*t
    t_square = t**2

    zeta = K + ((0.30188 - 0.000344*T) + 0.017998*t)*t_square # arc-seconds
    zeta = radians(zeta/3600.0)
    zappa = K + ((1.09468 + 0.000066*T) + 0.018203*t)*t_square # arc-seconds
    zappa = radians(zappa/3600.0)
    theta = (2004.3109 - (0.85330 + 0.000217*T)*T)*t # arc-seconds
    theta -= ((0.42665 + 0.000217*T) + 0.041833*t)*t_square # arc-seconds
    theta = radians(theta/3600.0)

    alpha = coord.a.rads + zeta
    delta0 = coord.b.rads

    A = cos(delta0)*sin(alpha)
    B = cos(theta)*cos(delta0)*cos(alpha) - sin(theta)*sin(delta0)
    C = sin(theta)*cos(delta0)*cos(alpha) + cos(theta)*sin(delta0)

    alpha = zappa + atan2(A,B)
    # If the object is close to the celestial pole
    if fabs(delta0 - pi/2) < radians(0.5):
        delta = acos(sqrt(A**2+B**2))
    else:
        delta = asin(C)

    return Equatorial(Longitude(alpha), Latitude(delta))

def equation_of_kepler_iterative(M_, e_, prec_):
    """Calculate the Eccentric anomaly from the Mean Anomaly using iteration.
    @param M_ The Mean anomaly as a Longitude object
    @param e_ The eccentricity of the orbit
    @param prec_ The precision at which to stop the iteration
    @return The Eccentric anomaly as a Longitude object
    """
    assert isinstance(M_, Longitude), 'M_ should be a Longitude'
    assert isinstance(e_, Number), 'e_ should be a Number'
    assert 0 <= e_ < 1, 'e_ should lie between 0 and 1'
    assert isinstance(prec_, Number), 'prec_ should be a Number'
    E = M_.rads
    while True:
        dE = (M_.rads + e_*sin(E) - E)/(1 - e_*cos(E))
        E += dE
        if fabs(dE) <= prec_: break
    return Longitude(E)

def equation_of_kepler_binarysearch(M_, e_, loop_=33):
    """Calculate the Eccentric anomaly from the Mean Anomaly using binary search.
    @param M_ The Mean anomaly as a Longitude object
    @param e_ The eccentricity of the orbit
    @param loop_ The number of times to
    @return The Eccentric anomaly as a Longitude object
    """
    assert isinstance(M_, Longitude), 'M_ should be a Longitude'
    assert isinstance(e_, Number), 'e_ should be a Number'
    assert 0 <= e_ < 1, 'e_ should lie between 0 and 1'
    assert type(loop_) is IntType, 'loop_ should be an integer'
    M_ = pfmod(M_.rads,2*pi)
    if M_ > pi:
        sign = -1
        M_ = 2*pi - M_
    else:
        sign = 1
    E0, D = pi/2, pi/4
    for idx in range(loop_):
        M1 = E0 - e*sin(E0)
        E0 += copysign(D, M_-M1)
        D /= 2
    E0 *= sign
    return Longitude(E0)


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
    print

    # Proper Motion tests
    # Sirius
    coord_0 = Equatorial(Longitude(radians(101.286962)),
                         Latitude(radians(-16.716108)))
    annual_pm = (
        Angle(-radians(0.03847/3600.0*15.0)),
        Angle(-radians(1.2053/3600.0))
    )
    epoch_0 = 2000.0
    epoch_targets = [1000.0, 0.0, -1000.0, -2000.0, -10000.0]
    for epoch in epoch_targets:
        epoch_yrs = epoch-epoch_0
        print proper_motion(coord_0, 2.64, -7.6/977792.0, annual_pm, epoch_yrs)
    print

    # theta Persei at epoch J2000
    epoch_start = epoch_j2000
    epoch_target = JulianDayNumber(Date(2028,11,13),Time(4,28,0)) # 2028 Nov 13.19 TD
    theta_persei_epoch_start = Equatorial(
        Longitude( radians((2.0 +44.0/60.0 + 11.986/3600.0)*15) ),
        Latitude( radians(49.0 +13.0/60.0 +42.48/3600.0) )
    )
    annual_pm = (
        Angle( radians(0.03425/3600.0*15) ),
        Angle( -radians(0.0895/3600.0) )
    )
    epoch_yrs = (epoch_target.jdn - epoch_start.jdn)/julian_year

    theta_persei_epoch_start_pm_corrected = \
        proper_motion_classical(theta_persei_epoch_start, annual_pm, epoch_yrs)
    print theta_persei_epoch_start_pm_corrected # 2h44m12.975s, +49deg 13'39.90"
    print

    # Precession test
    prec = Precession(epoch_start, epoch_target)
    theta_persei_epoch_target = \
        prec.apply_correction(theta_persei_epoch_start_pm_corrected)
    print theta_persei_epoch_target # 2h46m11.331s, +49deg20'54.54"
    print deprecated_precession(theta_persei_epoch_start_pm_corrected, epoch_start, epoch_target)
    print

    # Test for Kepler's Equation
    prec = 1e-6
    elist = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99 ]
    for e in elist:
        print degrees(equation_of_kepler_iterative(Longitude(radians(5)), e, prec).rads)
        print degrees(equation_of_kepler_binarysearch(Longitude(radians(5)), e, 55).rads)
    print
