from math import sin, cos, asin, acos, atan2, radians, pi, fabs, sqrt
from angle import *
from transforms import Equatorial
from calendar import Date, Time, JulianDayNumber

"""
Calculate the precession in the RA and Declination at the ending epoch of a
body given the mean coordinates at a starting epoch.
@param coord The mean coordinates referred to epoch_start as a Equatorial object
@param epoch_start The epoch to which the coord is referred as a
       Julian Day Number in TD
@param epoch_end The epoch for which the new coord is to be computed as a
       Julian Day Number in TD
@return The RA and Declination for epoch_end as a SphCoord object
@caution We assume that the coord is already corrected for the proper motion of
         the body over the epoch interval in question
"""
def precession_j2000(coord, epoch_start, epoch_end):
    assert isinstance(coord, Equatorial), 'coord must be a Equatorial'
    assert isinstance(epoch_start, JulianDayNumber)    \
           and isinstance(epoch_end, JulianDayNumber), \
           'epochs much be JulianDayNumbers'

    # In Julian Centuries since J2000
    T = (epoch_start.jdn - 2451545.0) / 36525.0
    t = (epoch_end.jdn - epoch_start.jdn) / 36525.0


    K = (2306.2181 + (1.39656 - 0.000139*T)*T)*t
    t_square = t**2
    zeta = K + ((0.30188 - 0.000344*T) + 0.017998*t)*t_square
    zeta = radians(zeta/3600)
    zappa = K + ((1.09468 + 0.000066*T) + 0.018203*t)*t_square
    zappa = radians(zappa/3600)
    theta = (2004.3109 - (0.85330 + 0.000217*T)*T)*t
    theta -= ((0.42665 + 0.000217*T) + 0.041833)*t_square
    theta = radians(theta/3600)

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


if __name__ == "__main__":
    # theta Persei (proper motion corrected)
    theta_persei_j2000 = Equatorial(Longitude(radians(41.054063)), \
                                    Latitude(radians(49.227750)))
    theta_persei_j20281113_19 = \
        precession_j2000(theta_persei_j2000,
                         JulianDayNumber(Date(2000,1,1),Time(12,0,0)),
                         JulianDayNumber(Date(2028,11,13),Time(4,28,0)))
    print theta_persei_j20281113_19 # 2h46m11.331s, +49deg20'54.54"
