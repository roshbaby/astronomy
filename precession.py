from math import sin, cos, asin, acos, atan2, radians, pi, fabs, sqrt
from angle import Longitude, Latitude
from coordinates import Equatorial
from calendar import Date, Time, JulianDayNumber
from constants import epoch_j2000, julian_century

class Precession:
    def __init__(self, epoch_start, epoch_end):
        assert isinstance(epoch_start, JulianDayNumber)    \
               and isinstance(epoch_end, JulianDayNumber), \
               'epochs much be JulianDayNumbers'
        # In Julian Centuries since J2000
        T = (epoch_start.jdn - epoch_j2000.jdn) / julian_century
        t = (epoch_end.jdn - epoch_start.jdn) / julian_century

        K = (2306.2181 + (1.39656 - 0.000139*T)*T)*t
        t_square = t**2

        self.zeta = K + ((0.30188 - 0.000344*T) + 0.017998*t)*t_square
        self.zeta = radians(self.zeta/3600.0)
        self.zappa = K + ((1.09468 + 0.000066*T) + 0.018203*t)*t_square
        self.zappa = radians(self.zappa/3600.0)
        self.theta = (2004.3109 - (0.85330 + 0.000217*T)*T)*t
        self.theta -= ((0.42665 + 0.000217*T) + 0.041833)*t_square
        self.theta = radians(self.theta/3600.0)


    def apply_correction(self, coord):
        assert isinstance(coord, Equatorial), 'coord must be a Equatorial'

        alpha  = coord.a.rads + self.zeta
        delta0 = coord.b.rads

        A = cos(delta0)*sin(alpha)
        B = cos(self.theta)*cos(delta0)*cos(alpha) - sin(self.theta)*sin(delta0)
        C = sin(self.theta)*cos(delta0)*cos(alpha) + cos(self.theta)*sin(delta0)

        alpha = self.zappa + atan2(A,B)
        # If the object is close to the celestial pole
        if fabs(delta0 - pi/2) < radians(0.5):
            delta = acos(sqrt(A**2+B**2))
        else:
            delta = asin(C)

        return Equatorial(Longitude(alpha), Latitude(delta))


if __name__ == "__main__":
    epoch_start = epoch_j2000
    epoch_end = JulianDayNumber(Date(2028,11,13),Time(4,28,0))

    prec = Precession(epoch_start, epoch_end)
    # theta Persei (proper motion corrected)
    theta_persei_epoch_start = Equatorial(
        Longitude(radians(41.054063)), Latitude(radians(49.227750))
    )
    print prec.apply_correction(theta_persei_epoch_start) # 2h46m11.331s, +49deg20'54.54"
