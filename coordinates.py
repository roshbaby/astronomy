import sys
from math import sin, cos, tan, asin, atan2, pi, radians
from angle import Angle, Longitude, Latitude
from constants import epsilon_j2000

DEFAULT_ENCODING = 'UTF-8'

"""
A Spherical Coordinate is a Longitude and Latitude pair.
It is a glorified tuple, but useful for type checking inputs to functions and
constructors.
"""
class SphCoord:
    """
    @param alpha_ The Longitude (can also be Right Ascension or Azimuth)
    @param delta_ The Latitude (can also be Declination or Altitude)
    """
    def __init__(self, a_, b_):
        assert isinstance(a_, Longitude), 'a_ must be a Longitude'
        assert isinstance(b_, Latitude), 'b_ must be a Latitude'
        self.a = a_ # Can be used as an R.A., longitude, azimuth
        self.b = b_ # can be used as a declination, latitude, or altitude


"""
Equatorial Coordinates
A Right Ascension and Declination pair
"""
class Equatorial(SphCoord):
    """
    Equatorial to Ecliptic conversion
    @param epsilon_ The obliquity of the Ecliptic as an Angle object
    @return Ecliptical object
    """
    def to_ecliptical(self, epsilon_):
        assert isinstance(epsilon_, Angle), 'epsilon_ must be Angle'
        alpha_rad, delta_rad = self.a.rads, self.b.rads
        epsilon_rad = epsilon_.rads
        lambda_rad = atan2(sin(alpha_rad)*cos(epsilon_rad)
                           + tan(delta_rad)*sin(epsilon_rad),
                           cos(alpha_rad))
        beta_rad   = asin(sin(delta_rad)*cos(epsilon_rad)
                          - cos(delta_rad)*sin(epsilon_rad)*sin(alpha_rad))
        return Ecliptical(Longitude(lambda_rad), Latitude(beta_rad))

    """
    Equatorial to AltAzimuthal conversion
    @param lati_ The observer's geographic latitude as a Latitude object
    @param hangle_ The hour angle as a Longitude object
    @return AltAzimuthal object
    """
    def to_altazimuthal(self, lati_, hangle_):
        assert isinstance(lati_, Latitude), 'lati_ must be a Latitude'
        assert isinstance(hangle_, Longitude), 'hangle_ must be a Longitude'
        lati_rad, decl_rad, hangle_rad = lati_.rads, self.b.rads, hangle_.rads
        azimuth_rad = atan2(sin(hangle_rad), cos(hangle_rad)*sin(lati_rad) \
                            - tan(decl_rad)*cos(lati_rad))
        altitude_rad = asin(sin(lati_rad)*sin(decl_rad) \
                            + cos(lati_rad)*cos(decl_rad)*cos(hangle_rad))
        return AltAzimuthal(Longitude(azimuth_rad),Latitude(altitude_rad))

    def __unicode__(self):
        return u'(\u03B1:' + self.a.hms() \
               + u', \u03B4:' + self.b.dms() + u')'

    def __str__(self):
        return unicode(self).encode(sys.stdout.encoding or DEFAULT_ENCODING,
                                    'replace')


"""
Ecliptical Coordinates
A Longitude and Latitude pair
"""
class Ecliptical(SphCoord):
    """
    Ecliptic to Equatorial conversion
    @param epsilon_ The obliquity of the Ecliptic as an Angle object
    @return Equatorial object
    """
    def to_equatorial(self, epsilon_):
        assert isinstance(epsilon_, Angle), 'epsilon_ must be Angle'
        lambda_rad, beta_rad = self.a.rads, self.b.rads
        epsilon_rad = epsilon_.rads
        alpha_rad = atan2(sin(lambda_rad)*cos(epsilon_rad)
                          - tan(beta_rad)*sin(epsilon_rad), cos(lambda_rad))
        delta_rad = asin(sin(beta_rad)*cos(epsilon_rad)
                         + cos(beta_rad)*sin(epsilon_rad)*sin(lambda_rad))
        return Equatorial(Longitude(alpha_rad), Latitude(delta_rad))

    def __unicode__(self):
        return u'(\u03BB:' + self.a.dms() \
               + u', \u03B2:' + self.b.dms() + u')'

    def __str__(self):
        return unicode(self).encode(sys.stdout.encoding or DEFAULT_ENCODING,
                                    'replace')

"""
Altitude-Azimuth (Local Horizontal) Coordinates
An Altitude and Azimuth pair
"""
class AltAzimuthal(SphCoord):
    def __unicode__(self):
        return u'(A:' + self.a.dms() + u', h:' + self.b.dms() + u')'
    def __str__(self):
        return unicode(self).encode(sys.stdout.encoding or DEFAULT_ENCODING,
                                    'replace')


if __name__ == "__main__":
    ## Spherical Coordinate Tests
    data = [
        (0, pi/2), # Not Longitudes or Latitudes
        (Latitude(0), Latitude(pi/2)), # alpha is not a Longitude
        (Longitude(0), Longitude(pi/2)), # delta is not a Latitude
        (Longitude(0), Angle(pi/2)) # delta is not a Latitude
    ]
    for tup in data:
        try:
            sph = SphCoord(tup[0], tup[1])
        except AssertionError as e:
            print 'Error:', e
    print

    equa = Equatorial(Longitude(0),Latitude(0))
    print equa
    ecli = equa.to_ecliptical(epsilon_j2000)
    print ecli
    print ecli.to_equatorial(epsilon_j2000)
    print

    equa = Equatorial(Longitude(pi/2),Latitude(0))
    print equa
    ecli = equa.to_ecliptical(epsilon_j2000)
    print ecli
    print ecli.to_equatorial(epsilon_j2000)
    print

    ecli = Ecliptical(Longitude(pi),Latitude(pi/4))
    print ecli
    equa = ecli.to_equatorial(epsilon_j2000)
    print equa
    print equa.to_ecliptical(epsilon_j2000)
    print

    ecli = Ecliptical(Longitude(0.75*pi),Latitude(-pi/4))
    print ecli
    equa = ecli.to_equatorial(epsilon_j2000)
    print equa
    print equa.to_ecliptical(epsilon_j2000)
    print

    equa = Equatorial(
        Longitude(0),
        Latitude(-radians(6 + 43/60.0 + 11.61/3600.0))
    )
    altazi = equa.to_altazimuthal(
        Latitude(radians(38 + 55/60.0 + 17/3600.0)),
        Longitude(radians(64.352133))
    )
    print altazi # A: 68deg 02'01.3", h: 15deg 07'29.6"
