import sys
from math import sin, cos, tan, asin, atan2, pi, radians
from angle import *

DEFAULT_ENCODING = 'UTF-8'

""" Equatorial Coordinates """
class Equatorial:
    """
    @param alpha_ Equatorial Right Ascension as an Angle object
    @param delta_ Equatorial Declination as an Angle object
    """
    def __init__(self, alpha_, delta_):
        assert isinstance(alpha_, Longitude), 'alpha_ must be a Longitude'
        assert isinstance(delta_, Latitude), 'delta_ must be a Longitude'
        self.alpha_angle = alpha_
        self.delta_angle = delta_

    """
    Equatorial to Ecliptic conversion
    Given an Equatorial object and Obliquity of the Ecliptic as an Angle object
    convert to an Ecliptical object
    @return Ecliptical object
    """
    def to_ecliptical(self, epsilon_):
        assert isinstance(epsilon_, Angle),'epsilon_ must be Angle'
        alpha_rad, delta_rad = self.alpha_angle.rads, self.delta_angle.rads
        epsilon_rad = epsilon_.rads
        lambda_rad = atan2(sin(alpha_rad)*cos(epsilon_rad)
                           + tan(delta_rad)*sin(epsilon_rad),
                           cos(alpha_rad))
        beta_rad   = asin(sin(delta_rad)*cos(epsilon_rad)
                          - cos(delta_rad)*sin(epsilon_rad)*sin(alpha_rad))
        return Ecliptical(Longitude(lambda_rad), Latitude(beta_rad))

    def __unicode__(self):
        return u'(\u03B1:' + self.alpha_angle.hms() \
               + u', \u03B4:' + self.delta_angle.dms() + u')'

    def __str__(self):
        return unicode(self).encode(sys.stdout.encoding or DEFAULT_ENCODING,
                                    'replace')


""" Ecliptical Coordinates """
class Ecliptical:
    """
    @param lambda_ Ecliptical Longitude in Radians
    @param beta_   Ecliptical Latitude in Radians
    """
    def __init__(self, lambda_, beta_):
        assert isinstance(lambda_,Angle), 'lambda_ must be a Longitude'
        assert isinstance(beta_, Angle), 'beta_ must be a Latitude'
        self.lambda_angle = lambda_
        self.beta_angle = beta_

    """
    Ecliptic to Equatorial conversion
    Given an Ecliptical object and Obliquity of the Ecliptic as an Angle object
    convert to an Equatorial object
    @return Equatorial object
    """
    def to_equatorial(self, epsilon_):
        assert isinstance(epsilon_, Angle),'epsilon_ must be Angle'
        lambda_rad, beta_rad = self.lambda_angle.rads, self.beta_angle.rads
        epsilon_rad = epsilon_.rads
        alpha_rad = atan2(sin(lambda_rad)*cos(epsilon_rad)
                          - tan(beta_rad)*sin(epsilon_rad), cos(lambda_rad))
        delta_rad = asin(sin(beta_rad)*cos(epsilon_rad)
                         + cos(beta_rad)*sin(epsilon_rad)*sin(lambda_rad))
        return Equatorial(Longitude(alpha_rad), Latitude(delta_rad))

    def __unicode__(self):
        return u'(\u03BB:' + self.lambda_angle.dms() \
               + u', \u03B2:' + self.beta_angle.dms() + u')'

    def __str__(self):
        return unicode(self).encode(sys.stdout.encoding or DEFAULT_ENCODING,
                                    'replace')

""" Altitude-Azimuth (Local Horizontal) Coordinates """
class AltAzimuthal:
    """
    @param lati_ The observer's geographic latitude as a Latitude object
    @param decl_ The declination of the celestial body as a Latitude object
    @param hangle_ The hour angle as a Longitude object
    """
    def __init__(self, lati_, decl_, hangle_):
        assert isinstance(lati_, Latitude), 'lati_ must be a Latitude'
        assert isinstance(decl_, Latitude), 'decl_ must be a Latitude'
        assert isinstance(hangle_, Longitude), 'hangle_ must be a Longitude'
        lati_rad, decl_rad, hangle_rad = lati_.rads, decl_.rads, hangle_.rads
        azimuth_rad = atan2(sin(hangle_rad), cos(hangle_rad)*sin(lati_rad) \
                            - tan(decl_rad)*cos(lati_rad))
        altitude_rad = asin(sin(lati_rad)*sin(decl_rad) \
                            + cos(lati_rad)*cos(decl_rad)*cos(hangle_rad))
        self.azimuth = Longitude(azimuth_rad)
        self.altitude = Latitude(altitude_rad)

    def __unicode__(self):
        return u'(A:' + self.azimuth.dms() + u', h:' + self.altitude.dms() + u')'

    def __str__(self):
        return unicode(self).encode(sys.stdout.encoding or DEFAULT_ENCODING,
                                    'replace')


if __name__ == "__main__":
    epsilon_j2000 = Angle(radians(28+(1+34.26/60.0)/60.0))
    equa1 = Equatorial(Longitude(0),Latitude(0))
    print equa1
    ecli1 = equa1.to_ecliptical(epsilon_j2000)
    print ecli1
    print ecli1.to_equatorial(epsilon_j2000)
    print

    equa2 = Equatorial(Longitude(pi/2),Latitude(0))
    print equa2
    ecli2 = equa2.to_ecliptical(epsilon_j2000)
    print ecli2
    print ecli2.to_equatorial(epsilon_j2000)
    print

    equa3 = Equatorial(Longitude(pi),Latitude(pi/4))
    print equa3
    ecli3 = equa3.to_ecliptical(epsilon_j2000)
    print ecli3
    print ecli3.to_equatorial(epsilon_j2000)
    print

    equa4 = Equatorial(Longitude(0.75*pi),Latitude(-pi/4))
    print equa4
    ecli4 = equa4.to_ecliptical(epsilon_j2000)
    print ecli4
    print ecli4.to_equatorial(epsilon_j2000)
    print

    altazi = AltAzimuthal(Latitude(radians(38 + 55/60.0 + 17/3600.0)),\
                          Latitude(-radians(6 + 43/60.0 + 11.61/3600.0)),\
                          Longitude(radians(64.352133)))
    print altazi
