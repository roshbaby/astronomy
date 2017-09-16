import sys
from math import sin, cos, tan, asin, atan2, pi, radians
from angle import Angle


DEFAULT_ENCODING = 'UTF-8'

""" Equatorial Coordinates """
class Equatorial:
    # alpha_rad - Equatorial Right Ascension in Radians
    # delta_rad - Equatorial Declination in Radians
    def __init__(self, alpha_rad_in, delta_rad_in):
        self.alpha_rad = alpha_rad_in
        self.delta_rad = delta_rad_in

    def __unicode__(self):
        return u'(\u03B1:' + Angle(self.alpha_rad).hms()   \
               + u', \u03B4:' + Angle(self.delta_rad).dms() + u')'

    def __str__(self):
        return unicode(self).encode(sys.stdout.encoding or DEFAULT_ENCODING,
                                    'replace')


""" Ecliptical Coordinates """
class Ecliptical:
    # lambda_rad - Ecliptical Longitude in Radians
    # beta_rad   - Ecliptical Latitude in Radians
    def __init__(self, lambda_rad_in, beta_rad_in):
        self.lambda_rad = lambda_rad_in
        self.beta_rad = beta_rad_in

    def __unicode__(self):
        return u'(\u03BB:' + Angle(self.lambda_rad).dms() \
               + u', \u03B2:' + Angle(self.beta_rad).dms() + u')'

    def __str__(self):
        return unicode(self).encode(sys.stdout.encoding or DEFAULT_ENCODING,
                                    'replace')


## Equatorial to Ecliptic conversion
## Given an Equatorial object and Obliquity of the Ecliptic in radians
#  convert to Ecliptical
# @return Ecliptical object
def equa2ecli(equa, epsilon_rad):
    alpha_rad, delta_rad = equa.alpha_rad, equa.delta_rad
    lambda_rad = atan2(sin(alpha_rad)*cos(epsilon_rad)
                       + tan(delta_rad)*sin(epsilon_rad),
                      cos(alpha_rad))
    beta_rad   = asin(sin(delta_rad)*cos(epsilon_rad)
                      - cos(delta_rad)*sin(epsilon_rad)*sin(alpha_rad))
    return Ecliptical(lambda_rad, beta_rad)

## Ecliptic to Equatorial conversion
# Given an Ecliptical object and Obliquity of the Ecliptic in radians
# convert to Equatorial
# @return Equatorial object
def ecli2equa(ecli, epsilon_rad):
    lambda_rad, beta_rad = ecli.lambda_rad, ecli.beta_rad
    alpha_rad = atan2(sin(lambda_rad)*cos(epsilon_rad)
                      - tan(beta_rad)*sin(epsilon_rad),
                      cos(lambda_rad))
    delta_rad = asin(sin(beta_rad)*cos(epsilon_rad)
                     + cos(beta_rad)*sin(epsilon_rad)*sin(lambda_rad))
    return Equatorial(alpha_rad, delta_rad)



if __name__ == "__main__":
    epsilon_j2000_rad = radians(28+(01+34.26/60.0)/60.0)
    equa1 = Equatorial(0,0)
    print equa1
    ecli1 = equa2ecli(equa1, epsilon_j2000_rad)
    print ecli1
    print ecli2equa(ecli1, epsilon_j2000_rad)
    print

    equa2 = Equatorial(pi/2,0)
    print equa2
    ecli2 = equa2ecli(equa2, epsilon_j2000_rad)
    print ecli2
    print ecli2equa(ecli2, epsilon_j2000_rad)
    print

    equa3 = Equatorial(pi,pi/4)
    print equa3
    ecli3 = equa2ecli(equa3, epsilon_j2000_rad)
    print ecli3
    print ecli2equa(ecli3, epsilon_j2000_rad)
    print

    equa4 = Equatorial(0.75*pi,-pi/4)
    print equa4
    ecli4 = equa2ecli(equa4, epsilon_j2000_rad)
    print ecli4
    print ecli2equa(ecli4, epsilon_j2000_rad)
