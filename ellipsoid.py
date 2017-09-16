import math
from angle import Angle
from numbers import Number

class Ellipsoid:
    """Ellipsoid with Functions to convert from Geocentric to Geographic Coordinates"""
    """
      @param a The Equatorial Radius for the Ellipsoid in meters
      @param f The flattening
    """
    def __init__(self, a, f):
        assert isinstance(a, Number), 'a should be number'
        assert isinstance(f, Number), 'f should be number'
        self.equ_radius = a
        self.flattening = f
        self.pol_radius = a*(1-f)
        self.eccentricity = math.sqrt((2-f)*f)

    """
      @param geog_lat_rad The Geographic Latitude in radians
      @param height_m The height above surface in meters
      @return Tuple of (rho*sin(phi'), rho*cos(phi'))
    """
    def geocentric(self, geog_lat_rad, height_m):
        factor = (1-self.flattening)
        height_rel = height_m/self.equ_radius
        u = math.atan(factor*math.tan(geog_lat_rad))
        rho_sin_phi = factor*math.sin(u) + height_rel*math.sin(geog_lat_rad)
        rho_cos_phi = math.cos(u) + height_rel*math.cos(geog_lat_rad)
        return (rho_sin_phi, rho_cos_phi)

    def geocentric_rho_phi(self, geog_lat_rad, height_m):
        rho_sin_phi, rho_cos_phi = self.geocentric(geog_lat_rad,height_m)
        rho_rel = math.sqrt(rho_sin_phi*rho_sin_phi + rho_cos_phi*rho_cos_phi)
        phi_rad = math.atan2(rho_sin_phi, rho_cos_phi)
        return (rho_rel,phi_rad)


if __name__ == "__main__":
    earth = Ellipsoid(6378140, 1/298.257)
    geog_lat_deg = round(33.0+(21.0+22/60.0)/60.0,4) # Round to nearest arc-second
    geog_lat_rad = math.radians(geog_lat_deg)
    height_m = 1706.0

    print 'Geog. Lat:\t', geog_lat_deg, 'degrees'

    rho_sin_phi, rho_cos_phi = earth.geocentric(geog_lat_rad, height_m)
    rho_rel, phi_rad = earth.geocentric_rho_phi(geog_lat_rad, height_m)

    print 'rho_sin_phi:\t', round(rho_sin_phi,6)
    print 'rho_cos_phi:\t', round(rho_cos_phi,6)
    print 'rho_m:\t\t', round(rho_rel*earth.equ_radius,2), 'metres'
    print 'phi_rad:\t', phi_rad, '(' + str(Angle(phi_rad)) + ')'
    print earth
