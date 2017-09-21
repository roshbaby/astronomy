from math import sin, cos, radians
from angle import Angle
from constants import epoch_j2000, julian_century
from calendar import Date, Time, JulianDayNumber

class Ecliptic:
    ########
    # @param jdn The Julian Day Number (Julian Ephemeris Day)
    ########
    def __init__(self, jdn_):
        assert isinstance(jdn_, JulianDayNumber), \
               'jdn_ should be a JulianDayNumber'
        T = (jdn_.jdn - epoch_j2000.jdn) / julian_century
        T_cube = T**3

        # @todo - Implement the IAU table accurate to 0.0003"

        ## Mean elongation of the Moon from Sun
        D = 297.85036 + (445267.111480 - 0.0019142*T)*T + T_cube/189474.0
        #print 'D:', D, D%360
        D %= 360
        D = radians(D)

        ## Mean anomaly of the Earth
        M_earth = 357.52772 + (35999.050340 - 0.0001603*T)*T - T_cube/3000000
        #print 'M_earth:', M_earth, M_earth%360
        M_earth %= 360
        M_earth = radians(M_earth)

        ## Mean anomaly of the Moon
        M_moon = 134.96298 + (477198.867398 + 0.0086972*T)*T + T_cube/56250
        #print 'M_moon:', M_moon, M_moon%360
        M_moon %= 360
        M_moon = radians(M_moon)

        ## Moon's argument of latitude
        F = 93.271191 + (483202.017538 - 0.0036825*T)*T + T_cube/327270
        #print 'F:', F, F%360
        F %= 360
        F = radians(F)

        ## Longitude of ascending node
        Omega = 125.04452 - (1934.136261 -  0.0020707*T)*T + T_cube/450000
        #print 'Omega:', Omega, Omega%360
        Omega %= 360
        Omega = radians(Omega)

        # Accurate to 0.5" in dPsi, 0.1" in dEpsilon
        L_sun = (280.4665 + 36000.7698*T) % 360
        L_moon = (218.3165 + 481267.8813*T) % 360
        #print 'L_sun:', L_sun, 'L_moon:', L_moon
        L_sun = radians(L_sun)
        L_moon = radians(L_moon)

        d_psi = -17.20*sin(Omega) - 1.32*sin(2*L_sun) \
                - 0.23*sin(2*L_moon) + 0.21*sin(2*Omega) # arc-sec
        self.d_psi_low = Angle(radians(d_psi/3600.0))

        d_epsilon = 9.20*cos(Omega) + 0.57*cos(2*L_sun) \
                    + 0.10*cos(2*L_moon) - 0.09*cos(2*Omega) #arc-sec
        self.d_epsilon_low = Angle(radians(d_epsilon/3600.0))

        ## Obliquity (valid only for |U|<1)
        U = T/100.0
        epsilon = (23+(26+21.448/60.0)/60.0)
        correction = (((((((((2.45*U + 5.79)*U +27.87)*U + 7.12)*U -39.05)*U - 249.67)*U -51.38)*U +1999.25)*U -1.55)*U -4680.93)*U # arc-sec
        epsilon += (correction/3600.0)
        self.epsilon_0 = Angle(radians(epsilon))

    ########
    # dPsi accurate to 0.5" as an Angle object
    ########
    def get_dpsi_low(self):
        """ dPsi accurate to 0.5" as an Angle object"""
        return self.d_psi_low

    ########
    # dEpsilon accurate to 0.1" as an Angle object
    ########
    def get_depsilon_low(self):
        """ dEpsilon accurate to 0.1" as an Angle object"""
        return self.d_epsilon_low

    ########
    # Obliquity NOT corrected for nutation
    ########
    def get_obliquity_uncorrected(self):
        """Obliquity NOT corrected for nutation as an Angle object"""
        return self.epsilon_0

    ########
    # Obliquity corrected for nutation
    ########
    def get_obliquity_corrected(self):
        """Obliquity corrected for nutation as an Angle object"""
        # @todo Use the more accurate value for d_epsilon
        return self.epsilon_0 + self.d_epsilon_low


if __name__ == "__main__":
    jdn = JulianDayNumber(Date(1987,4,10), Time(0,0,0))
    ecliptic = Ecliptic(jdn)

    print ecliptic.get_dpsi_low()     # -3.9"
    print ecliptic.get_depsilon_low() # +9.4"
    print ecliptic.get_obliquity_uncorrected() # 23deg 26'27.407"
    print ecliptic.get_obliquity_corrected()   # 23def 26'36.850"
    print

    ecliptic_j2k = Ecliptic(epoch_j2000)
    print 'Ecliptic J2000:', ecliptic_j2k.get_obliquity_uncorrected()
    print 'Ecliptic J2000 (corrected):', ecliptic_j2k.get_obliquity_corrected()
