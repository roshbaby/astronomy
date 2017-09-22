from math import pi, degrees, radians, sin, cos, atan2, asin, fabs
from mymath import pfmod
from angle import Angle, Longitude, Latitude
from coordinates import Ecliptical, Equatorial
from calendar import Date, Time, JulianDayNumber
from constants import epoch_j2000, julian_century
from ecliptic import Ecliptic

class Sun:
    def __init__(self, jdn_):
        assert isinstance(jdn_, JulianDayNumber), \
               'jdn_ should be a JulianDayNumber'

        two_pi = 2*pi

        T = (jdn_.jdn - epoch_j2000.jdn)/julian_century

        # Mean Longitude referred to the mean equinox
        L_0 = 280.46646 + (36000.76983 + 0.0003032*T)*T
        L_0 = pfmod(L_0, 360)
        L_0 = radians(L_0)

        # Mean anomaly
        M = 357.52911 + (35999.05029 - 0.0001537*T)*T
        M = pfmod(M, 360)
        M = radians(M)

        # Eccentricity of the Earth's orbit
        e = 0.016708634 - (0.000042037 + 0.0000001267*T)*T

        # Equation of the center
        C = (1.914602 - (0.004817 + 0.000014*T)*T)*sin(M) \
            + (0.019993 - 0.000101*T)*sin(2*M)            \
            + 0.000289*sin(3*M)
        C = pfmod(C, 360)
        C = radians(C)

        # True geometric Longitude referred to the mean equinox of the date
        self.L = pfmod(L_0 + C, two_pi)

        # True geometric latitude referred to the mean equinox
        self.beta = 0

        # True anomaly
        self.v = pfmod(M + C, two_pi)

        # Radius vector
        self.R = 1.000001018*(1-e**2)/(1 + e*cos(self.v))

        # Longitude of ascending node (see more accurate expression in ecliptic)
        Omega = radians(125.04 - 1934.136*T)

        # Apparent Longitude referred to the *true* equinox of the date
        self.lambda_apparent = degrees(self.L) - 0.00569 - 0.00478*sin(Omega)
        self.lambda_apparent = pfmod(self.lambda_apparent, 360)
        self.lambda_apparent = radians(self.lambda_apparent)

        # Apparent Latitude referred to the *true* equinox of the date
        self.beta_apparent = 0

        # Equatorial coordinates referred to the mean equinox of the date
        ecliptic = Ecliptic(jdn_)
        w_uncorrected = ecliptic.get_obliquity_uncorrected().rads

        self.alpha = atan2(cos(w_uncorrected)*sin(self.L), cos(self.L))
        self.alpha = pfmod(self.alpha, two_pi)

        self.delta = asin(sin(w_uncorrected)*sin(self.L))
        self.delta = pfmod(self.delta, two_pi)

        # Equatorial coordinates referred to the *true* equinox of the date
        w_corrected = w_uncorrected + radians(0.00256*cos(Omega))

        self.alpha_apparent = \
            atan2(cos(w_corrected)*sin(self.lambda_apparent), cos(self.L))
        self.alpha_apparent = pfmod(self.alpha_apparent, two_pi)

        self.delta_apparent = asin(sin(w_corrected)*sin(self.lambda_apparent))
        self.delta_apparent = pfmod(self.delta_apparent, two_pi)

        # Rectangular Coordinates referred to the mean equinox of the date
        self.X = self.R*cos(self.beta)*cos(self.L)
        self.Y = self.R*(cos(self.beta)*sin(self.L)*cos(w_uncorrected) \
                         - sin(self.beta)*sin(w_uncorrected))
        self.Z = self.R*(cos(self.beta)*sin(self.L)*sin(w_uncorrected) \
                         + sin(self.beta)*cos(w_uncorrected))


    def get_ecliptical(self):
        return Ecliptical(Longitude(self.L), Latitude(0))

    def get_ecliptical_apparent(self):
        return Ecliptical(Longitude(self.lambda_apparent), Latitude(0))

    def get_equatorial(self):
        return Equatorial(Longitude(self.alpha), Latitude(self.delta))

    def get_equatorial_apparent(self):
        return Equatorial(Longitude(self.alpha_apparent),
                          Latitude(self.delta_apparent))

    def get_rectangular(self):
        return (self.X, self.Y, self.Z)


########
# Table of equations for use with equinox and solstice computation
# for the years before +1000
# Note: Y = year/1000
########
spring_table27a = \
    lambda Y: (((-0.00071*Y +0.00111)*Y +0.06134)*Y +365242.13740)*Y \
              + 1721139.29189
summer_table27a = \
    lambda Y: (((0.00025*Y +0.00907)*Y -0.05323)*Y +365241.72562)*Y \
              + 1721233.25401
autumn_table27a = \
    lambda Y: (((0.00074*Y -0.00297)*Y -0.11677)*Y +365242.49558)*Y \
              + 1721325.70455
winter_table27a = \
    lambda Y: (((-0.00006*Y -0.00933)*Y -0.00769)*Y +365242.88257)*Y \
              + 1721414.39987
table27a = (
    spring_table27a, # Spring Equinox
    summer_table27a, # Summer Solstice
    autumn_table27a, # Autumn Equinox
    winter_table27a  # Winter Solstice
)

########
# Table of equations for use with equinox and solstice computation
# for the years after +1000
# Note: Y = (year-2000)/1000
########
spring_table27b = \
    lambda Y: (((-0.00057*Y -0.00411)*Y +0.05169)*Y +365242.37404)*Y \
              + 2451623.80984
summer_table27b = \
    lambda Y: (((-0.00030*Y +0.00888)*Y + 0.00325)*Y +365241.62603)*Y \
              + 2451716.56767
autumn_table27b = \
    lambda Y: (((0.00078*Y +0.00337)*Y -0.11575)*Y +365242.01767)*Y \
              + 2451810.21715
winter_table27b = \
    lambda Y: (((0.00032*Y -0.00823)*Y -0.06223)*Y +365242.74049)*Y \
              + 2451900.05952
table27b = (
    spring_table27b, # Spring Equinox
    summer_table27b, # Summer Solstice
    autumn_table27b, # Autumn Equinox
    winter_table27b  # Winter Solstice
)

########
# Table of coefficients for use with equinox and solstice computation
########
table27c = (
    (485, 324.96,   1934.136),
    (203, 337.23,  32964.467),
    (199, 342.08,     20.186),
    (182,  27.85, 445267.112),
    (156,  73.14,  45036.886),
    (136, 171.52,  22518.443),
    ( 77, 222.54,  65928.934),
    ( 74, 296.72,   3034.906),
    ( 70, 243.58,   9037.513),
    ( 58, 119.81,  33718.147),
    ( 52, 297.17,    150.678),
    ( 50,  21.02,   2281.226),
    ( 45, 247.54,  29929.562),
    ( 44, 325.15,  31555.956),
    ( 29,  60.93,   4443.417),
    ( 18, 155.12,  67555.328),
    ( 17, 288.79,   4562.452),
    ( 16, 198.04,  62894.029),
    ( 14, 199.76,  31436.921),
    ( 12,  95.39,  14577.848),
    ( 12, 287.11,  31931.756),
    ( 12, 320.81,  34777.259),
    (  9, 227.73,   1222.114),
    (  8,  15.45,  16859.074)
)

def mean_equinox_solstice(year_, idx_):
    if int(year_) <= 1000:
        Y = int(year_)/1000.0
        jde_0 = table27a[idx_](Y)
    else:
        Y = (int(year_)-2000)/1000.0
        jde_0 = table27b[idx_](Y)
    return jde_0

def equinox_solstice_helper(year_, idx_):
    jde_0 = mean_equinox_solstice(year_, idx_)

    #print 'Y:', Y
    #print 'jde_0:', jde_0
    #print 'jde2000', epoch_j2000.jdn

    T = (jde_0 - epoch_j2000.jdn)/julian_century
    W = radians(35999.373*T - 2.47)
    dlambda = 1 + 0.0334*cos(W) + 0.0007*cos(2*W)

    #print 'T:', T
    #print 'dlambda:', dlambda

    S = 0.0
    for (A,B,C) in table27c:
        S += A*cos(radians(B + C*T))

    #print 'S', S

    jde = jde_0 + (0.00001*S)/dlambda

    #print 'jde:', jde

    jdn_to_return = JulianDayNumber(Date(year_,1,1), Time(0,0,0))
    jdn_to_return.jdn = jde # Roundabout way of creating a JDN
    return jdn_to_return

########
# @todo - The iterative method is not accurate until we have the 'higher
# accuracy' calculations for the sun's circumstances
########
def equinox_solstice_iterative(year_, idx_, prec_):
    jde_0 = mean_equinox_solstice(year_, idx_)
    jdn_to_ret = JulianDayNumber(Date(year_,1,1), Time(0,0,0))
    jdn_to_ret.jdn = jde_0 # Hack
    while True:
        print 'jdn', jdn_to_ret.jdn
        sun = Sun(jdn_to_ret)
        print 'L:', degrees(sun.L)
        print 'l_app', degrees(sun.lambda_apparent)
        print 'R:', sun.R
        dpsi = Ecliptic(jdn_to_ret).get_dpsi_low().rads
        fk5_correction = -radians(0.09033/3600.0)
        aberration = -radians(20.4898/3600.0/sun.R)
        #print Angle(dpsi), Angle(fk5_correction), Angle(aberration)
        lambda_app = sun.L + dpsi + fk5_correction + aberration
        print 'lambda_app', degrees(lambda_app)%360
        dday = 58*sin(idx_*pi/2 - lambda_app)
        print 'dday', dday
        if fabs(dday) <= fabs(prec_): break
        jdn_to_ret.jdn += dday
    return jdn_to_ret

def spring_equinox(year_):
    """ Julian Day Number of the Spring Equinox in the given year"""
    return equinox_solstice_helper(year_,0)

def summer_solstice(year_):
    """ Julian Day Number of the Summer Solstice in the given year"""
    return equinox_solstice_helper(year_,1)

def autumn_equinox(year_):
    """ Julian Day Number of the Autumn Equinox in the given year"""
    return equinox_solstice_helper(year_,2)

def winter_solstice(year_):
    """ Julian Day Number of the Winter Solstice in the given year"""
    return equinox_solstice_helper(year_,3)


if __name__ == "__main__":
    # Sun tests
    sooraj = Sun(JulianDayNumber(Date(1992,10,13),Time(0,0,0)))
    print sooraj.get_ecliptical() # 199deg54'36", 0
    print sooraj.get_ecliptical_apparent() # 199deg54'32", 0
    print sooraj.get_equatorial()
    print sooraj.get_equatorial_apparent() # 13h13m31.4s, -7deg47'06"
    print sooraj.get_rectangular()
    print

    # Equinox/Solstice tests
    years = [ 1962, 1975, 1979, 2013, 2017, 2018, 2019, 2020 ]
    prec = 1e-6
    for year in years:
        print 'Year', year
        instant = spring_equinox(year)
        print 'Spring Equinox :', instant.get_date(), instant.get_time()
        instant = summer_solstice(year)
        print 'Summer Solstice:', instant.get_date(), instant.get_time()
        instant = autumn_equinox(year)
        print 'Autumn Equinox :', instant.get_date(), instant.get_time()
        instant = winter_solstice(year)
        print 'Winter Solstice:', instant.get_date(), instant.get_time()
        print
