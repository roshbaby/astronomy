from math import radians
from angle import Longitude
from calendar import Date, Time, JulianDayNumber
from constants import epoch_j2000, julian_century

"""
Mean Sidereal Time at Greenwich for a given UT
@return Greenwich MST as a Longitude object
"""
def MST(jdnumber):
    assert isinstance(jdnumber, JulianDayNumber), 'Invalid Julian Day Number'
    # Number of Julian centuries since J2000
    T = (jdnumber.jdn - epoch_j2000.jdn)/julian_century
    ret_deg = 280.46061837                                       \
              + 360.98564736629*(jdnumber.jdn - epoch_j2000.jdn) \
              + (0.000387933 - T/38710000.0)*T*T
    return Longitude(radians(ret_deg))

"""
Apparent Sidereal Time at Greenwich for a given UT
@return Greenwich AST as a Longitude object
"""
def AST(jdnumber):
    assert isinstance(jdnumber, JulianDayNumber), 'Invalid Julian Day Number'
    mst = MST(jdnumber)
    raise RuntimeError('@todo AST needs correction for nutation')


if __name__ == "__main__":
    jdn = JulianDayNumber(Date(2017,9,15),Time(22,20,11))
    theta_angle = MST(jdn)
    print theta_angle
    print MST(JulianDayNumber(Date(1987,4,10),Time(19,21,0))).hms()
    try:
        AST(jdn)
    except RuntimeError as e:
        print 'Error:', e
