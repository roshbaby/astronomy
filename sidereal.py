from angle import *
from calendar import *

"""
Mean Sidereal Time at Greenwich for a given UT
@return Greenwich MST as a Longitude object
"""
def theta_0(jdnumber):
    assert isinstance(jdnumber, JulianDayNumber), 'Invalid Julian Day Number'
    T = (jdnumber.jdn - 2451545.0)/36525.0 # Number of Julian centuries
    ret_deg = 280.46061837 + 360.98564736629*(jdnumber.jdn - 2451545.0) + \
              (0.000387933 - T/38710000.0)*T*T
    ret_deg = math.fmod(ret_deg, 360) # fmod preferred over % operator for floats
    if ret_deg < 0: ret_deg += 360
    return Longitude(math.radians(ret_deg))

if __name__ == "__main__":
    jdn = JulianDayNumber(Date(2017,9,15),Time(22,20,11))
    theta_angle = theta_0(jdn)
    print theta_angle
    print theta_0(JulianDayNumber(Date(1987,4,10),Time(19,21,0))).hms()
