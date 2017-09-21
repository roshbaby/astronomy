from math import radians
from angle import Angle
from calendar import Date, Time, JulianDayNumber

julian_year = 365.25
julian_century = 36525.0

epoch_j2000 = JulianDayNumber(Date(2000,1,1),Time(12,0,0))  # 2451545.0 TD

# Obliquity of the Ecliptic referred to J2000
# To use with apparent RA and declination, use epsilon_j2000_corrected which
# includes nutation
epsilon_j2000 = Angle( radians(23+(26+21.448/60.0)/60.0) )
epsilon_j2000_corrected = Angle(radians(23+(26+15.69/60.0)/60.0))

if __name__ == "__main__":
    print 'J2000:', epoch_j2000
    print 'epsilon:', epsilon_j2000
    print 'epsilon (corrected):', epsilon_j2000_corrected
