from math import radians
from angle import Angle
from calendar import Date, Time, JulianDayNumber

julian_year = 365.25
julian_century = 36525.0

epoch_j2000 = JulianDayNumber(Date(2000,1,1),Time(12,0,0))  # 2451545.0 TD
epsilon_j2000 = Angle(radians(28+(1+34.26/60.0)/60.0))

if __name__ == "__main__":
    print epoch_j2000
    print epsilon_j2000
