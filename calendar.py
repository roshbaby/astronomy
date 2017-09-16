import math
import sys
from types import *

DEFAULT_ENCODING = 'UTF-8'

weekdays = [u'Sunday', u'Monday', u'Tuesday', u'Wednesday',
            u'Thursday', u'Friday', u'Saturday']

months = [u'Jan', u'Feb', u'Mar', u'Apr', u'May', u'Jun',
          u'Jul', u'Aug', u'Sep', u'Oct', u'Nov', u'Dec']

class Time:
    def __init__(self, hr_,min_,sec_):
        assert type(hr_) is IntType, 'hour should be integer type'
        assert type(min_) is IntType, 'minute should be integer type'
        assert type(sec_) is IntType or type(sec_) is FloatType,\
               'seconds should be integer/float type'
        assert hr_ >= 0 and hr_ < 24,\
               'value {0:d} is invalid for hours'.format(hr_)
        assert min_ >= 0 and min_ < 60,\
               'value {0:d} is invalid for minutes'.format(min_)
        assert sec_ >= 0 and sec_ < 60.0,\
               'value {0:d} is invalid for seconds'.format(sec_)
        self.hr = hr_
        self.min = min_
        self.sec = sec_

    def __unicode__(self):
        return u'{0:02d}h'.format(self.hr) \
               + u'{0:02d}m'.format(self.min) \
               + u'{0:05.2f}s'.format(self.sec)

    def __str__(self):
        return unicode(self).encode(sys.stdout.encoding or DEFAULT_ENCODING,
                                    'replace')


class Date:
    def __init__(self, year_, month_, day_):
        assert type(year_) is IntType, 'year should be integer type'
        assert type(month_) is IntType, 'month should be integer type'
        assert type(day_) is IntType, 'day should be integer type'
        assert month_ > 0 and month_ < 13, \
               '{0:d} is not valid month'.format(month_)
        assert day_ >= 0, 'Day values should be zero/positive'
        self.year = year_
        self.month = month_
        self.day = day_

    """
    Check if the Date falls in the Gregorian calendar
    (i.e. 15 Oct 1582 or later)
    """
    def is_gregorian(self):
        return (self.year > 1582) or (self.year == 1582 and self.month > 10) \
               or (self.year == 1582 and self.month == 10 and self.day >= 15)

    """ Day of the week for the Date """
    def weekday(self):
        myjdn = JulianDayNumber(self, Time(0,0,0)).jdn # Force 0h UT for the Date
        myjdn = int(math.fmod(myjdn+1.5, 7))
        if myjdn < 0: myjdn += 7
        return weekdays[myjdn]

    def __unicode__(self):
        return u'{0: 04d} '.format(self.year) \
               + months[self.month-1]         \
               + u' {0:02d}'.format(self.day)

    def __str__(self):
        return unicode(self).encode(sys.stdout.encoding or DEFAULT_ENCODING,
                                    'replace')


class JulianDayNumber:
    """ Julian Day Number for given UT """
    def __init__(self, date, time):
        self.date = date
        self.time = time
        year, month, day = date.year, date.month, date.day
        hrs, mins, secs = time.hr, time.min, time.sec
        if month < 3:
            year -= 1
            month += 12
        B = 0 # Default for Julian Calendar
        if date.is_gregorian():
            A = int(math.floor(year/100.0))
            B = 2 - A + int(math.floor(A/4.0))
        # JDN at 0h UT
        self.jdn = int(math.floor(365.25*(year+4716))) \
                   + int(math.floor(30.6001*(month+1))) + day + B - 1524.5
        # JDN at given UT
        self.jdn += (hrs + (mins + secs/60.0)/60.0)/24.0

    """
    Convert the Julian Day Number to Calendar Date
    @caution Not valid for negative JDNs
    """
    def to_date(self):
        jdn = self.jdn + 0.5
        F,Z = math.modf(jdn) # Fractional, Integer parts
        A = Z # Initial value
        if Z >= 2299161: # Update A
            alpha = int(math.floor((Z-1867216.25)/36524.25))
            A += 1 + alpha - int(math.floor(alpha/4.0))
        B = A + 1524
        C = int(math.floor((B-122.1)/365.25))
        D = int(math.floor(365.25*C))
        E = int(math.floor((B-D)/30.6001))
        day = B - D - int(math.floor(30.6001*E)) + F
        month = E - 1 # default
        if E >= 14: month -= 12 # correction
        year = C - 4715 # default
        if month > 2: year -= 1 # correction
        return Date(year, month, int(day))

    """
    Compute the Modified JDN
    MJDN = 0.0 corresponds to 1858 Nov 17 0h UT
    """
    def get_mjdn(self):
        return self.jdn - 2400000.5


## Julian Day Number for an 'epoch' (i.e. Jan 0.0 UT of given year)
def jdn_epoch_ut(year):
    return JulianDayNumber(year-1,12,31,0,0,0)


if __name__ == "__main__":
    #t = Time(24,0,0)
    #t = Time(23,60,0)
    #t = Time(23,59,60)
    t = Time(23,59,59)
    d = Date(2017,9,15)
    print d
    print d.weekday()
    print Date(-4712,11,12)
    jdn = JulianDayNumber(Date(2015,9,15),Time(1,2,3))
    print jdn.to_date()
    print jdn.to_date().weekday()
    print
