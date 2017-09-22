from math import floor, fmod, modf, fabs
import sys
from types import IntType
from numbers import Number

DEFAULT_ENCODING = 'UTF-8'

weekdays = [u'Sunday', u'Monday', u'Tuesday', u'Wednesday',
            u'Thursday', u'Friday', u'Saturday']

months = [u'Jan', u'Feb', u'Mar', u'Apr', u'May', u'Jun',
          u'Jul', u'Aug', u'Sep', u'Oct', u'Nov', u'Dec']

days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


"""
Determine if the given date falls in the Julian or Gregorian calendar
"""
def is_gregorian(year_, month_, day_):
    return (year_ > 1582) or (year_ == 1582 and month_ > 10) \
           or (year_ == 1582 and month_ == 10 and day_ >= 15)

"""
Determine if the given date falls in a leap year
"""
def is_leap_year(year_, month_, day_):
    assert isinstance(year_, Number)      \
           and isinstance(month_, Number) \
           and isinstance(day_, Number),  \
           'year_, month_, and day_ should be Numbers'

    return (year_ % 4) == 0 and                        \
           ((not is_gregorian(year_, month_, day_)) or \
            (year_ % 100) != 0                      or \
            (year_ % 400) == 0)

class Time:
    def __init__(self, hr_,min_,sec_):
        assert type(hr_) is IntType, 'hour should be integer type'
        assert type(min_) is IntType, 'minute should be integer type'
        assert isinstance(sec_, Number), 'seconds should be a Number'
        self.hr, self.min, self.sec = hr_, min_, sec_
        self.__canonical()

    def __canonical(self):
        """ Reduce Time to a canonical form """
        time_in_hrs = self.hr + (self.min + self.sec/60.0)/60.0
        m, h = modf(time_in_hrs)
        s, m = modf(m*60)
        s   *= 60
        self.hr, self.min, self.sec = int(h), int(m), s

    def __add__(self,other):
        return Time(self.hr+other.hr, self.min+other.min, self.sec+other.sec)

    def __iadd__(self,other):
        self.hr  += other.hr
        self.min += other.min
        self.sec += other.sec
        self.__canonical()
        return self

    def __sub__(self,other):
        return Time(self.hr-other.hr, self.min-other.min, self.sec-other.sec)

    def __isub__(self,other):
        self.hr  -= other.hr
        self.min -= other.min
        self.sec -= other.sec
        self.__canonical()
        return self

    def __unicode__(self):
        return u'{0: 03d}h'.format(self.hr) \
               + u'{0:02d}m'.format(abs(self.min)) \
               + u'{0:05.2f}s'.format(fabs(self.sec))

    def __str__(self):
        return unicode(self).encode(sys.stdout.encoding or DEFAULT_ENCODING,
                                    'replace')


class Date:
    def __init__(self, year_, month_, day_):
        assert type(year_) is IntType, 'year should be integer type'
        assert type(month_) is IntType, 'month should be integer type'
        assert type(day_) is IntType, 'day should be integer type'
        assert month_ > 0 and month_ < 13, \
               '{0:d} is not a valid month'.format(month_)
        assert day_ >= 0, 'Day should be non negative'
        max_days = days[month_-1]
        # Leap Year Correction for February
        if is_leap_year(year_, month_, day_) and month_ == 2: max_days += 1
        assert day_ <= max_days, '{0:d} is invalid value for day'.format(day_)
        self.year = year_
        self.month = month_
        self.day = day_

    """
    Check if the Date falls in the Gregorian calendar
    (i.e. 15 Oct 1582 or later)
    """
    def is_gregorian(self):
        return is_gregorian(self.year, self.month, self.day)

    """ Day of the week for the Date """
    def weekday(self):
        myjdn = JulianDayNumber(self, Time(0,0,0)).jdn # Force 0h UT for the Date
        myjdn = int(fmod(myjdn+1.5, 7))
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
            A = int(floor(year/100.0))
            B = 2 - A + int(floor(A/4.0))
        # JDN at 0h UT
        self.jdn = int(floor(365.25*(year+4716))) \
                   + int(floor(30.6001*(month+1))) + day + B - 1524.5
        # JDN at given UT
        self.jdn += (hrs + (mins + secs/60.0)/60.0)/24.0

    """
    Get the Calendar Date corresponding to the Julian Day Number
    @caution Not valid for negative JDNs
    """
    def get_date(self):
        jdn = self.jdn + 0.5
        F,Z = modf(jdn) # Fractional, Integer parts
        A = Z # Initial value
        if Z >= 2299161: # Update A
            alpha = int(floor((Z-1867216.25)/36524.25))
            A += 1 + alpha - int(floor(alpha/4.0))
        B = A + 1524
        C = int(floor((B-122.1)/365.25))
        D = int(floor(365.25*C))
        E = int(floor((B-D)/30.6001))
        day = B - D - int(floor(30.6001*E)) + F
        month = E - 1 # default
        if E >= 14: month -= 12 # correction
        year = C - 4715 # default
        if month > 2: year -= 1 # correction
        return Date(year, month, int(day))

    """
    Get the Time of day corresponding to the Julian Day Number
    """
    def get_time(self):
        F,Z = modf(self.jdn-0.5) # Fractional, Integer parts
        M,H = modf(F*24)
        S,M = modf(M*60)
        S  *= 60
        return Time(int(H),int(M),S)

    """
    Compute the Modified JDN
    MJDN = 0.0 corresponds to 1858 Nov 17 0h UT
    """
    def get_mjdn(self):
        return self.jdn - 2400000.5

    def __unicode__(self):
        return u'{0:f}'.format(self.jdn)

    def __str__(self):
        return unicode(self).encode(sys.stdout.encoding or DEFAULT_ENCODING,
                                    'replace')


""" Julian Day Number for an 'epoch' (i.e. Jan 0.0 UT of given year) """
def jdn_epoch(year):
    return JulianDayNumber(year-1,12,31,0,0,0)


if __name__ == "__main__":
    t = Time(23,59,59)
    print t
    print t + Time(1,-59,-59)
    print t + Time(-1,59,61)
    t -= Time(-1,59,61)
    print t
    t += Time(1,-59,-61)
    print t
    t += t
    print t
    print Time(-23,-59,-59)

    try:
        d = Date(2017,9,31)
    except AssertionError as e:
        print 'Error:', e

    print Date(2016,2,29)
    print Date(2000,2,29)

    try:
        print Date(1900,2,29)
    except AssertionError as e:
        print 'Error:', e

    print Date(1580, 2, 29) #
    print Date(1500, 2, 29) # Leap Year in Julian calendar only
    print Date(1500, 2, 29).is_gregorian()

    d = Date(2017,9,15)
    print d
    print d.weekday()
    print d.is_gregorian() # Call class method
    print Date(-4712,11,12)
    jdn = JulianDayNumber(Date(1975,6,10),Time(8,18,00))
    print jdn.get_date()
    print jdn.get_date().weekday()
    print jdn.get_time()

    print JulianDayNumber(Date(2000,1,1), Time(12,0,0)) #JDE2000
