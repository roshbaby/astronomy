import math
import sys
from numbers import Number

DEFAULT_ENCODING = 'UTF-8'

"""
Angle class that can accumulate an angle value in radians
"""
class Angle:
    """ Create an Angle in radians """
    def __init__(self, val_in_rad):
        assert isinstance(val_in_rad, Number), 'Radians should be a Number'
        self.rads = val_in_rad
        self.canonical()

    def __add__(self, other):
        assert isinstance(other, Angle), 'Unsupported type for addition'
        return Angle(self.rads + other.rads)

    def __iadd__(self,other):
        assert isinstance(other, Angle), 'Unsupported type for i-addition'
        self.rads += other.rads
        self.canonical()
        return self

    def __sub__(self,other):
        assert isinstance(other, Angle), 'Unsupported type for subtraction'
        return Angle(self.rads - other.rads)

    def __isub__(self,other):
        assert isinstance(other, Angle), 'Unsupported type for i-subtraction'
        self.rads -= other.rads
        self.canonical()
        return self

    def __get_sign_str(self):
        sign = math.copysign(1,self.rads) # 1 or -1
        if (self.rads == 0): signstr = u' '
        elif (sign < 0):     signstr = u'-'
        else:                signstr = u'+'
        return signstr

    """ Convert the rads to a canonical form """
    def canonical(self):
        # Do nothing for the Angle class
        return self

    """ Return unicode string representation in DMS form """
    def dms(self):
        dgs = math.degrees(math.fabs(self.rads))
        mns,dgs = math.modf(dgs)
        mns *= 60
        scs, mns = math.modf(mns)
        scs *= 60
        return self.__get_sign_str()               \
               + u'{0:03d}\u00B0'.format(int(dgs)) \
               + u'{0:02d}\u2032'.format(int(mns)) \
               + u'{0:05.2f}\u2033'.format(scs)
               #+ u'{0:05.2f}'.format(scs).replace(u'.',u'\u2033')

    """ Return unicode string representation in HMS form """
    def hms(self):
        # 1hr equals 15 degs
        hrs = math.degrees(math.fabs(self.rads))/15.0
        mns, hrs = math.modf(hrs)
        mns *= 60
        scs, mns = math.modf(mns)
        scs *= 60
        return self.__get_sign_str()           \
               + u'{0:02d}h '.format(int(hrs)) \
               + u'{0:02d}m '.format(int(mns)) \
               + u'{0:05.2f}s'.format(scs)

    def __unicode__(self):
        return self.dms()

    def __str__(self):
        return unicode(self).encode(sys.stdout.encoding or DEFAULT_ENCODING,
                                    'replace')


"""
Angle subclass Longitude always keeps the radians value between 0 and 2pi
Useful for longitudes, Right Ascensions etc.
"""
class Longitude(Angle):
    """ Remove multiples of 2 pi from rads to get a value between 0 and 2 pi"""
    def canonical(self):
        two_pis = math.pi*2
        self.rads = math.fmod(self.rads,two_pis)
        if self.rads < 0:
            self.rads += two_pis
        return self


"""
Angle subclass Latitude always keeps the radians value between -pi and pi
Useful for latitudes, declinations, etc.
"""
class Latitude(Angle):
    """ Just check that the rads lie between -pi/2 and +pi/2 """
    def canonical(self):
        assert math.fabs(self.rads) <= math.pi/2, 'Latitude cannot exceed pi/2 radians either side of 0'
        return self

    """ Override hms to raise exception """
    def hms(self):
        raise RuntimeError('hms not supported for Latitude')


if __name__ == "__main__":
    myangle = Angle(-math.radians(math.pi))
    print str(myangle)
    print myangle
    print myangle.rads
    print str(Angle(math.pi))
    print Angle(math.pi/2)
    myangle = Angle(math.pi) + Angle(math.pi/2)
    print myangle
    myangle += Angle(math.pi/2)
    print myangle
    myangle = Angle(-2*math.pi-math.radians(30))
    print myangle, myangle.dms(), myangle.canonical()
    print myangle.hms()
    print Angle(0)
    print Angle(1)
    print Angle(1.0)
    print repr(Angle(math.pi))
    longitude = Longitude(-math.radians(math.pi))
    print longitude
    try:
        latitude = Latitude(-math.pi)
    except AssertionError as e:
        print 'Error:', e

    latitude = Latitude(-math.pi/2)
    print latitude
    try:
        latitude.hms()
    except RuntimeError as e:
        print 'Error:', e
