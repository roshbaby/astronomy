from math import pi, degrees, radians, fabs, modf, copysign
from mymath import pfmod
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
        sign = copysign(1,self.rads) # 1 or -1
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
        dgs = degrees(fabs(self.rads))
        mns, dgs = modf(dgs)
        scs, mns = modf(mns*60)
        scs     *= 60
        return self.__get_sign_str()               \
               + u'{0:03d}\u00B0'.format(int(dgs)) \
               + u'{0:02d}\u2032'.format(int(mns)) \
               + u'{0:05.2f}\u2033'.format(scs)
               #+ u'{0:05.2f}'.format(scs).replace(u'.',u'\u2033')

    """ Return unicode string representation in HMS form """
    def hms(self):
        # 1hr equals 15 degs
        hrs = degrees(fabs(self.rads))/15.0
        mns, hrs = modf(hrs)
        scs, mns = modf(mns*60)
        scs     *= 60
        return self.__get_sign_str()           \
               + u'{0:02d}h '.format(int(hrs)) \
               + u'{0:02d}m '.format(int(mns)) \
               + u'{0:06.3f}s'.format(scs)

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
    def __add__(self, other):
        assert isinstance(other, Longitude), 'Unsupported type for addition'
        return Longitude(self.rads + other.rads)

    def __sub__(self,other):
        assert isinstance(other, Longitude), 'Unsupported type for subtraction'
        return Longitude(self.rads - other.rads)

    """ Remove multiples of 2 pi from rads to get a value between 0 and 2 pi"""
    def canonical(self):
        two_pis = 2*pi
        self.rads = pfmod(self.rads,two_pis)
        return self


"""
Angle subclass Latitude always keeps the radians value between -pi/2 and pi/2
Useful for latitudes, declinations, etc.
"""
class Latitude(Angle):
    def __add__(self, other):
        assert isinstance(other, Latitude), 'Unsupported type for addition'
        return Latitude(self.rads + other.rads)

    def __sub__(self,other):
        assert isinstance(other, Latitude), 'Unsupported type for subtraction'
        return Latitude(self.rads - other.rads)

    """ Remap the rads to lie between -pi/2 and +pi/2 """
    def canonical(self):
        # First mod it to lie between 0 and 2*pi
        two_pis = 2*pi
        self.rads = pfmod(self.rads,two_pis)
        # Now map each quadrant so that the final result lies in -pi/2, pi/2
        # 0 - pi/2: no change
        # pi/2 - pi: map to range pi/2, 0
        # pi - 3pi/2: map to range 0, -pi/2
        # 3pi/2 - 2pi: map to range -pi/2, 0
        if 0 <= self.rads <= pi/2:
            pass
        elif self.rads <= pi:
            self.rads = pi - self.rads
        elif self.rads <= 3*pi/2:
            self.rads = -(self.rads - pi)
        else: # 3*pi/2 < self.rads < 2*pi
            self.rads -= two_pis

        # Check that we did every thing right
        assert fabs(self.rads) <= pi/2, 'Latitude should not exceed pi/2 radians either side of 0'
        return self

    """ Override hms to raise exception """
    def hms(self):
        raise RuntimeError('hms not supported for Latitude')


if __name__ == "__main__":
    myangle = Angle(-radians(pi))
    print str(myangle)
    print myangle
    print myangle.rads
    print str(Angle(pi))
    print Angle(pi/2)
    myangle = Angle(pi) + Angle(pi/2)
    print myangle
    myangle += Angle(pi/2)
    print myangle
    myangle = Angle(-2*pi-radians(30))
    print myangle, myangle.dms(), myangle.canonical()
    print myangle.hms()
    print Angle(0)
    print Angle(1)
    print Angle(1.0)
    print repr(Angle(pi))
    longitude = Longitude(-radians(pi))
    print longitude
    try:
        latitude = Latitude(-pi)
    except AssertionError as e:
        print 'Error:', e

    latitude = Latitude(-pi/2)
    print latitude
    try:
        latitude.hms()
    except RuntimeError as e:
        print 'Error:', e

    print latitude + Latitude(-pi/6) # -pi/3 (-60deg)

    try:
        latitude + longitude
    except AssertionError as e:
        print 'Error:', e

    try:
        longitude + latitude
    except AssertionError as e:
        print 'Error:', e

    print Angle(pi/4) + Latitude(-pi/3) # -15deg
    print Angle(pi/3) + Longitude(-pi/3) # +360
