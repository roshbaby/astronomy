import math
import sys
from types import *

DEFAULT_ENCODING = 'UTF-8'

class Angle:
    """ Create an Angle in radians """
    def __init__(self, val_in_rad):
        assert type(val_in_rad) is IntType or type(val_in_rad) is FloatType,\
               'Value in radians should be integer/float type'
        self.rads = val_in_rad

    def __add__(self, other):
        assert isinstance(other, Angle), 'Unsupported type for addition'
        return Angle(self.rads + other.rads)

    def __iadd__(self,other):
        assert isinstance(other, Angle), 'Unsupported type for i-addition'
        self.rads += other.rads
        return self

    def __sub__(self,other):
        assert isinstance(other, Angle), 'Unsupported type for subtraction'
        return Angle(self.rads - other.rads)

    def __isub__(self,other):
        assert isinstance(other, Angle), 'Unsupported type for i-subtraction'
        self.rads -= other.rads
        return self

    def __get_sign_str(self):
        sign = math.copysign(1,self.rads) # 1 or -1
        if (self.rads == 0): return u' '
        elif (sign < 0):     return u'-'
        else:                return u'+'

    """ Remove multiples of 2 pi to get a value between 0 and 2 pi"""
    def canonical(self):
        two_pis = math.pi*2
        self.rads = math.fmod(self.rads,two_pis)
        if self.rads < 0:
            self.rads += two_pis
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
               + u'{0:05.2f}'.format(scs).replace(u'.',u'\u2033')

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


if __name__ == "__main__":
    myangle = Angle(-math.radians(math.pi))
    print str(myangle)
    print myangle
    print myangle.rads
    print str(Angle(math.pi))
    print str(Angle(math.pi/2))
    myangle = Angle(math.pi) + Angle(math.pi/2)
    print myangle
    myangle += Angle(math.pi/2)
    print myangle
    myangle = Angle(-2*math.pi-math.radians(30))
    print myangle, myangle.dms(), myangle.canonical()
    print myangle.hms()
    print Angle(0)
    print repr(Angle(math.pi))
