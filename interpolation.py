from math import fabs, radians

## Helper functions first
"""
Check if all elements of the provided list are 2 element tuples
@precondition xyvalues_ is already verified to be a list
"""
def checktuples(xyvalues_):
    # Check that all elements are tuples
    idx = 0
    for tup in xyvalues_:
        assert isinstance(tup, tuple), \
               'Item at index {0:d} in xyvalues_ is not a tuple.' \
               ' xyvalues_ should contain only tuples.'.format(idx)
        assert len(tup) == 2, \
               'Item at index {0:d} in xyvalues_ is not a 2 element tuple.' \
               ' xyvalues_ should contain only 2 element tuples.'.format(idx)
        idx += 1


"""
Interpolate based on 3 (x,y) tuples where the x values are equidistant
"""
class Inter3polate:
    def __init__(self, xyvalues_):
        assert isinstance(xyvalues_, list), 'xyvalues_ should be a list'
        assert len(xyvalues_) == 3, 'xyvalues_ should have exactly 3 elements'
        checktuples(xyvalues_)
        # check that the x values are equidistant
        x_distance = xyvalues_[1][0] - xyvalues_[0][0]
        for i in range(1,len(xyvalues_)-1):
            assert xyvalues_[i+1][0] - xyvalues_[i][0] == x_distance, \
                   'xyvalues_ should have equidistant x values.'      \
                   ' Problem with index {0:d} and {1:d}'.format(i,i+1)
        x1, y1 = xyvalues_[0]
        self.x2, self.y2 = xyvalues_[1]
        x3, y3 = xyvalues_[2]
        self.a = float(self.y2 - y1)
        self.b = float(y3 - self.y2)
        self.c = self.b - self.a
        self.x_distance = x_distance

    """
    Interpolate for the given x value
    @param x_inter The value of the x-coordinate at which to interpolate for y
    @return The value of y corresponding to x_inter
    """
    def compute(self, x_inter):
        n = float(x_inter - self.x2)/self.x_distance
        return self.y2 + (self.a + self.b + n*self.c)*n/2

    """
    Determine the extrema (maxima or minima)
    @return (x,y) tuple corresponding to the extrema
    """
    def extrema(self):
        y_extrema = self.y2 - (self.a + self.b)**2/(8*self.c)
        n_extrema = -(self.a + self.b)/(2*self.c)
        x_extrema = self.x2 + n_extrema*self.x_distance
        return (x_extrema, y_extrema)

    """
    A better method to find the zero that works especially when other
    methods would diverge
    @param prec_ The value to use for the precision of the result
    @return The value of x for which the corresponding y is 0
    """
    def zero(self, prec_):
        n_0 = 0
        while True:
            n_0_delta = -(2*self.y2 + n_0*(self.a + self.b + n_0*self.c)) / \
                        (self.a + self.b + 2*self.c*n_0)
            if fabs(n_0_delta) <= fabs(prec_):
                break
            n_0 += n_0_delta
        return self.x2 + n_0*self.x_distance

    """
    Determine the zero using a method that works best for 'linear-like' data
    @param prec_ The value to use for the precision of the result
    @return The value of x for which the corresponding y is 0
    """
    def zero2(self, prec_):
        n_0 = 0
        while True:
            n_0_next = -2*self.y2/(self.a + self.b + n_0*self.c)
            if fabs(n_0_next - n_0) <= fabs(prec_):
                break
            n_0 = n_0_next
        return self.x2 + n_0*self.x_distance


class Inter5polate:
    def __init__(self, xyvalues_):
        assert isinstance(xyvalues_, list), 'xyvalues_ should be a list'
        assert len(xyvalues_) == 5, 'xyvalues_ should have exactly 5 elements'
        checktuples(xyvalues_)
        # check that the x values are equidistant
        x_distance = xyvalues_[1][0] - xyvalues_[0][0]
        for i in range(1,len(xyvalues_)-1):
            assert xyvalues_[i+1][0] - xyvalues_[i][0] == x_distance, \
                   'xyvalues_ should have equidistant x values.'      \
                   ' Problem with index {0:d} and {1:d}'.format(i,i+1)
        x1, y1 = xyvalues_[0]
        x2, y2 = xyvalues_[1]
        self.x3, self.y3 = xyvalues_[2]
        x4, y4 = xyvalues_[3]
        x5, y5 = xyvalues_[4]
        self.A = float(y2 - y1)
        self.B = float(self.y3 - y2)
        self.C = float(y4 - self.y3)
        self.D = float(y5 - y4)
        self.E = float(self.B - self.A)
        self.F = float(self.C - self.B)
        self.G = float(self.D - self.C)
        self.H = float(self.F - self.E)
        self.J = float(self.G - self.F)
        self.K = float(self.J - self.H)
        self.x_distance = x_distance

    """
    Interpolate for the given x value
    @param x_inter The value of the x-coordinate at which to interpolate for y
    @return The value of y corresponding to x_inter
    """
    def compute(self, x_inter):
        n = float(x_inter - self.x3)/self.x_distance
        C1 = n/2.0
        C2 = C1*n
        C3 = n*(n**2-1)/12.0
        C4 = C3*n/2.0
        return self.y3                \
               + (self.B+self.C)*C1   \
               + self.F*C2            \
               + (self.H + self.J)*C3 \
               + self.K*C4

    """
    Determine the extrema (maxima or minima)
    @param prec_ The value to use for the precision of the result
    @return (x,y) tuple corresponding to the extrema
    @todo Need a test for this
    """
    def extrema(self, prec_):
        n_extrema = 0
        while True:
            numer = 6*self.B + 6*self.C - self.H - self.J \
                    + 3*(n_extrema**2)*(self.H + self.J)  \
                    + 2*(n_extrema**3)*self.K
            denom = self.K - 12*self.F
            n_extrema_next = numer / denom
            if fabs(n_extrema_next - n_extrema) <= fabs(prec_):
                break
            n_extrema = n_extrema_next
        x_extrema = self.x3 + n_extrema*self.x_distance,
        y_extrema = self.compute(x_extrema)
        return (x_extrema, y_extrema)

    """
    Find the zero
    @param prec_ The value to use for the precision of the result
    @return The value of x for which the corresponding y is 0
    """
    def zero(self, prec_):
        M = self.K/24.0
        N = (self.H + self.J)/12.0
        P = self.F/2.0 - M
        Q = (self.B + self.C)/2.0 - N
        n_0 = 0
        while True:
            numer = (((M*n_0 + N)*n_0 + P)*n_0 + Q)*n_0 + self.y3
            denom = ((4.0*M*n_0 + 3.0*N)*n_0 + 2.0*P)*n_0 + Q
            n_0_delta = -numer/denom
            if fabs(n_0_delta) <= fabs(prec_):
                break
            n_0 += n_0_delta
        return self.x3 + n_0*self.x_distance


"""
Interpolate based on 'N' (x,y) tuples using Lagrange's formula
@param xyvalues List of 'N' (x,y) tuples
"""
class LagrangeInterpolate:
    """
    Iniitalize the LagrangeInterpolate Object
    @param xyvalues_ list of 'N' (x,y) tuples (N >= 3)
    """
    def __init__(self, xyvalues_):
        assert isinstance(xyvalues_, list), 'xyvalues_ should be a list'
        assert len(xyvalues_) > 2, 'xyvalues_ should have at least 3 elements'
        checktuples(xyvalues_)
        # Check that all x values are unique
        for i in range(len(xyvalues_)-1):
            for j in range(i+1, len(xyvalues_)):
                assert xyvalues_[i][0] != xyvalues_[j][0], \
                       'x values at indexes {0:d} and {1:d} are identical.'\
                       ' x values must be unique.'.format(i,j)
        self.xyvalues = xyvalues_


    """
    Interpolate for the given x value
    @param x_inter The value of the x-coordinate at which to interpolate for y
    """
    def compute(self, x_inter):
        # Generate the Li coefficients and consume them as we go along
        y_inter = 0.0
        for ival in self.xyvalues:
            L = 1.0
            for jval in self.xyvalues:
                if jval != ival:
                    L *= float(x_inter - jval[0]) / float(ival[0] - jval[0])
            y_inter += L*ival[1]
        return y_inter


if __name__ == "__main__":
    ## Inter3polate test
    data = [
        (0,1),                             # Not a list
        [ (0,1), (1,2) ],                  # 2 element list
        [ (3,4), (4,5), (6,7), (0,1,2) ],  # 4 element list
        [ (0,1), 1, (2,3) ],               # Not all tuples
        [ (3,4), (4,5), (0,1,2) ],         # One 3 element tuple
        [ (0,1), (1,2), (1.5,3) ],         # Not equidistant
        [ (0,1), (1.5,2), (2,3) ],         # Not equidistant
    ]
    for item in data:
        try:
            ipol = Inter3polate(item)
        except AssertionError as e:
            print 'Error:', e

    ipol = Inter3polate([(7,0.884226), (8,0.877366), (9,0.870531)])
    print ipol.compute(8.18125) # 0.876125

    ipol = Inter3polate([(12.0,1.3814294), (16.0,1.3812213), (20.0,1.3812453)])
    print ipol.extrema() # x,y = 17.58640,1.3812030

    ipol = Inter3polate([
        (26.0,-radians( (28+13.4/60.0)/60.0) ),
        (27.0, radians( (6 +46.3/60.0)/60.0) ),
        (28.0, radians( (38+23.2/60.0)/60.0) )
    ])
    print ipol.zero(1e-6) # 26.79873
    print ipol.zero2(1e-6) # 26.79873

    ipol = Inter3polate([(-1,-2), (0,3), (1,2)])
    print ipol.zero(1e-12) # -0.720759220056
    print


    ## Inter5polate test
    data = [
        (0,1),                             # Not a list
        [ (3,4), (4,5), (6,7), (0,1,2) ],  # 4 element list
        [ (0,1), (1,2), (2,3), (3,4), (4,5), (5,6) ], # 6 element list
        [ (0,1), 1, (2,3), (3,4), (4,5) ], # Not all tuples
        [ (3,4), (4,5), (6,7), (7,8), (0,1,2) ], # One 3 element tuple
        [ (0,1), (1,2), (2.5,3), (3,4), (4,5) ], # Not equidistant
        [ (0,1), (1,2), (2,3), (3.5,4), (4,5) ], # Not equidistant
    ]
    for item in data:
        try:
            ipol = Inter5polate(item)
        except AssertionError as e:
            print 'Error:', e

    data = [
        (27.0, 54+36.125/60.0),
        (27.5, 54+24.606/60.0),
        (28.0, 54+15.486/60.0),
        (28.5, 54+ 8.694/60.0),
        (29.0, 54+ 4.133/60.0)
    ]
    ipol = Inter5polate(data)
    print ipol.compute(28.0+(3+20.0/60.0)/24.0) # 54.2228

    data = [
        (25.0, -(1 +11/60.0 +21.20/3600.0)),
        (26.0, -(0 +28/60.0 +12.31/3600.0)),
        (27.0,  (0 +16/60.0 + 7.02/3600.0)),
        (28.0,  (1 + 1/60.0 +00.13/3600.0)),
        (29.0,  (1 +45/60.0 +46.33/3600.0))
    ]
    ipol = Inter5polate(data)
    print ipol.zero(1e-6)  # 26.638587
    print


    ## LagrangeInterpolate test
    data = [
        (0,1),                      # Not a list
        [ (0,1), (1,2) ],           # 2 element list
        [ (0,1), 1, (2,3) ],        # Not all tuples
        [ (3,4), (4,5), (6,7), (0,1,2) ], # One 3 element tuple
        [ (0,1), (0,2), (2,3) ],    # identical x values at elements 0, 1
        [ (0,1), (1,2), (1,3) ]     # identical x values at elements 1, 2
    ]
    for item in data:
        try:
            ipol = LagrangeInterpolate(item)
        except AssertionError as e:
            print 'Error:', e

    ipol = LagrangeInterpolate([(0,0),(1,1),(2,4)]) # y = x^2
    print ipol.compute(1.1)  # 1.21
    print ipol.compute(2.1)  # 4.41
    print ipol.compute(3.1)  # 9.61

    # Test for sin values of degrees
    isinpol = LagrangeInterpolate(
        [
            (29.43, 0.4913598528),
            (30.97, 0.5145891926),
            (27.69, 0.4646875083),
            (28.11, 0.4711658342),
            (31.58, 0.5236885653),
            (33.05, 0.5453707057)
        ]
    )
    print isinpol.compute(30.0)  # 0.5
    print isinpol.compute(0.0)   # 0.0 - really outside interpolation range
    print isinpol.compute(90.0)  # 1.0 - really outside interpolation range
