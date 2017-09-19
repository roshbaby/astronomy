"""
Interpolate based on 'N' (x,y) tuples using Lagrange's formula
@param xyvalues List of 'N' (x,y) tuples
@param
"""
class Interpolate:
    """
    Iniitalize the Interpolate Object
    @param xyvalues_ list of (x,y) tuples (at least 3 tuples in the list)
    """
    def __init__(self, xyvalues_,):
        assert isinstance(xyvalues_, list), 'xyvalues_ should be a list'
        assert len(xyvalues_) > 2, 'xyvalues_ should have at least 3 elements'
        self.xyvalues = xyvalues_

        # Check that all elements are tuples
        for i in range(len(xyvalues_)):
            assert isinstance(xyvalues_[i], tuple), \
                   'Item at index {0:d} in xyvalues_ is not a tuple.' \
                   ' xyvalues_ should contain only tuples.'.format(i)
            assert len(xyvalues_[i]) == 2, \
                   'Item at index {0:d} in xyvalues_ is not a 2 element tuple.' \
                   ' xyvalues_ should contain only 2 element tuples.'.format(i)

        # Check that all x values are unique
        for i in range(len(xyvalues_)-1):
            for j in range(i+1, len(xyvalues_)):
                assert xyvalues_[i][0] != xyvalues_[j][0], \
                       'x values at indexes {0:d} and {1:d} are identical.'\
                       ' x values must be unique.'.format(i,j)

    """
    Intepolate for the given x value
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
    data = [
        (0,1),                      # Not a list
        [ (0,1), (1,2) ],           # 2 element list
        [ (0,1), 1, (2,3) ],        # Not all tuples
        [ (3,4), (4,5), (6,7),  (0,1,2) ], # One 3 element tuple
        [ (0,1), (0,2), (2,3) ],    # identical x values at elements 0, 1
        [ (0,1), (1,2), (1,3) ]     # identical x values at elements 1, 2
    ]

    for item in data:
        try:
            ipol = Interpolate(item)
        except AssertionError as e:
            print 'Error:', e

    ipol = Interpolate([(0,0),(1,1),(2,4)]) # y = x^2
    print ipol.compute(1.1)  # 1.21
    print ipol.compute(2.1)  # 4.41
    print ipol.compute(3.1)  # 9.61

    # Test for sin values of degrees
    isinpol = Interpolate(
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
