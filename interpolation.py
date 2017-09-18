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
        assert type(xyvalues_) is list, 'xyvalues_ shouldbe a list'
        assert len(xyvalues_) > 2, 'xyvalues_ should at least be a 3 element list'
        self.xyvalues = xyvalues_

        # Check that all elements are tuples
        for i in range(len(xyvalues_)):
            assert type(xyvalues_[i]) is tuple, 'Item at index {0:d} in xyvalues_ is not a tuple. xyvalues_ should contain only tuples'.format(i)

        # Check that all x values are unique
        for i in range(len(xyvalues_)-1):
            for j in range(i+1, len(xyvalues_)):
                assert xyvalues_[i][0] != xyvalues_[j][0], 'x values at indexes {0:d} and {1:d} are identical. x values must be unique'.format(i,j)


    """
    Intepolate for the given x value
    @param x_inter The value of the x-coordinate at which to interpolate for y
    """
    def compute(self, x_inter):
        # Generate the Li coefficients and consume them as we go along
        y_inter = 0.0
        for i in range(len(self.xyvalues)):
            L = 1.0
            for j in range(len(self.xyvalues)):
                if j != i:
                    L *= float(x_inter - self.xyvalues[j][0]) /         \
                         float(self.xyvalues[i][0] - self.xyvalues[j][0])
            y_inter += L*self.xyvalues[i][1]
        return y_inter


if __name__ == "__main__":
    try:
        ipol = Interpolate((0,1))
    except AssertionError as e:
        print e

    try:
        ipol = Interpolate([(0,1),(1,2)])
    except AssertionError as e:
        print e

    try:
        ipol = Interpolate([(0,1),1,(2,3)])
    except AssertionError as e:
        print e

    try:
        ipol = Interpolate([(0,1),(0,2),(2,3)])
    except AssertionError as e:
        print e

    try:
        ipol = Interpolate([(0,1),(1,2),(1,3)])
    except AssertionError as e:
        print e

    ipol = Interpolate([(0,0),(1,1),(2,4)]) # y = x^2
    print ipol.compute(1.1)
    print ipol.compute(2.1)
    print ipol.compute(3.1)

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
    print isinpol.compute(30.0)
    print isinpol.compute(0.0)
    print isinpol.compute(90.0)
