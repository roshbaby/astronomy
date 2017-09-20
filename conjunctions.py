from math import radians
from calendar import Date, Time, JulianDayNumber
from angle import Latitude, Longitude
from coordinates import SphCoord
from interpolation import Inter5polate

# Helper functions first
"""
Helper function to verify that the input list contains only SphCoord objects
@param coords List to check contents for
"""
def check_sphcoords(coords_):
    for coord in coords_:
        assert isinstance(coord, SphCoord), \
               'coords_ should contain only SphCoord objects'

"""
@param dates_ list of 5 'dates' as JulianDayNumbers
@param coords1_ list of 5 SphCoords (of a celestial body) for the dates_
@param coords2_ list of 5 SphCoords (of another celestial body) for the dates_
@return (t, ddelta) tuple of the JDN of conjuction as a JulianDayNumber object
        and the difference in the 'declinations' (or 'latitudes') at the
        moment of conjunction as an Latitude object
"""
def conjunction(dates_, coords1_, coords2_, prec_):
    assert isinstance(dates_, list) and isinstance(coords1_, list) \
           and isinstance (coords2_, list), 'Arguments must be lists'
    assert len(dates_) == 5 and len(coords1_) == 5 and len(coords2_) == 5, \
           'Argument list must contain exactly 5 elements'
    for date in dates_:
        assert isinstance(date, JulianDayNumber), \
               'dates_ must contain only JulianDayNumber objects'
    check_sphcoords(coords1_)
    check_sphcoords(coords2_)

    # Construct a list of 'decimal' dates
    jdndates = []
    for date in dates_:
        jdndates.append(date.jdn)
    print

    # Construct a list of differences in radians for alpha and delta between
    # object1 and object 2 for the dates given
    dalphas = []
    ddeltas = []
    for coord1, coord2 in zip(coords1_,coords2_):
        dalpha = coord1.a.rads - coord2.a.rads
        dalphas.append(dalpha)
        ddelta = coord1.b.rads - coord2.b.rads
        ddeltas.append(ddelta)

    # Interpolate to find the zero for dalphas
    dalpha_data = zip(jdndates, dalphas)
    ipol = Inter5polate(dalpha_data)
    jdndate_conjunction = ipol.zero(prec_)

    # Interpolate to find the value of ddelta for the date of conjunction
    ddelta_data = zip(jdndates, ddeltas)
    ipol = Inter5polate(ddelta_data)
    ddelta_conjunction_rads = ipol.compute(jdndate_conjunction)

    # Create objects to return
    date_return = JulianDayNumber(Date(0,1,1),Time(0,0,0))
    date_return.jdn = jdndate_conjunction # Hack to build a JDN from a decimal value
    ddelta_return = Latitude(ddelta_conjunction_rads)

    return (date_return, ddelta_return)


"""
@param dates_ List of 5 'dates' as JulianDayNumber objects
@param star_ Coordinates of a star/fixed celestial body as a SphCoord object
@param coords_ List of 5 SphCoords (for a celestial body) for the dates_
@return (t, ddelta) tuple of the JDN of conjuction as a JulianDayNumber object
        and the difference in the 'declinations' (or 'latitudes') at the
        moment of conjunction as an Latitude object
"""
def conjunction_with_star(dates_, star_, coords_, prec_):
    assert isinstance(star_, SphCoord), 'star_ must be a SphCoord object'
    # We assume the star's coordinates stay constant over the list of dates_
    # provided. Create a dummy list of coords corresponding to the dates_
    star_coords = [star_ for date in dates]
    return conjunction(dates_, star_coords, coords_, prec_)


if __name__ == "__main__":
    dates = [
        JulianDayNumber(Date(1991,8,5),Time(0,0,0)),
        JulianDayNumber(Date(1991,8,6),Time(0,0,0)),
        JulianDayNumber(Date(1991,8,7),Time(0,0,0)),
        JulianDayNumber(Date(1991,8,8),Time(0,0,0)),
        JulianDayNumber(Date(1991,8,9),Time(0,0,0)),
    ]
    # Mercury
    coords1 = [
        SphCoord(
            Longitude( radians((10.0 +24.0/60.0 +30.125/3600.0)*15) ),
            Latitude( radians(6 +26.0/60.0 +32.05/3600.0) )
        ),
        SphCoord(
            Longitude( radians((10.0 +25.0/60.0 + 0.342/3600.0)*15) ),
            Latitude( radians(6 +10.0/60.0 +57.72/3600.0) )
        ),
        SphCoord(
            Longitude( radians((10.0 +25.0/60.0 +12.515/3600.0)*15) ),
            Latitude( radians(5 +57.0/60.0 +33.08/3600.0) )
        ),
        SphCoord(
            Longitude( radians((10.0 +25.0/60.0 + 6.235/3600.0)*15) ),
            Latitude( radians(5 +46.0/60.0 +27.07/3600.0) )
        ),
        SphCoord(
            Longitude( radians((10.0 +24.0/60.0 +41.185/3600.0)*15) ),
            Latitude( radians(5 +37.0/60.0 +48.45/3600.0) )
        )
    ]
    # Venus
    coords2 = [
        SphCoord(
            Longitude( radians((10.0 +27.0/60.0 +27.175/3600.0)*15) ),
            Latitude( radians(4 + 4.0/60.0 +41.83/3600.0) )
        ),
        SphCoord(
            Longitude( radians((10.0 +26.0/60.0 +32.410/3600.0)*15) ),
            Latitude( radians(3 +55.0/60.0 +54.66/3600.0) )
        ),
        SphCoord(
            Longitude( radians((10.0 +25.0/60.0 +29.042/3600.0)*15) ),
            Latitude( radians(3 +48.0/60.0 + 3.51/3600.0) )
        ),
        SphCoord(
            Longitude( radians((10.0 +24.0/60.0 +17.191/3600.0)*15) ),
            Latitude( radians(3 +41.0/60.0 +10.25/3600.0) )
        ),
        SphCoord(
            Longitude( radians((10.0 +22.0/60.0 +57.024/3600.0)*15) ),
            Latitude( radians(3 +35.0/60.0 +16.61/3600.0) )
        )
    ]

    jdn, ddelta = conjunction(dates, coords1, coords2, 1e-6)
    print 'Date/Time of conjunction:', jdn.get_date(), jdn.get_time()
    print 'Difference in declination/latitude:', ddelta

    dates = [
        JulianDayNumber(Date(1996,2, 7), Time(0,0,0)),
        JulianDayNumber(Date(1996,2,12), Time(0,0,0)),
        JulianDayNumber(Date(1996,2,17), Time(0,0,0)),
        JulianDayNumber(Date(1996,2,22), Time(0,0,0)),
        JulianDayNumber(Date(1996,2,27), Time(0,0,0))
    ]
    # Beta Librae
    star = SphCoord(
        Longitude( radians((15.0 +17.0/60.0 +0.446/3600.0)*15) ),
        Latitude( -radians(9.0 +22.0/60.0 +58.47/3600.00) )
    )
    # Planet 4 Vesta
    coords = [
        SphCoord(
            Longitude( radians((15.0 + 3.0/60.0 +51.937/3600.0)*15) ),
            Latitude( -radians( 8.0 +57.0/60.0 +34.51/3600.0) )
        ),
        SphCoord(
            Longitude( radians((15.0 + 9.0/60.0 +57.327/3600.0)*15) ),
            Latitude( -radians( 9.0 + 9.0/60.0 + 3.88/3600.0) )
        ),
        SphCoord(
            Longitude( radians((15.0 +15.0/60.0 +37.898/3600.0)*15) ),
            Latitude( -radians( 9.0 +17.0/60.0 +37.94/3600.0) )
        ),
        SphCoord(
            Longitude( radians((15.0 +20.0/60.0 +50.6332/3600.0)*15) ),
            Latitude( -radians( 9.0 +23.0/60.0 +16.25/3600.0) )
        ),
        SphCoord(
            Longitude( radians((15.0 + 25.0/60.0 +32.695/3600.0)*15) ),
            Latitude( -radians( 9.0 +26.0/60.0 + 1.01/3600.0) )
        ),
    ]
    jdn, ddelta = conjunction_with_star(dates, star, coords, 1e-6)
    print 'Date/Time of conjunction:', jdn.get_date(), jdn.get_time()
    print 'Difference in declination/latitude:', ddelta
