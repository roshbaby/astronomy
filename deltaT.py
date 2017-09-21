from calendar import Time
from interpolation import LagrangeInterpolate

"""
Table of (year, deltaT) tuples from 1620 to 2010 for interpolation
From "Astronomical Algorithms" by Jean Meeus, 2nd Ed., Table 10.A, pg. 79
deltaT is in seconds
"""
deltaT_table = [
    (1620, 121),
    (1622, 112),
    (1624, 103),
    (1626,  95),
    (1628,  88),

    (1630,  82),
    (1632,  77),
    (1634,  72),
    (1636,  68),
    (1638,  63),

    (1640,  60),
    (1642,  56),
    (1644,  53),
    (1646,  51),
    (1648,  48),

    (1650,  46),
    (1652,  44),
    (1654,  42),
    (1656,  40),
    (1658,  38),

    (1660,  35),
    (1662,  33),
    (1664,  31),
    (1666,  29),
    (1668,  26),

    (1670,  24),
    (1672,  22),
    (1674,  20),
    (1676,  18),
    (1678,  16),

    (1680,  14),
    (1682,  12),
    (1684,  11),
    (1686,  10),
    (1688,   9),

    (1690,   8),
    (1692,   7),
    (1694,   7),
    (1696,   7),
    (1698,   7),

    (1700,   7),
    (1702,   7),
    (1704,   8),
    (1706,   8),
    (1708,   9),

    (1710,   9),
    (1712,   9),
    (1714,   9),
    (1716,   9),
    (1718,  10),

    (1720,  10),
    (1722,  10),
    (1724,  10),
    (1726,  10),
    (1728,  10),

    (1730,  10),
    (1732,  10),
    (1734,  11),
    (1736,  11),
    (1738,  11),

    (1740,  11),
    (1742,  11),
    (1744,  12),
    (1746,  12),
    (1748,  12),

    (1750,  12),
    (1752,  13),
    (1754,  13),
    (1756,  13),
    (1758,  14),

    (1760,  14),
    (1762,  14),
    (1764,  14),
    (1766,  15),
    (1768,  15),

    (1770,  15),
    (1772,  15),
    (1774,  15),
    (1776,  16),
    (1778,  16),

    (1780,  16),
    (1782,  16),
    (1784,  16),
    (1786,  16),
    (1788,  16),

    (1790,  16),
    (1792,  15),
    (1794,  15),
    (1796,  14),
    (1798,  13),

    (1800,  13.1),
    (1802,  12.5),
    (1804,  12.2),
    (1806,  12.0),
    (1808,  12.0),

    (1810,  12.0),
    (1812,  12.0),
    (1814,  12.0),
    (1816,  12.0),
    (1818,  11.9),

    (1820,  11.6),
    (1822,  11.0),
    (1824,  10.2),
    (1826,   9.2),
    (1828,   8.2),

    (1830,   7.1),
    (1832,   6.2),
    (1834,   5.6),
    (1836,   5.4),
    (1838,   5.3),

    (1840,   5.4),
    (1842,   5.6),
    (1844,   5.9),
    (1846,   6.2),
    (1848,   6.5),

    (1850,   6.8),
    (1852,   7.1),
    (1854,   7.3),
    (1856,   7.5),
    (1858,   7.6),

    (1860,   7.7),
    (1862,   7.3),
    (1864,   6.2),
    (1866,   5.2),
    (1868,   2.7),

    (1870,   1.4),
    (1872,  -1.2),
    (1874,  -2.8),
    (1876,  -3.8),
    (1878,  -4.8),

    (1880,  -5.5),
    (1882,  -5.3),
    (1884,  -5.6),
    (1886,  -5.7),
    (1888,  -5.9),

    (1890,  -6.0),
    (1892,  -6.3),
    (1894,  -6.5),
    (1896,  -6.2),
    (1898,  -4.7),

    (1900,  -2.8),
    (1902,  -0.1),
    (1904,   2.6),
    (1906,   5.3),
    (1908,   7.7),

    (1910,  10.4),
    (1912,  13.3),
    (1914,  16.0),
    (1916,  18.2),
    (1918,  20.2),

    (1920,  21.1),
    (1922,  22.4),
    (1924,  23.5),
    (1926,  23.8),
    (1928,  24.3),

    (1930,  24.0),
    (1932,  23.9),
    (1934,  23.9),
    (1936,  23.7),
    (1938,  24.0),

    (1940,  24.3),
    (1942,  25.3),
    (1944,  26.2),
    (1946,  27.3),
    (1948,  28.2),

    (1950,  29.1),
    (1952,  30.0),
    (1954,  30.7),
    (1956,  31.4),
    (1958,  32.2),

    (1960,  33.1),
    (1962,  34.0),
    (1964,  35.0),
    (1966,  36.5),
    (1968,  38.3),

    (1970,  40.2),
    (1972,  42.2),
    (1974,  44.5),
    (1976,  46.5),
    (1978,  48.5),

    (1980,  50.5),
    (1982,  52.2),
    (1984,  53.8),
    (1986,  54.9),
    (1988,  55.8),

    (1990,  56.9),
    (1992,  58.3),
    (1994,  60.0),
    (1996,  61.6),
    (1998,  63.0),

    (2000,  63.8),
    (2002,  64.3),
    (2004,  64.6),
    (2006,  64.8),
    (2008,  65.5),

    (2010,  66.1)
]

deltaT_y2k_table = [
    (2000,  63.8),
    (2002,  64.3),
    (2004,  64.6),
    (2006,  64.8),
    (2008,  65.5),

    (2010,  66.1)
]


"""
\Delta T = Dynamical Time (TD) - Ephemeris Time (UT)
Chapront and Francou's *approximation*
@return Delta T as a Time object
"""
def deltaT(year):
    table_first_year = deltaT_table[0][0]
    table_last_year = deltaT_table[len(deltaT_table)-1][0]

    ret_sec = 0

    # If year is in our table use it instead
    if year >= table_first_year and year <= table_last_year:
        # First, find the index of the closest year
        my_index = 0
        for idx in range(len(deltaT_table)):
            if deltaT_table[idx][0] > year:
                my_index = idx
                break

        # Choose a slice of the table that is +/- 2 of my_index
        if my_index < 2:
            deltaT_table_slice = deltaT_table[0:5]
        elif len(deltaT_table) - my_index < 3:
            deltaT_table_slice = deltaT_table[len(deltaT_table)-5:]
        else:
            deltaT_table_slice = deltaT_table[my_index-2:my_index+3]

        ipol = LagrangeInterpolate(deltaT_table_slice)
        ret_sec = ipol.compute(year)

    else:
        # Table not available. Use Chapront & Francou's approximations
        t = (year - 2000.0)/100 # Centuries after epoch 2000.0

        if (year < 948):
            ret_sec = 2177 + (497+44.1*t)*t
        else:
            ret_sec = 102 + (102 + 25.3*t)*t
            if year >= table_last_year and year <= 2100:
                ret_sec += 0.37*(year-2100)

    return Time(0,0,ret_sec)


if __name__ == "__main__":
    years = [1619, 1621, 1625.5, 1627, 1821, 1869, 1871, 1871.04, 2009, 2005,
             2003, 2001, 2000, 1820, 1620]
    for yr in years:
        print 'year:', yr, '\tdeltaT:', deltaT(yr)
