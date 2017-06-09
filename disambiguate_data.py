"""
This example shows how to disambiguate the vector magnetic field data and construct the field vector in spherical coordinate components on the CCD grid.

This example uses data from the hmi.sharp_720s (Bobra et al. 2014; open-access at https://link.springer.com/article/10.1007%2Fs11207-014-0529-3), however; it works just the same for the full-disk series, hmi.B_720s (Hoeksema et al. 2014; open-access at https://link.springer.com/article/10.1007%2Fs11207-014-0516-8)
"""

import disambiguation

# fetch the data from JSOC by providing a recordset specification and a disambiguation method
query_info = disambiguation.Basic('hmi.sharp_720s[377][2011.02.15_00:00:00]', 2)
keys, azimuth, field, inclination, disambig = disambiguation.Basic.get_data(query_info)

# disambiguate the azimuthal component of the magnetic field
disambiguated_azimuth = disambiguation.Basic.perform_disambiguation(query_info, azimuth, disambig)

# construct the field vector in spherical coordinate components on the CCD grid
data_object = disambiguation.CoordinateTransform(disambiguated_azimuth, field, inclination, keys)
latlon, bptr = disambiguation.CoordinateTransform.ccd(data_object)
