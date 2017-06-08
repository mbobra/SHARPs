"""
This example shows how to disambiguate the vector magnetic field data and construct the field vector in spherical coordinate components on the CCD grid.

This example uses data from the hmi.sharp_720s (Bobra et al. 2014; open-access at https://link.springer.com/article/10.1007%2Fs11207-014-0529-3), however; it works just the same for the full-disk series, hmi.B_720s (Hoeksema et al. 2014; open-access at https://link.springer.com/article/10.1007%2Fs11207-014-0516-8)
"""

import disambiguation

# fetch the data from JSOC
keys, azimuth, field, inclination, disambig = disambiguation.basic.get_data('hmi.sharp_720s[377][2011.02.15_00:00:00]')

# disambiguate the azimuthal component of the magnetic field
disambiguated_azimuth = disambiguation.basic.perform_disambiguation(azimuth, disambig, 2)

# construct the field vector in spherical coordinate components on the CCD grid
latlon, bptr = disambiguation.coordinate_transform.ccd(disambiguated_azimuth, field, inclination, keys)
