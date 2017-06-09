"""
Purpose:   [1] To disambiguate the azimuthal component of the SDO/HMI vector magnetic field data (See Section 5 of Hoeksema et al. 2014 [open-access at https://link.springer.com/article/10.1007%2Fs11207-014-0516-8])
           [2] To decompose the three components of the vector magnetic field -- field strength, disambiguated azimuth, and inclination -- and project them into one of two coordinate systems: CCD or CEA (See Sun, X. 2013 [https://arxiv.org/abs/1309.2392] and Bobra et al. 2014 [open-access at https://link.springer.com/article/10.1007%2Fs11207-014-0529-3]). 

Usage:     This code depends on the numpy, scipy, astropy, sunpy, and drms libraries.
           The first four libraries are in the standard Anaconda distribution.
           The drms library can be obtained from pip: http://drms.readthedocs.io/en/stable/intro.html#installation.
           This code is compatible with python 3.5.x.

Examples:  See example disambiguate_data.py in the same repository.

Adapted:   From Xudong Sun's IDL code to do the same thing

Written:   Monica Bobra
           June 2017
"""

# import some modules
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from sunpy import wcs
import sunpy.coordinates
from sunpy.map import Map
from datetime import datetime as dt_obj
import drms

__all__ = ['basic','coordinate_transform']
__author__ = 'Monica Bobra'

class basic(object):
    """
    Class for the basic disambiguation functions.
    """
    #===========================================
    @staticmethod
    def get_data(single_record):

        """function: get_data
        This function reads the appropriate data and metadata.

        Parameters
        ----------
        single_record: string
            Single recordset specification, e.g. 'hmi.sharp_720s[377][2011.02.15_00:00:00]'
 
        Returns
        -------
        result: list
            List containing five items:
            [0] The relevant WCS keywords as a pandas dataframe
            [1] FITS file containing azimuthal component of magnetic field vector as an astropy HDUList object
            [2] FITS file containing field strength of magnetic field vector as an astropy HDUList object
            [3] FITS file containing inclination component of magnetic field vector as an astropy HDUList object
            [4] FITS file containing disambiguation information as an astropy HDUList object
        """

        c = drms.Client()

        try:
            keys, segments = c.query(single_record, key='T_REC, CRPIX1, CRPIX2, CRVAL1, CRVAL2, CDELT1, CRLN_OBS, CRLT_OBS, CROTA2, DSUN_REF, RSUN_REF', seg='inclination, azimuth, field, disambig')
        except:
            print("Invalid recordset specification")
            sys.exit(1)

        if (len(keys) > 1):
            print("Specify only one record")
            sys.exit(1)

        baseurl = 'http://jsoc.stanford.edu'
        azimuth = fits.open(baseurl+segments.azimuth[0])  
        field = fits.open(baseurl+segments.field[0]) 
        inclination = fits.open(baseurl+segments.inclination[0]) 
        disambig = fits.open(baseurl+segments.disambig[0]) 

        return [keys, azimuth, field, inclination, disambig]
 
    #===========================================
    @staticmethod
    def perform_disambiguation(azimuth, disambig, method):
    
        """function: perform_disambiguation
        This function performs the actual disambiguation.

        Parameters
        ----------
        azimuth: astropy HDUList object
            FITS file containing azimuthal component of magnetic field vector

        disambig: astropy HDUList object
            FITS file containing disambiguation information

        method: int
            Method used to disambiguate the data: 0 for potential acute, 1 for random, 2 for radial acute (suggested)
 
        Returns
        -------
        result: astropy HDUList object
            FITS file containing disambiguated azimuthal component of magnetic field vector
        """

        size_azimuth  = azimuth[1].shape
        size_disambig = disambig[1].shape
        if(size_disambig != size_azimuth):
            print("file_azimuth and file_disambig are not of the same dimensions")
            sys.exit(1)

        if(method < 0 or method > 2):
            method = 2
            print('Invalid disambiguation method, set to default method = 2')        
    
        disambig[1].data = (disambig[1].data / (np.power(2, method))).astype('uint8')
        disambig[1].data = disambig[1].data*180.
        azimuth[1].data = disambig[1].data + azimuth[1].data
        print("Disambiguated the data.")
     
        return azimuth

class coordinate_transform(object):
    """
    Class for coordinate transformations.
    """
    #===========================================
    @staticmethod
    def ccd(azimuth, field, inclination, keys):

        """function: ccd_conversion
        This function constructs HMI vector field from its native coordinates (field, inclination, and disambiguated azimuth) 
        to spherical coordinate components (B_phi, meridional B_theta, radial B_r) on the CCD grid.

        Parameters
        ----------
        azimuth: astropy HDUList object
            FITS file containing disambiguated azimuthal component of magnetic field vector

        field: astropy HDUList object
            FITS file containing field strength of magnetic field vector

        inclination: astropy HDUList object
            FITS file containing inclination component of magnetic field vector

        keys: pandas dataframe
            The relevant WCS keywords
 
        Returns
        -------
        latlon: astropy SkyCoord object (see: http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html)
            Astropy coordinate frame in heliographic stonyhurst coordinates

        bptr: numpy.ndarray
            Array containing the three spherical components of the vector magnetic field in 
            heliographic stonyhurst coordinates (Bp, Bt, Br) on the CCD grid
        """

        size_field  = field[1].shape
        size_inclination = inclination[1].shape
        size_azimuth = azimuth[1].shape
        if(size_field != size_inclination != size_azimuth):
            print("The three components are not of the same dimensions")
            print("Field",size_field,"inclination",size_inclination,"azimuth",size_azimuth)
            sys.exit(1)    

        # convert to radians
        gamma  =  inclination[1].data*(np.pi/180.)
        psi    =  azimuth[1].data*(np.pi/180.)

        # convert to b_xi, b_eta, and b_zeta as defined in Equation 1 of Sun et al. 2013
        b_xi   = -1*field[1].data * np.sin(gamma) * np.sin(psi)
        b_eta  = field[1].data * np.sin(gamma) * np.cos(psi)
        b_zeta = field[1].data * np.cos(gamma)

        # --- convert from helioprojective cartesian to stonyhurst heliographic coordinates ---

        # step 1: convert from pixels to helioprojective cartesian coordinates
        coord_x = np.ndarray([azimuth[1].header['NAXIS2'], azimuth[1].header['NAXIS1']])
        coord_y = np.ndarray([azimuth[1].header['NAXIS2'], azimuth[1].header['NAXIS1']])

        for pix_x in range(azimuth[1].header['NAXIS2']):                         # data[1].header gives the uncompressed header
            for pix_y in range(azimuth[1].header['NAXIS1']):
                coord_x[pix_x, pix_y] = pix_x + 1
                coord_y[pix_x, pix_y] = pix_y + 1

        x1 = (((coord_x-keys.CRPIX1[0])*np.cos(keys.CROTA2[0]*(np.pi/180.)) - (coord_y-keys.CRPIX2[0])*np.sin(keys.CROTA2[0]*(np.pi/180.)))*keys.CDELT1[0]+keys.CRVAL1[0])
        y1 = (((coord_y-keys.CRPIX2[0])*np.cos(keys.CROTA2[0]*(np.pi/180.)) + (coord_x-keys.CRPIX1[0])*np.sin(keys.CROTA2[0]*(np.pi/180.)))*keys.CDELT1[0]+keys.CRVAL2[0])

        # step 2: populate these values into an astropy coordinate frame
        # inputs to SkyCoord are defined here: 
        # https://github.com/sunpy/sunpy/blob/master/sunpy/coordinates/frames.py

        # the input dateobs takes a string! convert t_rec to a string of the required format
        def parse_tai_string(tstr,datetime=True):
            year   = int(tstr[:4])
            month  = int(tstr[5:7])
            day    = int(tstr[8:10])
            hour   = int(tstr[11:13])
            minute = int(tstr[14:16])
            if datetime: return dt_obj(year,month,day,hour,minute)
            else: return year,month,day,hour,minute

        out  = parse_tai_string(keys.T_REC[0])
        dateobs_out = out.strftime('%Y/%m/%dT%H:%M:%S')

        c = SkyCoord(x1 * u.arcsec, y1 * u.arcsec, frame='helioprojective', rsun = keys.RSUN_REF[0] * u.meter, L0=0 * u.degree, B0=keys.CRLT_OBS[0] * u.degree, D0=keys.DSUN_REF[0] * u.meter, dateobs=dateobs_out)

        # step 3: transform the skycoord frame from helioprojective cartesian to stonyhurst heliographic coordinates
        lonlat = c.transform_to('heliographic_stonyhurst')

        # construct transformation matrix for the heliographic components of the magnetic field
        # according to Eq (1) in Gary & Hagyard (1990)
        # see Eq (7)(8) in Sun (2013) for implementation
     
        b = keys.CRLT_OBS[0] * (np.pi/180.) # b-angle, disk center latitude, in radians
        p = - keys.CROTA2[0] * (np.pi/180.)   # p-angle, negative of CROTA2, in radians
        Phi = np.array(lonlat.lon) * (np.pi/180.)
        Lambda = np.array(lonlat.lat) * (np.pi/180.)

        sinb = np.sin(b)
        cosb = np.cos(b)
        sinp = np.sin(p)
        cosp = np.cos(p)
        sinphi = np.sin(Phi)
        cosphi = np.cos(Phi)
        sinlam = np.sin(Lambda)
        coslam = np.cos(Lambda)

        k11 = coslam * (sinb * sinp * cosphi + cosp * sinphi) - sinlam * cosb * sinp
        k12 = - coslam * (sinb * cosp * cosphi - sinp * sinphi) + sinlam * cosb * cosp
        k13 = coslam * cosb * cosphi + sinlam * sinb
        k21 = sinlam * (sinb * sinp * cosphi + cosp * sinphi) + coslam * cosb * sinp
        k22 = - sinlam * (sinb * cosp * cosphi - sinp * sinphi) - coslam * cosb * cosp
        k23 = sinlam * cosb * cosphi - coslam * sinb
        k31 = - sinb * sinp * sinphi + cosp * cosphi
        k32 = sinb * cosp * sinphi + sinp * cosphi
        k33 = - cosb * sinphi

        # create output arrays of (Bp,Bt,Br), which is identical to (Bxh, -Byh, Bzh)
        # in Gary & Hagyard (1990), see Appendix in Sun (2013)

        bptr = np.ndarray([azimuth[1].header['NAXIS2'], azimuth[1].header['NAXIS1'], 3])
        bptr[:,:,0] = k31 * b_xi + k32 * b_eta + k33 * b_zeta
        bptr[:,:,1] = k21 * b_xi + k22 * b_eta + k23 * b_zeta
        bptr[:,:,2] = k11 * b_xi + k12 * b_eta + k13 * b_zeta

        return [lonlat, bptr]
