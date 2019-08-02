"""Calculate conditions using LSST OpSim packages."""
import os

import numpy as np
import astropy.time
from astropy.time import Time, TimeDelta
import healpy

# import lsst.sims.ocs.kernel.time_handler
# import lsst.sims.speedObservatory
# import lsst.sims.seeingModel
# import lsst.sims.cloudModel
# import lsst.sims.skybrightness_pre
import lsst.sims.utils

from lsst.ts.observatory.model import ObservatoryModel, Target
from lsst.ts.observatory.model import ObservatoryState
from lsst.ts.observatory.model import ObservatoryPosition
from lsst.sims.seeingModel import SeeingData, SeeingModel
from lsst.sims.cloudModel import CloudData
from lsst.sims.skybrightness_pre import SkyModelPre


# constants

SKY_DATA_PATH = '/data/des70.a/data/neilsen/sims_skybrightness_pre/data'
DEFAULT_START_TIME = astropy.time.Time('2022-01-01 00:00:00Z', format='iso')

# exception classes

# interface functions

calc_five_sigma_depth = lsst.sims.utils.m5_flat_sed

# classes


class SeeingSource(object):
    """Callable object to calculate seeing values."""

    def __init__(self, survey_start_time=DEFAULT_START_TIME):
        """Initialize source of simulated seeing values.

        Args:
           - survey_start_time :: an astropy.time.Time object
             designating the start of the survey.
        """
        self.start_time = survey_start_time
        self.seeing_data = SeeingData(self.start_time)
        self.seeing_model = SeeingModel()

        self.seeing_filter_index = {}
        for filter_index, filter_name in enumerate(self.seeing_model.filter_list):
            self.seeing_filter_index[filter_name] = filter_index

        
    def __call__(self, seconds_into_survey, band, airmass):
        """Return simulated seeing parameters.

        Args:
           - seconds_into_survey :: the time into the survey (in seconds)
             for which to return the seeing parameters.
           - band :: the filter for which to return the seeing values.
           - airmass :: the airmass of the pointing for which seeing
             should be calculated.

        Returns:
           a tuple with fwhm500, fwhmGeom, and fwhmEff
        """
        fwhm500 = self.seeing_data(
            self.start_time + TimeDelta(seconds_into_survey, format='sec'))
        seeings = self.seeing_model(fwhm500, airmass)
        fwhm_geom = seeings['fwhmGeom'][self.seeing_filter_index[band]]
        fwhm_eff = seeings['fwhmEff'][self.seeing_filter_index[band]]
        return fwhm500, fwhm_geom, fwhm_eff


class CloudSource(object):
    """Callable object to return simulated cloud levels."""

    def __init__(self, survey_start_time=DEFAULT_START_TIME):
        """Initialize source of simulated cloud levels.

        Args:
           - survey_start_time :: an astropy.time.Time object
             designating the start of the survey.
        """
        self.start_time = survey_start_time        
        self.cloud_data = CloudData(DEFAULT_START_TIME)
        

    def __call__(self, seconds_into_survey):
        """Return the simulated cloud level.

        Args:
           - seconds_into_survey :: the time into the survey (in seconds)
             for which to return the cloud level.

        Returns:
           the fraction cloud coverage
        """
        cloud_level = self.cloud_data(
            self.start_time + TimeDelta(seconds_into_survey, format='sec'))

        return cloud_level


class SlewTimeSource(object):
    """Callable object to calculate slew times."""

    pre_position = ObservatoryPosition()
    post_position = ObservatoryPosition()
    
    def __init__(self, survey_start_time=DEFAULT_START_TIME):
        """Initialize slew time calculator.
        """
        self.start_time = survey_start_time
        self.observatory = ObservatoryModel()
        self.observatory.configure_from_module()
        self.observatory.params.rotator_followsky = True

    def __call__(self, seconds_into_survey,
                 pre_alt, pre_az, pre_band,
                 post_alt, post_az, post_band):
        """Calculate the slew time.

        Args:
           - seconds_into_survey :: the time into the survey (in seconds)
           - pre_alt :: the alt of the starting pointing, in degrees
           - pre_az :: the az of the starting pointing, in degrees
           - pre_band :: the filter at the starting pointing
           - post_alt :: the alt of the ending pointing, in degrees
           - post_az :: the az of the ending pointing, in degrees
           - post_band :: the filter at the starting pointing
        """        
        current_time = self.start_time + TimeDelta(seconds_into_survey, format='sec')
        self.pre_position.time = current_time
        self.pre_position.tracking = True
        self.pre_position.alt_rad = np.radians(pre_alt)
        self.pre_position.az_rad = np.radians(pre_az)
        self.pre_position.rot_rad = self.observatory.current_state.rot_rad
        self.pre_position.filter = pre_band
        pre_state = self.observatory.get_closest_state(self.pre_position)

        self.post_position.time = current_time
        self.post_position.tracking = True
        self.post_position.alt_rad = np.radians(post_alt)
        self.post_position.az_rad = np.radians(post_az)
        self.post_position.rot_rad = self.observatory.current_state.rot_rad
        self.post_position.filter = post_band
        post_state = self.observatory.get_closest_state(self.post_position)

        slew_time = self.observatory.get_slew_delay_for_state(
            pre_state, post_state, False)

        return slew_time


class SkyBrightnessSource(object):
    """Interface to pre-computed sky brightnesses."""

    def __init__(self, data_path=SKY_DATA_PATH):
        """Initialize the interface to pre-computed sky brightness.

        Args:
           - data_path :: the full path to the parent directory
             with pre-computed sky brightness values.
        """
        self.data_path = data_path
        self._sky_sim_hpx = None
        self._sky_sim_field = None
        self._nside = None

    @property
    def sky_sim_hpx(self):
        """An lsst.sims.skybrightness_pre.SkyModelPre for healpix values."""
        if self._sky_sim_hpx is None:
            self._sky_sim_hpx = lsst.sims.skybrightness_pre.SkyModelPre(
                data_path=os.path.join(self.data_path, 'healpix'),
                opsimFields=False,
                verbose=False)
        return self._sky_sim_hpx

    @property
    def sky_sim_field(self):
        """An lsst.sims.skybrightness_pre.SkyModelPre for LSST field values."""
        if self._sky_sim_field is None:
            self._sky_sim_field = lsst.sims.skybrightness_pre.SkyModelPre(
                data_path=os.path.join(self.data_path, 'opsimFields'),
                opsimFields=True,
                verbose=False)
        return self._sky_sim_field

    @property
    def nside(self):
        """The healpix nside for the healpix pre-computed sky values."""
        if self._nside is None:
            full = self.sky_sim_hpx.returnMags(
                59853.02, airmass_mask=False, planet_mask=False,
                moon_mask=False, zenith_mask=False)
            self._nside = healpy.npix2nside(full["r"].size)
        return self._nside

    def field_brightness(self, mjd, idx, band):
        """Return the pre-computed sky brightness of an LSST field.

        Args:
           - mjd :: the MJD for which to return the sky brightness
           - idx :: the LSST field id
           - filter :: the filter for which to return the sky brightness

        Return:
           the sky brightness in mag/asec^2
        """
        # sky indexes are 0 indexed, field indexes are 1 indexed
        sky_idx = idx-1
        sky_mag = self.sky_sim_field.returnMags(mjd, [sky_idx])[band]
        return sky_mag

    def hpx_brightness(self, mjd, idx, band):
        """Return the pre-computed sky brightness at a healpix on the sky.

        Args:
           - mjd :: the MJD for which to return the sky brightness
           - idx :: the healpix id of the desired location on the sky
           - filter :: the filter for which to return the sky brightness

        Return:
           the sky brightness in mag/asec^2
        """
        sky_mag = self.sky_sim_hpx.returnMags(mjd, [idx])[band]
        return sky_mag

    def coord_brightness(self, mjd, ra, decl, band):
        """Return the pre-computed sky brightness at a set of coordinates.

        Args:
           - mjd :: the MJD for which to return the sky brightness
           - ra :: the Right Ascension, in decimal degrees
           - decl :: the declination, in decimal degrees
           - filter :: the filter for which to return the sky brightness

        Return:
           the sky brightness in mag/asec^2
        """
        theta = np.radians(90-decl)
        phi = np.radians(ra)
        hpix = healpy.ang2pix(self.nside, theta, phi)
        sky_mag = self.hpx_brightness(mjd, hpix, band)
        return sky_mag

# internal functions & classes
