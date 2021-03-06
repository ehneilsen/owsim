"""Calculate circumstances for a set of visits."""
import argparse
from configparser import ConfigParser
import code
import logging
import sys
import os
from logging import debug, info, warning, error, critical
from collections import OrderedDict

import numpy as np
import pandas as pd

import astropy.utils.iers

# Fix broken URL for IERS_A
#IERS_A_URL = 'http://toshi.nofs.navy.mil/ser7/finals2000A.all'
#astropy.utils.iers.IERS_A_URL = IERS_A_URL
#astropy.utils.iers.iers.IERS_A_URL = IERS_A_URL
#astropy.utils.iers.iers.conf.iers_auto_url = IERS_A_URL


import astropy
import astropy.time
import astropy.coordinates
from astropy.coordinates import SkyCoord
from astropy.coordinates import EarthLocation
import astropy.units as u


from lsst_conditions import SeeingSource
from lsst_conditions import SlewTimeSource
from lsst_conditions import CloudSource
from lsst_conditions import SkyBrightnessSource
from lsst_conditions import calc_five_sigma_depth
from lsst_conditions import SKY_DATA_PATH
from lsst_conditions import DEFAULT_START_TIME

try:
    import skybright
except ModuleNotFoundError:
    pass

# constants

SUMMARYALLPROPS_TYPES = {
    "observationId": int,
    "night": int,
    "observationStartTime": float,
    "observationStartMJD": float,
    "observationStartLST": float,
    "filter": (np.str_, 1),
    "proposals": (np.str_, 150),
    "fieldId": int,
    "fieldRA": float,
    "fieldDec": float,
    "altitude": float,
    "azimuth": float,
    "numExposures": int,
    "visitTime": float,
    "visitExposureTime": float,
    "airmass": float,
    "skyBrightness": float,
    "cloud": float,
    "seeingFwhm500": float,
    "seeingFwhmGeom": float,
    "seeingFwhmEff": float,
    "fiveSigmaDepth": float,
    "moonRA": float,
    "moonDec": float,
    "moonAlt": float,
    "moonAz": float,
    "moonDistance": float,
    "moonPhase": float,
    "sunRA": float,
    "sunDec": float,
    "sunAlt": float,
    "sunAz": float,
    "solarElong": float,
    "slewTime": float,
    "slewDistance": float,
    "paraAngle": float,
    "rotTelPos": float,
    "rotSkyPos": float
}

# exception classes

# interface functions


def visit_circumstances(visits,
                        site=EarthLocation.of_site('Cerro Pachon'),
                        survey_start_time=DEFAULT_START_TIME,
                        sky_data_path=SKY_DATA_PATH,
                        skybright_config=None):
    """Calculate conditions for a set of visits.

    Args:
       - visits :: a pandas.DataFrame of visits
       - site :: on astropy.coordinates.EarthLocation
                 for the observing site
       - survey_start_time :: the astropy.time.Time of the start
                 of the survey
       - sky_data_path :: the path to precomputed sky brightness data

    Returns:
       a pandas.DataFrame of visits with observing conditions

    """
    info('Calculating astropy time')
    t = astropy.time.Time(
        visits.mjd, format='mjd', scale='utc', location=site)
    info('Calculating POSIX time')
    clocktime = t.unix
    info('Calculating seconds into survey')
    seconds_into_survey = clocktime - survey_start_time.unix
    info('Calculating LST')
    lst = t.sidereal_time('mean')    
    info('Calculating night')
    night = calc_night(visits.mjd, site)

    info('Creating astropy SkyCoords')
    coords = SkyCoord(ra=visits.ra.values*u.degree,
                      dec=visits.decl.values*u.degree,
                      frame='icrs')
    info('Creating astropy SkyCoords in alt, az')
    alt_az = astropy.coordinates.AltAz(obstime=t, location=site)
    alt_az_coords = coords.transform_to(alt_az)
    info('Calculating airmass')
    airmass = calc_airmass(alt_az_coords)

    info('Calculating sun coordinates')
    sun_coords = astropy.coordinates.get_sun(t)
    info('Calculating sun alt and az')
    sun_alt_az_coords = sun_coords.transform_to(alt_az)
    info('Calculating sun separation')
    sun_distance = coords.separation(sun_coords)

    info('Calculating moon coordinates')
    moon_coords = astropy.coordinates.get_moon(t, site)
    info('Calculating moon alt and az')
    moon_alt_az_coords = moon_coords.transform_to(alt_az)
    info('Calculating moon separation')
    moon_distance = coords.separation(moon_coords)
    info('Calculating moon phase')
    moon_phase = calc_moon_phase(sun_coords, moon_coords)

    info('Calculating seeing')
    seeing_source = np.vectorize(SeeingSource(survey_start_time))
    seeing_array = seeing_source(
        seconds_into_survey, visits['filter'], airmass)
    seeing = pd.DataFrame.from_items(
        zip(['fwhm500', 'fwhmGeom', 'fwhmEff'], seeing_array))

    if skybright_config is not None:
        info('Use skybright to calculate skybrightness')
        calc_coord_brightness = skybright.MoonSkyModel(skybright_config)
    else:
        info('Looking up sky brightness')
        calc_coord_brightness = np.vectorize(
            SkyBrightnessSource(sky_data_path).coord_brightness)

    sky_brightness = calc_coord_brightness(
        visits.mjd, visits.ra, visits.decl, visits['filter'])

    info('Looking up cloud level')
    cloud_source = np.vectorize(CloudSource(survey_start_time))
    clouds = cloud_source(seconds_into_survey)

    info('Calculationg 5-sigma depth')
    five_sigma_depth = np.zeros(len(visits), dtype=float)
    for band in visits['filter'].unique():
        band_idxs = np.where(visits['filter'].values == band)
        five_sigma_depth[band_idxs] = calc_five_sigma_depth(
            band,
            sky_brightness[band_idxs],
            seeing.fwhmEff.iloc[band_idxs],
            visits.exptime.iloc[band_idxs],
            airmass[band_idxs])

    info('Calculating slew time')
    pre_alt = alt_az_coords.alt.deg[:-1]
    pre_az = alt_az_coords.az.deg[:-1]
    pre_band = visits['filter'][:-1]
    post_alt = alt_az_coords.alt.deg[1:]
    post_az = alt_az_coords.az.deg[1:]
    post_band = visits['filter'][1:]
    slew_time = np.zeros(len(t), dtype=float)
    calc_slew_time = np.vectorize(SlewTimeSource())
    slew_time[1:] = calc_slew_time(seconds_into_survey[1:],
        pre_alt, pre_az, pre_band, post_alt, post_az, post_band)

    info('Calculating slew distance')
    slew_distance = np.zeros(len(t), dtype=float)
    slew_distance[1:] = coords[:-1].separation(coords[1:])

    circ_dict = OrderedDict((
        ('observationId', np.arange(len(visits))),
        ('night', night),
        ('observationStartTime', clocktime),
        ('observationStartMJD', visits.mjd),
        ('observationStartLST', lst.deg),
        ('filter', visits['filter']),
        ('proposals', visits.proposals),
        ('fieldId', visits.target.astype(int)),
        ('fieldRA', visits.ra),
        ('fieldDec', visits.decl),
        ('altitude', alt_az_coords.alt.deg),
        ('azimuth', alt_az_coords.az.deg),
        ('numExposures', visits.nexp.astype(int)),
        ('visitTime', visits.duration),
        ('visitExposureTime', visits.exptime),
        ('airmass', airmass),
        ('skyBrightness', sky_brightness),
        ('cloud', clouds),
        ('seeingFwhm500', seeing.fwhm500.values),
        ('seeingFwhmGeom', seeing.fwhmGeom.values),
        ('seeingFwhmEff', seeing.fwhmEff.values),
        ('fiveSigmaDepth', five_sigma_depth),
        ('moonRA', moon_coords.ra.deg),
        ('moonDec', moon_coords.dec.deg),
        ('moonAlt', moon_alt_az_coords.alt.deg),
        ('moonAz', moon_alt_az_coords.az.deg),
        ('moonDistance', moon_distance.deg),
        ('moonPhase', moon_phase),
        ('sunRA', sun_coords.ra.deg),
        ('sunDec', sun_coords.dec.deg),
        ('sunAlt', sun_alt_az_coords.alt.deg),
        ('sunAz', sun_alt_az_coords.az.deg),
        ("solarElong", sun_distance),
        ("slewTime", slew_time),
        ("slewDistance", slew_distance),
        ("paraAngle", None),
        ("rotTelPos", None),
        ("rotSkyPos", None)
        ))

    circumstances = pd.DataFrame(circ_dict, index=visits.index.values)
    circumstances = circumstances.astype(SUMMARYALLPROPS_TYPES)
    return circumstances

# classes

# internal functions & classes


def calc_night(mjd, site):
    """Determint the night MJD for a given floating mjd at a given site.

    Args:
       - mjd :: the floating point UT MJD.
       - site :: an astropy.coordinates.EarthLocation
                 for the observing site

    Returns:
       an intereger MJD for the night on which mjd falls at
       the requested site
    """
    lon_offset = site.lon.deg/360
    first_noon = np.floor(mjd.min() - lon_offset) + lon_offset
    nights = 1 + np.floor(mjd-first_noon).astype(int)
    return nights


def calc_airmass(alt_az_coords):
    """Calculate the airmass from horizon coordinates.

    This function uses Hiltner's formula for airmass.

    Args:
       - alt_az_coords :: the astropy.coordinates.SkyCoord
           for which to calculate airmass

    Returns:
       the airmass

    """
    # Hiltner's formula for airmass
    x0 = 1/np.sin(alt_az_coords.alt.radian)
    dx = x0 - 1
    x = x0 - (0.0018167*dx) - (0.002875*(dx**2)) - (0.0008083*(dx**3))
    return x


def calc_moon_phase(sun_coords, moon_coords):
    """Calculate the moon phase.

    Args:
       - sun_coords :: the astropy.coordinates.SkyCoords of the sun
       - moon_coords :: the astropy.coordinates.SkyCoords of the moon

    Returns:
       the moon phase
    """
    elongation = sun_coords.separation(moon_coords)
    moon_phase = np.arctan2(
        sun_coords.distance * np.sin(elongation),
        moon_coords.distance - sun_coords.distance*np.cos(elongation))
    return moon_phase


def _main():
    parser = argparse.ArgumentParser()
    parser.add_argument('config', help='configuration file name')
    parser.add_argument('visits', help='file from which to read visits')
    parser.add_argument('circumstances',
                        help='file to which to write visit circumstances')
    if 'skybright' in sys.modules:
        parser.add_argument('-s', '--skybright',  default='',
                            help='Use the skybright package with this config file')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='Print less information')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print more information')
    parser.add_argument('-i', '--interactive', action='store_true',
                        help='Drop into an interactive shell after test.')

    args = parser.parse_args()
    if not args.quiet:
        logging.basicConfig(format='%(asctime)s %(message)s',
                            level=logging.INFO)

    if args.verbose:
        logging.basicConfig(format='%(asctime)s %(message)s',
                            level=logging.DEBUG)

    if 'skybright' in sys.modules and len(args.skybright)>0:
        skybright_config = ConfigParser()
        skybright_config.read(args.skybright)
    else:
        skybright_config = None

    interactive = args.interactive
    visit_fname = args.visits
    circumstances_fname = args.circumstances

    config = ConfigParser()
    config.read(args.config)

    site = EarthLocation.of_site(config['site']['name'])

    info("Loading visits")
    visits = pd.read_table(args.visits, sep=" ")

    if interactive:
        expsamp = visits.iloc[:5]
        print("visit_circumstances(expsamp, site)")

        vars = globals().copy()
        vars.update(locals())
        shell = code.InteractiveConsole(vars)
        shell.interact()
    else:
        info("Scheduling followup visits")
        circumstances = visit_circumstances(
            visits, site, skybright_config=skybright_config)

        info("Writing followup visits")
        _, extension = os.path.splitext(circumstances_fname)
        circumstances.to_csv(circumstances_fname,
                             sep=" ", index=False, header=True)

    return 0


if __name__ == '__main__':
    status = _main()
    sys.exit(status)
