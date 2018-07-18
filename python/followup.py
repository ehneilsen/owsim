"""Create tables of followup visits."""
import argparse
from configparser import ConfigParser
import readline
import code
import logging
import sys
import csv
from logging import debug, info, warning, error, critical
from collections import OrderedDict

import numpy as np
import pandas as pd

import astropy
import astropy.coordinates
from astropy.coordinates import SkyCoord
import astropy.units as u

import astroplan

import apsupp
from apsupp import DirectScheduler

from lsst.sims.speedObservatory import Telescope


# constants

# exception classes

class StrategyFailed(Exception):
    pass

# interface functions


def schedule_followup(events, fields, config):
    """Schedule follow-up for a given set of events.

    Args:
       - events :: a pandas.DataFrame of events
       - fields :: a pandas.DataFrame of fields for followup
       - config :: a ConfigParser with follow-up configuration data.

    Returns:
       a pandas.DataFrame with the schedule of follow-up visits.
    """
    readout_time = float(config['instrument']['readout_time']) * u.second
    shutter_time = float(config['instrument']['shutter_time']) * u.second
    exptime = float(config['obs_block']['exptime']) * u.second
    nexp = int(config['obs_block']['nexp'])
    bands = config['search']['bands'].split()

    #
    # Create an astroplan.Observer
    #
    observer = astroplan.Observer.at_site(config['site']['name'])

    #
    # Create a Transitioner
    #
    telescope = Telescope()
    transitioner = apsupp.OpsimTransitioner(telescope)

    #
    # Build astroplan constraints
    #
    alt_constraint = astroplan.AltitudeConstraint(
        float(config['constraints']['min_alt'])*u.deg,
        float(config['constraints']['max_alt'])*u.deg)
    airmass_constraint = astroplan.AirmassConstraint(
        float(config['constraints']['max_airmass']))

    constraints = [alt_constraint, airmass_constraint]

    twilight_type = config['constraints']['twilight']
    if twilight_type == 'nautical':
        constraints.append(astroplan.AtNightConstraint.twilight_nautical())
    elif twilight_type == 'civil':
        constraints.append(astroplan.AtNightConstraint.twilight_civil())
    elif twilight_type == 'astronomical':
        constraints.append(astroplan.AtNightConstraint.twilight_astronomical())
    else:
        raise NotImplementedError

    #
    # Create an astroplan.Scheduler
    #
    time_resolution = float(config['scheduler']['time_resolution']) * u.second
    scheduler = DirectScheduler(
        constraints=constraints, observer=observer,
        transitioner=transitioner, time_resolution=1*u.hour)

    # helper function to turn targets into blocks
    def obsblock(target, band):
        block = astroplan.ObservingBlock.from_exposures(
            target, 1, exptime, nexp, readout_time,
            configuration={'filter': band})
        block.shutter_time = shutter_time
        block.duration = block.duration \
                         + block.number_exposures * block.shutter_time
        return block

    # Schedule follow-up for each event
    latest_visit_end = astropy.time.Time(0, format='mjd')
    search_area = float(config['search']['area'])
    scan_schedules = []
    for event, targets in followup_targets(events, fields, search_area):
        event_time = astropy.time.Time(np.min(event.mjd.min()), format='mjd')
        trigger_time = event_time \
                       + float(config['search']['trigger_delay'])*u.second
        info(f'Scheduling event {event.event_id} at {event.ra}, {event.decl} on {event_time.iso}')

        # Only work on one event at a time
        if trigger_time < latest_visit_end:
            continue

        # Attempt the first strategy -- two scans the area in one
        # band, one hour apart, within 12 hours of the event.
        # If this isn't possible, two scans with larger separation,
        # in two bands each.
        try:
            obsblocks = [obsblock(t, bands[0]) for t in targets]
            schedule = two_scan_schedule(scheduler, obsblocks,
                                         trigger_time, event_time,
                                         obs_window=12*u.hour,
                                         scan_separation=1*u.hour)
        except StrategyFailed:
            obsblocks = [obsblock(t, b) for b in bands for t in targets]
            schedule = two_scan_schedule(scheduler, obsblocks,
                                         trigger_time, event_time,
                                         obs_window=24*u.hour,
                                         scan_separation=1*u.hour,
                                         ignore_failure=True)

        scan_schedules.append(schedule)

    schedule = pd.concat([s for s in scan_schedules],
                         axis=0)

    schedule['proposals'] = config['scheduler']['proposals']
    return schedule

# classes

# internal functions & classes


# 1.75 degrees is the radius of the LSST field of view in degrees
def search_radius(area, fov_radius=1.75):
    """Return the radius within which to pick pointings to cover area on the sky

    Args:
       - area :: the area, in square degrees
       - fov_radius :: the field of view of one pointing
    """
    area_rad = np.radians(np.radians(area))
    radius_rad = np.arccos(1-area_rad/(2*np.pi))
    radius = np.degrees(radius_rad) + fov_radius
    return radius


def followup_targets(events, fields, area):
    """A generator of event/astroplan target pairs given events and fields

    All ra and dec values should be in decimal degrees.

    Args:
       - events :: a pandas.DataFrame with event coordinates
                   in ra and decl columns
       - fields :: a pandas.DataFrame with field coordinates
                   in fieldRA and fieldDec columns
       - area :: the search area in square degrees

    Returns:
       a succession of event/astroplan.FixedTarget pairs
    """

    event_fields = find_event_fields(events, fields, area)
    targets = OrderedDict()

    # there should be a faster (vectorized) way to do this.
    for idx, event in events.iterrows():
        targets = []
        for field_idx, field in event_fields.query(
                f'event_id=={event.event_id}').iterrows():
            field_coords = SkyCoord(ra=field.fieldRA*u.deg,
                                    dec=field.fieldDec*u.deg)
            field_name = "%d" % field.fieldID
            target = astroplan.FixedTarget(coord=field_coords, name=field_name)
            targets.append(target)

        yield event, targets


def find_event_fields(events, fields, area):
    """Create a pandas.DataFrame of fields covering an area around events.

    All ra and dec values should be in decimal degrees.

    Args:
       - events :: a pandas.DataFrame with event coordinates
                   in ra and decl columns
       - fields :: a pandas.DataFrame with field coordinates
                   in fieldRA and fieldDec columns
       - area :: the search area in square degrees

    Returns:
       a pandas.DataFrame with field/event pairs and separations.
    """

    event_coords = SkyCoord(events.ra.values*u.deg,
                            events.decl.values*u.deg)
    field_coords = SkyCoord(fields.fieldRA.values*u.deg,
                            fields.fieldDec.values*u.deg)
    field_id, event_idx, d2d, d3d = event_coords.search_around_sky(
        field_coords, search_radius(area)*u.deg)

    event_fields = pd.DataFrame(
        {'event_id': events.event_id.iloc[event_idx].values,
         'mjd': events.mjd.iloc[event_idx].values,
         'ra': events.ra.iloc[event_idx].values,
         'decl': events.decl.iloc[event_idx].values,
         'distance': d2d.deg,
         'fieldID': fields.fieldID.iloc[field_id].values,
         'fieldRA': fields.fieldRA.iloc[field_id].values,
         'fieldDec': fields.fieldDec.iloc[field_id].values
         })

    event_fields.sort_values(['event_id', 'fieldRA'], inplace=True)
    event_fields = event_fields[
        ['event_id', 'mjd', 'ra', 'decl',
         'fieldID', 'fieldRA', 'fieldDec', 'distance']]
    event_fields.set_index(['event_id'], inplace=True)

    return event_fields


def two_scan_schedule(scheduler, obsblocks, trigger_time, event_time,
                      obs_window=12*u.hour, scan_separation=1*u.hour,
                      ignore_failure=False):
    """Schedule two scans of a set of observing blocks

    Args:
       - scheduler :: and astroplan.scheduler to schedule each scan
       - obsblocks  :: a sequence of astroplan.blocks to schedule
       - trigger_time :: reference time for the observing window
       - obs_window :: duration of window in which to schedule the first scan
       - scan_separation :: the duration of the gap between the starts of scans

    trigger_time is an astropy.time.Time object. obs_window and
    scan_separation are in astropy time units.

    Returns:
       A pandas.DataFrame holding a table of visits

    """

    scan_time = [trigger_time, 0]
    scan_schedules = []

    info(f"Scheduling scan 1")
    scan_time[0] = trigger_time
    scan_schedule = apsupp.Schedule(
        scan_time[0], scan_time[0] + obs_window)
    scheduler(obsblocks, scan_schedule)
    info(f"Scheduled {len(scan_schedule.observing_blocks)} out of {len(obsblocks)}")
    scan_schedules.append(scan_schedule)

    # Only reobserve blocks we observed the first time
    obsblocks = scan_schedule.observing_blocks

    # Start an hour after we stated the previous scan, if possible
    for slot in scan_schedule.slots:
        if slot.occupied:
            scan_time[1] = slot.start + scan_separation
            break

    info(f"Scheduling scan 2")
    scan_schedule = apsupp.Schedule(
        scan_time[1], scan_time[1] + obs_window)
    scheduler(obsblocks, scan_schedule)

    if not ignore_failure:
        if len(scan_schedule.observing_blocks) < len(obsblocks):
            raise StrategyFailed

    scan_schedules.append(scan_schedule)

    schedule = pd.concat([s.to_df() for s in scan_schedules],
                         axis=0)

    return schedule


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('config', help='configuration file name')
    parser.add_argument('events', help='file from which to read events')
    parser.add_argument('fields', help='file from which to read fields')
    parser.add_argument('visits', help='file to which to write visits')
    parser.add_argument('-e', '--episode', type=int, default=0,
                        help='the episode index')
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

    interactive = args.interactive
    visit_fname = args.visits
    episode = args.episode

    config = ConfigParser()
    config.read(args.config)

    info("Loading fields")
    fields = pd.read_table(args.fields)

    info("Loading events")
    events = pd.read_table(args.events).query(f'episode == {episode}')

    if interactive:
        vars = globals().copy()
        vars.update(locals())
        shell = code.InteractiveConsole(vars)
        shell.interact()
    else:
        info("Scheduling followup visits")
        visits = schedule_followup(events, fields, config)

        info("Writing followup visits")
        visits.to_csv(visit_fname,
                      sep=" ",
                      quoting=csv.QUOTE_MINIMAL,
                      index=False,
                      header=True)

    return 0


if __name__ == '__main__':
    status = main()
    sys.exit(status)
