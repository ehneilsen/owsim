"""tests for followup."""
import unittest
from configparser import ConfigParser

import numpy as np
import pandas as pd

from astropy.coordinates import SkyCoord
import astropy.units as u
import astroplan

import followup

rand_seed = 6463


class TestFollowup(unittest.TestCase):
    def setUp(self):
        self.most_matches = 30
        self.least_matches = 10
        np.random.seed(rand_seed)
        self.config = ConfigParser()
        self.config.read('etc/followup.conf')
        self.fields = pd.read_csv('data/fieldID.txt', sep='\s+')
        self.events = pd.read_csv('data/events.txt', sep='\s+').query('episode==0')

    def test_search_radius(self):
        sr = followup.search_radius(60)
        self.assertAlmostEqual(sr, 6.12125378)

    def test_find_event_fields(self):
        ef = followup.find_event_fields(self.events, self.fields, 60)

        event_coords = SkyCoord(ef.ra.values*u.deg, ef.decl.values*u.deg)
        field_coords = SkyCoord(ef.fieldRA.values*u.deg,
                                ef.fieldDec.values*u.deg)
        distances = event_coords.separation(field_coords)
        self.assertLess(np.max(distances), 6.2*u.deg)

        num_matches_by_event = ef.groupby(level=0).mjd.count()
        self.assertLess(num_matches_by_event.max(), self.most_matches)
        self.assertGreater(num_matches_by_event.min(), self.least_matches)
        self.assertEqual(len(num_matches_by_event), len(self.events))

    def test_followup_targets(self):
        sample_events = self.events.sample(3)
        event_targets = followup.followup_targets(
            sample_events, self.fields, 60)
        for event, targets in event_targets:
            self.assertIsInstance(event, pd.Series)
            for target in targets:
                self.assertIsInstance(target, astroplan.FixedTarget)

    def test_schedule_followup(self):
        num_events = 3
        checked_events = 0
        while checked_events < num_events:
            sample_events = self.events.sample()
            visits = followup.schedule_followup(
                sample_events, self.fields, self.config)

            # If the sampled event is not observable, find another
            if len(visits) < 1:
                continue
            else:
                checked_events += 1

            # Time since event, in days
            visits['dt'] = (visits.mjd - sample_events.iloc[0].mjd)

            # Time since first visit, in hours
            visits['dt0'] = (visits.mjd - visits.iloc[0].mjd)*24

            # Check scan 1
            scan1 = visits.query('dt0 < 0.99999')
            self.assertGreater(scan1.dt.min(), 0, "Scan starts before event")
            self.assertLess(scan1.dt.min(), 1, "First scan started too late")
            self.assertLess(scan1.dt0.max(), 1, "First scan too long")
            nbands = len(scan1['filter'].unique())
            self.assertLessEqual(nbands, 2,
                                 "First scan had wrong number of filters")
            if nbands == 1:
                # One night strategy
                self.assertLess(scan1.dt.min(), 12,
                                "One night strategy first scan started too late")
                scan2 = visits.query('0.99999 <= dt0 < 2')
                self.assertEqual(len(scan1), len(scan2),
                                 "Different scans in one night strategy not equal")

            if nbands == 2:
                # Two night strategy
                scan2 = visits.query('dt >= 1')
                self.assertLessEqual(len(scan2.dt), len(scan1.dt))


if __name__ == '__main__':
    unittest.main()
