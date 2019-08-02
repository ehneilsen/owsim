"""tests for apsupp."""
import unittest

import astropy
from astropy.coordinates import SkyCoord
import astropy.units as u

import astroplan

from lsst.ts.observatory.model import ObservatoryModel

# Hack for when opsim is not really installed, but you
# have a checked out version of sims_speedObservatory
#
# import os
# import importlib.util
# opsim_telescope_python = os.environ['OPSIM_TELESCOPE_PYTHON']
# spec = importlib.util.spec_from_file_location(
#     "Telescope", opsim_telescope_python)
# telescope = importlib.util.module_from_spec(spec)
# spec.loader.exec_module(telescope)
# Telescope = telescope.Telescope

import apsupp
from lsst_conditions import SlewTimeSource

class TestDirectScheduler(unittest.TestCase):

#    @unittest.skip("skipping test_schedule")
    def test_schedule(self):
        pachon = astroplan.Observer.at_site("Cerro Pachon")
        start_time = astropy.time.Time(59598.66945625747, format='mjd')
        readout_time = 2 * u.second
        exptime = 15 * u.second
        nexp = 2
        slew_rate = 0.5*(0.66+0.57) * u.deg/u.second
        filter_change_time = {'filter': {'default': 120 * u.second}}

        transitioner = astroplan.Transitioner(slew_rate, filter_change_time)

        targets = [
            astroplan.FixedTarget(
                coord=SkyCoord(ra=81.1339111328125*u.deg,
                               dec=-17.532396316528303*u.deg),
                name="1857"),
            astroplan.FixedTarget(
                coord=SkyCoord(ra=81.1758499145508*u.deg,
                               dec=-12.2103328704834*u.deg),
                name="2100")]

        constraints = [astroplan.AltitudeConstraint(30*u.deg, 89*u.deg),
                       astroplan.AirmassConstraint(2.2),
                       astroplan.AtNightConstraint.twilight_astronomical()]
        obsblocks = [
            astroplan.ObservingBlock.from_exposures(
                target, 1, exptime, nexp, readout_time,
                configuration={'filter': band},
                constraints=constraints)
            for band in ('g', 'i') for target in targets]

        scheduler = apsupp.DirectScheduler(constraints=constraints,
                                           observer=pachon,
                                           transitioner=transitioner,
                                           time_resolution=1*u.hour)
        schedule = astroplan.scheduling.Schedule(
            start_time, start_time + 1*u.day)
        scheduler(obsblocks, schedule)

        self.assertIsInstance(schedule.slots[1].block.target,
                              astroplan.FixedTarget)
        self.assertEqual(schedule.slots[1].block.target.name,
                         '1857')
        self.assertIsInstance(schedule.slots[2].block,
                              astroplan.TransitionBlock)
        self.assertIsInstance(schedule.slots[3].block.target,
                              astroplan.FixedTarget)
        self.assertEqual(schedule.slots[3].block.target.name,
                         '2100')
        self.assertIsInstance(schedule.slots[4].block,
                              astroplan.TransitionBlock)
        self.assertIsInstance(schedule.slots[5].block.target,
                              astroplan.FixedTarget)
        self.assertEqual(schedule.slots[5].block.target.name,
                         '1857')
        self.assertIsInstance(schedule.slots[6].block,
                              astroplan.TransitionBlock)
        self.assertIsInstance(schedule.slots[7].block.target,
                              astroplan.FixedTarget)
        self.assertEqual(schedule.slots[7].block.target.name, '2100')


class TestLSSTTransitioner(unittest.TestCase):

    def test_transitioner(self):
        pachon = astroplan.Observer.at_site("Cerro Pachon")
        start_time = astropy.time.Time(59598.66945625747, format='mjd')
        readout_time = 2 * u.second
        exptime = 15 * u.second
        nexp = 2

        targets = [
            astroplan.FixedTarget(
                coord=SkyCoord(ra=81.1339111328125*u.deg,
                               dec=-17.532396316528303*u.deg),
                name="1857"),
            astroplan.FixedTarget(
                coord=SkyCoord(ra=81.1758499145508*u.deg,
                               dec=-12.2103328704834*u.deg),
                name="2100")]

        constraints = [astroplan.AltitudeConstraint(30*u.deg, 89*u.deg),
                       astroplan.AirmassConstraint(2.2),
                       astroplan.AtNightConstraint.twilight_astronomical()]
        obs_blocks = [
            astroplan.ObservingBlock.from_exposures(
                target, 1, exptime, nexp, readout_time,
                configuration={'filter': band},
                constraints=constraints)
            for band in ('g', 'i') for target in targets]

        transitioner = apsupp.LSSTTransitioner()
        for old_block, new_block in zip(obs_blocks[:-1], obs_blocks[1:]):
            block = transitioner(old_block, new_block, start_time, pachon)
            duration_seconds = (block.duration/u.second).value

            old_filter = old_block.configuration['filter']
            new_filter = new_block.configuration['filter']
            if old_filter == new_filter:
                min_expected_duration = 2
            else:
                min_expected_duration = 120

            self.assertGreaterEqual(duration_seconds, min_expected_duration)
            self.assertLess(duration_seconds, 600)


class TestOpsimTransitioner(unittest.TestCase):

    def test_transitioner(self):
        observatory = ObservatoryModel()
        observatory.configure_from_module()
        pachon = astroplan.Observer.at_site("Cerro Pachon")
        start_time = astropy.time.Time(59598.66945625747, format='mjd')
        readout_time = observatory.params.readouttime * u.second
        exptime = 15 * u.second
        nexp = 2

        targets = [
            astroplan.FixedTarget(
                coord=SkyCoord(ra=81.1339111328125*u.deg,
                               dec=-17.532396316528303*u.deg),
                name="1857"),
            astroplan.FixedTarget(
                coord=SkyCoord(ra=81.1758499145508*u.deg,
                               dec=-12.2103328704834*u.deg),
                name="2100")]

        constraints = [astroplan.AltitudeConstraint(30*u.deg, 89*u.deg),
                       astroplan.AirmassConstraint(2.2),
                       astroplan.AtNightConstraint.twilight_astronomical()]
        obs_blocks = [
            astroplan.ObservingBlock.from_exposures(
                target, 1, exptime, nexp, readout_time,
                configuration={'filter': band},
                constraints=constraints)
            for band in ('g', 'i') for target in targets]

        self.assertGreater(observatory.params.readouttime, 1.5)
        self.assertGreater(observatory.params.filter_changetime, 15.0)

        slew_src = SlewTimeSource()
        transitioner = apsupp.OpsimTransitioner(slew_src)
        whitepaper_transitioner = apsupp.LSSTTransitioner()
        for old_block, new_block in zip(obs_blocks[:-1], obs_blocks[1:]):
            block = transitioner(old_block, new_block, start_time, pachon)
            wp_block = whitepaper_transitioner(
                old_block, new_block, start_time, pachon)
            duration_seconds = (block.duration/u.second).value
            wp_duration_seconds = (wp_block.duration/u.second).value

            self.assertAlmostEqual(
                duration_seconds, wp_duration_seconds, delta=5)

            old_filter = old_block.configuration['filter']
            new_filter = new_block.configuration['filter']
            if old_filter == new_filter:
                min_expected_duration = observatory.params.readouttime
            else:
                min_expected_duration = max(observatory.params.readouttime,
                                            observatory.params.filter_changetime)

            self.assertGreaterEqual(duration_seconds, min_expected_duration)
            self.assertLess(duration_seconds, 600)


if __name__ == '__main__':
    unittest.main()
