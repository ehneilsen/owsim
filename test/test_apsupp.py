"""tests for apsupp"""
import unittest

import astropy
from astropy.coordinates import SkyCoord
import astropy.units as u

import astroplan
from astroplan import Scheduler
from astroplan import AltitudeConstraint

import apsupp

class TestDirectScheduler(unittest.TestCase):
    def test_schedule(self):
        pachon = astroplan.Observer.at_site("Cerro Pachon")
        start_time = astropy.time.Time(59598.66945625747, format='mjd')
        readout_time = 3 * u.second
        exptime = 15 * u.second
        nexp = 2
        slew_rate = 1.0 * u.deg/u.second
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
                target, 1, exptime, nexp, readout_time, configuration={'filter': band},
                constraints=constraints)
            for band in ('g', 'i') for target in targets]

        scheduler = apsupp.DirectScheduler(constraints=constraints,
                                           observer=pachon,
                                           transitioner=transitioner,
                                           time_resolution=1*u.hour)
        schedule = astroplan.scheduling.Schedule(
            start_time, start_time + 1*u.day)
        scheduler(obsblocks, schedule)

        self.assertIsInstance(schedule.slots[1].block.target, astroplan.FixedTarget)
        self.assertEqual(schedule.slots[1].block.target.name, '1857')
        self.assertIsInstance(schedule.slots[2].block, astroplan.TransitionBlock)
        self.assertIsInstance(schedule.slots[3].block.target, astroplan.FixedTarget)
        self.assertEqual(schedule.slots[3].block.target.name, '2100')
        self.assertIsInstance(schedule.slots[4].block, astroplan.TransitionBlock)
        self.assertIsInstance(schedule.slots[5].block.target, astroplan.FixedTarget)
        self.assertEqual(schedule.slots[5].block.target.name, '1857')
        self.assertIsInstance(schedule.slots[6].block, astroplan.TransitionBlock)
        self.assertIsInstance(schedule.slots[7].block.target, astroplan.FixedTarget)
        self.assertEqual(schedule.slots[7].block.target.name, '2100')
        
if __name__ == '__main__':
    unittest.main()
