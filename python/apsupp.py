"""Classes to supplement astroplan with new functionality"""

from copy import deepcopy
import unittest
from collections import OrderedDict

import pandas as pd

import astropy.units as u
import astroplan
from astroplan import Scheduler
import logging
from logging import debug, info, warning, error, critical

#logging.basicConfig(format='%(asctime)s %(message)s', level=logging.DEBUG)

class DirectScheduler(Scheduler):
    """astroplan.Scheduler that schedules in the order provided

    Implements and astroplan scheduler that preserves the order of
    exposures, just shifting the block of exposures until it can
    schedule as many consecutively as it can, and inserting
    appropriate transition blocks between supplied target blocks.

    """
    
    def _make_schedule(self, blocks):

        global_constraints = [] if self.constraints is None else self.constraints
        def block_allowed(block, t):
            block_constraints = [] if block.constraints is None else block.constraints
            constraints = set(global_constraints + block_constraints)
            for constraint in constraints:
                if not constraint(self.observer, block.target, t):
                    return False
            return True

        def try_schedule(start_time, end_time):
            debug(f'Testing in slot {start_time.iso} to {end_time.iso}') 
            schedule = Schedule(start_time, end_time)
            if end_time - start_time < 3*u.second:
                return 0, schedule
            
            current_time = start_time + 0.01*u.second
            number_scheduled = 0
            for block in blocks:
                block_start = current_time
                if len(schedule.scheduled_blocks) > 0:
                    transition_block = self.transitioner(
                        schedule.scheduled_blocks[-1], block,
                        current_time, self.observer)
                    block_start += transition_block.duration
                else:
                    transition_block = None

                # If we cannot finish this block, go on to the next
                block_end = block_start + block.duration
                if block_end > end_time:
                    continue
                    
                if block_allowed(block, block_start):
                    # Add a transition block, if needed
                    if transition_block is not None:
                        debug(f'Testing transition insert at {current_time.iso} for {transition_block.duration}')
                        schedule.insert_slot(current_time, transition_block)
                        current_time += transition_block.duration

                    debug(f'Testing insert at {current_time.iso} for {block.duration}') 
                    schedule.insert_slot(current_time, block)
                    number_scheduled += 1
                    current_time += block.duration
                    
            return number_scheduled, schedule


        # Work through open slots and start times until we find one at
        # which we can observe the whole set.  If there are none,
        # return the earliest instance where we get as many as
        # possible.
        number_scheduled = 0
        best_schedule = None
        for slot in self.schedule.open_slots:
            start_time = slot.start
            current_time = start_time

            while number_scheduled < len(blocks) and start_time < slot.end:
                new_number_scheduled, schedule = try_schedule(start_time, slot.end)
                if new_number_scheduled > number_scheduled:
                    number_scheduled = new_number_scheduled
                    best_schedule = schedule

                start_time += self.time_resolution

            if number_scheduled == len(blocks):
                break

        if best_schedule is None:
            return
            
        # Find the last occupied slot before the first we are trying to
        # schedule to add a transition, if necessary
        prior_slot = None
        first_slot = best_schedule.slots[0]
        for slot in self.schedule.slots:
            if not slot.occupied:
                continue
            if slot.start >= first_slot.start:
                break
            prior_slot = slot

        sched_shift = 0 * u.second
        if prior_slot is not None:
            prior_block = prior_slot.block
            transition_block = self.transitioner(
                prior_block, first_slot.block,
                prior_slot.end, self.observer)
            transition_end = proir_slot.end + transition_block.duration
            if transition_end > first_slot.start:
                sched_shift = transition_end - first_slot.start
                self.schedule.insert_slot(prior_slot.end, transition_block)            
        
        # Add blocks from the best schedule to the final schedule
        for slot in best_schedule.slots:
            if slot.block is not None:
                self.schedule.insert_slot(slot.start+sched_shift, slot.block)


class Schedule(astroplan.Schedule):
    def to_df(self):
        start_time = []
        mjd = []
        target_name = []
        ra = []
        decl = []
        exptime = []
        band = []
        for slot in self.slots:
            if hasattr(slot.block, 'target'):
                block = slot.block
                target = block.target
                for exp_num in range(block.number_exposures):
                    frac_exposures = exp_num/block.number_exposures
                    exp_start_time = slot.start \
                                     + exp_num * block.duration * frac_exposures
                    start_time.append(exp_start_time.iso.replace(' ', 'T'))
                    mjd.append(exp_start_time.mjd)
                    target_name.append(target.name)
                    ra.append(target.ra.deg)
                    decl.append(target.dec.deg)
                    band.append(block.configuration['filter'])
                    exptime.append(block.time_per_exposure.to(u.second).value)

        exposures = OrderedDict()
        exposures['start_time'] = start_time
        exposures['mjd'] = mjd
        exposures['target'] = target_name
        exposures['ra'] = ra
        exposures['decl'] = decl
        exposures['filter'] = band
        exposures['exptime'] = exptime

        df = pd.DataFrame(exposures)
        return df
