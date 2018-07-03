"""Classes to supplement astroplan with new functionality"""

from copy import deepcopy
import unittest

import astroplan
from astroplan import Scheduler

class DirectScheduler(Scheduler):
    """astroplan.Scheduler that schedules in the order provided

    Implements and astroplan scheduler that preserves the order of
    exposures, just shifting the block of exposures until it can
    schedule as many consecutively as it can, and inserting
    appropriate transition blocks between supplied target blocks.

    """
    
    def _make_schedule(self, blocks):

        if len(self.schedule.scheduled_blocks) > 0:
            raise NotImplementedError
        
        global_constraints = [] if self.constraints is None else self.constraints
        def block_allowed(block, t):
            block_constraints = [] if block.constraints is None else block.constraints
            constraints = set(global_constraints + block_constraints)
            for constraint in constraints:
                if not constraint(self.observer, block.target, t):
                    return False
            return True

        def try_schedule(start_time):
            schedule = deepcopy(self.schedule)
            current_time = start_time
            number_scheduled = 0
            for block in blocks:
                if block_allowed(block, current_time):
                    # Add a transition block, if needed
                    if len(schedule.scheduled_blocks) > 0:
                        transition_block = self.transitioner(
                            schedule.scheduled_blocks[-1], block,
                            current_time, self.observer)
                        schedule.insert_slot(current_time, transition_block)
                        current_time += transition_block.duration

                    schedule.insert_slot(current_time, block)
                    current_time += block.duration
                    number_scheduled += 1
                    
            return number_scheduled, schedule

        start_time = self.schedule.start_time
        number_scheduled = 0
        
        # Work through start times until we find one at which we can
        # observe the whole set.  If there are none, return the
        # earliest instance where we get as many as possible.
        
        while number_scheduled < len(blocks) and start_time < self.schedule.end_time:
            new_number_scheduled, schedule = try_schedule(start_time)
            if new_number_scheduled > number_scheduled:
                number_scheduled = new_number_scheduled
                self.schedule.slots = schedule.slots

            start_time += self.time_resolution


