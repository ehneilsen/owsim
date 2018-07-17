"""tests for apsupp"""
import unittest

import numpy as np
import pandas as pd
import owsched

test_db_file_name = 'data/baseline2018a.db'
test_replacement_exposures_fname = 'data/conditions.txt'

# import code
# vars = globals().copy()
# vars.update(locals())
# shell = code.InteractiveConsole(vars)
# shell.interact()        

class TestOwSched(unittest.TestCase):

    @unittest.skip('')
    def test_schedule(self):
        dfs = owsched.query_summary(test_db_file_name)
        self.assertIn('propId', dfs.Proposal.columns)
        self.assertIn('propName', dfs.Proposal.columns)
        self.assertIn('observationStartTime', dfs.SummaryAllProps.columns)
        self.assertIn('fieldRA', dfs.SummaryAllProps.columns)
        self.assertIn('fieldDec', dfs.SummaryAllProps.columns)
        self.assertIn('visitTime', dfs.SummaryAllProps.columns)
        self.assertIn('slewTime', dfs.SummaryAllProps.columns)
        self.assertIn('filter', dfs.SummaryAllProps.columns)

    @unittest.skip('')
    def test_split_into_sequences(self):
        exposures = pd.read_table(test_replacement_exposures_fname)
        sequences = [s[1] for s in owsched.split_into_sequences(exposures)]

        for i, seq in enumerate(sequences[1:]):
            prev_seq = sequences[i]
            gap_days = seq.observationStartMJD.min() \
                       - prev_seq.observationStartMJD.max()
            gap_seconds = 24*60*60*gap_days
            self.assertGreater(gap_seconds, 120)

    def test_overwrite_schedule(self):
        reference = owsched.query_summary(test_db_file_name).SummaryAllProps

        # Force all replacement exposures to have a common, new proposalId
        # so we can easily extract them after the overwrite.
        test_proposal_id = reference.proposalId.max() + 1
        replacements = pd.read_table(test_replacement_exposures_fname)
        replacements.proposalId = test_proposal_id
        replacements.observationId = replacements.observationId \
                                     + reference.observationId.max() \
                                     - replacements.observationId.min() \
                                     + 1
        schedule = owsched.overwrite_schedule(reference, replacements)
        schedule.sort_values('observationStartTime', inplace=True)
        
        # Check that we added the right number of visits
        added_visits = schedule.query(f'proposalId == {test_proposal_id}')
        self.assertEqual(len(added_visits), len(replacements))

        # Check that there are no more reference visits than we started with
        preserved_rows = schedule.query(f'proposalId != {test_proposal_id}')
        self.assertLessEqual(len(preserved_rows), len(reference))

        # If a visit is associated with multiple programs, it may appear once for
        # each program.
        unique_schedule = schedule.drop_duplicates(subset=['observationId'])

        previous_end = unique_schedule.eval('observationStartTime + visitTime')[:-1]
        next_start = unique_schedule.eval('observationStartTime - slewTime')[1:]
        self.assertLess(np.min(next_start.values - previous_end.values), -1.0)

        
if __name__ == '__main__':
    unittest.main()
