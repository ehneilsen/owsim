"""tests for apsupp"""
import unittest

import numpy as np
import pandas as pd
import owsched

test_db_file_name = 'data/baseline2018a.db'
test_replacement_exposures_fname = 'data/conditions.txt'


class TestOwSched(unittest.TestCase):

    def test_query_summary(self):
        dfs = owsched.query_summary(test_db_file_name)
        self.assertIn('propId', dfs.Proposal.columns)
        self.assertIn('propName', dfs.Proposal.columns)
        self.assertIn('observationStartTime', dfs.SummaryAllProps.columns)
        self.assertIn('fieldRA', dfs.SummaryAllProps.columns)
        self.assertIn('fieldDec', dfs.SummaryAllProps.columns)
        self.assertIn('visitTime', dfs.SummaryAllProps.columns)
        self.assertIn('slewTime', dfs.SummaryAllProps.columns)
        self.assertIn('filter', dfs.SummaryAllProps.columns)

    def test_split_into_sequences(self):
        exposures = pd.read_table(test_replacement_exposures_fname, sep=' ')
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
        replacements = pd.read_table(test_replacement_exposures_fname, sep=' ')
        replacements['proposalId'] = test_proposal_id
        replacements.observationId = replacements.observationId \
                                     + reference.observationId.max() \
                                     - replacements.observationId.min() \
                                     + 1
        replacement_sequences = owsched.split_into_sequences(replacements)
        schedule = owsched.overwrite_schedule(
            reference, replacement_sequences)
        schedule.sort_values('observationStartTime', inplace=True)

        # Check that we added the right number of visits
        added_visits = schedule.query(f'proposalId == {test_proposal_id}')
        self.assertEqual(len(added_visits), len(replacements))

        # If a visit is associated with multiple programs, it may
        # appear once for each program.
        unique_schedule = schedule.drop_duplicates(subset=['observationId'])

        # Check to make sure there are no overlapping exposures
        previous_end = (unique_schedule
                        .eval('observationStartTime + visitTime')[:-1])
        next_start = (unique_schedule
                      .eval('observationStartTime - slewTime')[1:])
        self.assertLess(np.min(next_start.values - previous_end.values), -1.0)

    def test_expand_by_proposal(self):
        proposals = owsched.query_summary(test_db_file_name).Proposal
        replacements = pd.read_table(test_replacement_exposures_fname, sep=' ')
        prop_visits, new_proposals = owsched.expand_by_proposal(
            replacements, proposals)

        expected_num = replacements.proposals.str.split().apply(len).sum()
        self.assertEqual(len(prop_visits), expected_num)


if __name__ == '__main__':
    unittest.main()
