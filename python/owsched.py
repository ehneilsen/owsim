"""Create tables of followup visits"""
import os
import sys
import argparse
import logging
from collections import namedtuple
from logging import debug, info, warning, error, critical
import contextlib
import sqlite3
import datetime

import pandas as pd

# constants

# exception classes

# interface functions


def overwrite_schedule(reference_schedule, replacement_sequences, gap=120):
    """Create a new schedule, replacing visits in a reference schedule.

    Args:
       - reference_schedule :: a pandas.DataFrame with the original visits
       - replacements :: pandas.GroupBy of replacement visit sequences
       - gap :: gap duration (in seconds) between reference
                and overwrite visits
    """
    prior_clock = 0
    subsequences = []
    for _, replacement_sequence in replacement_sequences:
        seq_start_clock = replacement_sequence.observationStartTime.min()
        info("Overwriting with replacement sequence starting at "
             + datetime.datetime.fromtimestamp(seq_start_clock).isoformat())
        # Find the start time of the next replacement sequence
        replacement_start_clock = (replacement_sequence
                                   .observationStartTime
                                   .min())

        # Find visits in the reference schedule between the end of the previous
        # replacement sequence and the start of the next one, and add them
        # to the list of sequences if there are any.
        ref_subset = reference_schedule.query(
            f'(observationStartTime > {prior_clock+gap}) and ((observationStartTime+visitTime+{gap})<{replacement_start_clock})')
        if len(ref_subset) > 0:
            subsequences.append(ref_subset)

        # Actually add the replacement sequence
        subsequences.append(replacement_sequence)

        # Record the end for use in determining the next window.
        prior_clock = (subsequences[-1]
                       .eval('observationStartTime+visitTime')
                       .max())

    ref_subset = reference_schedule.query(
        f'(observationStartTime > {prior_clock+gap})')
    if len(ref_subset) > 0:
        subsequences.append(ref_subset)

    info("Combining %d subsequences" % len(subsequences))

    if len(subsequences) > 10:
        group_size = int(len(subsequences)/48)
    else:
        group_size = 10

    if len(subsequences) > group_size:
        all_visits = subsequences[0]
        subsequences = subsequences[1:]
        while len(subsequences) > group_size:
            new_subseqs = subsequences[:group_size]
            subsequences = subsequences[group_size:]
            all_visits = pd.concat([all_visits] + new_subseqs)
        info("all_visits now %d long", len(all_visits))
        if len(subsequences) > 0:
            all_visits = pd.concat([all_visits] + subsequences)
    else:
        all_visits = pd.concat(subsequences, sort=True)
        
    return all_visits


def query_summary(fname):
    """Query and opsim database for visits

    Args:
       - fname :: the name of the sqlite3 database file

    Returns:
       a named tuple with the summary and proposal tables
    """
    table_contents = []
    with contextlib.closing(sqlite3.connect(fname)) as conn:
        for table_name in OpsimDatabaseDataFrames._fields:
            info("Querying " + table_name)
            table_contents.append(
                pd.read_sql_query('SELECT * FROM ' + table_name, conn))

    dbframes = OpsimDatabaseDataFrames(*table_contents)
    return dbframes


def split_into_sequences(visits, min_gap=120):
    """Split a set of visits into sequences at gaps

    Args:
       - visits :: the pandas.DataFrame of visits to split
       - min_gap :: the minimum gap between sequences, in seconds

    Returns:
       a pandas.GroupBy object grouped by sequence
    """
    start_gaps = visits.observationStartTime.diff()
    durations = visits.visitTime + visits.slewTime
    idle_time = start_gaps - durations
    first_in_sequence = idle_time > min_gap
    sequence_ids = first_in_sequence.cumsum()
    sequences = visits.groupby(sequence_ids)
    return sequences


def expand_by_proposal(visits, proposals):
    """Split visits into separate rows by proposal

    Args:
       - visits :: a pandas.DataFrame of visits to split by proposal
       - proposals :: a pandas.DataFrame of proposals

    Returns:
       a pandas.DataFrame of split visits
    """
    new_proposals = proposals.copy()
    visits_props = []
    for _, visit in visits.iterrows():
        for proposal_name in visit.proposals.split():
            visit = visit.copy()
            try:
                matching_proposals = (new_proposals
                                      .set_index('propName')
                                      .loc[proposal_name]
                                      .propId)
            except KeyError:
                next_propid = proposals.propId.max() + 1
                new_proposals = new_proposals.append(
                    {'propId': next_propid, 'propName': proposal_name},
                    ignore_index=True)
                matching_proposals = (new_proposals
                                      .set_index('propName')
                                      .loc[proposal_name]
                                      .propId)
            visit['proposalId'] = matching_proposals
            visits_props.append(visit)

    split_visits = pd.DataFrame.from_records(visits_props)
    split_visits.drop('proposals', axis=1, inplace=True)
    return split_visits, new_proposals


# classes

OpsimDatabaseDataFrames = namedtuple('OpsimDatabaseDataFrames',
                                     ('SummaryAllProps', 'Proposal'))


# internal functions & classes

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('reference_db',
                        help='database file for the reference simulation')
    parser.add_argument('replacements',
                        help='table of visits to add to reference database')
    parser.add_argument('overwrote_db',
                        help='database file with visits overwritten')
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
    overwrote_db = args.overwrote_db
    if os.path.isfile(overwrote_db):
        error(f'File {overwrote_db} already exists; not overwriting')
        return 1

    info("Loading reference db")
    reference_sim = query_summary(args.reference_db)

    info("Loading replacement visits")
    replacements = pd.read_table(args.replacements, sep=' ')

    replacements.sort_values('observationStartTime', inplace=True)
    max_ref_id = reference_sim.SummaryAllProps['observationId'].max()
    replacements['observationId'] = replacements['observationId'] \
                                    - replacements['observationId'].min() \
                                    + max_ref_id + 1
    
    if interactive:
        vars = globals().copy()
        vars.update(locals())
        shell = code.InteractiveConsole(vars)
        shell.interact()
    else:
        info("Replicating visits by proposal")
        visits_by_proposal, proposals = expand_by_proposal(
            replacements, reference_sim.Proposal)
        info("Splitting replacements into sequences")
        replacement_sequences = split_into_sequences(visits_by_proposal)
        info("Creating overwritten schedule")
        schedule = overwrite_schedule(
            reference_sim.SummaryAllProps, replacement_sequences)

        with contextlib.closing(sqlite3.connect(overwrote_db)) as conn:
            info("Writing revised proposals")
            proposals.to_sql('Proposal', conn)

            info("Writing the revised schedule")
            schedule.to_sql('SummaryAllProps', conn)

    return 0


if __name__ == '__main__':
    status = main()
    sys.exit(status)
