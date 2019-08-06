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

# constants

# exception classes

# interface functions

# classes

# internal functions & classes

def _main():
    parser = argparse.ArgumentParser()
    parser.add_argument('events', help='input event file')
    parser.add_argument('num_batches', type=int, help='number of batches')
    parser.add_argument('min_gap', type=float, help='minimum gap size (days)')
    parser.add_argument('out_fname_template', help='template for output filenames')
    args = parser.parse_args()

    events = pd.read_csv(args.events, sep='\s+')
    in_cols = events.columns
    num_batches = args.num_batches
    min_gap = args.min_gap

    sorted_events = events.sort_values('mjd')
    gaps = sorted_events.mjd.diff()
    event_groups = np.where(gaps>min_gap, 1, 0).cumsum()
    unbatched_group = np.arange(len(sorted_events))//(len(sorted_events)/num_batches)
    group_batch = pd.DataFrame({'event_group': event_groups, 'unbatched_group': unbatched_group}).groupby('event_group')['unbatched_group'].min()
    sorted_events['event_group'] = event_groups
    sorted_events.set_index('event_group', inplace=True)
    sorted_events['batch'] = group_batch + 1

    for batch_id, batch_events in sorted_events.groupby('batch'):
        fname = args.out_fname_template % batch_id
        batch_events[in_cols].to_csv(fname, sep="\t", index=False)
        
    return 0


if __name__ == '__main__':
    status = _main()
    sys.exit(status)
