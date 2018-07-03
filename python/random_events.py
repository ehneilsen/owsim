import argparse
from sys import stdout
from math import sqrt
from configparser import ConfigParser
import readline
import code
import logging
import numpy as np
import pandas as pd

max_iter_npoints = 1000000
min_iter_npoints = 100
info = logging.info

def num_events(rate, volume, uptime, nyears):
    """Randomly determine a number of everts.

    Args:
       rate - events per Gpc^3 per year
       volume - detectable volume, in Gpc^3 
       uptime - fraction of detector uptime
       nyears - the number of years
    """
    expected_events = rate*volume*uptime*nyears
    n_events = np.random.poisson(expected_events)
    return n_events

    
def sample_sphere(n, truncate = True):
    """Randomly generate pointings on a sphere

    Use a conceptually straightforward (but inefficient) way to
    generate random points on a sphere: generate random poinds in a
    cube, throw away points outside a sphere contained entirely in the
    cube, and project onto the sphere.

    Args: 
       n - the number of points to generate 
       truncate - if the algorithm ends up with more points than
           requested, return only the number requested, not all
           generated points.

    Returns:
       a pandas.DataFrame of the randomly generated points

    """
    done = False
    point_dfs = []
    accumulated_samples = 0
    while accumulated_samples < n:
        # (2*r)^3 / (4/3 pi r^3) = 6/pi
        iter_npoints = min(int(np.round((n-accumulated_samples)*6/np.pi)), max_iter_npoints)
        # do 3-sigma more
        iter_npoints = iter_npoints + np.int(3*np.sqrt(iter_npoints))
        iter_npoints = max(iter_npoints, min_iter_npoints)
        
        x = np.random.uniform(-1, 1, iter_npoints)
        y = np.random.uniform(-1, 1, iter_npoints)
        z = np.random.uniform(-1, 1, iter_npoints)

        r = np.sqrt(x*x+y*y+z*z)
        in_sphere = r < 1.0
        
        r = r[in_sphere]
        x = x[in_sphere]/r
        y = y[in_sphere]/r
        z = z[in_sphere]/r
        
        theta = np.arccos(z)
        phi = np.arctan2(y, x)
        ra = (np.degrees(phi) + 360 ) % 360
        decl = 90.0-np.degrees(theta)

        new_df = pd.DataFrame({'ra': ra, 'decl': decl})
        new_df = new_df[['ra', 'decl']]
        
        point_dfs.append(new_df)
        new_samples = ra.shape[0]
        accumulated_samples += new_samples
        info('completed %d samples' % accumulated_samples)
        
    points = pd.concat(point_dfs)
    if truncate:
        points = points[:n]

    points.reset_index(drop=True, inplace=True)
    points.index.rename('sample_idx', inplace=True)
    
    return points

def episode_triggers(rate, volume, uptime, start_mjd, end_mjd):
    """Generate a table of random triggers

    Args:
       rate - events per Gpc^3 per year
       volume - detectable volume, in Gpc^3 
       uptime - fraction of detector uptime
       start_mjd - the MJD at the start of the period
       end_mjd - the MJD at the end of the period

    Returns: 
       a pandas DataFrame with coordinates and the MJD date of events

    """

    nyears = (end_mjd - start_mjd)/365.2422
    n = num_events(rate, volume, uptime, nyears)
    events = sample_sphere(n, truncate=True)
    events['mjd'] = np.random.uniform(start_mjd, end_mjd, n)
    events = events[['mjd', 'ra', 'decl']]
    events.sort_values('mjd', inplace=True)
    return events

def triggers(rate, volume, uptime, start_mjd, end_mjd, episodes):
    """Generate a table of random triggers

    Args:
       rate - events per Gpc^3 per year
       volume - detectable volume, in Gpc^3 
       uptime - fraction of detector uptime
       start_mjd - the MJD at the start of the period
       end_mjd - the MJD at the end of the period
       episodes - times to simulate the survey

    Returns: 
       a pandas DataFrame with coordinates and the MJD date of events

    """
    episode_events = []
    for episode in range(episodes):
        events = episode_triggers(rate, volume, uptime, start_mjd, end_mjd)
        events['episode'] = episode
        episode_events.append(events[['episode', 'mjd', 'ra', 'decl']])

    events = pd.concat(episode_events, axis=0)

    return events



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config', help='configuration file name')
    parser.add_argument('fname', help='output file name')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='Print less information')
    parser.add_argument('-i', '--interactive', action='store_true',
                        help='Drop into an interactive shell after test.')

    args = parser.parse_args()
    config = ConfigParser()
    config.read(args.config)

    rate = config.getfloat('astrophysical', 'rate')
    volume = config.getfloat('survey', 'volume')
    uptime = config.getfloat('survey', 'uptime')
    min_mjd = config.getfloat('survey', 'min_mjd')
    max_mjd = config.getfloat('survey', 'max_mjd')
    np.random.seed = config.getint('simulation', 'npseed')
    episodes = config.getint('simulation', 'episodes')
    
    interactive = args.interactive
    fname = args.fname

    if not args.quiet:
        logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

    info("Generating triggers")
    events = triggers(rate, volume, uptime, min_mjd, max_mjd, episodes)
    info("Writing " + fname)
    events.to_csv(fname, sep="\t", index=False, header=True)
    
    if interactive:
        vars = globals().copy()
        vars.update(locals())
        shell = code.InteractiveConsole(vars)
        shell.interact()
