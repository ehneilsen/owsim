# Preparing an environment for owsim

``owsim`` relies on packages in the LSST software environment. The
easiest way to prepare such an environment on the DES cluster at
Fermilab is through cvmfs:

```bash
bash --noprofile --norc
source /cvmfs/sw.lsst.eu/linux-x86_64/lsst_sims/sims_w_2019_29/loadLSST.bash
setup sims_featureScheduler
```

Currently, ``owsim`` is not packaged as an eups package, so it needs
to be added to the PYTHONPATH manually:

```bash
OWSIM_DIR=/data/des70.a/data/$(whoami)/owsim ;# or wherever it is
export PYTHONPATH=${OWSIM_DIR}/python:$PYTHONPATH
```

``owsim`` uses the ``astropy.astroplan`` package, which is not part of
the LSST environment. You can use ``pip`` to install it somewhere
convenient *outside the LSST environment*:

```bash
EXTRA_PYTHON_PACKAGES=/data/des70.a/data/$(whoami)/lsst_supplement
mkdir -p ${EXTRA_PYTHON_PACKAGES}
pip install -t ${EXTRA_PYTHON_PACKAGES} astroplan

```

This can be added to the path as necessary:
```bash
EXTRA_PYTHON_PACKAGES=/data/des70.a/data/$(whoami)/lsst_supplement
export PYTHONPATH=$PYTHONPATH:${EXTRA_PYTHON_PACKAGES}
```

To test ``owsim``:

```bash
  cd ${OWSIM_DIR}
  make test
```

If the tests pass, you are now in an environment in which you can run
``owsim``.

# Using `skybright` to calculate sky brightness

An alternate method of calculating sky brightness (developed for DES)
is also supported. To use it, ``skybright`` must be added to the
``PYTHONPATH``, for example:

```bash
export PYTHONPATH=/data/des70.a/data/neilsen/skybright:$PYTHONPATH
```

`skybright` can be obtained here: [https://github.com/ehneilsen/skybright](https://github.com/ehneilsen/skybright)
