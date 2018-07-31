# Preparing an environment for owsim

``owsim`` requires ``python`` 3.6 or later, and the LSST opsim4
stack. The easiest way to get this on the DES cluster at Fermilab is
through a ``singularity`` container. Such a container is similar to a
virtual machine, in that it appears (from the users point of view) to
be a whole different system from the host, but really it shares the
same kernel.

The only host on the DES cluster to have singularity is
des70.fnal.gov; so, use that host.

Make your own copy of the container thus:

```bash
  mkdir /data/des70.a/data/$(whoami)/singularity
  cd /data/des70.a/data/$(whoami)/singularity
  tar -xf /data/des70.a/data/neilsen/singularity/tarballs/opsim4_fbs_py3-2018-07-11_mod.tgz
```

This will yield many ``implausibly old time stamp`` complaints; these
can be ignored.

Start the container:

```bash
  cd /data/des70.a/data/$(whoami)/singularity
  env -i singularity exec \
    --home $HOME:/hosthome \
    --bind /data:/data \
    opsim4_fbs_py3-2018-07-11 \
    bash --login
```

This should give you a ``Singularity> `` prompt. Set up your
environment in the container thus::

```bash
  OWSIM_DIR=/data/des70.a/data/$(whoami)/owsim ;# or wherever it is
  source /home/opsim/.opsim4_profile_fbs
  export PYTHONPATH=${OWSIM_DIR}/python:$PYTHONPATH
```

Test ``owsim`` in the container::

```bash
  cd ${OWSIM_DIR}
  make test
```

If the tests pass, you are now in an environment in which you can run
``owsim``.