.. _how-to-setup-nightly-builds:

Setup Nightly Builds
--------------------

Subversion Repository
~~~~~~~~~~~~~~~~~~~~~

There is a Subversion repository for the nightly build scripts at http://svn.mi.fu-berlin.de/seqan-nightly/trunk/.

Caveats
~~~~~~~

32bit: bzlib and zlib
^^^^^^^^^^^^^^^^^^^^^

On Debian, you need the packages *libz-dev* and *libbz2-dev*.
If you want to do 32 bit builds, you have to install *lib32bz2-dev* and *lib32z1-dev*.
