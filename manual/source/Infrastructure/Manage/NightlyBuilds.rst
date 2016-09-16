.. sidebar:: ToC

    .. contents::

.. _infra-manage-nightly:

Nightly Builds
==============

Every night the master and develop branches of SeqAn are fetched and built on a variety of platforms. The results can be seen at the `SeqAn CDash site <http://www.seqan.de/cdash/index.php?project=SeqAn>`__.

The scripts that facilitate this are hosted `here <http://svn.mi.fu-berlin.de/seqan-nightly/trunk/>`__. Please note that the ``linux`` and ``macosx`` directories are outdated, both are now handled by the ``unix`` directory.

Unix Script variables
---------------------

+---------------------+--------------------------------------------------------------------+
| Variable            | Description                                                        |
+=====================+====================================================================+
| ``BITS``            | ``32`` or ``64`` (64 by default)                                   |
+---------------------+--------------------------------------------------------------------+
| ``GIT_BRANCH``      | ``master``, ``develop`` or a valid branch name (develop by default)|
+---------------------+--------------------------------------------------------------------+
| ``COMPILERS``       | list of compiler-binaries to use                                   |
+---------------------+--------------------------------------------------------------------+
| ``COMPILER_FLAGS``  | flags to append to the compiler calls                              |
+---------------------+--------------------------------------------------------------------+
| ``WITH_MEMCHECK``   | if set to anything but 0 CTEST will perform memchecks              |
+---------------------+--------------------------------------------------------------------+
| ``WITH_COVERAGE``   | if set to anything but 0 CTEST will perform coverage checks        |
+---------------------+--------------------------------------------------------------------+
| ``MODEL``           | ``Nightly``, ``Experimental`` or ``Continuous`` (defaults to       |
|                     | Experimental);                                                     |
|                     | this only influences the section where it is printed in CDash      |
+---------------------+--------------------------------------------------------------------+
| ``TMPDIR``          | place to store temporary files of run (will be pruned after        |
|                     | run; defaults to ``/tmp``)                                         |
+---------------------+--------------------------------------------------------------------+
| ``TESTROOT``        | The place checkouts and builds take place (if unset defaults       |
|                     | to ``TMPDIR``, which means it will be pruned; otherwise it will    |
|                     | be reused on next run)                                             |
+---------------------+--------------------------------------------------------------------+
| ``THREADS``         | number of threads to use (defaults to 1)                           |
+---------------------+--------------------------------------------------------------------+

Please see the up-to-date variables `here <http://svn.mi.fu-berlin.de/seqan-nightly/trunk/unix/bin/misc.sh>`__.

Unix cron jobs
--------------

To setup the build, checkout the subversion directory mentioned above and decide on the variables you wish to set. Remember to give ``TMPDIR`` and ``TESTROOT`` enough space.

Then open your crontab with

.. code-block:: console

    crontab -e

And add the jobs that you wish to have executed. It could look like this:

.. code-block:: console

    # Shell variable for cron
    SHELL=/bin/sh
    # PATH variable for cron
    PATH=/usr/local/libexec/ccache:/usr/local/bin:/usr/local/sbin:/sbin:/usr/sbin:/bin:/usr/bin
    #m h d m w
    5  1 * * * MODEL=Nightly TMPDIR=/tmp TESTROOT=${HOME}/nightly-builds/testroot GIT_BRANCH=master  BITS=32 COMPILERS="clang++35 clang++36 clang++37 clang++38 clang++-devel" THREADS=4 nice -n 10 ${HOME}/nightly-builds/unix/bin/run.sh >/dev/null
    5  1 * * * MODEL=Nightly TMPDIR=/tmp TESTROOT=${HOME}/nightly-builds/testroot GIT_BRANCH=develop BITS=32 COMPILERS="clang++35 clang++36 clang++37 clang++38 clang++-devel" THREADS=4 nice -n 10 ${HOME}/nightly-builds/unix/bin/run.sh >/dev/null

    5  3 * * * MODEL=Nightly TMPDIR=/tmp TESTROOT=${HOME}/nightly-builds/testroot GIT_BRANCH=master  BITS=64 COMPILERS="clang++35 clang++36 clang++37 clang++38 clang++-devel g++49 g++5 g++6" THREADS=4 nice -n 10 ${HOME}/nightly-builds/unix/bin/run.sh >/dev/null
    5  3 * * * MODEL=Nightly TMPDIR=/tmp TESTROOT=${HOME}/nightly-builds/testroot GIT_BRANCH=develop BITS=64 COMPILERS="clang++35 clang++36 clang++37 clang++38 clang++-devel g++49 g++5 g++6" THREADS=4 nice -n 10 ${HOME}/nightly-builds/unix/bin/run.sh >/dev/null


The first columns mean that on the 5th minute of the 1st/3rd hour (1:05am/3:05am) of every day, of every month and every week the subsequent command is executed. The variables are described above.

``nice -n 10`` ensures that the cron jobs get a low priority so the system remains responsive. The path after that is the path to your svn checkout. The redirection in the end prevents output from spamming your mail account.

Remember to add all folders to the PATH variable that contain binaries which you have added to ``COMPILERS=``.

Windows
-------

TODO double check this.


Now, get the build scripts:

.. code-block:: console

    copy seqan-src\misc\ctest\run_nightly.sh .
    copy seqan-src\misc\ctest\Seqan_Nightly.cmake.example Seqan_Nightly.cmake
    copy seqan-src\util\cmake\CTestConfig.cmake seqan-src\

Adjust the build name and site name in ``Seqan_Nightly.cmake``.
Now, test the setup by running:

.. code-block:: console

    run_nightly.bat

Add ``run_nightly.bat`` to nightly Scheduled Tasks of Windows (analogously to the `CTest Tutorial <http://www.vtk.org/Wiki/CMake_Scripting_Of_CTest#On_Windows_.2F_Cygwin_.2F_MinGW>`_):

   #.   Open ``Scheduled Tasks`` from Control Panel.
   #.   Select ``Add Scheduled Task``.
   #.   Select ``Next`` to select command.
   #.   Click **Browse...** and select ``run_nightly.bat``.
   #.   Click **Next** and select name and repetition date. Repetition date for Nightly dashboards should be ``Daily``.
   #.   Click **Next** and select time to start the dashboard.
   #.   Click **Next** and select ``Open advanced properties...`` to fine tune the scheduled task.
   #.   Select **Next** and type password of the user.``
   #.   Task is created. The Advanced Properties dialog should open.
   #.   In advanced properties, specify full command name. It is very important that you use double quotes in case you have spaces in your path.
   #.   Select ``Ok``, which will ask for password again.
   #.   The new task should be created.
