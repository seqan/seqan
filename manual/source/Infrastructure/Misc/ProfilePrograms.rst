.. sidebar:: ToC

    .. contents::

.. _infra-misc-profile-programs:

Profiling Programs
==================

Linux Perf Tools (Linux)
------------------------

*  https://perf.wiki.kernel.org/
*  Requires ``echo '-1' > /proc/sys/kernel/perf_event_paranoid`` as root.

Useful commands:

*  ``perf top`` - display ``top``-like display but on function granularity
*  ``perf record PROGRAM`` - execute PROGRAM with profiling
*  ``perf report PROGRAM`` - display report for PROGRAM

Google Perftools (Linux, Mac Os X)
----------------------------------

*  Download and install http://code.google.com/p/gperftools/ (also available through Ubuntu/Debian packages)
*  Compile your program with debug symbols (you probably want to enable optimization as well).

.. code-block:: console

    # Tell the profiler where to write its output.
    export CPUPROFILE=${OUT}
    LD_PRELOAD="/usr/lib/libprofiler.so.0" ${PROGRAM} ${COMMAND_LINE}
    google-pprof ${PROGRAM} ${OUT}

Interesting commands:

*  ``gv``/``web`` - display weighted call graph in gv or in your browser
*  ``top``/``topX`` - display top 10/X hitters
*  ``disasm NAME`` - disassemble functions matching NAME

Instruments (Mac Os X)
----------------------

.. todo:: Write me!

