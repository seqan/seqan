.. sidebar:: ToC

   .. contents::


.. _infra-use-install:

TODO


Delete the following?

User Library Installation
^^^^^^^^^^^^^^^^^^^^^^^^^

The user could also want to install the library headers only.
The checkout step is the same as above, but he has to create a new build directory and execute CMake with different parameters.
The library will be installed to ``~/local/seqan``.

.. code-block:: console

    ~ # git clone https://github.com/seqan/seqan seqan-src
    ~ # mkdir -p seqan-build/library_only
    ~ # cd seqan-build/library_only
    library_only # cmake ../../seqan-src -DCMAKE_INSTALL_PREFIX=~/local/seqan \
                     -DSEQAN_BUILD_SYSTEM=SEQAN_RELEASE_LIBRARY
    library_only # make dox
    library_only # make install

The user can now find the SeqAn library in ``~/local/seqan/include``:

.. code-block:: console

    library_only # tree ~/local/seqan/
    /home/${USER}/local/seqan/
    ├── include
    │   └── seqan
    │       ├── align
    │       │   ├── align_base.h
    │       │   ├── align_cols.h
    │       │   ├── align_config.h
    │       │   ├── align_iterator_base.h
    │       │   ├── alignment_algorithm_interface.h
    │       │   ├── alignment_algorithm_tags.h
    │       │   ├── alignment_operations.h
    │       │   ├── align_metafunctions.h
    │       │   ├── align_traceback.h
    │       │   ├── gap_anchor.h
    ...
    │       ├── system.h
    │       └── version.h
    └── share
        └── doc
            └── seqan
                ├── LICENSE
                └── README
