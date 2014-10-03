.. sidebar:: ToC

   .. contents::


.. _how-to-use-the-code-generator:

Using The Code Generator
------------------------

SeqAn comes with a code generator to create new sandboxes, library modules, apps, tests, and demos.
You must have Python (>= 2.5) installed to use this generator.

The generator can be used to create:

* library modules,
* tests for library modules,
* application (app) directories,
* demo programs,
* headers for applications, and
* library headers.

For more information, see the articles :ref:`infrastructure-repository-structure` and :ref:`infrastructure-build-system`.

The top level help screen of the ``skel.py`` program looks as follows:

.. code-block:: console

    Usage: skel.py [options] repository NAME
           skel.py [options] [module|test|app|demo|header|lheader] NAME LOCATION

    The SeqAn code generator.  The first version ("repository") is to be be called
    to create your new entries below the directory sandbox.  The second version is
    to be called to create new library modules, tests, apps, and demos inside a
    sandbox.

    Options:
      -h, --help            show this help message and exit
      -s SKEL_ROOT, --skel-root=SKEL_ROOT
                            Set path to the directory where the skeletons live in.
                            Taken from environment variable SEQAN_SKELS if
                            available.
      -a AUTHOR, --author=AUTHOR
                            Set author to use.  Should have the format USER
                            <EMAIL>.  Taken from environment variable SEQAN_AUTHOR
                            if it exists.
      -d, --dry-run         Do not change anything, just simulate.
      -c, --cmakelists-only
                            Only create CMakeLists.txt files
      -i, --infos-only      Only create INFO files

Creating A Sandbox
~~~~~~~~~~~~~~~~~~

Creating a sandbox is easy.
The following command will create a new sandbox called ``my_sandbox`` within the directory ``sandbox`` of the repository.

.. code-block:: console

    seqan # ./util/bin/skel.py repository sandbox/my_sandbox
    ...
    seqan # $ tree sandbox/my_sandbox/
    sandbox/my_sandbox/
    ├── CMakeLists.txt
    ├── apps
    │   └── CMakeLists.txt
    ├── demos
    │   └── CMakeLists.txt
    ├── include
    │   └── seqan
    └── tests
        └── CMakeLists.txt

Creating An App
~~~~~~~~~~~~~~~

Inside the sandboxes, you can create new apps using the *skel.py* tool, too.
We create a new app in our sandbox called "my\_app".

.. code-block:: console

    seqan # ./util/bin/skel.py app my_app sandbox/my_sandbox
    ...
    seqan # tree sandbox/my_sandbox/apps
    sandbox/my_sandbox/apps
    ├── CMakeLists.txt
    └── my_app
        ├── CMakeLists.txt
        ├── INFO
        └── my_app.cpp

Creating A Library Module
~~~~~~~~~~~~~~~~~~~~~~~~~

Library modules can be created likewise.

.. code-block:: console

    seqan # ./util/bin/skel.py module my_module sandbox/my_sandbox
    ...
    seqan # tree sandbox/my_sandbox/include/seqan
    sandbox/my_sandbox/include/seqan
    ├── my_module
    │   ├── INFO
    │   └── my_module_base.h
    └── my_module.h

Creating A Test
~~~~~~~~~~~~~~~

Now, we can also create a test for our module.

.. code-block:: console

    seqan # ./util/bin/skel.py test my_module sandbox/my_sandbox
    ...
    sandbox/my_sandbox/tests/my_module/
    ├── CMakeLists.txt
    ├── test_my_module.cpp
    └── test_my_module.h

Create A Demo
~~~~~~~~~~~~~

Demos can be created in a similar way:

.. code-block:: console

    seqan # ./util/bin/skel.py demo my_demo sandbox/my_sandbox
    ...
    seqan # tree ./sandbox/my_sandbox/demos
    ./sandbox/my_sandbox/demos
    ├── CMakeLists.txt
    └── my_demo.cpp

Create A Header
~~~~~~~~~~~~~~~

To create a header in your application, use the following command:

.. code-block:: console

    seqan # ./util/bin/skel.py header header_name.h sandbox/my_sandbox/apps/my_app

Create A Library Header
~~~~~~~~~~~~~~~~~~~~~~~

To create a library header (one with a ``namespace seqan { ... }``
construct), use the following:

.. code-block:: console

    seqan # ./util/bin/skel.py lheader lheader_name.h sandbox/my_sandbox/include/seqan/my_module

Setting The Author String
~~~~~~~~~~~~~~~~~~~~~~~~~

By default, SeqAn uses ``Your Name <your.email@example.net>`` as the author string in the generated files.
There are two ways to change this:

First, you can set the environment variable ``SEQAN_AUTHOR``:

.. code-block:: console

    seqan # export SEQAN_AUTHOR='Me <me@example.com>'
    seqan # ./util/bin/skel.py demo my_demo sandbox/my_sandbox

Second, you can use the ``--author`` parameter when calling *skel.py*.
This will override the environment variable ``SEQAN_AUTHOR`` if set.

.. code-block:: console

    seqan # ./util/bin/skel.py --author 'Me <me@example.com>' demo my_demo sandbox/my_sandbox

Creating App Tests
~~~~~~~~~~~~~~~~~~

See the article :ref:`how-to-write-app-tests` for this.

Using Your Own Templates
~~~~~~~~~~~~~~~~~~~~~~~~

You might want to use your own templates, e.g. to replace the license comment header at each top of the file.
To do this, you have to create a copy of **/util/skel**, e.g. to ``~/.seqan/templates``.
Then, you can edit the files, but you have to keep the file names intact.
You can specify the location of your template by specifying ``--skel-root`` or having set the environment variable **SEQAN_SKELS** when calling ``skel.py``.

In your templates, the following placeholders will be replaced.
Note the trailing **s**, leaving this out is a common source of error.

::

      %(AUTHOR)s  will be replaced by the author's name, either given on command
                  line or taken from environment variable SEQAN_AUTHOR.

      %(NAME)s    will be replaced by the name of the generated code.
      %(TITLE)s   will be replaced by the name of the generated, but centered in
                  74 characters, to be used in the file header comment.

      %(YEAR)d    will be replaced by the current year.
      %(DATE)s    will be replaced by the current date.
      %(TIME)s    will be replaced by the current time.

      %(HEADER_GUARD)s  will be replaced by the UPPER_CASE_PATH_H_ to the file.

      %(CMAKE_PROJECT_PATH)s  will be replaced by lower_case_path to the target
                              directory.
