.. sidebar:: ToC

   .. contents::

.. _how-to-automate-tests-with-ctest:

Automating Tests With CTest
---------------------------

The dashboard lives at `the SeqAn CDash site <http://www.seqan.de/cdash/index.php?project=SeqAn>`_.

For Linux and Mac OS X
~~~~~~~~~~~~~~~~~~~~~~

Create ``~/Nightly`` where everything will take place and check out the trunk:

.. code-block:: console

    cd ~
    mkdir Nightly
    cd Nightly
    svn co http://svn.seqan.de/seqan/trunk seqan-trunk

Now, get the build scripts:

.. code-block:: console

    cp seqan-trunk/misc/ctest/run_nightly.sh .
    cp seqan-trunk/misc/ctest/Seqan_Nightly.cmake.example Seqan_Nightly.cmake
    cp seqan-trunk/util/cmake/CTestConfig.cmake seqan-trunk/

Adjust the build name and site name in ``Seqan_Nightly.cmake``.
Now, test the setup by running:

.. code-block:: console

    chmod u+x run_nightly.sh
    ./run_nightly.sh

Add ``run_nightly.sh`` to your nightly *cron* jobs:

.. code-block:: console

    crontab -e

Example crontab file:

.. code-block:: console

    #min    hour    mday    month   wday    command
    05      1       *       *       *       sh -l ${HOME}/Nightly/run_nightly.sh > /dev/null

For Windows
~~~~~~~~~~~

Create ``Nightly`` in your home directory where everything will take place and check out the trunk:

.. code-block:: console

    cd /D %HOME%
    mkdir Nightly
    cd Nightly
    svn co http://svn.seqan.de/seqan/trunk seqan-trunk

Now, get the build scripts:

.. code-block:: console

    copy seqan-trunk\misc\ctest\run_nightly.sh .
    copy seqan-trunk\misc\ctest\Seqan_Nightly.cmake.example Seqan_Nightly.cmake
    copy seqan-trunk\util\cmake\CTestConfig.cmake seqan-trunk\

Adjust the build name and site name in ``Seqan_Nightly.cmake``.
Now, test the setup by running:

.. code-block:: console

    run_nightly.bat

Add ``run_nightly.bat`` to nightly Scheduled Tasks of Windows (analogously to the `CTest Tutorial <http://www.vtk.org/Wiki/CMake_Scripting_Of_CTest#On_Windows_.2F_Cygwin_.2F_MinGW>`_):

   #.   Open "Scheduled Tasks" from Control Panel.
   #.   Select Add Scheduled Task``
   #.   Select Next to select command.``
   #.   Click Browse... and select ``\ ``run_nightly.bat``\ ``.``
   #.   Click Next and select name and repetition date. Repetition date for Nightly dashboards should be Daily.``
   #.   Click Next and select time to start the dashboard.``
   #.   Click Next and select Open advanced properties... to fine tune the scheduled task.``
   #.   Select Next and type password of the user.``
   #.   Task is created. The Advanced Properties dialog should open.``
   #.   In advanced properties, specify full command name. This is very important that you use double quotes in case you have space in your path.``
   #.   Select 'Ok, which will ask for password again.``
   #.   The new task should be created.``

Sparse Checkouts
~~~~~~~~~~~~~~~~

This is only necessary/interesting if you are a developer with read permissions to more than *core*, *extras*, and your own sandbox.
You can checkout only a subset of the directories in the repository using a Subversion feature called *sparse directories.*

Also consult the `Subversion reference on sparse checkouts <http://svnbook.red-bean.com/en/1.5/svn.advanced.sparsedirs.html>`_.

.. code-block:: console

    #sh
    svn co --depth immediates https://svn.seqan.de/seqan/trunk seqan-trunk-sparse
    cd seqan-trunk-sparse
    svn update --set-depth infinity build core docs extras misc util
    svn update --set-depth files sandbox

