.. sidebar:: ToC

   .. contents::


.. _how-to-use-seqan-nodes-in-knime:

Creating Workflows with KNIME
=============================

`KNIME <http://www.knime.org>`_ is a well established data analysis framework which supports the generation of workflows for data analysis.
In the following, we describe how to use SeqAn applications in KNIME.

Install SeqAn in KNIME
----------------------

The Installation of the SeqAn NGS Toolbox in KNIME is very easy.
Download the latest KNIME release from the KNIME website.
In KNIME click on ``Help > Install new Software``.

.. figure:: install-knime-1.png

In the opening dialog choose ``Add...``.

.. figure:: install-knime-2.png

In the opening dialog fill in the following Information:

``Name``
  ``KNIME Nightly Unstable``
``Location``
  ``http://update.knime.org/community-contributions/trunk/``

.. figure:: install-knime-3.png

After pressing OK, KNIME will show you all the contents of the added Update Site, containing also the SeqAn nodes.

.. figure:: install-knime-4.png

Select the SeqAn NGS Toolbox and click Next.
Follow the instructions.
After a restart of KNIME the SeqAn nodes will be available under ``Community Nodes``.

Add your own application to KNIME
---------------------------------

Using the CTD and a node generator program, all SeqAn applications that use the :dox:`ArgumentParser` can be made available to run in KNIME.
This is done automatically and nightly for all applications in the master branch on `github <https://github.com/seqan/seqan/tree/master>`_ that are listed in the CMAKE variable ``SEQAN_CTD_EXECUTABLES``.
The auto-generated KNIME nodes of these apps are then uploaded to the KNIME community node server and can easily be used by all KNIME users.

The following two steps are required to make your application KNIME-ready.

Adapt your applications to use the argument parser
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Follow the :ref:`tutorial-parsing-command-line-arguments` tutorial and adapt your application to use only the :dox:`ArgumentParser` to parse command line arguments.
Especially, take care to:

#. Declare your input and output file names as such via ``ArgParseArgument::INPUT_FILE`` and ``ArgParseArgument::OUTPUT_FILE``.
#. Detect the file format from the file extension (and not from a dedicated file format option).
   This can be done, for example, with :dox:`guessFormatFromFilename guessFormatFromFilename()` on an :dox:`AutoSeqFormat` object to detect a particular sequence format (e.g. FASTA) in a predefined set of formats.
#. For input/output files define a list of possible extensions via :dox:`ArgumentParser#setValidValues setValidValues()` (e.g. "fa fasta"). This list of possible extensions can be generated with :dox:`ArgumentParser#getFileExtensions getFileExtensions()` for a :dox:`TagSelector` of predefined file formats (e.g. AutoSeqFormat).
#. Avoid mutual exclusive options or other constraints that cannot be not represented by the ArgumentParser, simply ignore one of them (depending on a behavioral option).
   See the ArgumentParser tutorial if you need to define a numerical interval of possible values or a finite set of argument options.
#. Give default values.

Register your application to be considered by the node generator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Add the following section to the ``CMakeLists.txt`` file in your application folder (replace ``razers`` by your executable name):

.. code-block:: cmake

    # ----------------------------------------------------------------------------
    # Setup Common Tool Description for Generic Workflow Nodes
    # ----------------------------------------------------------------------------

    # Include executable razers in CTD structure.
    set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} razers CACHE INTERNAL "")

Use existing and contribute new workflows
-----------------------------------------

With the steps described above you will be able to set up your own workflows in KNIME.
If you want to contribute a workflow to the SeqAn community or use workflows from others you can do that on
https://github.com/seqan/knime_seqan_workflows

To contribute your own workflow, simply clone the workflow git repository into your own github repository and add a new folder ``WORKFLOWNAME_workflow``.
In KNIME export your workflow without the data files as a ``.zip`` file into that folder.
Provide a README, a screenshot and some examples as well.
Just have a look into existing workflow folders to get a notion.

After everything is ready, add and commit the new folder into your github repository and make a github pull request to the original workflow repository (https://github.com/seqan/knime\_seqan\_workflows) and - voila - it will be shared with the community.
