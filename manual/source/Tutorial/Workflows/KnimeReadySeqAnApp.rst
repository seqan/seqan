.. sidebar:: ToC

    .. contents::

.. _tutorial-workflows-knime-ready-seqan-app:

Make Your SeqAn App KNIME Ready
===============================

Learning Objective
  You will learn how to use the the seqan:::dox:`ArgumentParser` and the SeqAn build system 
  so that, at the end, a new SeqAn application can be integrated in KNIME easily. 
  After completing this tutorial, you will be able write a new SeqAn application 
  that can be imported into a KNIME Eclipse plugin with a couple of commands.

Difficulty
  Basic

Duration
  1.5 h

Prerequisites
  :ref:`tutorial-getting-started-first-steps-in-seqan`, :ref:`tutorial-getting-started-parsing-command-line-arguments`

In this tutorial you will learn how to write a SeqAn app, which can be, easily converted into a KNIME node.

The first part consists of preparing a dummy app such that it can be used in a KNIME workflow and in the second part you are asked to adapt the app such that it becomes a simple quality control tool.

Using the seqan::ArgumentParser
-------------------------------

When we add options to the parser using :dox:`ArgumentParser#addOption`, 
we pass an :dox:`ArgParseOption` object together with the parser. The :dox:`ArgParseArgument::ArgumentType` of this :dox:`ArgParseOption` object is highly correlated to how the node generated from our application will look like.

The ArgumentType can be one of the following 

::

  *STRING:*  Argument is a string.
  *INTEGER:* Argument is a signed 32 bit integer.
  *INT64:* Argument is a signed 64 bit integer.
  *DOUBLE:* Argument is a floating point number stored as double.
  *INPUT_FILE:*  Argument is an input file.
  *OUTPUT_FILE:* Argument is an output file.
  *INPUT_PREFIX:*  Argument is a prefix to input file(s).
  *OUTPUT_PREFIX:* Argument is a prefix to output file(s).

Consider the following example application.

.. includefrags:: demos/tutorial/workflows/knime_node.cpp
   :fragment: all

While adding an :dox:`ArgParseOption` to your :dox:`ArgumentParser` you should consider the following points.

- KNIME needs to know the input and output ports of a node. Therefore we must specify them using ``ArgParseArgument::INPUT_FILE`` or ``ArgParseArgument::OUTPUT_FILE`` as can be seen in the example above.
- In addition, KNIME needs to know the valid file endings, which you can specify with :dox:`ArgParseArgument#setValidValues`, which is also shown in the example.

.. tip::

  Later, when building workflows, you can only connect an output-port of a node to the input-port of the next one if they have a compatible file endings.

- There are special types of input/output ports which are prefixes to a list of files. Such ports are specified using ``ArgParseArgument::INPUT_PREFIX`` or ``ArgParseArgument::OUTPUT_PREFIX``. You can only connect an output prefix port to an input prefix port and vise-versa.

Using the SeqAn build system to generate KNIME nodes 
----------------------------------------------------

If you are using the SeqAn build system you can generate a workflow plugin directory for all the SeqAn apps including your new one using the target ``prepare_workflow_plugin``.

In order for your application to turn into a KNIME node, you should register your app ``e.g. my_app``, by simply adding the line:

.. code-block:: cmake

    set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} <my_app> CACHE INTERNAL "")

to the end of the *CMakeList.txt* file of your application. All applications with this line in their  *CMakeList.txt* file will be included in the generated plugin when building the target``prepare_workflow_plugin``.

.. tip::

  *If You are not using the SeqAn build system for some reason*, but you used the seqan::ArgumentParser as recommended above, you still can generate a CTD file of your application.
  After building your application and go to the directory containing the executable of your application and run the following.
  
  .. code-block:: console

    ./seqan_app_name -write-ctd seqan_app_name.ctd

  This will give you the CTD file of your command-line tool. Then you can follow :ref:`tutorial-workflows-generating-knime-nodes-overview` section of the tutorial Generating KNIME Nodes to prepare a plugin directory of your application. 

