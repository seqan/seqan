.. _tutorial-workflows-index :

Workflows
=========
.. toctree::
    :hidden:
    :titlesonly:

    UseSeqAnNodesInKnime
    GenerateKnimeNodes
    KnimeReadySeqAnApp
    GenerateSeqAnKnimeNodes

Introduction
------------

In bioinformatics Workflows contain an interconnected and orchestrated series of computational or data manipulation steps. In order to compose and execute such workflows one needs workflow management systems. Workflow management systems can represent how computation proceeds from one step to the next one in a form of directed graph in which nodes represent tasks to be executed and edges represent either the flow of data or dependencies between different tasks. Among the many workflow management systems suited for bioinformatics workflows we will consider two common frameworks that are used to create, manage and execute workflows. namely:

- KNIME and 
- Galaxy

KNIME Workflows
---------------
In KNIME workflows are composed of (KNIME) nodes connected to each other by edges. The nodes are a representation of an application/algorithm that takes an input, process it and produces a desired output. The edges represent a flow of data from one specific application to the next one.

Generic KNIME nodes 
^^^^^^^^^^^^^^^^^^^
KNIME nodes are usually shipped as eclipse plugins. The term Generic KNIME node refers to KNIME node (eclipse plugin) generated from any command line tool. This is done via **GenericKnimeNodes** (GWN) package which provides an infrastructure to automatically generate such nodes from the description of their command line.

.. important::

    - If you only want to use existing SeqAn apps in KNIME follow :ref:`tutorial-workflows-use-seqan-nodes-in-knime`.

    - If you want to learn how to convert any command-line-tool into a KNIME node read the tutorial :ref:`tutorial-workflows-generating-knime-nodes`

    - If you are a SeqAn application developer and you want to make your application KNIME ready follow the tutorial :ref:`tutorial-workflows-knime-ready-seqan-app`

    - If you want to learn how to generate a KNIME node out of a SeqAn application  follow the tutorial :ref:`tutorial-workflows-generating-seqan-knime-nodes`


Galaxy Workflows
----------------
