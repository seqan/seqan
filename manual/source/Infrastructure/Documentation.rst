.. sidebar:: ToC

   .. contents::


.. _infrastructure-documentation:

Documentation Infrastructure
============================

The documentation of SeqAn has two parts: (1) the API documentation and (2) the manual that you are reading right now.

SeqAn API Documentation
-----------------------

The SeqAn API documentation is created using a customly-written system called *dox*.
You can find out more about the syntax in :ref:`style-guide-dox-api-docs`.

You can build the documentation in the `dox` folder:

.. code-block:: console

   dox # ./dox_only.sh

SeqAn Manual
------------

The SeqAn manual is created using the `Sphinx <http://sphinx-doc.org/>`_ documentation system.
You can build it from the `manual` folder:

.. code-block:: console

   manual # make html

Note that you have to first build the dox documentation since plugins for generating the ``:dox:`` links rely on the generated search index for checks.
