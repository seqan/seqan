.. _infra-contribute:

Contributer Guide
=================

SeqAn is on GitHub: http://github.com/seqan/seqan

Use the GitHub page to fork SeqAn, create tickets and/or pull requests.
You can also follow and like the project there!

Contributing Code or Documentation
----------------------------------

If you are unfamiliar with git, you need to learn about it first.
See the `Atlassian Git Tutoial <https://www.atlassian.com/git/tutorials/>`_
for an introduction to Git.

Next learn about the specific Git Workflow that we use and how we mark commits:

  * :ref:`Git Workflow <infra-contribute-git>`
  * :ref:`Writing Commit Messages <infra-contribute-git-commits>`

If you are just changing something small, try to follow the style of whatever you are changing. If you contribute more code, please take the time to read:

  * :ref:`C++ Code Style <infra-contribute-style-cpp>`
  * :ref:`Other Code Styles <infra-contribute-style-other>`

SeqAn's documentation system, called **dox**, is similar to doxygen, but not identical. Read about it here if you want to contribute documentation (all code should be documented!):

  * :ref:`API Documentation System (dox) <infra-contribute-dox>`


.. toctree::
   :glob:
   :titlesonly:
   :hidden:

   GitWorkflow
   WriteCommitMessages
   StyleCpp
   StyleOther
   DoxApiDocs
