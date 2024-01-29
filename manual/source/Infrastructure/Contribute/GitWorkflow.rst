.. sidebar:: ToC

    .. contents::

.. _infra-contribute-git:

Git Workflow
============

Getting Started
---------------

Install the command line client, download a GUI and have a look at the basic Atlassian tutorial.

GUI
^^^

* `SourceTree <https://www.sourcetreeapp.com>`_ - Windows or MacOS X.
* `Gitg <https://wiki.gnome.org/Gitg>`_ - Linux/GNOME.

Documentation
^^^^^^^^^^^^^

* `Atlassian`__ git tutorials - easy and recommended.
* `Official <https://git-scm.com/doc>`_ git manual - complete but more advanced.

.. __: https://www.atlassian.com/git/tutorial/git-basics

Clone the SeqAn repository
^^^^^^^^^^^^^^^^^^^^^^^^^^

SeqAn is hosted on `GitHub <https://github.com/seqan/>`_.
Execute the following command to get the last sources:

.. code-block:: console

    ~ # git clone https://github.com/seqan/seqan.git seqan


SeqAn Workflow
--------------

The SeqAn workflow is based on the `Gitflow <https://www.atlassian.com/git/tutorials/comparing-workflows>`_ workflow by `Atlassian`__.

.. __: https://www.atlassian.com

Develop a feature or fix a bug
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Follow the `steps <https://www.atlassian.com/git/workflows#workflow-gitflow>`_ in “Mary and John begin new features” and “Mary finishes her feature”.

* Create a new `branch <https://www.atlassian.com/git/tutorial/git-branches#branch>`_ based on `main <https://github.com/seqan/seqan/tree/main>`_.
* Perform your changes and `commit <https://www.atlassian.com/git/tutorial/git-basics#commit>`_ them onto your feature branch.
* Keep your commit history concise (see below) and `write proper commit messages <infra-contribute-git-commits>`_.
* When the development is complete, push the feature branch to your repository on GithHub.
* Make sure that you have `signed the Contributor License Agreement <https://www.clahub.com/agreements/seqan/seqan>`_
* `Create a GitHub pull request <https://github.com/seqan/seqan/compare/main>`_ to `main <https://github.com/seqan/seqan/tree/main>`_.
* Delete your branch once it has been merged.

What to include in a commit
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The git history of your contribution should be concise. Please follow the following hints:

* A single commit should be a logical unit; don't split a logical change over multiple commits and don't address different issues in one commit.
* Do not include revisions to your changes in your history, i.e. if you receive comments on your PR, change your previous commits via ``git commit --amend`` or ``git rebase``, don't just push more changes onto the history.
* Always split functional changes and style changes, including whitespace changes, into separate commits.
* Follow our style for `commit messages <infra-contribute-git-commits>`_.
* If you don't follow these rules your contribution will be squashed into a single commit by the project member doing the merge.

An example of a good git log:

.. code-block:: console

  [FIX-#666] fix bug in sequence i/o module
  [INTERNAL] remove empty lines
  [FIX] repair apps that depended on broken behaviour
  [TEST] add test that triggers #666

An example of a bad git log:

.. code-block:: console

  [FIX] fix bug in sequence i/o module
  [INTERNAL] remove empty line
  [FIX] forgot to change foo/bar.h
  revert previous changes
  revert "revert previous changes"
  [FIX] correctly this time
  [INTERNAL] remove another empty line
  [FIX,TEST] fix apps and add test
