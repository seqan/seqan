How To: Use Git
---------------

TOC

Getting Started
~~~~~~~~~~~~~~~

Install the command line client, download a GUI and have a look at the
basic Atlassian tutorial.

GUI
^^^

-

   -  `SourceTree <http://www.sourcetreeapp.com>`__ - Windows or MacOS
      X.
   -  `Gitg <http://wiki.gnome.org/Gitg>`__ - Linux/GNOME.

Documentation
^^^^^^^^^^^^^

-

   -  `Atlassian <https://www.atlassian.com/git/tutorial/git-basics>`__
      git tutorials - easy and recommended.
   -  `Official <http://git-scm.com/doc>`__ git manual - complete but
      more advanced.

Clone the SeqAn repository
^^^^^^^^^^^^^^^^^^^^^^^^^^

SeqAn is hosted on `GitHub <http://github.com/seqan/>`__. Execute the
following command to get the last sources:

::

    #html
    <pre class="wiki" style="background-color:black;color:lightgray">
    $ git clone https://github.com/seqan/seqan.git SeqAn

SeqAn Workflow
~~~~~~~~~~~~~~

The SeqAn workflow is based on the
`Gitflow <https://www.atlassian.com/git/workflows#workflow-gitflow>`__
workflow by `Atlassian <https://www.atlassian.com>`__. The workflow is
based on two persistent branches:
`master <https://github.com/seqan/seqan/tree/master>`__ and
`develop <https://github.com/seqan/seqan/tree/develop>`__. Development
of new library and app features usually occurs on develop. The master
branch receives only new library and app releases, in addition to
hot-fixes to previous releases. Thus, the master branch is always stable
and safe to use, and the develop branch contains the last development
but might occasionally break overnight. The most frequent development
use cases are documented below.

Develop a feature in a core module or app
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Follow the
`steps <https://www.atlassian.com/git/workflows#workflow-gitflow>`__ in
“Mary and John begin new features” and “Mary finishes her feature”.

-

   -  Create a new feature
      `branch <https://www.atlassian.com/git/tutorial/git-branches#branch>`__
      based on
      `develop <https://github.com/seqan/seqan/tree/develop>`__.
   -  Perform your changes and
      `commit <https://www.atlassian.com/git/tutorial/git-basics#commit>`__
      them onto your feature branch.
   -  When the development is complete, push the feature branch to your
      repository on GithHub.
   -  `Create a GitHub pull
      request <https://github.com/seqan/seqan/compare/develop>`__ to
      `develop <https://github.com/seqan/seqan/tree/develop>`__.
   -  Delete your feature branch once it has been merged.

Fix an existing bug in a core module or app
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Follow the
`steps <https://www.atlassian.com/git/workflows#workflow-gitflow>`__ in
“End-user discovers a bug”.

-

   -  Create a new hotfix
      `branch <https://www.atlassian.com/git/tutorial/git-branches#branch>`__
      based on `master <https://github.com/seqan/seqan/tree/master>`__.
   -  Perform your changes and
      `commit <https://www.atlassian.com/git/tutorial/git-basics#commit>`__
      them onto your hotfix branch.
   -  When the fix is read, push your hotfix branch to repository on
      GitHub. Then:

#.

   #. `Create a GitHub pull
      request <https://github.com/seqan/seqan/compare/master>`__ to
      `master <https://github.com/seqan/seqan/tree/master>`__.

| ``    2. ``\ ```Create`` ``a`` ``GitHub`` ``pull``
``request`` <https://github.com/seqan/seqan/compare/develop>`__\ `` to ``\ ```develop`` <https://github.com/seqan/seqan/tree/develop>`__\ ``.``
| ``    3. The pull requests should contain only the commits from your hotfix branch.``

-

   -  Delete your hotfix branch once it has been merged through the pull
      request.

Develop new modules and apps
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a new module or app
`branch <https://www.atlassian.com/git/tutorial/git-branches#branch>`__
where to develop your new module or application. The branch should be
based on master if your module or application doesn’t rely on any
recently developed features. If a new feature becomes necessary later
on, the branch can be
`rebased <https://www.atlassian.com/git/tutorial/rewriting-git-history#rebase>`__
onto develop. When the development is complete, the branch can be merged
back into the corresponding base branch - either master or develop.

Rules
~~~~~

-

   -  Never push feature branches to the SeqAn repository.
   -  Sumbit code reviews through GitHub.

.. raw:: mediawiki

   {{TracNotice|{{PAGENAME}}}}
