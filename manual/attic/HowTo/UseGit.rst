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
workflow by `Atlassian <https://www.atlassian.com>`__.

Develop a feature in a core module or app
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Follow the
`steps <https://www.atlassian.com/git/workflows#workflow-gitflow>`__ in
“Mary and John begin new features” and “Mary finishes her feature”.

-

   -  Create a new feature
      `branch <https://www.atlassian.com/git/tutorial/git-branches#branch>`__
      based on
      `main <https://github.com/seqan/seqan/tree/main>`__.
   -  Perform your changes and
      `commit <https://www.atlassian.com/git/tutorial/git-basics#commit>`__
      them onto your feature branch.
   -  When the development is complete, push the feature branch to your
      repository on GithHub.
   -  `Create a GitHub pull
      request <https://github.com/seqan/seqan/compare/main>`__ to
      `main <https://github.com/seqan/seqan/tree/main>`__.
   -  Delete your feature branch once it has been merged.

Fix an existing bug in a core module or app
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Follow the
`steps <https://www.atlassian.com/git/workflows#workflow-gitflow>`__ in
“End-user discovers a bug”.

-

   -  Create a new hotfix
      `branch <https://www.atlassian.com/git/tutorial/git-branches#branch>`__
      based on `main <https://github.com/seqan/seqan/tree/main>`__.
   -  Perform your changes and
      `commit <https://www.atlassian.com/git/tutorial/git-basics#commit>`__
      them onto your hotfix branch.
   -  When the fix is read, push your hotfix branch to repository on
      GitHub. Then:

#.

   #. `Create a GitHub pull
      request <https://github.com/seqan/seqan/compare/main>`__ to
      `main <https://github.com/seqan/seqan/tree/main>`__.

| ``    2. ``\ ```Create`` ``a`` ``GitHub`` ``pull``
``request`` <https://github.com/seqan/seqan/compare/main>`__\ `` to ``\ ```main`` <https://github.com/seqan/seqan/tree/main>`__\ ``.``
| ``    3. The pull requests should contain only the commits from your hotfix branch.``

-

   -  Delete your hotfix branch once it has been merged through the pull
      request.

Develop new modules and apps
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a new module or app
`branch <https://www.atlassian.com/git/tutorial/git-branches#branch>`__
where to develop your new module or application. The branch should be
based on main.

Rules
~~~~~

-

   -  Never push feature branches to the SeqAn repository.
   -  Sumbit code reviews through GitHub.

.. raw:: mediawiki

   {{TracNotice|{{PAGENAME}}}}
