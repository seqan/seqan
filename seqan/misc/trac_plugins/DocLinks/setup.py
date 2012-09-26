from setuptools import setup

setup(
    name='DocLinksMacro',
    version='0.2',
    packages=['doc_links'],
    package_data={'doc_links' : ['htdocs/css/*.css',
                                 'htdocs/doc-link.png',
                                 ]},
    entry_points = {'trac.plugins': ['doclinksmacro = doc_links.macro']},
    )
