.. sidebar:: ToC

    .. contents::

.. _how-to-recipes-custom-file-endings:

Custom File Endings
===================

SeqAn's File-I/O uses file endings to determine the format of the file, e.g. whether to read in a sequence file
in **FASTA** format or an alignment file in **BAM** format.
Without this, it would be very cumbersome to guess the correct file format in order to use the correct parsing
algorithm.
However, in some use cases application developer might want to extend the fixed endings by their own endings.
This might be useful when writing applications for workflow systems, which might add arbitrary endings to
the respective files. Using SeqAn to read in a file with unsupported file ending would raise an exception,
which helps to spot user errors much easier.

In the following we will describe how to extend the file formats for existing parsers in SeqAn.
The central data structure when reading formatted files is the :dox:`FormattedFile` class.
The formatted file is used to select the file type to be parsed and whether it should open an input (for reading) or an
output stream (for writing).
For example you can use the :dox:`SeqFileIn` to read sequence files.
We will use the ``SeqFileIn`` to show how to extend the known file extensions with a custom one, namely a ``.fa.dat``
extension.

In the beginning we include all necessary headers:

.. includefrags:: demos/howto/custom_file_endings.cpp
   :fragment: includes

The first step to our own format is the definition of our own sequence input file which we will call ``MySeqFileIn``.
To do so, we define a new tag with a unique name, like ``MyFastaAdaptor``.
We use this tag to specialize the :dox:`SeqFileIn`, which can be done by using the third template parameter.
The following code snippet shows the defintion.

.. includefrags:: demos/howto/custom_file_endings.cpp
   :fragment: custom_file

Next, we define our custom format.
Again we define a new tag, which uniquely represents our new format.
In this example we call it ``MySeqFormat``. Given this new format tag, we extend the already existing format TagList
of the input sequence file by defining a new TagList.

.. hint::

   The :dox:`TagList` allows us to create a list of tags, which can be recursively iterated by a tag-apply function, to
   map a runtime value to it's corresponding tag and by thus employing tag-dispatching to the corresponding function.
   This might induce a certain compile time overhead, but it does not infer any runtime polymorphism.

Two more steps are required.
First, we have to overload a metafunction called ``FileFormat`` for our newly defined ``MySeqFileIn`` type, which we use
to declare a :dox:`TagSelector` type for our extended format TagList called ``MySeqInFormats``.
This meta-function will be used internally to test if the provided file extension format is contained in the format list
by the principle explained in the hint box above.
To finish the format definition we need to tell how the magic header looks like.
A magic header is used to determine the correct file format if the extension cannot be known.
This is for example the case, if the data is read from an input stream rather than a file.

The following code snippet will demonstrate the adaption we need to make to our code:

.. includefrags:: demos/howto/custom_file_endings.cpp
   :fragment: custom_format

After we added our custom file and custom format, we now have to specify the actual extension.
This is shown here:

.. includefrags:: demos/howto/custom_file_endings.cpp
   :fragment: custom_extension

Therefore, we overload the ``FileExtensions`` value meta-function with our defined ``MySeqFormat`` tag, which defines
in an array of ``char *``` with one element in it, namely our `.fa.dat` extension.

The last step before we can use our extended format tag is to tell SeqAn what to do, if a file is read with our custom
file ending.
This can be simply achieved by overloading the respective :dox:`FormattedFileIn#readRecord` function by using our
``MySeqFormat`` tag:

.. includefrags:: demos/howto/custom_file_endings.cpp
   :fragment: custom_read_record

Now we are ready to use our new file extension in a real application, which would otherwise cause an :dox:`IOError`.

.. includefrags:: demos/howto/custom_file_endings.cpp
   :fragment: main

The output of this example would be:

.. includefrags:: demos/howto/custom_file_endings.cpp.stdout
