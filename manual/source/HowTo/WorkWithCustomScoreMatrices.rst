.. sidebar:: ToC

   .. contents::


.. _how-to-work-with-custom-score-matrices:

Working With Custom Score Matrices
==================================

This How To describes how to create new scoring matrices for Amino Acids and DNA alphabets and how to load score matrices from files.

Creating A New Built-In Score Matrix
------------------------------------

The following program demonstrate how to implement a new built-in score matrix.

.. includefrags:: core/demos/howto/scores/init_score.cpp
   :fragment: includes

Then, we perform the necessary definitions for the matrix.
This consists of three steps:

* defining a tag struct
* specializing the class ``ScoringMatrixData_`` with your tag

Note how we use enum values to compute the matrix size which itself is retrieved from the :dox:`FiniteOrderedAlphabetConcept#ValueSize` metafunction.

.. includefrags:: core/demos/howto/scores/init_score.cpp
   :fragment: user-defined-matrix

We define a function ``showScoringMatrix`` for displaying a matrix.

.. includefrags:: core/demos/howto/scores/init_score.cpp
   :fragment: show-scoring-matrix

Finally, the function ``main`` function demostrates some of the things you can do with scores:

* Construct empty score matrix object (2.)
* Programatically fill the matrix with a built-in matrix values (3.1)
* Programmatically fill the score matrix in a loop (3.2)
* Programatically fill the matrix with the user-defined matrix values (3.3)
* Directly create a score matrix with the user-defined matrix values (4)

.. includefrags:: core/demos/howto/scores/init_score.cpp
   :fragment: main

Here is the output of the program:

.. code-block:: console

    $ make tutorial_init_score
    $ ./demos/tutorial_init_score
    BLOSUM 30
        A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
    A   4   -1  0   0   -3  1   0   0   -2  0   -1  0   1   -2  -1  1   1   -5  -4  1   0   0   0   -7
    R   -1  8   -2  -1  -2  3   -1  -2  -1  -3  -2  1   0   -1  -1  -1  -3  0   0   -1  -2  0   -1  -7
    N   0   -2  8   1   -1  -1  -1  0   -1  0   -2  0   0   -1  -3  0   1   -7  -4  -2  4   -1  0   -7
    D   0   -1  1   9   -3  -1  1   -1  -2  -4  -1  0   -3  -5  -1  0   -1  -4  -1  -2  5   0   -1  -7
    C   -3  -2  -1  -3  17  -2  1   -4  -5  -2  0   -3  -2  -3  -3  -2  -2  -2  -6  -2  -2  0   -2  -7
    Q   1   3   -1  -1  -2  8   2   -2  0   -2  -2  0   -1  -3  0   -1  0   -1  -1  -3  -1  4   0   -7
    E   0   -1  -1  1   1   2   6   -2  0   -3  -1  2   -1  -4  1   0   -2  -1  -2  -3  0   5   -1  -7
    G   0   -2  0   -1  -4  -2  -2  8   -3  -1  -2  -1  -2  -3  -1  0   -2  1   -3  -3  0   -2  -1  -7
    H   -2  -1  -1  -2  -5  0   0   -3  14  -2  -1  -2  2   -3  1   -1  -2  -5  0   -3  -2  0   -1  -7
    I   0   -3  0   -4  -2  -2  -3  -1  -2  6   2   -2  1   0   -3  -1  0   -3  -1  4   -2  -3  0   -7
    L   -1  -2  -2  -1  0   -2  -1  -2  -1  2   4   -2  2   2   -3  -2  0   -2  3   1   -1  -1  0   -7
    K   0   1   0   0   -3  0   2   -1  -2  -2  -2  4   2   -1  1   0   -1  -2  -1  -2  0   1   0   -7
    M   1   0   0   -3  -2  -1  -1  -2  2   1   2   2   6   -2  -4  -2  0   -3  -1  0   -2  -1  0   -7
    F   -2  -1  -1  -5  -3  -3  -4  -3  -3  0   2   -1  -2  10  -4  -1  -2  1   3   1   -3  -4  -1  -7
    P   -1  -1  -3  -1  -3  0   1   -1  1   -3  -3  1   -4  -4  11  -1  0   -3  -2  -4  -2  0   -1  -7
    S   1   -1  0   0   -2  -1  0   0   -1  -1  -2  0   -2  -1  -1  4   2   -3  -2  -1  0   -1  0   -7
    T   1   -3  1   -1  -2  0   -2  -2  -2  0   0   -1  0   -2  0   2   5   -5  -1  1   0   -1  0   -7
    W   -5  0   -7  -4  -2  -1  -1  1   -5  -3  -2  -2  -3  1   -3  -3  -5  20  5   -3  -5  -1  -2  -7
    Y   -4  0   -4  -1  -6  -1  -2  -3  0   -1  3   -1  -1  3   -2  -2  -1  5   9   1   -3  -2  -1  -7
    V   1   -1  -2  -2  -2  -3  -3  -3  -3  4   1   -2  0   1   -4  -1  1   -3  1   5   -2  -3  0   -7
    B   0   -2  4   5   -2  -1  0   0   -2  -2  -1  0   -2  -3  -2  0   0   -5  -3  -2  5   0   -1  -7
    Z   0   0   -1  0   0   4   5   -2  0   -3  -1  1   -1  -4  0   -1  -1  -1  -2  -3  0   4   0   -7
    X   0   -1  0   -1  -2  0   -1  -1  -1  0   0   0   0   -1  -1  0   0   -2  -1  0   -1  0   -1  -7
    *   -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  1

    Coordinate Products
        A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
    A   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    R   0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21  22  23
    N   0   2   4   6   8   10  12  14  16  18  20  22  24  26  28  30  32  34  36  38  40  42  44  46
    D   0   3   6   9   12  15  18  21  24  27  30  33  36  39  42  45  48  51  54  57  60  63  66  69
    C   0   4   8   12  16  20  24  28  32  36  40  44  48  52  56  60  64  68  72  76  80  84  88  92
    Q   0   5   10  15  20  25  30  35  40  45  50  55  60  65  70  75  80  85  90  95  100 105 110 115
    E   0   6   12  18  24  30  36  42  48  54  60  66  72  78  84  90  96  102 108 114 120 126 132 138
    G   0   7   14  21  28  35  42  49  56  63  70  77  84  91  98  105 112 119 126 133 140 147 154 161
    H   0   8   16  24  32  40  48  56  64  72  80  88  96  104 112 120 128 136 144 152 160 168 176 184
    I   0   9   18  27  36  45  54  63  72  81  90  99  108 117 126 135 144 153 162 171 180 189 198 207
    L   0   10  20  30  40  50  60  70  80  90  100 110 120 130 140 150 160 170 180 190 200 210 220 230
    K   0   11  22  33  44  55  66  77  88  99  110 121 132 143 154 165 176 187 198 209 220 231 242 253
    M   0   12  24  36  48  60  72  84  96  108 120 132 144 156 168 180 192 204 216 228 240 252 264 276
    F   0   13  26  39  52  65  78  91  104 117 130 143 156 169 182 195 208 221 234 247 260 273 286 299
    P   0   14  28  42  56  70  84  98  112 126 140 154 168 182 196 210 224 238 252 266 280 294 308 322
    S   0   15  30  45  60  75  90  105 120 135 150 165 180 195 210 225 240 255 270 285 300 315 330 345
    T   0   16  32  48  64  80  96  112 128 144 160 176 192 208 224 240 256 272 288 304 320 336 352 368
    W   0   17  34  51  68  85  102 119 136 153 170 187 204 221 238 255 272 289 306 323 340 357 374 391
    Y   0   18  36  54  72  90  108 126 144 162 180 198 216 234 252 270 288 306 324 342 360 378 396 414
    V   0   19  38  57  76  95  114 133 152 171 190 209 228 247 266 285 304 323 342 361 380 399 418 437
    B   0   20  40  60  80  100 120 140 160 180 200 220 240 260 280 300 320 340 360 380 400 420 440 460
    Z   0   21  42  63  84  105 126 147 168 189 210 231 252 273 294 315 336 357 378 399 420 441 462 483
    X   0   22  44  66  88  110 132 154 176 198 220 242 264 286 308 330 352 374 396 418 440 462 484 506
    *   0   23  46  69  92  115 138 161 184 207 230 253 276 299 322 345 368 391 414 437 460 483 506 529

    User defined matrix (also BLOSUM 30)...
        A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
    A   4   -1  0   0   -3  1   0   0   -2  0   -1  0   1   -2  -1  1   1   -5  -4  1   0   0   0   -7
    R   -1  8   -2  -1  -2  3   -1  -2  -1  -3  -2  1   0   -1  -1  -1  -3  0   0   -1  -2  0   -1  -7
    N   0   -2  8   1   -1  -1  -1  0   -1  0   -2  0   0   -1  -3  0   1   -7  -4  -2  4   -1  0   -7
    D   0   -1  1   9   -3  -1  1   -1  -2  -4  -1  0   -3  -5  -1  0   -1  -4  -1  -2  5   0   -1  -7
    C   -3  -2  -1  -3  17  -2  1   -4  -5  -2  0   -3  -2  -3  -3  -2  -2  -2  -6  -2  -2  0   -2  -7
    Q   1   3   -1  -1  -2  8   2   -2  0   -2  -2  0   -1  -3  0   -1  0   -1  -1  -3  -1  4   0   -7
    E   0   -1  -1  1   1   2   6   -2  0   -3  -1  2   -1  -4  1   0   -2  -1  -2  -3  0   5   -1  -7
    G   0   -2  0   -1  -4  -2  -2  8   -3  -1  -2  -1  -2  -3  -1  0   -2  1   -3  -3  0   -2  -1  -7
    H   -2  -1  -1  -2  -5  0   0   -3  14  -2  -1  -2  2   -3  1   -1  -2  -5  0   -3  -2  0   -1  -7
    I   0   -3  0   -4  -2  -2  -3  -1  -2  6   2   -2  1   0   -3  -1  0   -3  -1  4   -2  -3  0   -7
    L   -1  -2  -2  -1  0   -2  -1  -2  -1  2   4   -2  2   2   -3  -2  0   -2  3   1   -1  -1  0   -7
    K   0   1   0   0   -3  0   2   -1  -2  -2  -2  4   2   -1  1   0   -1  -2  -1  -2  0   1   0   -7
    M   1   0   0   -3  -2  -1  -1  -2  2   1   2   2   6   -2  -4  -2  0   -3  -1  0   -2  -1  0   -7
    F   -2  -1  -1  -5  -3  -3  -4  -3  -3  0   2   -1  -2  10  -4  -1  -2  1   3   1   -3  -4  -1  -7
    P   -1  -1  -3  -1  -3  0   1   -1  1   -3  -3  1   -4  -4  11  -1  0   -3  -2  -4  -2  0   -1  -7
    S   1   -1  0   0   -2  -1  0   0   -1  -1  -2  0   -2  -1  -1  4   2   -3  -2  -1  0   -1  0   -7
    T   1   -3  1   -1  -2  0   -2  -2  -2  0   0   -1  0   -2  0   2   5   -5  -1  1   0   -1  0   -7
    W   -5  0   -7  -4  -2  -1  -1  1   -5  -3  -2  -2  -3  1   -3  -3  -5  20  5   -3  -5  -1  -2  -7
    Y   -4  0   -4  -1  -6  -1  -2  -3  0   -1  3   -1  -1  3   -2  -2  -1  5   9   1   -3  -2  -1  -7
    V   1   -1  -2  -2  -2  -3  -3  -3  -3  4   1   -2  0   1   -4  -1  1   -3  1   5   -2  -3  0   -7
    B   0   -2  4   5   -2  -1  0   0   -2  -2  -1  0   -2  -3  -2  0   0   -5  -3  -2  5   0   -1  -7
    Z   0   0   -1  0   0   4   5   -2  0   -3  -1  1   -1  -4  0   -1  -1  -1  -2  -3  0   4   0   -7
    X   0   -1  0   -1  -2  0   -1  -1  -1  0   0   0   0   -1  -1  0   0   -2  -1  0   -1  0   -1  -7
    *   -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  1

        A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
    A   4   -1  0   0   -3  1   0   0   -2  0   -1  0   1   -2  -1  1   1   -5  -4  1   0   0   0   -7
    R   -1  8   -2  -1  -2  3   -1  -2  -1  -3  -2  1   0   -1  -1  -1  -3  0   0   -1  -2  0   -1  -7
    N   0   -2  8   1   -1  -1  -1  0   -1  0   -2  0   0   -1  -3  0   1   -7  -4  -2  4   -1  0   -7
    D   0   -1  1   9   -3  -1  1   -1  -2  -4  -1  0   -3  -5  -1  0   -1  -4  -1  -2  5   0   -1  -7
    C   -3  -2  -1  -3  17  -2  1   -4  -5  -2  0   -3  -2  -3  -3  -2  -2  -2  -6  -2  -2  0   -2  -7
    Q   1   3   -1  -1  -2  8   2   -2  0   -2  -2  0   -1  -3  0   -1  0   -1  -1  -3  -1  4   0   -7
    E   0   -1  -1  1   1   2   6   -2  0   -3  -1  2   -1  -4  1   0   -2  -1  -2  -3  0   5   -1  -7
    G   0   -2  0   -1  -4  -2  -2  8   -3  -1  -2  -1  -2  -3  -1  0   -2  1   -3  -3  0   -2  -1  -7
    H   -2  -1  -1  -2  -5  0   0   -3  14  -2  -1  -2  2   -3  1   -1  -2  -5  0   -3  -2  0   -1  -7
    I   0   -3  0   -4  -2  -2  -3  -1  -2  6   2   -2  1   0   -3  -1  0   -3  -1  4   -2  -3  0   -7
    L   -1  -2  -2  -1  0   -2  -1  -2  -1  2   4   -2  2   2   -3  -2  0   -2  3   1   -1  -1  0   -7
    K   0   1   0   0   -3  0   2   -1  -2  -2  -2  4   2   -1  1   0   -1  -2  -1  -2  0   1   0   -7
    M   1   0   0   -3  -2  -1  -1  -2  2   1   2   2   6   -2  -4  -2  0   -3  -1  0   -2  -1  0   -7
    F   -2  -1  -1  -5  -3  -3  -4  -3  -3  0   2   -1  -2  10  -4  -1  -2  1   3   1   -3  -4  -1  -7
    P   -1  -1  -3  -1  -3  0   1   -1  1   -3  -3  1   -4  -4  11  -1  0   -3  -2  -4  -2  0   -1  -7
    S   1   -1  0   0   -2  -1  0   0   -1  -1  -2  0   -2  -1  -1  4   2   -3  -2  -1  0   -1  0   -7
    T   1   -3  1   -1  -2  0   -2  -2  -2  0   0   -1  0   -2  0   2   5   -5  -1  1   0   -1  0   -7
    W   -5  0   -7  -4  -2  -1  -1  1   -5  -3  -2  -2  -3  1   -3  -3  -5  20  5   -3  -5  -1  -2  -7
    Y   -4  0   -4  -1  -6  -1  -2  -3  0   -1  3   -1  -1  3   -2  -2  -1  5   9   1   -3  -2  -1  -7
    V   1   -1  -2  -2  -2  -3  -3  -3  -3  4   1   -2  0   1   -4  -1  1   -3  1   5   -2  -3  0   -7
    B   0   -2  4   5   -2  -1  0   0   -2  -2  -1  0   -2  -3  -2  0   0   -5  -3  -2  5   0   -1  -7
    Z   0   0   -1  0   0   4   5   -2  0   -3  -1  1   -1  -4  0   -1  -1  -1  -2  -3  0   4   0   -7
    X   0   -1  0   -1  -2  0   -1  -1  -1  0   0   0   0   -1  -1  0   0   -2  -1  0   -1  0   -1  -7
    *   -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  -7  1

    User DNA scoring scheme...
        A   C   G   T   N
    A   1   0   0   0   0
    C   0   1   0   0   0
    G   0   0   1   0   0
    T   0   0   0   1   0
    N   0   0   0   0   0

Loading Score Matrices From File
------------------------------------

This small demo program shows how to load a score matrix from a file.
Examples for score file are ``core/demos/howto/scores/dna_example.txt`` for DNA alphabets and ``core/tests/score/PAM250`` for amino acids.

Include the necessary headers.

.. includefrags:: core/demos/howto/scores/load_score.cpp
   :fragment: includes

We define a function that can show a scoring matrix.

.. includefrags:: core/demos/howto/scores/load_score.cpp
   :fragment: show-scoring-matrix

Finally, the main program loads the scoring matrix from the file given on the command line and then shows it.

.. includefrags:: core/demos/howto/scores/load_score.cpp
   :fragment: main

Here's the program output.

.. code-block:: console

   $ make tutorial_load_score
   $ ./demos/tutorial_load_score ../../demos/howto/scores/dna_example.txt
       A   C   G   T
   A   1   -1  -1  -1
   C   -1  1   -1  -1
   G   -1  -1  1   -1
   T   -1  -1  -1  1

