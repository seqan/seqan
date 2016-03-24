.. sidebar:: ToC

    .. contents::

.. _tutorial-algorithms-graph-algorithms:

Graph Algorithms
================

Learning Objective
  This tutorial shows how to use some graph algorithms in SeqAn. In particular we will use the dijkstra algorithm to find shortest path and viterbi Algorithm to compute Viterbi path of a sequence.

Difficulty
  Average

Duration
  1 h

Prerequisites
  :ref:`tutorial-datastructures-graphs`

Overview
--------

The following graph algorithms are currently available in SeqAn:

Elementary Graph Algorithms
  * Breadth-First Search (:dox:`breadthFirstSearch`)
  * Depth-First Search (:dox:`depthFirstSearch`)
  * Topological Sort (:dox:`topologicalSort`)
  * Strongly Connected Components (:dox:`stronglyConnectedComponents`)

Minimum Spanning Tree
  * Prim's Algorithm  (:dox:`primsAlgorithm`)
  * Kruskal's Algorithm (:dox:`kruskalsAlgorithm`)

Single-Source Shortest Path
  * DAG Shortest Path (:dox:`dagShortestPath`)
  * Bellman-Ford (:dox:`bellmanFordAlgorithm`)
  * Dijkstra (:dox:`dijkstra`)

All-Pairs Shortest Path
 * All-Pairs Shortest Path (:dox:`allPairsShortestPath`)
 * Floyd Warshall (:dox:`floydWarshallAlgorithm`)

Maximum Flow
 * Ford-Fulkerson (:dox:`fordFulkersonAlgorithm`)

Transitive Closure
 * Transitive Closure (:dox:`transitiveClosure`)

Bioinformatics Algorithms
 * Needleman-Wunsch (:dox:`globalAlignment`)
 * Gotoh (:dox:`globalAlignment`)
 * Hirschberg with Gotoh (:dox:`globalAlignment`)
 * Smith-Waterman (:dox:`localAlignment`)
 * Multiple Sequence Alignment (:dox:`globalMsaAlignment`)
 * UPGMA (:dox:`upgmaTree`)
 * Neighbor Joining (:dox:`njTree`)

The biological algorithms use heavily the alignment graph.
Most of them are covered in the tutorial :ref:`tutorial-datastructures-alignment`.
All others use the appropriate standard graph.
All algorithms require some kind of additional input, e.g., the Dijkstra algorithm requires a distance property map, alignment algorithms sequences and a score type and the network flow algorithm capacities on the edges.

Generally, only a single function call is sufficient to carry out all the calculations of a graph algorithm.
In most cases you will have to define containers that store the algorithms results prior to the function call.

In our example, we apply the shortest-path algorithm of Dijkstra. It is implemented in the function :dox:`dijkstra`.

Let's have a look at the input parameters.
The first parameter is of course the graph, ``g``.
Second, you will have to specify a vertex descriptor.
The function will compute the distance from this vertex to all vertices in the graph.
The last input parameter is an edge map containing the distances between the vertices.
One may think that the distance map is already contained in the graph.
Indeed this is the case for our graph type but it is not in general.
The cargo of a graph might as well be a string of characters or any other type.
So, we first have to find out how to access our internal edge map.
We do not need to copy the information to a new map.
Instead we can define an object of the type :dox:`InternalPropertyMap` of our type ``TCargo``.
It will automatically find the edge labels in the graph when the function :dox:`PropertyMapConcept#property` or :dox:`PropertyMapConcept#getProperty` is called on it with the corresponding edge descriptor.

The output containers of the shortest-path algorithm are two property maps, ``predMap`` and ``distMap``.
The ``predMap`` is a vertex map that determines a shortest-paths-tree by mapping the predecessor to each vertex.
Even though we are not interested in this information, we have to define it and pass it to the function.
The ``distMap`` indicates the length of the shortest path to each vertex.

.. includefrags:: demos/tutorial/graph/graph_dijkstra.cpp
   :fragment: dijkstra-containers

Having defined all these property maps, we can then call the function :dox:`dijkstra`:

.. includefrags:: demos/tutorial/graph/graph_dijkstra.cpp
   :fragment: dijkstra

Finally, we have to output the result.
Therefore, we define a second vertex iterator ``itV2`` and access the distances just like the city names with the function :dox:`PropertyMapConcept#property` on the corresponding property map.

.. includefrags:: demos/tutorial/graph/graph_dijkstra.cpp
   :fragment: dijkstra-output

Assignment 1
^^^^^^^^^^^^

.. container:: assignment

   Type
     Application

   Objective
     Write a program which calculates the connected components of the graph defined in :ref:`tutorial-datastructures-graphs-assignment-2` of the Graphs tutorial and Output the connected component for each vertex.

   Solution
     .. container:: foldable

        SeqAn provides the function :dox:`stronglyConnectedComponents` to compute the connected components of a directed graph.
        The first parameter of this function is of course the graph.
        The second parameter is an output parameter.
        It is a vertex map that will map a component id to each vertex. Vertices that share the same id are in the same component.

        .. includefrags:: demos/tutorial/graph/graph_algo_scc.cpp
            :fragment: connected-components

        Now, the only thing left to do is to walk through our graph and ouput each vertex and the corresponding component using the function :dox:`PropertyMapConcept#getProperty`.
        One way of doing so is to define a :dox:`VertexIterator`.

        .. includefrags:: demos/tutorial/graph/graph_algo_scc.cpp
            :fragment: output-connected-components

        .. includefrags:: demos/tutorial/graph/graph_algo_scc.cpp
            :fragment: return

        The output for the graph defined in the `Assignment 1`_ looks as follows:

        .. includefrags:: demos/tutorial/graph/graph_algo_scc.cpp.stdout
            :fragment: output-connected-components

        The graph consists of four components.
        The first contains vertex ``a``, ``b``, and ``e``, the second contains vertex ``c`` and ``d``, the third
        contains vertex ``f`` and ``g`` and the last contains only vertex ``h``.


Assignment 2
^^^^^^^^^^^^

.. container:: assignment

   Type
     Application

   Objective
      Extend the program from the :ref:`tutorial-datastructures-graphs-assignment-3` of the Graphs tutorial.
      Given the sequence ``s = "CTTCATGTGAAAGCAGACGTAAGTCA"``.

      #. calculate the Viterbi path of ``s`` and output the path as well as the probability of the path and
      #. calculate the probability that the HMM generated ``s`` with the forward and backward algorithm.

   Solution
     .. container:: foldable

        In :ref:`tutorial-datastructures-graphs-assignment-3` of the Graphs tutorial we defined an HMM with three states: exon, splice, and intron.

        The Viterbi path is the sequence of states that is most likely to produce a given output.
        In SeqAn, it can be calculated with the function :dox:`HmmAlgorithms#viterbiAlgorithm`.
        The produced output for this assignment is the DNA sequence ``s``.

        The first parameter of the function :dox:`HmmAlgorithms#viterbiAlgorithm` is of course the HMM, and the second parameter is the sequence ``s``.
        The third parameter is an output parameter that will be filled by the function.
        Since we want to compute a sequence of states, this third parameter is a :dox:`String` of :dox:`VertexDescriptor VertexDescriptors` which assigns a state to each character of the sequence ``s``.

        The return value of the function :dox:`HmmAlgorithms#viterbiAlgorithm` is the overall probability of this sequence of states, the Viterbi path.

        The only thing left is to output the path.
        The path is usually longer than the given sequence.
        This is because the HMM may have silent states, e.g. the begin and end state.
        To check if a state is silent SeqAn provides the function :dox:`HmmGraph#isSilent`.

        .. includefrags:: demos/tutorial/graph/graph_hmm.cpp
            :fragment: viterbi

        The output of the above piece of code is:

        .. includefrags:: demos/tutorial/graph/graph_hmm.cpp.stdout
            :fragment: viterbi

        It is even simpler to use the forward algorithm in SeqAn since it needs only the HMM and the sequence as parameters and returns a single probability.
        This is the probability of the HMM to generate the given sequence. The corresponding function is named :dox:`HmmAlgorithms#forwardAlgorithm`.

        .. includefrags:: demos/tutorial/graph/graph_hmm.cpp
            :fragment: forward-algorithm

        Analogously, the function :dox:`HmmAlgorithms#backwardAlgorithm` implements the backward algorithm in SeqAn.

        .. includefrags:: demos/tutorial/graph/graph_hmm.cpp
            :fragment: backward-algorithm

        The output of these two code fragments is:

        .. includefrags:: demos/tutorial/graph/graph_hmm.cpp.stdout
            :fragment: forward-backward
