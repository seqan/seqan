// ==========================================================================
//                                  ANISE
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

// TODO(holtgrew): In case of overlaps, we could create an exact link length!

#include "gpm.h"

#include <algorithm>
#include <map>
#include <vector>

#include <lemon/smart_graph.h>
#include <lemon/path.h>

#include <seqan/basic.h>

#include "asm/overlapper.h"

#include "scaffolder/gpm_options.h"
#include "scaffolder/internal_shared.h"
#include "scaffolder/mate_link.h"
#include "scaffolder/scaffolding_result.h"
#include "scaffolder/utils.h"

namespace scaffolder {

namespace {  // anonymous namespace

// Minimal bandwidth to use.
int const MIN_BAND = 20;

template<typename T>
struct identity { using type = T; };

// An interval that is fuzzily placed but has a fixed length.
template <typename T>
struct FuzzyInterval
{
    FuzzyInterval() = default;
    FuzzyInterval(T value, T sd, T length, int mult = 3) : value(value), sd(sd), length(length), mult(mult)
    {}

    T left() const { return value; }
    T minLeft() const { return value - mult * sd; }
    T maxLeft() const { return value + mult * sd; }

    T right() const { return value + length; }
    T minRight() const { return value - mult * sd + length; }
    T maxRight() const { return value + mult * sd + length; }

    T fuzziness() const { return sd * mult; }

    T value { 0 }, sd { 0 }, length { 0 };
    int mult { 3 };
};

template <typename TStream, typename T>
TStream & operator<<(TStream & out, FuzzyInterval<T> const & val)
{
    return out << "FuzzyInterval(value=" << val.value << ", sd=" << val.sd << ", length=" << val.length
               << ", mult=" << val.mult << ")";
}

template <typename T>
T overlapLen(FuzzyInterval<T> const & lhs, FuzzyInterval<T> const & rhs)
{
    typedef typename std::conditional<std::is_integral<T>::value,
                                      std::make_signed<T>,
                                      identity<T> >::type::type S;
    return std::max(S(0), std::min(S(lhs.maxRight()), S(rhs.maxRight())) - std::max(S(lhs.minLeft()), S(rhs.minLeft())));
}

template <typename T>
bool overlap(FuzzyInterval<T> const & lhs, FuzzyInterval<T> const & rhs)
{
    return (rhs.minLeft() < lhs.maxRight()) && (lhs.minLeft() < rhs.maxRight());
}

template <typename T>
bool overlap(FuzzyInterval<T> const & lhs, FuzzyInterval<T> const & rhs, int minLen)
{
    return overlap(lhs, rhs) && (overlapLen(lhs, rhs) > (T)minLen);
}

// A value that also has a standard deviation and a multiplier for that.
template <typename T>
struct FuzzyValue
{
    FuzzyValue() = default;
    FuzzyValue(T value, T sd, int mult = 3) : value(value), sd(sd), mult(mult)
    {}

    T fuzziness() const { return sd * mult; }

    T min() const { return value - mult * sd; }
    T max() const { return value + mult * sd; }

    FuzzyValue & operator+=(FuzzyValue const & rhs)
    {
        value += rhs.value;
        sd += rhs.sd;
        mult = std::max(mult, rhs.mult);
        return *this;
    }

    FuzzyValue operator+(FuzzyValue const & rhs) const
    {
        FuzzyValue tmp = *this;
        tmp += rhs;
        return tmp;
    }

    FuzzyValue & operator+=(T const & rhs)
    {
        value += rhs;
        return *this;
    }

    FuzzyValue operator+(T const & rhs) const
    {
        FuzzyValue tmp = *this;
        tmp += rhs;
        return tmp;
    }

    T value { 0 }, sd { 0 };
    int mult { 3 };
};

template <typename TStream, typename T>
TStream & operator<<(TStream & out, FuzzyValue<T> const & val)
{
    return out << "FuzzyValue(value=" << val.value << ", sd=" << val.sd << ", mult=" << val.mult << ")";
}

template <typename T, typename T2>
FuzzyInterval<T> intervalFromValue(FuzzyValue<T> const & val, T2 length)
{
    return FuzzyInterval<T>(val.value, val.sd, length, val.mult);
}

// State of the PathMerger.

struct PathMergingAlgoState
{
    // The underlying graph.
    lemon::SmartGraph graph;
    // Mapping vertex ids to graph nodes.
    std::vector<lemon::SmartGraph::Node> nodes;
    // Mapping from sortedLinks (in Pathmerger::run()) to edges in graph.
    std::vector<lemon::SmartGraph::Arc> edges;
    // Maps nodes to the id of the adjacent path.
    lemon::SmartGraph::NodeMap<int> pathForNode;
    // Whether or not an edge is selected.
    lemon::SmartGraph::EdgeMap<bool> selected;
    // Mapping of edge to the length information.
    lemon::SmartGraph::EdgeMap<ReductionGraphLabel> labels;
    // The paths, indexed.
    std::map<int, lemon::SimplePath<lemon::SmartGraph> > paths;
    // Functor that returns whether the overlap between seq0 and seq1 (seq0 left of seq1) exists with length ovlLen,
    // allowing for "band". Parameters are (seq0, seq1, ovlLen, band).
    TOverlapExistsFunc overlapExists;
    // Functor that returns the overlap between (seq0, seq1, ovlLen, band) if any.  score == -1 if there is no such overlap.
    TComputeOverlapFunc computeOverlap;
    // The next path id.
    int nextPathID;

    PathMergingAlgoState(TOverlapExistsFunc overlapExists,
                         TComputeOverlapFunc computeOverlap) :
            pathForNode(graph), selected(graph, false), labels(graph), overlapExists(overlapExists),
            computeOverlap(computeOverlap), nextPathID(0)
    {}
};


void printPath(std::ostream & out, lemon::SimplePath<lemon::SmartGraph> const & path, lemon::SmartGraph const & graph)
{
    out << "path {";
    for (int i = 0; i < path.length(); ++i)
        out << "(" << graph.id(graph.source(path.nth(i))) << ", "
            << graph.id(graph.target(path.nth(i))) << ")";
    out << "}";
}

void printPath(std::ostream & out, std::list<lemon::SmartGraph::Arc> const & path, lemon::SmartGraph const & graph)
{
    out << "path {";
    for (auto arc : path)
        out << "(" << graph.id(graph.source(arc)) << ", " << graph.id(graph.target(arc)) << ")";
    out << "}";
}

// Copy path to list of arcs.

void copyPath(std::list<lemon::SmartGraph::Arc> & dest,
              lemon::SimplePath<lemon::SmartGraph> const & src)
{
    for (int i = 0; i < src.length(); ++i)
        dest.push_back(src.nth(i));
}

void copyPath(lemon::SimplePath<lemon::SmartGraph> & dest,
              std::list<lemon::SmartGraph::Arc> const & src)
{
    for (auto const & arc : src)
        dest.addBack(arc);
}

// Helper functions for Zipper* classes.

class ZipperMixin
{
protected:
    // The state with the graph, labels etc.
    PathMergingAlgoState & state;
    // The input paths as lists, path2 will be copied to output after merging.  We use lists here since we need the
    // splicing but also reverse iteration in case of zippering to the left.
    typedef std::list<lemon::SmartGraph::Arc> TArcList;
    typedef TArcList::iterator TArcListIter;
    TArcList & path1, & path2;
    // Options.
    PathMergingOptions const & options;

public:

    ZipperMixin(PathMergingAlgoState & state,
                TArcList & path1,
                TArcList & path2,
                PathMergingOptions const & options) :
            state(state), path1(path1), path2(path2), options(options)
    {}

protected:

    std::string arcStr(lemon::SmartGraph::Arc e)
    {
        std::stringstream ss;
        ss << "(" << state.graph.id(state.graph.source(e)) << ", "
           << state.graph.id(state.graph.target(e)) << ")";
        return ss.str();
    }

    template <typename TEdge>
    ReductionGraphLabel & label(TEdge const & edge)
    {
        return state.labels[edge];
    }

    template <typename TEdge>
    ReductionGraphLabel const & label(TEdge const & edge) const
    {
        return state.labels[edge];
    }

    // Helper function that filters a vector via copy if and resizing.
    template <typename TValue, typename TFunc>
    void filterTo(std::vector<TValue> & values, TFunc func)
    {
        auto itEnd = std::copy_if(values.begin(), values.end(), values.begin(), func);
        values.resize(itEnd - values.begin());
    }
};

// Helper class for zippering right.

class ZipperRightHelper : public ZipperMixin
{

public:
    ZipperRightHelper(PathMergingAlgoState & state,
                      TArcList & path1,
                      TArcList & path2,
                      PathMergingOptions const & options) :
            ZipperMixin(state, path1, path2, options)
    {}

    // Zipper path1 into path2 with guiding arc arc starting at start1 in path1 and start2 in path2.
    bool run(lemon::SmartGraph::Arc arc, TArcListIter start1, TArcListIter start2)
    {
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
        {
            std::cerr << "    STARTING ZIPPERING RIGHT\n"
                      << "      PATH1\t";
            printPath(std::cerr, path1, state.graph);
            std::cerr << "\n        LENS";
            for (auto arc : path1)
                std::cerr << "\t" << label(arc).lengthMean;
            std::cerr << "\n      PATH2\t";
            printPath(std::cerr, path2, state.graph);
            std::cerr << "\n        LENS";
            for (auto arc : path2)
                std::cerr << "\t" << label(arc).lengthMean;
            std::cerr << "\n"
                      << "      GUIDING EDGE\t" << arcStr(arc) << "\n"
                      << "      START1\t" << arcStr(*start1) << "\n"
                      << "      START2\t" << arcStr(*start2) << "\n";
        }
        SEQAN_CHECK(!label(arc).isContig(), "Must not be contig.");

        SEQAN_CHECK(label(*start1).isContig(), "Must start with contig edge.");
        SEQAN_CHECK(label(*start2).isContig(), "Must start with contig edge.");

        // If we already are at the end of path2 then we can simply append the rest of path1.  NB: We cannot be at the
        // end of path1, otherwise there would be no guiding edge.
        if (std::next(start2) == path2.end())
        {
            // Verify indicated overlap.
            FuzzyValue<double> arcLen(label(arc).lengthMean, label(arc).lengthStdDev);
            int band = arcLen.fuzziness();
            band = std::max(band, 5);
            if (arcLen.max() < 0)  // overlap required
            {
                if (!state.overlapExists(label(*start2).idx, label(*start1).idx, -arcLen.value,
                                         std::max(band, MIN_BAND)))
                {
                    if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                        std::cerr << "      NO OVERLAP, ZIPPERING FAILED!\n";
                    return false;
                }
            }

            // In the case that there is an overlap in the standard deviation, update the arc length to the lenght of
            // the overlap.  If the overlap is significant then set standard deviation to 0.
            if (arcLen.min() < 0)
            {
                // TODO(holtgrew): Make undo-able.
                auto ovl = state.computeOverlap(label(*start2).idx, label(*start1).idx, -arcLen.value, band);
                if (ovl.errors != assembler::Overlap::INVALID)
                {
                    SEQAN_ASSERT_GT(ovl.len0, ovl.begin1);
                    label(arc).lengthMean = -((int)ovl.len0 - (int)ovl.begin1);
                    if ((int)ovl.len0 - (int)ovl.begin1 > options.significantOverlap)
                        label(arc).lengthStdDev = 0;
                }
                else  // no overlap, indicate through lengthMean
                {
                    label(arc).lengthMean = 0;
                }
            }

            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                std::cerr << "      START2 AT END OF PATH2, TAKING SHORTCUT\n";
            path2.push_back(arc);
            std::copy(start1, path1.end(), std::back_inserter(path2));

            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            {
                std::cerr << "    DONE ZIPPERING RIGHT\n"
                          << "      RESULTING PATH\t";
                printPath(std::cerr, path2, state.graph);
                std::cerr << "\n      LENS";
                for (auto arc : path2)
                    std::cerr << "\t" << label(arc).lengthMean;
                std::cerr << "\n";
            }
            return true;
        }

        // We first have to find the best position for the active contig in path2.  We can do this by scanning over
        // path2 and computing the inferred distances.
        auto candidates = computeCandidates(arc, start1, start2);

        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
        {
            std::cerr << "    CANDIDATES\n";
            for (auto it : candidates)
                if (it != path2.end())
                    std::cerr << "      " << arcStr(*it) << "\n";
                else
                    std::cerr << "      <end>\n";
        }

        // In the simplest case, we have no candidates and zippering right failed.
        if (candidates.empty())
        {
            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                std::cerr << "    NO CANDIDATES.\n"
                          << "    ZIPPERING RIGHT FAILED!\n";
            return false;
        }

        // In the next simplest case, we have to fit the active contig (and the remainder of path1) to the end of path2.
        if (candidates.front() == path2.end() || candidates.front() == std::prev(path2.end()))
        {
            // Compute path length and standard deviation.
            double pathLen = 0, pathSD = 0;
            for (auto it = std::next(start2); it != path2.end(); ++it)
            {
                pathLen += label(*it).lengthMean;
                pathSD += label(*it).lengthStdDev;
            }

            double arcLen = label(arc).lengthMean - pathLen;  // of new arc
            double arcStdDev = label(arc).lengthStdDev + pathSD;
            // We cannot allow an overlap if after last, overlap with last contig would have occured in filtered candidates.
            if (candidates.front() != std::prev(path2.end()) && arcLen < 0)
                arcLen = 0;

            // In the case that there is an overlap in the standard deviation, update the arc length to the length of
            // the overlap.  If the overlap is significant then set standard deviation to 0.
            if (arcLen < 0)
            {
                // TODO(holtgrew): Make undo-able.
                int band = options.mult * arcStdDev;
                band = std::max(band, MIN_BAND);
                auto ovl = state.computeOverlap(label(path2.back()).idx, label(*start1).idx, -arcLen, band);
                if (ovl.errors != assembler::Overlap::INVALID)
                {
                    SEQAN_ASSERT_LT(ovl.begin1, ovl.len0);
                    arcLen = ovl.len0 - ovl.begin1;
                    if ((int)ovl.len0 - (int)ovl.begin1 > options.significantOverlap)
                        arcStdDev = 10;
                }
            }

            // Create new arc.
            lemon::SmartGraph::Node v = state.graph.target(path2.back());
            lemon::SmartGraph::Node w = state.graph.target(arc);
            lemon::SmartGraph::Edge f = state.graph.addEdge(v, w);
            label(f) = ReductionGraphLabel(ReductionGraphLabel::INFERRED, arcLen, arcStdDev, seqan::maxValue<unsigned>(),
                                           label(arc).count, label(arc).weight);
            path2.push_back(state.graph.direct(f, true));
            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                std::cerr << "      new arc: " << arcStr(path2.back()) << " len = " << arcLen << ", SD = " << arcStdDev << "\n";
            // Add remainder of path1.
            std::copy(start1, path1.end(), std::back_inserter(path2));

            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                std::cerr << "    PATH1 fit behind PATH2.\n"
                          << "    ZIPPERING RIGHT SUCCEEDED!\n";
            return true;
        }

        // In all other cases, we have found a good target position.  We simply take the first candidate and fit the
        // active contig here.
        bool success = false;
        if (label(*candidates.front()).isContig())
            success = placeAfterContig(candidates.front(), arc, start1, start2);
        else
            success = placeAfterContig(std::prev(candidates.front()), arc, start1, start2);
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
        {
            std::cerr << "    RESULTING PATH AFTER ZIPPERING STEP RIGHT RECURSION (SUCCESS == " << success << ")\n      ";
            printPath(std::cerr, path2, state.graph);
            std::cerr << "\n";
        }
        return success;
    }

    // Place the contig arc *start1 after contig arc *dest and continue zippering right.  Guiding arc is arc.
    bool placeAfterContig(TArcListIter dest, lemon::SmartGraph::Arc arc, TArcListIter start1, TArcListIter start2)
    {
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            std::cerr << "placeAfterContig(" << arcStr(*dest) << ", " << arcStr(arc) << ", " << arcStr(*start1)
                      << ", " << arcStr(*start2) << "\n";

        // Compute path length and stardard deviation until (including) dest but excluding start contig length.
        double pathLen = 0, pathStdDev = 0;
        for (TArcListIter it2 = std::next(start2); it2 != std::next(dest); ++it2)
        {
            pathLen += label(*it2).lengthMean;
            pathStdDev += label(*it2).lengthStdDev;
        }

        // Subtract length of path (and add std deviation) from guiding arc to get estimate for new edge length.
        double newArcLen = label(arc).lengthMean - pathLen;
        double newArcSD = pathStdDev + label(arc).lengthStdDev;
        lemon::SmartGraph::Node v = state.graph.target(*dest);
        lemon::SmartGraph::Node w = state.graph.source(*start1);
        lemon::SmartGraph::Edge f = state.graph.addEdge(v, w);  // create new link edge
        label(f) = ReductionGraphLabel(ReductionGraphLabel::INFERRED, newArcLen, newArcSD,
                                       seqan::maxValue<unsigned>(), label(arc).count,
                                       label(arc).weight);
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            std::cerr << "      new arc: " << arcStr(state.graph.direct(f, true)) << " len = " << label(f).lengthMean
                      << ", SD = " << label(f).lengthStdDev << "\n";
        // Insert new link edge into path 2.
        dest = path2.insert(++dest, state.graph.direct(f, true));
        // Insert edge *start into path2.
        dest = path2.insert(++dest, *start1);
        // Create new link edge.
        lemon::SmartGraph::Node x = state.graph.target(*start1);
        lemon::SmartGraph::Node y = state.graph.target(*++dest);
        lemon::SmartGraph::Edge g = state.graph.addEdge(x, y);
        label(g) = label(*dest);
        label(g).lengthMean -= newArcLen + label(*start1).lengthMean;;
        *dest = state.graph.direct(g, true);
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            std::cerr << "      new arc: " << arcStr(state.graph.direct(g, true)) << " len = " << label(g).lengthMean
                      << ", SD = " << label(g).lengthStdDev << "\n";

        // We might be done already.
        if (std::next(start1) == path1.end())
        {
            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            {
                std::cerr << "    NO MORE IN PATH1\n";
                std::cerr << "      RESULTING PATH AFTER ZIPPERING RIGHT\n      ";
                printPath(std::cerr, path2, state.graph);
                std::cerr << "\n";
            }
            return true;
        }

        // Zipper right recursion.
        SEQAN_ASSERT_NOT(std::next(start1) == path1.end());
        SEQAN_ASSERT_NOT(std::next(std::next(start1)) == path1.end());
        return run(*std::next(start1), std::next(std::next(start1)), std::prev(dest));
    }

    // Compute candidates when fitting start1 into path2 (starting at start2) using guiding edge arc.
    std::vector<TArcListIter> computeCandidates(lemon::SmartGraph::Arc arc, TArcListIter start1, TArcListIter start2)
    {
        std::vector<TArcListIter> candidates;

        // We will look for placements that are within the deviation of the guiding arc.
        FuzzyValue<double> arcLen(label(arc).lengthMean, label(arc).lengthStdDev, options.mult);

        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            std::cerr << "      minLen=" << arcLen.min() << "\tmaxLen=" << arcLen.max() << "\n";

        // Collect edges on path2 that overlap with [minLen, maxLen).  We add path2.end() if it would fit right of
        // path2.  These candidates are later filtered based on the overlap list.
        FuzzyValue<double> pathLen(0, 0);
        for (TArcListIter it2 = start2; it2 != path2.end(); ++it2)
        {
            FuzzyValue<double> it2Len(label(*it2).lengthMean, label(*it2).lengthStdDev, options.mult);

            // If *it2 is a contig edge then we check whether the active contig overlaps with it.  If *it2 is a link
            // edge then the active contig has to fit completely into the link edge.
            if (label(*it2).isContig())
            {
                auto it2Interval = intervalFromValue(pathLen, it2Len.value);
                auto start1Interval = intervalFromValue(arcLen + it2Len.value, label(*start1).lengthMean);
                if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                    std::cerr << "        overlap significantly? " << options.significantOverlap << "\n"
                              << "            " << label(*it2).idx << ": " << it2Interval << "\n"
                              << "            " << label(*start1).idx << ": " << start1Interval << "\n"
                              << "            => " << overlap(it2Interval, start1Interval, options.significantOverlap) << "\n";
                if (overlap(it2Interval, start1Interval, options.significantOverlap))
                {
                    unsigned seq0 = label(*it2).idx;
                    unsigned seq1 = label(*start1).idx;
                    int ovlLen = it2Interval.right() - start1Interval.left();
                    int band = it2Interval.fuzziness() + start1Interval.fuzziness();
                    if (it2Interval.minLeft() > start1Interval.minLeft())
                        std::swap(seq0, seq1);
                    if (state.overlapExists(seq0, seq1, ovlLen, std::max(band, MIN_BAND)))
                    {
                        candidates.push_back(it2);
                        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                            std::cerr << "YES, OVERLAP EXISTS: seq0" << seq0 << ", seq1" << seq1
                                      << ", ovlLen=" << ovlLen << ", band=" << band << "\n";
                    }
                }
            }
            else
            {
                // First condition: Contig must fit into link edge length-wise (slack > 0)
                // Second condition: Must be able to start far enough to the left (within slack).
                double slack = it2Len.max() - label(*start1).lengthMean;
                double d = pathLen.fuzziness();
                if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                    std::cerr << "        d == " << d << ", slack == " << slack << "\n"
                              << "            overlaps \n"
                              << "              " << intervalFromValue(arcLen, 0) << "\n"
                              << "              " << intervalFromValue(pathLen, it2Len.value) << "\n"
                              << "              => " << overlap(intervalFromValue(arcLen, 0),
                                                                intervalFromValue(pathLen, it2Len.value)) << "\n";
                if (slack > 0 && overlap(intervalFromValue(arcLen, 0), intervalFromValue(pathLen, it2Len.value)))
                    candidates.push_back(it2);
            }

            // Traverse edge *it2 and increase length/standard deviation.
            if (it2 != start2)
                pathLen += it2Len;
        }

        // Reverse candidates to obtain rightmost feasible first.  This protects us against "pushing" contigs along the
        // path if they fully overlap with their left neighbour.
        reverse(candidates.begin(), candidates.end());

        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            std::cerr << "maxLen == " << arcLen.max() << ", pathLen == " << pathLen << "\n";
        if (arcLen.max() > pathLen.min())
            candidates.push_back(path2.end());

        return candidates;
    }
};

class ZipperLeftHelper : public ZipperMixin
{
public:
    ZipperLeftHelper(PathMergingAlgoState & state,
                     TArcList & path1,
                     TArcList & path2,
                     PathMergingOptions const & options) :
            ZipperMixin(state, path1, path2, options)
    {}


    // Wrapper for zipperLeft() that does translation of start1 and start2.
    bool run(lemon::SmartGraph::Arc arc, TArcListIter start1, TArcListIter start2)
    {
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
        {
            std::cerr << "    STARTING ZIPPERING LEFT WRAPPER\n"
                      << "      PATH1\t";
            printPath(std::cerr, path1, state.graph);
            std::cerr << "\n        LENS";
            for (auto arc : path1)
                std::cerr << "\t" << label(arc).lengthMean;
            std::cerr << "\n      PATH2\t";
            printPath(std::cerr, path2, state.graph);
            std::cerr << "\n        LENS";
            for (auto arc : path2)
                std::cerr << "\t" << label(arc).lengthMean;
            std::cerr << "\n"
                      << "      GUIDING EDGE\t" << arcStr(arc) << "\n"
                      << "      START1\t" << arcStr(*start1) << "\n"
                      << "      START2\t" << arcStr(*start2) << "\n";
        }
        SEQAN_CHECK(!label(arc).isContig(), "Must not be contig.");

        // The guiding arc is the right one for zippering right but when zippering left then we have to take the edge
        // left of start1.  We are done when start1 is the begin of path1 then path2 is already complete.

        // Check for shortcut.
        if (start1 == path1.begin())
        {
            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            {
                std::cerr << "    start1 AT BEGINNING OF path1, DONE\n"
                          << "    DONE ZIPPERING LEFT\n"
                          << "      RESULTING PATH\t";
                printPath(std::cerr, path2, state.graph);
                std::cerr << "\n";
            }
            return true;
        }

        // The real start edge in path1 is two edges left of start1.
        --start1;
        TArcListIter leftArc = start1;  // guiding edge for zippering left
        SEQAN_CHECK(start1 != path1.begin(), "Cannot be at beginning yet.");
        --start1;
        // The real start edge in path2 is two edges right of start2.
        ++start2;
        SEQAN_CHECK(start2 != path2.end(), "Cannot be at end yet.");
        ++start2;
        SEQAN_CHECK(start2 != path2.end(), "Cannot be at end yet.");

        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
        {
            std::cerr << "    REAL start1 == " << arcStr(*start1) << "\n"
                      << "    REAL start2 == " << arcStr(*start2) << "\n"
                      << "    REAL arc    == " << arcStr(*leftArc) << "\n";
        }

        bool success = zipperLeft(*leftArc, start1, start2);

        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
        {
            std::cerr << "    DONE ZIPPERING LEFT WRAPPER (SUCCESS == " << success << ")\n"
                      << "      RESULTING PATH\t";
            printPath(std::cerr, path2, state.graph);
            std::cerr << "\n";
        }

        return success;
    }

private:

    // Zipper path1 into path2 with guiding arc arc starting at start1 in path1 and start2 in path2.
    bool zipperLeft(lemon::SmartGraph::Arc arc, TArcListIter start1, TArcListIter start2)
    {
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
        {
            std::cerr << "    STARTING ZIPPERING LEFT\n"
                      << "      PATH1\t";
            printPath(std::cerr, path1, state.graph);
            std::cerr << "\n        LENS";
            for (auto arc : path1)
                std::cerr << "\t" << label(arc).lengthMean;
            std::cerr << "\n      PATH2\t";
            printPath(std::cerr, path2, state.graph);
            std::cerr << "\n        LENS";
            for (auto arc : path2)
                std::cerr << "\t" << label(arc).lengthMean;
            std::cerr << "\n"
                      << "      GUIDING EDGE\t" << arcStr(arc) << "\n"
                      << "      START1\t" << arcStr(*start1) << "\n"
                      << "      START2\t" << arcStr(*start2) << "\n";
        }
        SEQAN_CHECK(!label(arc).isContig(), "Must not be contig.");

        SEQAN_CHECK(label(*start1).isContig(), "Must start with contig edge.");
        SEQAN_CHECK(label(*start2).isContig(), "Must start with contig edge.");

        // If we end up at the beginning of path2 then we can simply concatenate the remainder of path1 to path2.
        if (start2 == path2.begin())
        {
            // Verify indicated overlap.
            double arcLen = label(arc).lengthMean;
            double arcSD = label(arc).lengthStdDev;
            if (arcLen + options.mult * arcSD < 0)  // indicating overlap
            {
                int band = options.mult * arcSD;
                band = std::max(band, MIN_BAND);
                if (!state.overlapExists(label(*start1).idx, label(*start2).idx, -arcLen, band))
                {
                    if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                        std::cerr << "      NO OVERLAP, ZIPPERING FAILED!\n";
                    return false;
                }
            }

            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                std::cerr << "      START2 AT BEGINNING OF PATH2, TAKING SHORTCUT\n";
            path2.push_front(arc);
            for (TArcListIter it1 = start1; it1 != path1.begin(); --it1)
                path2.push_front(*it1);
            path2.push_front(path1.front());

            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            {
                std::cerr << "    DONE ZIPPERING LEFT\n"
                          << "      RESULTING PATH\t";
                printPath(std::cerr, path2, state.graph);
                std::cerr << "\n";
            }
            return true;
        }

        // We could take no shortcuts and the hard work begins here.

        // Look for the position of the active contig in path2.  We can do this by scanning over path2 to the left and
        // computing inferred distances.

        // We first have to find the best position for the active contig in path2.  We can do this by scanning over
        // path2 and computing the inferred distances.
        auto candidates = computeCandidates(arc, start1, start2);

        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
        {
            std::cerr << "    CANDIDATES\n";
            for (auto it : candidates)
                if (it != path2.end())
                    std::cerr << "      " << arcStr(*it) << "\n";
                else
                    std::cerr << "      <front>\n";
        }

        // In the simplest case we have no candidates and zippering left failed.
        if (candidates.empty())
        {
            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                std::cerr << "    NO CANDIDATES.\n"
                          << "    ZIPPERING LEFT FAILED!\n";
            return false;
        }

        // In the next simplest case, we have to fit the active contig (and the remainder of path1) to the front of path2.
        if (candidates.front() == path2.end() || candidates.front() == path2.begin())
        {
            // Compute path length and standard deviation.
            double pathLen = 0, pathSD = 0;
            for (auto it = path2.begin(); it != start2; ++it)
            {
                pathLen += label(*it).lengthMean;
                pathSD += label(*it).lengthStdDev;
            }

            double arcLen = label(arc).lengthMean - pathLen;  // of new arc
            double arcStdDev = label(arc).lengthStdDev + pathSD;
            // We cannot allow an overlap if before first, overlap with last contig would have occured in filtered
            // candidates.
            if (candidates.front() != path2.begin() && arcLen < 0)
                arcLen = 10;
            // Create new arc.
            lemon::SmartGraph::Node v = state.graph.source(arc);
            lemon::SmartGraph::Node w = state.graph.source(path2.front());
            lemon::SmartGraph::Edge f = state.graph.addEdge(v, w);
            label(f) = ReductionGraphLabel(ReductionGraphLabel::INFERRED, arcLen, arcStdDev, seqan::maxValue<unsigned>(),
                                           label(arc).count, label(arc).weight);
            path2.push_front(state.graph.direct(f, true));
            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                std::cerr << "      new arc: " << arcStr(path2.front()) << " len = " << arcLen << ", SD = " << arcStdDev << "\n";
            // Add remainder of path1.
            auto itDst = path2.begin();
            for (auto it = path1.begin(); it != std::next(start1); ++it)
                itDst = std::next(path2.insert(itDst, *it));

            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                std::cerr << "    PATH1 fit before PATH2.\n"
                          << "    ZIPPERING LEFT SUCCEEDED!\n";
            return true;
        }

        // In all other cases, we have found a good target position.  We simply take the first candidate and fit the
        // active contig here.
        bool success = false;
        if (label(*candidates.front()).isContig())
            success = placeBeforeContig(candidates.front(), arc, start1, start2);
        else
            success = placeBeforeContig(std::next(candidates.front()), arc, start1, start2);
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
        {
            std::cerr << "    RESULTING PATH AFTER ZIPPERING STEP LEFT RECURSION (SUCCESS == " << success << ")\n      ";
            printPath(std::cerr, path2, state.graph);
            std::cerr << "\n";
        }
        return success;
    }

   // Place the contig arc *start1 before contig arc *dest and continue zippering left.  Guiding arc is arc.
    bool placeBeforeContig(TArcListIter dest, lemon::SmartGraph::Arc arc, TArcListIter start1, TArcListIter start2)
    {
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            std::cerr << "placeBeforeContig(" << arcStr(*dest) << ", " << arcStr(arc) << ", " << arcStr(*start1)
                      << ", " << arcStr(*start2) << "\n";

        // Compute path length and stardard deviation until (including) dest but excluding start contig length.
        double pathLen = 0, pathStdDev = 0;
        for (TArcListIter it2 = dest; it2 != start2; ++it2)
        {
            pathLen += label(*it2).lengthMean;
            pathStdDev += label(*it2).lengthStdDev;
        }

        // Subtract length of path (and add std deviation) from guiding arc to get estimate for new edge length.
        double newArcLen = label(arc).lengthMean - pathLen;
        double newArcSD = pathStdDev + label(arc).lengthStdDev;
        lemon::SmartGraph::Node v = state.graph.target(*start1);
        lemon::SmartGraph::Node w = state.graph.source(*dest);
        lemon::SmartGraph::Edge f = state.graph.addEdge(v, w);  // create new link edge
        label(f) = ReductionGraphLabel(ReductionGraphLabel::INFERRED, newArcLen, newArcSD,
                                       seqan::maxValue<unsigned>(), label(arc).count,
                                       label(arc).weight);
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            std::cerr << "      new arc: " << arcStr(state.graph.direct(f, true)) << " len = " << label(f).lengthMean
                      << ", SD = " << label(f).lengthStdDev << "\n";
        // Insert new link edge into path 2.
        dest = path2.insert(dest, state.graph.direct(f, true));
        // Insert edge *start1 into path2.
        dest = path2.insert(dest, *start1);
        // Create new link edge.
        lemon::SmartGraph::Node x = state.graph.source(*std::prev(dest));
        lemon::SmartGraph::Node y = state.graph.source(*start1);
        lemon::SmartGraph::Edge g = state.graph.addEdge(x, y);
        label(g) = label(*--dest);
        label(g).lengthMean -= newArcLen + label(*start1).lengthMean;;
        *dest = state.graph.direct(g, true);
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            std::cerr << "      new arc: " << arcStr(state.graph.direct(g, true)) << " len = " << label(g).lengthMean
                      << ", SD = " << label(g).lengthStdDev << "\n";

        
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
        {
            std::cerr << "    INTERMEDIATE RESULT\n      ";
            printPath(std::cerr, path2, state.graph);
            std::cerr << "\n";
        }

        // We might be done already.
        if (start1 == path1.begin())
        {
            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            {
                std::cerr << "    NO MORE IN PATH1\n";
                std::cerr << "      RESULTING PATH AFTER ZIPPERING LEFT\n      ";
                printPath(std::cerr, path2, state.graph);
                std::cerr << "\n";
            }
            return true;
        }

        // Zipper right recursion.
        SEQAN_ASSERT_NOT(std::prev(start1) == path1.begin());
        return zipperLeft(*std::prev(start1), std::prev(start1, 2), std::next(dest));
    }

    // Compute candidates when fitting start1 into path2 (starting at start2) using guiding edge arc.  Case for
    // zippering left.
    std::vector<TArcListIter> computeCandidates(lemon::SmartGraph::Arc arc, TArcListIter start1,
                                                TArcListIter start2)
    {
        std::vector<TArcListIter> candidates;

        // We will look for placements that are within the deviation of the guiding arc.
        double const mean = label(arc).lengthMean, stdDev = label(arc).lengthStdDev;
        double const minLen = mean - options.mult * stdDev;
        double const maxLen = mean + options.mult * stdDev;

        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            std::cerr << "      minLen=" << minLen << "\tmaxLen=" << maxLen << "\n";

        // Helper function that returns whether two intervals overlap and overlap significantly.
        auto overlaps = [](double b1, double e1, double b2, double e2) { return b2 < e1 && b1 < e2; };
        auto ovlSign = [&](double b1, double e1, double b2, double e2) {
            if (e2 - b2 <= options.significantOverlap)
                return false;
            return (overlaps(b1, e1, b2 + options.significantOverlap, e2) &&
                    overlaps(b1, e1, b2, e2 - options.significantOverlap));
        };

        // Collect edges on path2 that overlap with [minLen, maxLen).  We add path2.end() if it would fit right of
        // path2.  These candidates are later filtered based on the overlap list.
        double pathLen = 0, pathStdDev = 0;
        for (TArcListIter it2 = start2; /*see below*/; --it2)
        {
            double it2Len = label(*it2).lengthMean, it2SD = label(*it2).lengthStdDev;

            // If *it2 is a contig edge then we check whether the active contig overlaps with it.  If *it2 is a link
            // edge then the active contig has to fit completely into the link edge.
            if (label(*it2).isContig())
            {
                double d = options.mult * pathStdDev;
                if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                    std::cerr << "        ovlSign(" << pathLen - d << ", "
                              << pathLen + d + label(*it2).lengthMean << ", "
                              << minLen << ", " << maxLen + label(*start1).lengthMean << ")\n";
                if (ovlSign(-(pathLen + d), -pathLen + d + it2Len,  // *it2's interval
                            -minLen - label(*start1).lengthMean, -maxLen))  // *start1's interval
                {
                    unsigned seq0 = label(*it2).idx;
                    unsigned seq1 = label(*start1).idx;
                    std::pair<int, int> range0(-pathLen, -pathLen + it2Len);
                    std::pair<int, int> range1(-mean - label(*start1).lengthMean, -mean);
                    int band = d + options.mult * stdDev;
                    if (range0.first > range1.first)
                    {
                        std::swap(seq0, seq1);
                        std::swap(range0, range1);
                    }
                    int overlapLen = std::min(range0.second, range1.second) - std::max(range0.first, range1.first);
                    SEQAN_ASSERT_GT_MSG(overlapLen + band, 0, "Claimed to be significant previously!");
                    if (state.overlapExists(seq0, seq1, overlapLen, std::max(band, 5)))
                    {
                        candidates.push_back(it2);
                        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                            std::cerr << "YES, OVERLAP EXISTS: seq0" << seq0 << ", seq1" << seq1
                                      << ", overlapLen=" << overlapLen << ", band=" << d + options.mult * stdDev << "\n";
                    }
                }
            }
            else
            {
                // First condition: Contig must fit into link edge length-wise (slack > 0)
                // Second condition: Must be able to start far enough to the left (within slack).
                double slack = it2Len + options.mult * it2SD - label(*start1).lengthMean;
                double d = options.mult * pathStdDev;
                if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                    std::cerr << "        d == " << d << ", slack == " << slack << "\n"
                              << "          overlaps(" << minLen << ", " << maxLen << ", "
                              << pathLen - d << ", "
                              << pathLen + d + it2Len << ")\n";
                if (slack > 0 && overlaps(minLen, maxLen, pathLen - d, pathLen + d + it2Len))
                    candidates.push_back(it2);
            }

            // Traverse edge *it2 and increase length/standard deviation.
            if (it2 != start2)
            {
                pathLen += label(*it2).lengthMean;
                pathStdDev += label(*it2).lengthStdDev;
            }

            if (it2 == path2.begin())
                break;  // done
        }

        // Reverse candidates to obtain leftmost feasible first.  This protects us against "pushing" contigs along the
        // path if they fully overlap with their right neighbour.
        reverse(candidates.begin(), candidates.end());

        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            std::cerr << "maxLen == " << maxLen << ", pathLen == " << pathLen
                      << ", options.mult * pathStdDev == " << options.mult * pathStdDev << "\n";
        if (maxLen > pathLen - options.mult * pathStdDev)
            candidates.push_back(path2.end());

        return candidates;
    }
};

// Algorithm for merging two paths.

class PathMerger
{
    // The state with the graph, labels etc.
    PathMergingAlgoState & state;
    // The guiding edge.
    lemon::SmartGraph::Arc guidingArc;
    // The input paths as lists, path2 will be copied to output after merging.  We use lists here since we need the
    // splicing but also reverse iteration in case of zippering to the left.
    typedef std::list<lemon::SmartGraph::Arc> TArcList;
    typedef TArcList::iterator TArcListIter;
    TArcList path1, path2;
    // Options.
    PathMergingOptions const & options;

public:

    std::string arcStr(lemon::SmartGraph::Arc e)
    {
        std::stringstream ss;
        ss << "(" << state.graph.id(state.graph.source(e)) << ", "
           << state.graph.id(state.graph.target(e)) << ")";
        return ss.str();
    }

    PathMerger(PathMergingAlgoState & state,
               lemon::SmartGraph::Arc guidingArc,
               lemon::SimplePath<lemon::SmartGraph> const & origPath1,
               lemon::SimplePath<lemon::SmartGraph> const & origPath2,
               PathMergingOptions const & options) :
            state(state), guidingArc(guidingArc), options(options)
    {
        copyPath(path1, origPath1);
        copyPath(path2, origPath2);
    }

    // Target of guidingArc must be in origPath1.
    bool run(lemon::SimplePath<lemon::SmartGraph> & result)
    {
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
        {
            std::cerr << "BEGIN PATH MERGING\n"
                      << "    PATH1\t";
            printPath(std::cerr, path1, state.graph);
            std::cerr << "\n"
                      << "    PATH2\t";
            printPath(std::cerr, path2, state.graph);
            std::cerr << "\n"
                      << "    GUIDING EDGE\t" << arcStr(guidingArc) << "\n";
        }
        // Get starting point in both paths.  Both the index into origPath1 and the iterator into path point to a contig
        // edge.  In the case of path1, it is the contig that is to be placed.  In the case of path, it is the contig
        // that the path1 contig has to be placed after.
        std::pair<TArcListIter, TArcListIter> start = findStartPoints(guidingArc);

        // Zipper left starting at the given start position.
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            std::cerr << "  ZIPPERING RIGHT\n";
        ZipperRightHelper zipperRight(state, path1, path2, options);
        if (!zipperRight.run(guidingArc, start.first, start.second))
        {
            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                std::cerr << "  => ZIPPERING RIGHT FAILED!\n";
            return false;
        }

        // Zipper right.  We use the same start/end position as for the right but the zipperLeftWrapper() function will
        // interpret it differently.
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            std::cerr << "  ZIPPERING LEFT...\n";
        ZipperLeftHelper zipperLeft(state, path1, path2, options);
        if (!zipperLeft.run(guidingArc, start.first, start.second))
        {
            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                std::cerr << "  => ZIPPERING LEFT FAILED!\n";
            return false;
        }

        // Finally, copy out the merged path2 to the result.
        copyPath(result, path2);
        return true;
    }

private:

    template <typename TEdge>
    ReductionGraphLabel & label(TEdge const & edge)
    {
        return state.labels[edge];
    }

    template <typename TEdge>
    ReductionGraphLabel const & label(TEdge const & edge) const
    {
        return state.labels[edge];
    }

    // Search for start point in each graph.  The start point must exist and represent two contig edges.  In path2, we
    // look for an edge whose target is equal to the source of arc, in path1, we look for an arc whose source is equal
    // to the target of arc.
    std::pair<TArcListIter, TArcListIter> findStartPoints(lemon::SmartGraph::Arc arc)
    {
        // Look in path1.
        TArcListIter it1 = path1.begin();
        for (; it1 != path1.end(); ++it1)
            if (state.graph.source(*it1) == state.graph.target(arc))
                break;
        SEQAN_CHECK(it1 != path1.end(), "source(arc) not found in path1!");
        SEQAN_CHECK(state.labels[*it1].isContig(), "Must start off contig edge.");

        // Look in path2.
        TArcListIter it2 = path2.begin();
        for (; it2 != path2.end(); ++it2)
            if (state.graph.target(*it2) == state.graph.source(arc))
                break;
        SEQAN_CHECK(it2 != path2.end(), "source(arc) not found in path2!");
        SEQAN_CHECK(state.labels[*it2].isContig(), "Must start off contig edge.");

        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            std::cerr << "  STARTING OFF\n"
                      << "    IN PATH1: " << arcStr(*it1) << "\n"
                      << "    IN PATH2: " << arcStr(*it2) << "\n";

        return std::make_pair(it1, it2);
    }
};


// Algorithm implementation for the actual path merging.

class PathMergingAlgo
{
    // The input links.
    std::vector<MateLink> const & links;
    // The information about the contigs.
    std::vector<ContigEdgeLabel> const & contigInfos;
    // Configuration for the path merging.
    PathMergingOptions const & options;

    // The state of the path merger.
    PathMergingAlgoState state;

public:
    PathMergingAlgo(std::vector<MateLink> const & links,
                    std::vector<ContigEdgeLabel> const & contigInfos,
                    TOverlapExistsFunc overlapExists,
                    TComputeOverlapFunc computeOverlap,
                    PathMergingOptions const & options)
            : links(links), contigInfos(contigInfos), options(options), state(overlapExists, computeOverlap)
    {}

    // Execute the greedy path merging algorithm.
    void run(ScaffoldingResult & result)
    {
        // Sort links by weight.
        std::vector<MateLink> sortedLinks(links);
        std::stable_sort(sortedLinks.begin(), sortedLinks.end(),
                         [](MateLink lhs, MateLink rhs) {
                             return lhs.label.count > rhs.label.count;
                         });
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
        {
            std::cerr << "SORTED LINKS\n"
                      << "============\n\n";
            for (auto link : sortedLinks)
                std::cerr << link << "\n";
            std::cerr << "\n";

            std::cerr << "CONTIG INFOS\n"
                      << "============\n\n";
            for (auto link : contigInfos)
                std::cerr << link << "\n";
            std::cerr << "\n";
        }

        // Build underlying graph and initial path set, mark contig edges as selected.
        init(sortedLinks);

        // Run main loop.
        mainLoop(sortedLinks);

        // Write out the result.
        collectResult(result);
    }

private:

    void collectResult(ScaffoldingResult & result)
    {
        for (auto const & pathWithID : state.paths)
        {
            ScaffoldingResult::TScaffold scaffold;
            auto const & path = pathWithID.second;
            int pos = 0, posSD = 0;
            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            {
                std::cerr << "RESULT PATH\t";
                printPath(std::cerr, path, state.graph);
                std::cerr << "\n";
                for (int i = 0; i < path.length(); ++i)
                    std::cerr << "  " << state.labels[path.nth(i)].lengthMean << " "
                              << state.labels[path.nth(i)].lengthStdDev << "\n";
            }
            for (int i = 0; i < path.length(); i += 2)
            {
                SEQAN_ASSERT_MSG(state.labels[path.nth(i)].isContig(), "i == %d", i);
                if (i + 1 < path.length())
                    SEQAN_ASSERT_MSG(!state.labels[path.nth(i + 1)].isContig(), "i == %d", i);
                scaffold.push_back(PositionedContig(pos, posSD, state.labels[path.nth(i)].idx,
                                                    state.labels[path.nth(i)].lengthMean));
                if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                    std::cerr << scaffold.back() << "\n";
                pos += state.labels[path.nth(i)].lengthMean;
                if (i + 1 < path.length())
                {
                    pos += state.labels[path.nth(i + 1)].lengthMean;
                    posSD = state.labels[path.nth(i + 1)].lengthStdDev;  // NO plus!
                }
            }

            std::sort(scaffold.begin(), scaffold.end(), [](PositionedContig lhs, PositionedContig rhs) {
                        return lhs.pos < rhs.pos;
                      });

            result.scaffolds.push_back(scaffold);
        }
    }

    void mainLoop(std::vector<MateLink> const & sortedLinks)
    {
        for (auto indexed : enumerate(sortedLinks))
        {
            auto const & link = indexed.second;

            // Let v, w denote the two nodes connected by link.
            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                std::cerr << "Edge (" << (2 * link.source + 1) << ", " << (2 * link.target) << ")\n";
            lemon::SmartGraph::Node v = state.nodes[2 * link.source + 1];
            lemon::SmartGraph::Node w = state.nodes[2 * link.target];
            int pathIdxV = state.pathForNode[v];
            int pathIdxW = state.pathForNode[w];
            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                std::cerr << "v == " << pathIdxV << ", w == " << pathIdxW << "\n";
            SEQAN_ASSERT(state.paths.count(pathIdxV));
            SEQAN_ASSERT(state.paths.count(pathIdxW));

            // Ignore if adjacent to the same path.
            if (pathIdxV == pathIdxW)
                continue;

            // Obtain shortcuts to paths.
            lemon::SimplePath<lemon::SmartGraph> & pathV = state.paths[pathIdxV];
            lemon::SimplePath<lemon::SmartGraph> & pathW = state.paths[pathIdxW];

            // Try to merge the two paths.
            lemon::SimplePath<lemon::SmartGraph> path;
            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                std::cerr << "MERGING PATHS\t" << pathIdxV << "\tand\t" << pathIdxW << "\tguided by\t"
                          << "edge #" << indexed.first << " (len= " << indexed.second.label.lengthMean << ", SD=" << indexed.second.label.lengthStdDev << ")\n";
            if (mergePaths(path, state.edges[indexed.first], pathW, pathV))
            {
                for (int i = 0; i + 1 < path.length(); ++i)
                    SEQAN_CHECK(state.graph.target(path.nth(i)) == state.graph.source(path.nth(i + 1)), "");
                // Check for invalid overlaps generated by zippering that were not detected during zippering.
                bool overlapsOK = !hasConflictingOverlap(path);

                // If the increase of happiness is higher than the increase of unhappiness then replace pathV and pathW
                // by path.
                if (overlapsOK && (happiness(path) - happiness(pathV) - happiness(pathW) > 0))
                {
                    int pathIdx = state.nextPathID++;
                    state.paths[pathIdx] = path;
                    replacePaths(pathV, pathW, pathIdx);
                    state.paths.erase(pathIdxV);
                    state.paths.erase(pathIdxW);
                }
            }
        }
    }

    // Returns whether there are no conflicts in the path indicated by the overlaps implied by the path.
    bool hasConflictingOverlap(lemon::SimplePath<lemon::SmartGraph> const & path) const
    {
        if (path.empty())
            return false;

        // Build list of positioned contigs.
        struct PosContig
        {
            PosContig() = default;
            PosContig(int contigID, int length, FuzzyValue<double> pos) :
                    contigID(contigID), length(length), pos(pos) {}
            PosContig(int contigID, int length, double p, double posSD) :
                    contigID(contigID), length(length), pos(p, posSD) {}

            int contigID { 0 };
            int length { 0 };
            FuzzyValue<double> pos;
        };
        std::vector<PosContig> positionedContigs;

        auto toPos = [](ReductionGraphLabel const & label) {
            return FuzzyValue<double>(label.lengthMean, label.lengthStdDev); };

        auto label = [&](lemon::SmartGraph::Arc const & arc) { return state.labels[arc]; };

        for (lemon::SimplePath<lemon::SmartGraph>::ArcIt first(path), prev(first), arc = first;
             arc != lemon::INVALID; ++arc)
        {
            if (arc == first)
            {
                positionedContigs.push_back(PosContig(label(arc).idx, label(arc).lengthMean, 0, 0));
            }
            else
            {
                auto prevArc = arc;
                ++arc;
                positionedContigs.push_back(PosContig(label(arc).idx, label(arc).lengthMean, toPos(label(prevArc))));
                positionedContigs.back().pos += toPos(label(prev));
            }
            prev = arc;
        }

        // Check whether contigs overlap.
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            std::cerr << "Checking whether overlaps exist...\n";
        for (auto it = positionedContigs.begin(); it != positionedContigs.end(); ++it)
            for (auto it2 = std::next(it); it2 != positionedContigs.end(); ++it2)
            {
                if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                    std::cerr << "it->pos.value=" << it->pos.value << "\n"
                              << "it->length=" << it->length << "\n"
                              << "it2->pos.value=" << it2->pos.value << "\n"
                              << "it2->length=" << it2->length << "\n";
                int ovlLen = (it->pos.value + it->length) - it2->pos.value;
                if (ovlLen > options.significantOverlap)
                    if (!state.overlapExists(it->contigID, it2->contigID, ovlLen,
                                             std::max(MIN_BAND, (int)it2->pos.fuzziness())))
                        return true;  // overlap supposed to be significant but nonexistant
            }

        return false;  // no conflicting contig found
    }

    // Return whether path1 and path2 could be merged into path, guided by e.
    bool mergePaths(lemon::SimplePath<lemon::SmartGraph> & path,
                    lemon::SmartGraph::Arc e,
                    lemon::SimplePath<lemon::SmartGraph> const & path1,
                    lemon::SimplePath<lemon::SmartGraph> const & path2)
    {
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
        {
            std::cerr << "MERGE PATHS GUIDED BY ("
                      << state.graph.id(state.graph.source(e)) << ", "
                      << state.graph.id(state.graph.target(e)) << ")\n"
                      << "  path1 == ";
            printPath(std::cerr, path1, state.graph);
            std::cerr << " INTO\n"
                      << "  path2 == ";
            printPath(std::cerr, path2, state.graph);
            std::cerr << "\n";
        }

        PathMerger merger(state, e, path1, path2, options);
        bool result = merger.run(path);

        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            std::cerr << "MERGING SUCCESS? " << result << "\n";
        return result;
    }

    // Build state and initial graph.
    void init(std::vector<MateLink> const & sortedLinks)
    {
        // Add nodes for contigs.
        for (auto indexedInfo : enumerate(contigInfos))
        {
            auto const & info = indexedInfo.second;
            // Add nodes for the contig.
            lemon::SmartGraph::Node u = state.graph.addNode(), v = state.graph.addNode();
            state.nodes.push_back(u);
            state.nodes.push_back(v);
            // Add edge for the contig.
            lemon::SmartGraph::Arc e = state.graph.direct(state.graph.addEdge(u, v), true);
            state.labels[e] = ReductionGraphLabel(ReductionGraphLabel::CONTIG, info.length, 0, indexedInfo.first);
            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                std::cerr << "LABEL FOR EDGE\t"
                          << state.graph.id(u) << " -> " << state.graph.id(v)
                          << "\t" << state.graph.id(e) << "\t IS \t" << state.labels[e] << "\n";
            // Update edge labels.
            state.selected[e] = true;
            // Create path.x
            unsigned pathID = state.nextPathID++;
            state.paths[pathID].addBack(e);
            state.pathForNode[u] = pathID;
            state.pathForNode[v] = pathID;
        }

        // Add edges for links.
        for (auto indexedLink : enumerate(sortedLinks))
        {
            auto const & link = indexedLink.second;
            unsigned idxU = 2 * link.source + 1;
            unsigned idxV = 2 * link.target;
            lemon::SmartGraph::Arc e = state.graph.direct(state.graph.addEdge(state.nodes[idxU], state.nodes[idxV]), true);
            state.labels[e] = ReductionGraphLabel(ReductionGraphLabel::MATE, link.label.lengthMean, link.label.lengthStdDev,
                                                  indexedLink.first, link.label.count, link.label.weight);
            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                std::cerr << "    LABEL FOR EDGE\t" << idxU << " -> " << idxV << "\t"
                          << state.graph.id(e) << "\t IS \t" << state.labels[e] << "\n";
            state.edges.push_back(e);
        }
    }

    // Return amount of happiness with respect to the given path.
    int happiness(lemon::SimplePath<lemon::SmartGraph> const & path) const
    {
        if (path.empty())
            return 0;

        // Build position map for the vertices.
        std::map<int, int> id2pos;
        id2pos[state.graph.id(state.graph.source(path.front()))] = 0;
        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
            std::cerr << "      id2pos[" << state.graph.id(state.graph.source(path.front())) << "] = " << 0 << "\n";
        int currentPos = 0;
        for (int i = 0; i < path.length(); ++i)
        {
            currentPos += state.labels[path.nth(i)].lengthMean;
            id2pos[state.graph.id(state.graph.target(path.nth(i)))] = currentPos;
            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                std::cerr << "      id2pos[" << state.graph.id(state.graph.target(path.nth(i))) << "] = " << currentPos << "\n";
        }

        // Check happiness with respect to links from the input.
        int result = 0;
        for (auto link : links)
        {
            if (!id2pos.count(2 * link.source)  || !id2pos.count(2 * link.target))
                continue;  // Skip, not both on path.
            bool happy = true;
            if (id2pos[2 * link.source] > id2pos[2 * link.target])  // correct orientation
                happy = false;
            int distance = id2pos[2 * link.target] - id2pos[2 * link.source + 1];
            int band = options.mult * std::max(link.label.lengthStdDev, 1.0);
            band = std::max(MIN_BAND, band);
            int minLen = link.label.lengthMean - band;
            int maxLen = link.label.lengthMean + band;
            happy = happy && (distance >= minLen) && (distance <= maxLen);
            if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
                std::cerr << "link.source == " << link.source << ", link.target == " << link.target << "\n"
                          << "distance == " << distance << ", minLen = " << minLen << ", maxLen = " << maxLen << "\n"
                          << "id2pos[2 * link.source + 1] == " << id2pos[2 * link.source + 1] << "\t"
                          << "id2pos[2 * link.target] == " << id2pos[2 * link.target] << "\n"
                          << "happy = " << happy << "\n";
            result += (happy ? 1 : -1) * (int)link.label.count;
        }

        if (options.verbosity >= PathMergingOptions::VERY_VERBOSE)
        {
            std::cerr << "HAPPYNESS OF FOLLOWING PATH IS " << result << "\n  ";
            printPath(std::cerr, path, state.graph);
            std::cerr << "\n\nLENGTHS";
            for (int i = 0; i < path.length(); ++i)
                std::cerr << "\t" << state.labels[path.nth(i)].lengthMean;
            std::cerr << "\n";
        }

        return result;
    }

    // Update book-keeping such that the paths are replaced by the path with idx newPathIdx.
    void replacePaths(lemon::SimplePath<lemon::SmartGraph> const & path1,
                      lemon::SimplePath<lemon::SmartGraph> const & path2,
                      int newPathIdx)
    {
        // Update selection status.
        selectPath(path1, false);
        selectPath(path2, false);
        SEQAN_CHECK(state.paths.count(newPathIdx), "Must be known!");
        selectPath(state.paths[newPathIdx], true);

        // Update vertices on path1 and path2 to point to newPathIdx.
        for (int i = 0; i < path1.length(); ++i)
        {
            state.pathForNode[state.graph.source(path1.nth(i))] = newPathIdx;
            state.pathForNode[state.graph.target(path1.nth(i))] = newPathIdx;
        }
        for (int i = 0; i < path2.length(); ++i)
        {
            state.pathForNode[state.graph.source(path2.nth(i))] = newPathIdx;
            state.pathForNode[state.graph.target(path2.nth(i))] = newPathIdx;
        }
    }

    // Mark all edges from path as de/selected.
    void selectPath(lemon::SimplePath<lemon::SmartGraph> const & path, bool flag)
    {
        for (int i = 0; i < path.length(); ++i)
            state.selected[path.nth(i)] = flag;
    }
};

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Function greedyPathMerging()
// ----------------------------------------------------------------------------

void greedyPathMerging(ScaffoldingResult & result,
                       std::vector<MateLink> const & links,
                       std::vector<ContigEdgeLabel> const & contigInfos,
                       TOverlapExistsFunc overlapExists,
                       TComputeOverlapFunc computeOverlap,
                       PathMergingOptions const & options)
{
    PathMergingAlgo algo(links, contigInfos, overlapExists, computeOverlap, options);
    algo.run(result);
    result.shiftScaffolds();
}

}  // namespace scaffolder
