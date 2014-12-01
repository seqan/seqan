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

#include "gpm.h"

#include <list>
#include <vector>

#include "scaffolder/mate_link.h"
#include "scaffolder/gpm_options.h"
#include "scaffolder/utils.h"

namespace scaffolder {

namespace {  // anonymous namespace

// --------------------------------------------------------------------------
// Class LinkBundler
// --------------------------------------------------------------------------

class LinkBundler
{
    std::vector<MateLink> links;
    PathMergingOptions const & options;

public:
    LinkBundler(std::vector<MateLink> const & links,
                PathMergingOptions const & options) :
            links(links), options(options)
    {}

    void run(std::vector<MateLink> & bundled)
    {
        auto lt = [](MateLink const & lhs, MateLink const & rhs) {
            return std::make_pair(lhs.source, lhs.target) < std::make_pair(rhs.source, rhs.target);
        };
        std::sort(links.begin(), links.end(), lt);

        // Enumerate links for each pair of vertices.
        auto cmp = [](MateLink lhs, MateLink rhs) { return lhs.source == rhs.source && lhs.target == rhs.target; };
        auto r = er(links.begin(), links.end(), cmp), e = r.make_end();
        for (; r != e; ++r)
        {
            std::list<MateLink> localLinks(r.begin(), r.end());
            while (!localLinks.empty())
            {
                MateLink firstLink = localLinks.front();
                mergeLinks(firstLink, localLinks, std::back_inserter(bundled));
            }
        }
    }

private:

    // Merge remaining links in link with link.
    template <typename TIt>
    void mergeLinks(MateLink link, std::list<MateLink> & links, TIt outIt) const
    {
        int count = 0;
        double weight = 0;
        double p = 0, q = 0;
        // Iterate over links and remove the entries that match link.
        for (auto it = links.begin(); it != links.end(); /* see below */)
        {
            if (abs(link.label.lengthMean - it->label.lengthMean) > options.mult * link.label.lengthStdDev)
            {
                // no bundling
                ++it;
            }
            else  // bundling
            {
                count += it->label.count;
                weight += it->label.weight;
                p += count * it->label.lengthMean / it->label.lengthStdDev / it->label.lengthStdDev;
                q += count / it->label.lengthStdDev / it->label.lengthStdDev;
                auto it2 = it;
                ++it;
                links.erase(it2);
            }
        }

        // We modify link and write it to outIt.  We keep source and target but overwrite the label.
        link.label = MateEdgeLabel(p / q, 1 / sqrt(q), count, weight);
        *outIt++ = link;
    }
};

}  // anonymous namespace


// --------------------------------------------------------------------------
// Function bundleLinks()
// --------------------------------------------------------------------------

void bundleLinks(std::vector<MateLink> & bundled,
                 std::vector<MateLink> const & links,
                 PathMergingOptions const & options)
{
    LinkBundler bundler(links, options);
    bundler.run(bundled);
}

std::vector<MateLink> bundleLinks(std::vector<MateLink> const & links,
                                  PathMergingOptions const & options)
{
    std::vector<MateLink> bundled;
    bundleLinks(bundled, links, options);
    return bundled;
}

}  // namespace scaffolder
