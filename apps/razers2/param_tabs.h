/*==========================================================================
             RazerS - Fast Read Mapping with Controlled Loss Rate
                   http://www.seqan.de/projects/razers.html

 ============================================================================
  Copyright (C) 2008 by David Weese

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your options) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ==========================================================================
  Data tables with precomputed parameters.  Generated from gapped_params
  folder.
 ==========================================================================*/

#ifndef APPS_RAZERS_PARAM_TABS_H_
#define APPS_RAZERS_PARAM_TABS_H_

// NOTE: DDDoc style comments are used here but not activated on purpose. Only there for RazerS developers.

#include <seqan/sequence.h>

struct GappedParamsRecord
{
    // Length of the read, "N??" in gapped_params file name.
    unsigned readLength;
    // Distance type, "H" or "L" in gapped_params file name.
    char type;

    // Number of errors, first column in gapped params file.
    unsigned errors;
    // Shape, second column in gapped params file.
    char const * shape;
    // Threshold, third column in gapped params file.
    unsigned t;
    // Loss rate, fourth column in gapped params file.
    double lossRate;
    // Potential match count on simulated data from which the parameters were
    // computed.
    unsigned measure;
};

bool getGappedParamsRecords(seqan::String<GappedParamsRecord> & records,
                            unsigned n,
                            char errorModel);

#endif  // #ifndef APPS_RAZERS_PARAM_TABS_H_
