/*==========================================================================
             RazerS - Fast Read Mapping with Controlled Loss Rate
                   http://www.seqan.de/projects/razers.html

 ============================================================================
  Copyright (C) 2012 by David Weese

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
  ==========================================================================*/

#include <iostream>

#include "param_tabs.h"

// The 68k lines of parameters are included from another file.
static GappedParamsRecord RECORDS[] =
{
     #include "param_tabs.inc"
};

bool getGappedParamsRecords(seqan::String<GappedParamsRecord> & records,
                            unsigned n,
                            char errorModel)
{
    if (n < 15u || n > 75u)
        return false;  // We do not have parameter for these settings.

    if (errorModel != 'L' && errorModel != 'H')
        return false;  // Invalid error model.

    // We can iterate until readLength == 0 since the generator script creates
    // a terminator record with this property and read length 0 does not make
    // sense otherwise.
    for (unsigned i = 0; RECORDS[i].readLength != 0; ++i)
        if (RECORDS[i].readLength == n && RECORDS[i].type == errorModel)
            appendValue(records, RECORDS[i]);

    return true;
}
