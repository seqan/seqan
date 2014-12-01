/*=========================================================================
  Copyright (C) 2009 by Stephan Aiche

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
  $Id$
 ==========================================================================*/

#ifndef REPSEP_HEADER_RGRAPH_SH
#define REPSEP_HEADER_RGRAPH_SH

//////////////////////////////////////////////////////////////////////////////

struct GraphScoring{
    typedef double TScoreValue;

    TScoreValue matchScore;
    TScoreValue mismatchScore;
    TScoreValue matePairScore;

    GraphScoring(){
        matchScore = -0.3f;
        mismatchScore = 1.5f;
        matePairScore = -3.5f;    
    }
};

//////////////////////////////////////////////////////////////////////////////

GraphScoring::TScoreValue & matchScore(GraphScoring & me)
{
    return me.matchScore;
}

GraphScoring::TScoreValue const & matchScore(GraphScoring const & me)
{
    return me.matchScore;
}

//////////////////////////////////////////////////////////////////////////////


GraphScoring::TScoreValue & mismatchScore(GraphScoring & me)
{
    return me.mismatchScore;
}


GraphScoring::TScoreValue const & mismatchScore(GraphScoring const & me)
{
    return me.mismatchScore;
}

//////////////////////////////////////////////////////////////////////////////


GraphScoring::TScoreValue & matePairScore(GraphScoring & me)
{
    return me.matePairScore;
}


GraphScoring::TScoreValue const & matePairScore(GraphScoring const & me)
{
    return me.matePairScore;
}

//////////////////////////////////////////////////////////////////////////////

#endif
