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

#ifndef SEQAN_APPS_RR_UTILS_H
#define SEQAN_APPS_RR_UTILS_H

// includes 
#include <cmath>
#include <seqan/sequence.h>

using namespace seqan;

/*
    approximated cumulative normal distribution
    http://www.sitmo.com/doc/Calculating_the_Cumulative_Normal_Distribution
    Abromowitz and Stegun approximation 
*/

double cumulated_normal(const double x)
{
  const double b1 =  0.319381530;
  const double b2 = -0.356563782;
  const double b3 =  1.781477937;
  const double b4 = -1.821255978;
  const double b5 =  1.330274429;
  const double p  =  0.2316419;
  const double c  =  0.39894228;

  if(x >= 0.0) {
      double t = 1.0 / ( 1.0 + p * x );
      return (1.0 - c * exp( -x * x / 2.0 ) * t *
      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
  }
  else {
      double t = 1.0 / ( 1.0 - p * x );
      return ( c * exp( -x * x / 2.0 ) * t *
      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    }
}

double cumulated_normal(const double x, const double mu, const double sd){
    return cumulated_normal(x/sd -  mu/sd);
}

int factorial(const int n)
{
    if(n < 2) return 1;
    else
    {
        int f = 2;
        for(int s = 3; s <= n;++s) f *= s;
        return f;
    }
}

double poisson(const int k, const double lambda)
{    return (1.0/static_cast<double>(factorial(k))) * pow(lambda,static_cast<double>(k)) * exp((-1.0) * lambda);    }

double p_value(const double x,const double mu, const double sd){
    return (1.0 - cumulated_normal(x-1,mu,sd));
}

double p_value(const unsigned n_u, const double lambda)
{
    double p=0.0;
    for(unsigned i = 0; i < n_u;++i)
    { p += poisson(i,lambda); }
    
    return (1 - p);
}

/////////////////////////////////////////////////////////////////////////////
// some accessor stuff for the common repsep structures
// 1. Triple< TAlphabet, .. gapped sequence char
//            TId,       .. readId
//            TReadPos > .. position of the char in ungapped read (gap is length(read) + 1)
/////////////////////////////////////////////////////////////////////////////
template<typename TAlphabet, typename TId, typename TReadPos>
TAlphabet _sequenceCharacter(Triple<TAlphabet, TId, TReadPos> arc) 
{
    return arc.i1;
}

template<typename TAlphabet, typename TId, typename TReadPos>
TId _readId(Triple<TAlphabet, TId, TReadPos> arc) 
{
//IOREV _notio_
    return arc.i2;
}

template<typename TAlphabet, typename TId, typename TReadPos>
TReadPos _positionInRead(Triple<TAlphabet, TId, TReadPos> arc) 
{
    return arc.i3;
}

///////////////////////////////////

template<typename TAlphabet, typename TId, typename TReadPos>
void _setSequenceCharacter(Triple<TAlphabet, TId, TReadPos> & arc, TAlphabet val) 
{
    arc.i1 = val;
}

template<typename TAlphabet, typename TId, typename TReadPos>
void _setReadId(Triple<TAlphabet, TId, TReadPos> & arc, TId id) 
{
    arc.i2 = id;
}

template<typename TAlphabet, typename TId, typename TReadPos>
void _setPositionInRead(Triple<TAlphabet, TId, TReadPos> & arc, TReadPos pos) 
{
    arc.i3 = pos;
}

/////////////////////////////////////////////////////////////////////////////
// 1. Triple< TAlphabet, .. gapped sequence char
//            TReadPos > .. position of the char in ungapped read (gap is length(read) + 1)
/////////////////////////////////////////////////////////////////////////////
template<typename TAlphabet, typename TReadPos>
TAlphabet _sequenceCharacter(Pair<TAlphabet, TReadPos> arc) 
{
    return arc.i1;
}

template<typename TAlphabet, typename TReadPos>
TReadPos _positionInRead(Pair<TAlphabet, TReadPos> arc) 
{
    return arc.i2;
}

///////////////////////////////////

template<typename TAlphabet, typename TReadPos>
void _setSequenceCharacter(Pair<TAlphabet, TReadPos> & arc, TAlphabet val) 
{
    arc.i1 = val;
}

template<typename TAlphabet, typename TReadPos>
void _setPositionInRead(Pair<TAlphabet, TReadPos> & arc, TReadPos pos) 
{
    arc.i2 = pos;
}

/////////////////////////////////////////////////////////////////////////////
// (currently unused parser)
// TODO: check removal
template<typename TReadId>
Pair<TReadId,TReadId>
_parseIds(String<char> & fastaTag){
//IOREV _delcandidate_ marked for removal in code comment
    // format is >257,288[id=1,mateId=431]
    Pair<TReadId,TReadId> ret;

    typedef typename Iterator< String<char> >::Type TIter;
    TIter iter = begin(fastaTag);

    // first find start of interesting part [
    while(value(iter) != '[') ++iter;
    // now step to =
    while(value(iter) != '=') ++iter;
    ++iter;
    ret.i1 = value(iter);
    ++iter;
    while(value(iter) != ','){
        append(ret.i1,value(iter));
        ++iter;
    }
    // now step to =
    while(value(iter) != '=') ++iter;
    ++iter;
    ret.i2 = value(iter);
    ++iter;
    while(value(iter) != ','){
        append(ret.i2,value(iter));
        ++iter;
    }

    // after this there will only be the libraryId and so

    return ret;
}

template <typename TString>
int _parseGroupId(TString & tag)
{
//IOREV _notio_
    String<char> groupId;

    for (unsigned i = 0; i < length(tag); ++i)
    {
        if(static_cast<char>(value(tag,i)) == '-')
        {
            ++i;
            while(i < length(tag))
            {
                append(groupId, static_cast<char>(value(tag,i)));
                ++i;
            }
        }
    }
    return atoi(toCString(groupId));
}

/////////////////////////////////////////////////////////////////////////////

#endif
