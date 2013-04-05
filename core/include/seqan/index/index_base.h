// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_INDEX_BASE_H
#define SEQAN_HEADER_INDEX_BASE_H

//#define SEQAN_TEST_INDEX

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// needful forward declarations

	// suffix array construction specs
	struct Skew3;
	struct Skew7;
	struct LarssonSadakane;
	struct ManberMyers;
	struct SAQSort;
	struct QGramAlg;

	// lcp table construction algorithms
	struct Kasai;
	struct KasaiOriginal;	// original, but more space-consuming algorithm

	// enhanced suffix array construction algorithms
	struct Childtab;
	struct Bwt;

/**
.Tag.Index Find Algorithm
..summary:Tag to specify the index search algorithm.
..remarks:These tag can be used to specify the @Function.find@ algorithm 
for @Class.Index@ based substring searches.
..cat:Index

..tag.EsaFindMlr:Binary search with mlr-heuristic.
...remarks:Exact string matching using a suffix array binary search with the mlr-heuristic.

..tag.EsaFindLcpe:Binary search using lcp values.
...remarks:Exact string matching using a suffix array binary search and a lcp-interval tree.

..tag.FinderSTree:Suffix tree search.
...remarks:Exact string matching using a suffix tree.

..see:Class.Finder
..see:Spec.IndexEsa
..see:Spec.IndexQGram
..include:seqan/index.h
*/

	// finder tags
    struct FinderMlr_;     // simple Suffix Array finder with mlr-heuristic
    struct FinderLcpe_;    // Suffix Array finder using an enhanced LCP-Table
    struct FinderSTree_;    // Suffix Array finder using an enhanced LCP-Table

    typedef Tag<FinderMlr_> const EsaFindMlr;
    typedef Tag<FinderLcpe_> const EsaFindLcpe;
    typedef Tag<FinderSTree_> const FinderSTree;

	template <typename TSpec = void>
	struct IndexEsa {};


//////////////////////////////////////////////////////////////////////////////
/**
.Metafunction.DefaultIndexSpec:
..cat:Index
..summary:Default @Class.Index@ specialization type.
..signature:DefaultIndexSpec<TText>::Type
..class:Class.Index
..param.TText:The given text type.
..returns:Can be @Spec.IndexEsa@ or $IndexQGram$, etc.
..remarks:Currently @Spec.IndexEsa@ is default if $TText$ is a @Class.String@.
..include:seqan/index.h
*/
    template < typename TObject >
    struct DefaultIndexSpec {
        typedef IndexEsa<> Type;
    };

//////////////////////////////////////////////////////////////////////////////
/**
.Metafunction.DefaultIndexStringSpec:
..cat:Index
..summary:Default @Class.String@ specialization type of the @Metafunction.Fibre@ of an @Class.Index@. 
..signature:DefaultIndexStringSpec<TIndex>::Type
..class:Class.Index
..param.TIndex:An @Class.Index@ Type.
..returns:If the underlying text is a @Class.String@ or a set of Strings (see @Class.StringSet@) the String's spec. type is returned.
..remarks:Most of the @Class.Index@ fibres are strings. The @Class.String@ specialization type is chosen by this meta-function.
..include:seqan/index.h
*/
    template < typename TIndex >
    struct DefaultIndexStringSpec {
        typedef Alloc<> Type;
    };

    template < typename TValue, typename TSpec >
    struct DefaultIndexStringSpec< String<TValue, External<TSpec> > > {
        typedef External<TSpec> Type;
    };

	template < typename TString, typename TSpec >
	struct DefaultIndexStringSpec< StringSet<TString, TSpec> >:
		DefaultIndexStringSpec<TString> {};


//////////////////////////////////////////////////////////////////////////////
/**
.Class.Index:
..summary:Contains preprocessing data of a fixed text. Allows fast dictionary look-up and advanced computations.
..cat:Index
..signature:Index<TText[, TSpec]>
..param.TText:The text type.
...type:Class.String
...type:Class.StringSet
...metafunction:Metafunction.Host
..param.TSpec:The index type.
...default:The result of @Metafunction.DefaultIndexSpec@
...metafunction:Metafunction.Spec
..remarks:An index contains various arrays or objects, also called fibres (see @Metafunction.Fibre@).
..remarks:These fibres are created on demand depending on the requirements of an algorithm.
..include:seqan/index.h
*/

///.Function.setHaystack.param.haystack.type:Class.Index

	// index as a haystack
	template < 
        typename TObject, 
        typename TSpec = typename DefaultIndexSpec<TObject>::Type > 
	class Index;

	template <typename TObject, typename TSpec>
	struct Host< Index<TObject, TSpec> > {
		typedef TObject Type;
	};

	template <typename TObject, typename TSpec>
	struct Spec< Index<TObject, TSpec> > {
		typedef TSpec Type;
	};

// TODO(singer): include @ sign before WaveletTree and WaveletTreeStructure when incoorporated.
//////////////////////////////////////////////////////////////////////////////
/**
.Metafunction.Fibre:
..summary:Type of a specific container member (fibre).
..signature:Fibre<TContainer, TSpec>::Type
..class:Class.Index
..cat:Index
..param.TContainer:The container type.
...type:Class.Index
..param.TSpec:Tag to specify the fibre.
..returns:Fibre type.
..remarks:Some containers, such as @Class.Index@, can be seen as a bundle consisting of various fibres. Because not every table is a fibre we did not call them tables, however, in many cases one can think of fibres as tables. The fibre interface was designed to unify the access to the members of the different fibres.
To get a reference or the type of a specific fibre use @Function.getFibre@ or @Metafunction.Fibre@.		
..remarks:A @Metafunction.Fibre@ does not need to be a real container. It can also be a view (see @Tag.ESA Index Fibres.EsaRawText@).
..include:seqan/index.h
*/

// In most cases this type is $String<Size<TIndex>::Type>$.

	// meta function to get the type of a bundle fibre
	template < typename TIndex, typename TSpec >
	struct Fibre {
		typedef String< typename Size<TIndex>::Type > Type;
	};

	template < typename TIndex, typename TSpec >
	struct Fibre<TIndex const, TSpec> {
		typedef typename Fibre<TIndex, TSpec>::Type const Type;
	};

	struct FibreRecord {
		unsigned	id;
		void*		ptr;
		bool		owner;
	};

	// less function to search in sorted list for fibre id
	struct FibreLess: public ::std::binary_function<FibreRecord, unsigned, bool>
	{	// functor for operator>
		inline bool operator()(FibreRecord const & _Left, unsigned const Right_) const
		{	// apply operator> to operands
			return (_Left.id < Right_);
		}
	};

//////////////////////////////////////////////////////////////////////////////
/**
.Metafunction.DefaultIndexCreator:
..cat:Index
..summary:Default algorithm to create a demanded and not yet existing @Metafunction.Fibre@.
..signature:DefaultIndexCreator<TIndex, TFibre>::Type
..class:Class.Index
..param.TIndex:An @Class.Index@ Type.
..param.TFibre:A tag specifying the fibre (e.g. @Tag.ESA Index Fibres.EsaSA@).
..returns:A tag specifying the default algorithm to create the fibre with.
..include:seqan/index.h
*/
    // standard algorithm for indices creation
    template < typename TIndex, typename TFibre >
	struct DefaultIndexCreator {
		typedef Default Type;
	};


//////////////////////////////////////////////////////////////////////////////

	template < 
		typename TSA,
		typename TText,
        typename TAlgSpec >
    struct SACreatorRandomAccess_
    {
        typedef typename AllowsFastRandomAccess<TSA>::Type   TRandomSA;
        typedef typename AllowsFastRandomAccess<TText>::Type TRandomText;
        typedef typename And<TRandomText,TRandomSA>::Type Type;
    };

	template < 
        typename TLCP,
		typename TText,
		typename TSA,
        typename TAlgSpec >
    struct LcpCreatorRandomAccess_
    {
        typedef typename AllowsFastRandomAccess<TText>::Type TRandomText;
        typedef typename AllowsFastRandomAccess<TLCP>::Type  TRandomLCP;
        typedef typename AllowsFastRandomAccess<TSA>::Type   TRandomSA;
        typedef typename And<TRandomLCP, typename And<TRandomText,TRandomSA>::Type>::Type Type;
    };


//////////////////////////////////////////////////////////////////////////////
/*
	.Class.Bundle:
	..summary:General purpose container of various members.
	..signature:Bundle<TValue, TSize>
	..param.TValue:The value type, that is the type of the items/characters stored in the string.
	...remarks:Use @Metafunction.Value@ to get the value type for a given class.
	..param.TSpec:The specializing type.
	...default:$Alloc<>$, see @Spec.Alloc String@.
*/
/*
	template < typename TSpec = void >
	truct Bundle {
		typedef ::std::vector<FibreRecord>	TFibreRecords;
		TFibreRecords						fibres;
	};

	template < typename TBundleSpec, typename TFibreSpec >
	inline FibreRecord& getRecord(Bundle<TBundleSpec> &bundle, TFibreSpec const) {
		unsigned id = (unsigned)ClassIdentifier_<TFibreSpec>::getID();

		typename Bundle<TBundleSpec>::TFibreRecords::iterator first = lower_bound(bundle.fibres.begin(), bundle.fibres.end(), id, FibreLess());
		if (!first->id != id) {
			FibreRecord rec;
			rec.id = id;
			rec.ptr = NULL;
			rec.owner = true;
			bundle.fibres.insert(first, rec);
		} else
			return *first;
	}

	template < typename TBundleSpec, typename TFibreSpec >
	inline typename Fibre<Bundle<TBundleSpec>, TFibreSpec>::Type & getFibre(Bundle<TBundleSpec> &bundle, TFibreSpec const) {
		typedef typename Fibre<Bundle<TBundleSpec>, TFibreSpec>::Type Type;
		unsigned id = (unsigned)ClassIdentifier_<TFibreSpec>::getID();

		FibreRecord &rec = getRecord(bundle, TFibreSpec());
		if (!rec.ptr)
			rec.ptr = new Type();
		return *reinterpret_cast<Type*>(rec.ptr);
	}

	template < typename TBundleSpec, typename TFibreSpec >
	inline typename Fibre<Bundle<TBundleSpec>, TFibreSpec>::Type const & getFibre(Bundle<TBundleSpec> const &bundle, TFibreSpec const) {
		typedef typename Fibre<Bundle<TBundleSpec>, TFibreSpec>::Type Type;
		unsigned id = (unsigned)ClassIdentifier_<TFibreSpec>::getID();

		FibreRecord &rec = getRecord(bundle, TFibreSpec());
		return *reinterpret_cast<Type*>(rec.ptr);
	}
*/

//////////////////////////////////////////////////////////////////////////////
// various fibre specs for enhanced suffix arrays

	struct FibreText_;		// Original text. Can be a String or a StringSet
	struct FibreRawText_;	// Concatenation of the strings above
	struct FibreSA_;		// suffix array (of raw text with virtual $-delimiters) with Pair entries
	struct FibreRawSA_;	// suffix array with integer entries
	struct FibreSae_;		// suffix array reordered in a b-tree
	struct FibreLcp_;		// lcp table of raw text
	struct FibreLcpe_;		// lcp interval tree
	struct FibreChildtab_;	// childtab (Kurtz et al.) of raw text
	struct FibreBwt_;		// burrows wheeler table of raw text

	typedef Tag<FibreText_> const		FibreText;
	typedef Tag<FibreRawText_> const	FibreRawText;
	typedef Tag<FibreSA_> const         FibreSA;
	typedef Tag<FibreRawSA_> const		FibreRawSA;
	typedef Tag<FibreSae_> const		FibreSae;
	typedef Tag<FibreLcp_> const		FibreLcp;
	typedef Tag<FibreLcpe_> const		FibreLcpe;
	typedef Tag<FibreChildtab_> const	FibreChildtab;
	typedef Tag<FibreBwt_> const		FibreBwt;

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.SAValue:
..cat:Index
..summary:The default alphabet type of a suffix array, i.e. the type to store a position of a string or string set.
..signature:SAValue<TObject>::Type
..class:Class.Index
..param.TObject:A string, string set, or index type.
...type:Class.String
...type:Class.StringSet
...type:Class.Index
..returns:A type to store a position. 
...text:If $TObject$ is a @Class.String@, it is a single integer value. By default this is the @Metafunction.Size@ type of $TObject$.
...text:If $TObject$ is a @Class.StringSet@, it could be a single integer too (called global position, see @Spec.ConcatDirect@) or a @Class.Pair@ (called local position, see @Spec.Owner@).
Currently SeqAn defaults to a local position for @Class.StringSet@ classes (index_base.h).
..example.code:template < typename TString, typename TSpec >
struct SAValue< StringSet<TString, TSpec> > {
	typedef Pair<
		typename Size< StringSet<TString, TSpec> >::Type,
		typename SAValue<TString>::Type,
		Pack
	> Type;
};
..remarks.note:SAValue is the return type of various function, e.g. @Function.position@ for the @Class.Index@ @Class.Finder@ class, @Function.getOccurrence@, @Function.getOccurrences@ etc.
You should always use the type of this meta-function to store the return values.
If you want to write algorithms for both variants (local and global positions) you 
should use the functions @Function.posLocalize@, @Function.posGlobalize@, @Function.getSeqNo@ and @Function.getSeqOffset@.
..remarks.note:If $TObject$ is an @Class.Index@, @Metafunction.Position@ returns the same value as $SAValue$. You can change the position type of an index by overloading $SAValue$, not @Metafunction.Position@.
..include:seqan/index.h
*/
	template <typename TObject>
	struct SAValue:
		Size<TObject> {};
	
	template <typename TObject>
	struct SAValue<TObject const>:
		SAValue<TObject> {};
	
	// to speed up sequence number computation
	// we use a pair of seqNo and localPosition
	template < typename TString, typename TSpec >
	struct SAValue< StringSet<TString, TSpec> > {
		typedef Pair<
			typename Size< StringSet<TString, TSpec> >::Type,
			typename SAValue<TString>::Type,
			Pack
		> Type;
	};

/*
	template < typename TString, typename TSpec >
	struct SAValue< StringSet<TString, TSpec> > {
		typedef Pair<
			typename Size< StringSet<TString, TSpec> >::Type,
			typename SAValue<TString>::Type,
			BitPacked<2,30>						    // max. 4 sequences 
		> Type;										// max. 2^30 characters each
	};
*/
	template < typename TText, typename TSpec >
	struct SAValue< Index<TText, TSpec> >:
		SAValue<TText> {};

	template < typename TObject, typename TSpec >
	struct DefaultIndexStringSpec< Index<TObject, TSpec> >:
		DefaultIndexStringSpec<TObject> {};

//////////////////////////////////////////////////////////////////////////////
// value and size type of an index

	template < typename TText, typename TSpec >
    struct Value< Index<TText, TSpec> > {
		typedef typename Value<
			typename Fibre< Index<TText, TSpec>, FibreRawText>::Type 
		>::Type Type;
    };

	template < typename TText, typename TSpec >
    struct Size< Index<TText, TSpec> > {
		typedef typename Size<
			typename Fibre< Index<TText, TSpec>, FibreRawText>::Type 
		>::Type Type;
    };

	template < typename TText, typename TSpec >
	struct Position< Index<TText, TSpec> >:
		SAValue< Index<TText, TSpec> > {};

//////////////////////////////////////////////////////////////////////////////
// infix of an index

	template < typename TText, typename TSpec >
    struct Infix< Index<TText, TSpec> >:
		public Infix<TText> {};

	template < typename TText, typename TSpec >
    struct Infix< Index<TText, TSpec> const >:
		public Infix<TText> {};

//////////////////////////////////////////////////////////////////////////////
// default table type

	template < typename TObject, typename TSpec, typename TFibre >
	struct Fibre< Index<TObject, TSpec>, Tag<TFibre> const > {
		typedef String< 
			typename Size< Index<TObject, TSpec> >::Type,
			typename DefaultIndexStringSpec< Index<TObject, TSpec> >::Type 
		> Type;
	};

//////////////////////////////////////////////////////////////////////////////
// original text

	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, TSpec>, FibreText> {
		typedef TText Type;
	};

//////////////////////////////////////////////////////////////////////////////
// concatenated text

	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, TSpec>, FibreRawText> {
		typedef typename Concatenator<TText>::Type Type;
	};

//////////////////////////////////////////////////////////////////////////////
// suffix array type

	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, TSpec>, FibreSA> {
		typedef String<
			typename SAValue< Index<TText, TSpec> >::Type,
			typename DefaultIndexStringSpec< Index<TText, TSpec> >::Type 
		> Type;
	};

//////////////////////////////////////////////////////////////////////////////
// globalize functor

	template <typename InType, typename TLimitsString, typename Result = typename Value<TLimitsString>::Type>
	struct FunctorGlobalize : public ::std::unary_function<InType,Result>
	{
		TLimitsString const *limits;
		FunctorGlobalize() {}
		FunctorGlobalize(TLimitsString const &_limits) : limits(&_limits) {}

		inline Result operator()(InType const &x) const
		{
			return posGlobalize(x, *limits);
		}
    };

	template <typename InType, typename Result>
	struct FunctorGlobalize<InType, Nothing, Result> : public ::std::unary_function<InType,InType>
	{
		FunctorGlobalize() {}
		FunctorGlobalize(Nothing const &) {}

        inline InType operator()(InType const &x) const
        {
			return x;
		}
    };

//////////////////////////////////////////////////////////////////////////////
// raw suffix array contains integer offsets relative to raw text
/*
	template < typename TString, typename TSpec >
	struct Fibre< Index<TString, TSpec>, FibreRawSA>:
		public Fibre< Index<TString, TSpec> const, FibreSA> {};

	template < typename TString, typename TSSetSpec, typename TSpec >
	struct Fibre< Index<StringSet<TString, TSSetSpec>, TSpec>, FibreRawSA> 
	{
		typedef Index< StringSet<TString, TSSetSpec>, TSpec> TIndex;
		typedef ModifiedString<
			typename Fibre<TIndex, FibreSA>::Type,
			ModView< FunctorGlobalize< 
				typename Value< typename Fibre<TIndex, FibreSA>::Type >::Type,
				typename StringSetLimits<TString>::Type >
			>
		> Type;
	};
*/
	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, TSpec>, FibreRawSA> 
	{
		typedef Index<TText, TSpec> TIndex;
		typedef ModifiedString<
			typename Fibre<TIndex, FibreSA>::Type,
			ModView< FunctorGlobalize< 
				typename Value< typename Fibre<TIndex, FibreSA>::Type >::Type,
				typename StringSetLimits<TText>::Type >
			>
		> Type;
	};

//////////////////////////////////////////////////////////////////////////////
// default burrows-wheeler table

	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, TSpec>, FibreBwt> {
		typedef String <
			typename Value< Index<TText, TSpec> >::Type,
			typename DefaultIndexStringSpec< Index<TText, TSpec> >::Type
		> Type;
	};


//////////////////////////////////////////////////////////////////////////////
// default fibre creators

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, TSpec>, FibreSA> {
        typedef Skew7 Type;							// standard suffix array creator is skew7
    };

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, TSpec>, FibreLcp> {
        typedef Kasai Type;
    };

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, TSpec>, FibreBwt> {
        typedef Bwt Type;
    };

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, TSpec>, FibreChildtab> {
        typedef Childtab Type;
    };


//////////////////////////////////////////////////////////////////////////////
// fibre interface to access the enhanced suffix array tables

/**
.Function.getFibre:
..summary:Returns a specific fibre of a container.
..signature:getFibre(container, fibreTag)
..class:Class.Index
..cat:Index
..param.container:The container holding the fibre.
...type:Class.Index
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.ESA Index Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
..example.code:
Index< String<char> > index_esa("tobeornottobe");

String<char> & text = getFibre(indexEsa, EsaText());
*/

	template <typename TText, typename TSpec>
	inline Holder<TText> & _dataHost(Index<TText, TSpec> &index) {
		return index.text;
	}
	template <typename TText, typename TSpec>
	inline Holder<TText> const & _dataHost(Index<TText, TSpec> const &index) {
		return index.text;
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreText>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreText) {
		return value(index.text);
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreText>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreText) {
		return value(index.text);
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreRawText>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreRawText) {
		return concat(value(index.text));
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreRawText>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreRawText) {
		return concat(value(index.text));
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreSA>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreSA) {
		return index.sa;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreSA>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreSA) {
		return index.sa;
	}

//////////////////////////////////////////////////////////////////////////////
/*
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreRawSA>::Type const & 
	getFibre(Index<TText, TSpec> &index, FibreRawSA) {
		return indexSA(index);
	}

	template <typename TString, typename TSSetSpec, typename TSpec>
	inline typename Fibre<Index<StringSet<TString, TSSetSpec>, TSpec>, FibreRawSA>::Type
	getFibre(Index<StringSet<TString, TSSetSpec>, TSpec> &index, FibreRawSA) 
	{
		typedef Index< StringSet<TString, TSSetSpec>, TSpec> TIndex;
		
		typedef FunctorGlobalize<
			typename Value< typename Fibre<TIndex, FibreSA>::Type >::Type,
			typename StringSetLimits<StringSet<TString, TSSetSpec> >::Type
		> TFunctor;
		
		typedef ModifiedString<
			typename Fibre<Index<StringSet<TString, TSSetSpec>, TSpec>, FibreSA>::Type,
			ModView< TFunctor >
		> ModString;

		return ModString(indexSA(index), TFunctor(stringSetLimits(indexText(index))));
	}
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreRawSA>::Type
	getFibre(Index<TText, TSpec> &index, FibreRawSA) 
	{
		typedef Index<TText, TSpec> TIndex;
		
		typedef FunctorGlobalize<
			typename Value< typename Fibre<TIndex, FibreSA>::Type >::Type,
			typename StringSetLimits<TText>::Type
		> TFunctor;
		
		typedef ModifiedString<
			typename Fibre<Index<TText, TSpec>, FibreSA>::Type,
			ModView< TFunctor >
		> ModString;

		return ModString(indexSA(index), TFunctor(stringSetLimits(indexText(index))));
	}
//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreLcp>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreLcp) {
		return index.lcp;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreLcp>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreLcp) {
		return index.lcp;
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreLcpe>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreLcpe) {
		return index.lcpe;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreLcpe>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreLcpe) {
		return index.lcpe;
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreChildtab>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreChildtab) {
		return index.childtab;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreChildtab>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreChildtab) {
		return index.childtab;
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreBwt>::Type & 
	getFibre(Index<TText, TSpec> &index, FibreBwt) {
		return index.bwt;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreBwt>::Type & 
	getFibre(Index<TText, TSpec> const &index, FibreBwt) {
		return index.bwt;
	}

//////////////////////////////////////////////////////////////////////////////

///.Function.length.param.object.type:Class.Index
///.Function.length.remarks:If $object$ is of type @Class.Index@, the number of characters in the raw underlying text of the index is returned.

	template <typename TText, typename TSpec>
	inline typename Size<Index<TText, TSpec> >::Type 
	length(Index<TText, TSpec> const &index) {
		return length(indexRawText(index));
	}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.countSequences
..cat:Index
..summary:Return the number of sequences in an index' underlying text.
..signature:countSequences(index)
..class:Class.Index
..param.index:The index to return the number of sequences of.
...type:Class.Index
..returns:The number of sequences in the index' underlying text.
...metafunction:Metafunction.Size
..include:seqan/index.h
 */

	template <typename TText, typename TSpec>
	inline typename Size<TText>::Type 
	countSequences(Index<TText, TSpec> const &index) {
		return countSequences(indexText(index));
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	struct GetSequenceByNo< Index<TText, TSpec> >
	{
		typedef typename GetSequenceByNo<TText>::Type Type;
	};

	template <typename TText, typename TSpec>
	struct GetSequenceByNo< Index<TText, TSpec> const>
	{
		typedef typename GetSequenceByNo<TText const>::Type Type;
	};

//////////////////////////////////////////////////////////////////////////////

	template <typename TSeqNo, typename TText, typename TSpec>
	inline typename GetSequenceByNo< Index<TText, TSpec> >::Type
	getSequenceByNo(TSeqNo seqNo, Index<TText, TSpec> &index)
	{
		return getSequenceByNo(seqNo, indexText(index));
	}

	template <typename TSeqNo, typename TText, typename TSpec>
	inline typename GetSequenceByNo< Index<TText, TSpec> const>::Type
	getSequenceByNo(TSeqNo seqNo, Index<TText, TSpec> const &index)
	{
		return getSequenceByNo(seqNo, indexText(index));
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TSeqNo, typename TText, typename TSpec>
	inline typename Size<Index<TText, TSpec> >::Type 
	sequenceLength(TSeqNo seqNo, Index<TText, TSpec> const &index) {
		return sequenceLength(seqNo, indexText(index));
	}

//////////////////////////////////////////////////////////////////////////////

	template <typename TPos, typename TText, typename TSpec>
	inline typename Size<Index<TText, TSpec> >::Type 
	suffixLength(TPos pos, Index<TText, TSpec> const &index)
    {
		return length(indexText(index)) - pos;
	}

	template <typename TPos, typename TString, typename TSSetSpec, typename TSpec>
	inline typename Size<Index<StringSet<TString, TSSetSpec>, TSpec> >::Type 
	suffixLength(TPos pos, Index<StringSet<TString, TSSetSpec>, TSpec> const &index)
    {
        typename StringSetLimits<StringSet<TString, TSSetSpec> >::Type const &limits = stringSetLimits(index);
		return sequenceLength(getSeqNo(pos, limits), index) - getSeqOffset(pos, limits);
	}


//////////////////////////////////////////////////////////////////////////////
/**
.Function.textAt:
..summary:Shortcut for $value(indexText(..), ..)$.
..cat:Index
..signature:textAt(position, index)
..class:Class.Index
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
*/

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, FibreRawText>::Type>::Type 
	textAt(TPos i, TIndex &index) {
		return value(getFibre(index, FibreRawText()), i);
	}
	template <typename TPos, typename TString, typename TSSetSpec, typename TSpec>
	inline typename Reference<typename Fibre< Index< StringSet<TString, TSSetSpec>, TSpec>, FibreRawText>::Type>::Type 
	textAt(TPos i, Index< StringSet<TString, TSSetSpec>, TSpec> &index) {
		return value(getFibre(index, FibreRawText()), posGlobalize(i, stringSetLimits(index)));
	}
	template <typename TPos, typename TString, typename TSSetSpec, typename TSpec>
	inline typename Reference<typename Fibre< Index< StringSet<TString, TSSetSpec>, TSpec> const, FibreRawText>::Type>::Type 
	textAt(TPos i, Index< StringSet<TString, TSSetSpec>, TSpec> const &index) {
		return value(getFibre(index, FibreRawText()), posGlobalize(i, stringSetLimits(index)));
	}
	template <typename TPos, typename TString, typename TSpec>
	inline typename Reference<typename Fibre< Index< StringSet<TString, Owner<Default> >, TSpec>, FibreRawText>::Type>::Type 
	textAt(TPos i, Index< StringSet<TString, Owner<Default> >, TSpec> &index) {
		Pair <
			typename Size< StringSet<TString, Owner<Default> > >::Type,
			typename Size< TString >::Type > locPos;
		posLocalize(locPos, i, stringSetLimits(index));
		return value(value(getFibre(index, FibreText()), getValueI1(locPos)), getValueI2(locPos));
	}
	template <typename TPos, typename TString, typename TSpec>
	inline typename Reference<typename Fibre< Index< StringSet<TString, Owner<Default> >, TSpec> const, FibreRawText>::Type>::Type 
	textAt(TPos i, Index< StringSet<TString, Owner<Default> >, TSpec> const &index) {
		Pair <
        typename Size< StringSet<TString, Owner<Default> > >::Type,
        typename Size< TString >::Type > locPos;
		posLocalize(locPos, i, stringSetLimits(index));
		return value(value(getFibre(index, FibreText()), getValueI1(locPos)), getValueI2(locPos));
	}

//////////////////////////////////////////////////////////////////////////////
// infix

	template <typename TText, typename TSpec, typename TPosBegin, typename TPosEnd>
	inline typename Infix<TText>::Type
	infix(Index<TText, TSpec> &index, TPosBegin pos_begin, TPosEnd pos_end)
	{
		return infix(indexText(index), pos_begin, pos_end);
	}

	template <typename TText, typename TSpec, typename TPosBegin, typename TPosEnd>
	inline typename Infix<TText>::Type
	infix(Index<TText, TSpec> const &index, TPosBegin pos_begin, TPosEnd pos_end)
	{
		return infix(indexText(index), pos_begin, pos_end);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.rawtextAt:
..summary:Shortcut for $value(indexRawText(..), ..)$.
..cat:Index
..signature:rawtextAt(position, index)
..class:Class.Index
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
*/

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, FibreRawText>::Type>::Type rawtextAt(TPos i, TIndex &index) {
		return value(getFibre(index, FibreRawText()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, FibreRawText>::Type>::Type rawtextAt(TPos i, TIndex const &index) {
		return value(getFibre(index, FibreRawText()), i);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.saAt:
..summary:Shortcut for $value(indexSA(..), ..)$.
..cat:Index
..signature:saAt(position, index)
..class:Class.Index
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
*/

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, FibreSA>::Type>::Type saAt(TPos i, TIndex &index) {
		return value(getFibre(index, FibreSA()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, FibreSA>::Type>::Type saAt(TPos i, TIndex const &index) {
		return value(getFibre(index, FibreSA()), i);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.rawsaAt:
..summary:Shortcut for $value(indexRawSA(..), ..)$.
..cat:Index
..signature:saAt(position, index)
..class:Class.Index
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
*/

	template <typename TPos, typename TIndex>
	inline typename Value<typename Fibre<TIndex const, FibreRawSA>::Type>::Type rawsaAt(TPos i, TIndex const &index) {
		return posGlobalize(saAt(i, index), stringSetLimits(indexText(index)));
	}


//////////////////////////////////////////////////////////////////////////////
/**
.Function.lcpAt:
..summary:Shortcut for $value(indexLcp(..), ..)$.
..cat:Index
..signature:lcpAt(position, index)
..class:Class.Index
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
*/

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, FibreLcp>::Type>::Type lcpAt(TPos i, TIndex &index) {
		return value(getFibre(index, FibreLcp()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, FibreLcp>::Type>::Type lcpAt(TPos i, TIndex const &index) {
		return value(getFibre(index, FibreLcp()), i);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.lcpeAt:
..summary:Shortcut for $value(indexLcpe(..), ..)$.
..cat:Index
..signature:lcpeAt(position, index)
..class:Class.Index
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
*/

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, FibreLcpe>::Type>::Type lcpeAt(TPos i, TIndex &index) {
		return value(getFibre(index, FibreLcpe()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, FibreLcpe>::Type>::Type lcpeAt(TPos i, TIndex const &index) {
		return value(getFibre(index, FibreLcpe()), i);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.childAt:
..summary:Shortcut for $value(indexChildtab(..), ..)$.
..cat:Index
..signature:childAt(position, index)
..class:Class.Index
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
*/

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, FibreChildtab>::Type>::Type childAt(TPos i, TIndex &index) {
		return value(getFibre(index, FibreChildtab()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, FibreChildtab>::Type>::Type childAt(TPos i, TIndex const &index) {
		return value(getFibre(index, FibreChildtab()), i);
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.bwtAt:
..summary:Shortcut for $value(indexBwt(..), ..)$.
..cat:Index
..signature:bwtAt(position, index)
..class:Class.Index
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference or proxy to the value.
..include:seqan/index.h
*/

	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, FibreBwt>::Type>::Type bwtAt(TPos i, TIndex &index) {
		return value(getFibre(index, FibreBwt()), i);
	}
	template <typename TPos, typename TIndex>
	inline typename Reference<typename Fibre<TIndex const, FibreBwt>::Type>::Type bwtAt(TPos i, TIndex const &index) {
		return value(getFibre(index, FibreBwt()), i);
	}

//////////////////////////////////////////////////////////////////////////////

    template <typename TIndex, typename TPos, typename TSize>
	inline typename SAValue<TIndex>::Type toSuffixPosition(TIndex &, TPos i, TSize) {
        return i;
	}
    template <typename TIndex, typename TPos, typename TSize>
	inline typename SAValue<TIndex const>::Type toSuffixPosition(TIndex const &, TPos i, TSize) {
        return i;
	}

//////////////////////////////////////////////////////////////////////////////
// interface for infinity/invalid values

	template <typename TValue>
	inline void _setSizeInval(TValue &v) {
		v = MaxValue<TValue>::VALUE;
	}

	template <typename TValue>
	inline bool _isSizeInval(TValue const &v) {
//IOREV _notio_
		return v == MaxValue<TValue>::VALUE;
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexText:
..summary:Shortcut for $getFibre(.., EsaText)$.
..cat:Index
..signature:indexText(index)
..class:Class.Index
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference to the @Tag.ESA Index Fibres.EsaText@ fibre (original text).
..include:seqan/index.h
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreText>::Type & indexText(Index<TText, TSpec> &index) { return getFibre(index, FibreText()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreText>::Type & indexText(Index<TText, TSpec> const &index) { return getFibre(index, FibreText()); }

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename StringSetLimits<TText const>::Type
	stringSetLimits(Index<TText, TSpec> &) { 
		return Nothing(); 
	}

	template <typename TText, typename TSpec>
	inline typename StringSetLimits<TText const>::Type
	stringSetLimits(Index<TText, TSpec> const &) { 
		return Nothing(); 
	}

	template <typename TString, typename TSSetSpec, typename TSpec>
	inline typename StringSetLimits< StringSet<TString, TSSetSpec> const >::Type & 
	stringSetLimits(Index<StringSet<TString, TSSetSpec>, TSpec> &index) {
		return stringSetLimits(indexText(index)); 
	}

	template <typename TString, typename TSSetSpec, typename TSpec>
	inline typename StringSetLimits< StringSet<TString, TSSetSpec> const >::Type & 
	stringSetLimits(Index<StringSet<TString, TSSetSpec>, TSpec> const &index) {
		return stringSetLimits(indexText(index)); 
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexRawText:
..summary:Shortcut for $getFibre(.., EsaRawText)$.
..cat:Index
..signature:indexRawText(index)
..class:Class.Index
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference to the @Tag.ESA Index Fibres.EsaRawText@ fibre (concatenated input text).
..include:seqan/index.h
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreRawText>::Type & indexRawText(Index<TText, TSpec> &index) { return getFibre(index, FibreRawText()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreRawText>::Type & indexRawText(Index<TText, TSpec> const &index) { return getFibre(index, FibreRawText()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexSA:
..summary:Shortcut for $getFibre(.., EsaSA)$.
..cat:Index
..signature:indexSA(index)
..class:Class.Index
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference to the @Tag.ESA Index Fibres.EsaSA@ fibre (suffix array).
..include:seqan/index.h
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreSA>::Type & indexSA(Index<TText, TSpec> &index) { return getFibre(index, FibreSA()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreSA>::Type & indexSA(Index<TText, TSpec> const &index) { return getFibre(index, FibreSA()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexRawSA:
..summary:Shortcut for $getFibre(.., EsaRawSA)$.
..cat:Index
..signature:indexRawSA(index)
..class:Class.Index
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference to the @Tag.ESA Index Fibres.EsaRawSA@ fibre (suffix array).
..include:seqan/index.h
*/
/*
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreRawSA>::Type const & indexRawSA(Index<TText, TSpec> &index) { return getFibre(index, FibreRawSA()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreRawSA>::Type const & indexRawSA(Index<TText, TSpec> const &index) { return getFibre(index, FibreRawSA()); }
*/
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreRawSA>::Type indexRawSA(Index<TText, TSpec> &index) { return getFibre(index, FibreRawSA()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreRawSA>::Type indexRawSA(Index<TText, TSpec> const &index) { return getFibre(index, FibreRawSA()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexLcp:
..summary:Shortcut for $getFibre(.., EsaLcp)$.
..cat:Index
..signature:indexLcp(index)
..class:Class.Index
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference to the @Tag.ESA Index Fibres.EsaLcp@ fibre (lcp table).
..include:seqan/index.h
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreLcp>::Type & indexLcp(Index<TText, TSpec> &index) { return getFibre(index, FibreLcp()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreLcp>::Type & indexLcp(Index<TText, TSpec> const &index) { return getFibre(index, FibreLcp()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexLcpe:
..summary:Shortcut for $getFibre(.., EsaLcpe)$.
..cat:Index
..signature:indexLcpe(index)
..class:Class.Index
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference to the @Tag.ESA Index Fibres.EsaLcpe@ fibre (enhanced lcp table).
..include:seqan/index.h
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreLcpe>::Type & indexLcpe(Index<TText, TSpec> &index) { return getFibre(index, FibreLcpe()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreLcpe>::Type & indexLcpe(Index<TText, TSpec> const &index) { return getFibre(index, FibreLcpe()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexBwt:
..summary:Shortcut for $getFibre(.., EsaBwt)$.
..cat:Index
..signature:indexBwt(index)
..class:Class.Index
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference to the @Tag.ESA Index Fibres.EsaBwt@ fibre (Burrows-Wheeler table).
..include:seqan/index.h
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreBwt>::Type & indexBwt(Index<TText, TSpec> &index) { return getFibre(index, FibreBwt()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreBwt>::Type & indexBwt(Index<TText, TSpec> const &index) { return getFibre(index, FibreBwt()); }

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexChildtab:
..summary:Shortcut for $getFibre(.., EsaChildtab)$.
..cat:Index
..signature:indexChildtab(index)
..class:Class.Index
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.IndexEsa
..returns:A reference to the @Tag.ESA Index Fibres.EsaChildtab@ fibre (child table).
..include:seqan/index.h
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, FibreChildtab>::Type & indexChildtab(Index<TText, TSpec> &index) { return getFibre(index, FibreChildtab()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, FibreChildtab>::Type & indexChildtab(Index<TText, TSpec> const &index) { return getFibre(index, FibreChildtab()); }

}

#endif

