// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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

#ifndef SEQAN_HEADER_FILE_PAGE_RAID0_H
#define SEQAN_HEADER_FILE_PAGE_RAID0_H

/* IOREV
 * _nottested_
 * _nodoc_
 *
 * not tested by any test or app
 * no documentation for the functions
 *
 * hard to say how/if this works, since there is no doc
 * 
 */


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{


    //////////////////////////////////////////////////////////////////////////////
    // page based read/write for striped files

	template < typename TValue, unsigned FILE_COUNT, typename TFile, typename TSpec >
	inline bool 
	readPage(
		int pageNo, 
		PageFrame<TValue, File< Striped<FILE_COUNT, TFile> >, TSpec> &pf, 
		File< Striped<FILE_COUNT, TFile> > &file)
	{
//IOREV _nodoc_
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "readPage:  " << ::std::hex << (TValue*)pf.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << ::std::endl;
		#endif
		pf.dirty = false;
		pf.status = pf.READING;
		return asyncReadAt(
			file[pageNo % FILE_COUNT], 
			(TValue*)pf.begin, 
			size(pf), 
			(pos_t)(pageNo / FILE_COUNT) * (pos_t)pageSize(pf), 
			pf.request);
	}

	template < typename TValue, unsigned FILE_COUNT, typename TFile, typename TSpec >
	inline bool 
	writePage(
		PageFrame<TValue, File< Striped<FILE_COUNT, TFile> >, TSpec> &pf, 
		int pageNo, 
		File< Striped<FILE_COUNT, TFile> > &file)
	{
//IOREV _nodoc_ different signature from read function
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "writePage: " << ::std::hex << (TValue*)pf.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << ::std::endl;
		#endif
		pf.status = pf.WRITING;
		return asyncWriteAt(
			file[pageNo % FILE_COUNT], 
			(TValue*)pf.begin, 
			size(pf), 
			(pos_t)(pageNo / FILE_COUNT) * (pos_t)pageSize(pf), 
			pf.request);
	}

	template < typename TValue, unsigned FILE_COUNT, typename TFile, typename TSpec, typename TSize >
	inline bool 
	readLastPage(
		int pageNo, 
		PageFrame<TValue, File< Striped<FILE_COUNT, TFile> >, TSpec> &pf, 
		File< Striped<FILE_COUNT, TFile> > &file,
		TSize size)
	{
//IOREV _nodoc_
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "readPage:  " << ::std::hex << (TValue*)pf.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << " size " << size << ::std::endl;
		#endif
		pf.dirty = false;
		pf.status = pf.READY;
		return readAt(
			file[pageNo % FILE_COUNT], 
			(TValue*)pf.begin, 
			size, 
			(pos_t)(pageNo / FILE_COUNT) * (pos_t)pageSize(pf));
	}

	template < typename TValue, unsigned FILE_COUNT, typename TFile, typename TSpec, typename TSize >
	inline bool 
	writeLastPage(
		PageFrame<TValue, File< Striped<FILE_COUNT, TFile> >, TSpec> &pf, 
		int pageNo, 
		File< Striped<FILE_COUNT, TFile> > &file,
		TSize size)
	{
//IOREV _nodoc_ different signature from read function
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "writePage: " << ::std::hex << (TValue*)pf.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << " size " << size << ::std::endl;
		#endif
		pf.dirty = false;
		pf.status = pf.READY;
//        resize(pf, size);
		return writeAt(
			file[pageNo % FILE_COUNT], 
			(TValue*)pf.begin, 
			size, 
			(pos_t)(pageNo / FILE_COUNT) * (pos_t)pageSize(pf));
	}


	//////////////////////////////////////////////////////////////////////////////
	// bucket based read/write methods for striped files

	template < typename TValue, unsigned FILE_COUNT, typename TFile >
	inline unsigned 
	readBucket(
		PageBucket<TValue> &b, 
		int pageNo, 
		unsigned pageSize, 
		unsigned dataSize, 
		File< Striped<FILE_COUNT, TFile> > &file) 
	{
//IOREV _nodoc_
		typedef typename Position<TFile>::Type pos_t;
        unsigned readSize = _min(dataSize - b.pageOfs, (unsigned)(b.end - b.begin));
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "readBucket:  " << ::std::hex << b.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << " at " << (pos_t)(pageNo / FILE_COUNT) * (pos_t)pageSize + b.pageOfs;
			::std::cerr << " size " << readSize << ::std::endl;
		#endif
        if (readSize && readAt(file[pageNo % FILE_COUNT], b.begin, readSize, (pos_t)(pageNo / FILE_COUNT) * (pos_t)pageSize + b.pageOfs)) {
            b.pageOfs += readSize;
            b.cur = b.begin;
            b.end = b.begin + readSize;
            return readSize;
        } else
            return 0;
	}

	template < typename TValue, unsigned FILE_COUNT, typename TFile >
	inline bool 
	writeBucket(
		PageBucket<TValue> &b,
		int pageNo, 
		unsigned pageSize, 
		File< Striped<FILE_COUNT, TFile> > &file) 
	{
//IOREV _nodoc_
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "writeBucket: " << ::std::hex << b.begin;
			::std::cerr << " from page " << ::std::dec << pageNo << " at " << (pos_t)(pageNo / FILE_COUNT) * (pos_t)pageSize + b.pageOfs;
			::std::cerr << " size " << b.cur - b.begin << ::std::endl;
		#endif
        if ((b.cur == b.begin) || writeAt(file[pageNo % FILE_COUNT], b.begin, b.cur - b.begin, (pos_t)(pageNo / FILE_COUNT) * (pos_t)pageSize + b.pageOfs)) {
            b.pageOfs += b.cur - b.begin;
            b.cur = b.begin;
            return true;
        } else
            return false;
	}

	template < typename TValue, unsigned FILE_COUNT, typename TFile, typename TSpec >
	inline bool 
	writeBucket(
		PageFrame<TValue, File< Striped<FILE_COUNT, TFile> >, Dynamic<TSpec> > &pf, 
		unsigned &pageOfs, 
		File< Striped<FILE_COUNT, TFile> > &file) 
	{
//IOREV _nodoc_
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			::std::cerr << "writeBucket: " << ::std::hex << pf.begin;
			::std::cerr << " from page " << ::std::dec << pf.pageNo << " at " << (pos_t)(pf.pageNo / FILE_COUNT) * (pos_t)pageSize(pf) + pageOfs;
			::std::cerr << " size " << size(pf) << ::std::endl;
		#endif
        if (pf.end == pf.begin) return true;
        if (asyncWriteAt(file[pf.pageNo % FILE_COUNT], pf.begin, size(pf), (pos_t)(pf.pageNo / FILE_COUNT) * (pos_t)pageSize(pf) + pageOfs, pf.request)) {
            pf.status = pf.WRITING;
            pageOfs += size(pf);
            return true;
        } else
            return false;
	}

}

#endif
