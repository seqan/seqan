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

#ifndef SEQAN_HEADER_FILE_CGVIZ_H
#define SEQAN_HEADER_FILE_CGVIZ_H

/* IOREV
 * _tested_
 * _nodoc_
 *
 * tested in tests/file/test_file.h
 * tag mentionen in doc, but no further documentation, no link to spec
 * 
 */


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File Formats - CGViz
//////////////////////////////////////////////////////////////////////////////


/**
.Tag.File Format.tag.CGViz:
	CGViz file format for sequences. Only output.
..include:seqan/file.h
*/
struct TagCGViz_;
//IOREV
typedef Tag<TagCGViz_> const CGViz; //IOREV

/////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
// goNext
//////////////////////////////////////////////////////////////////////////////

template <typename TFile>
void goNext(TFile & file, CGViz) {
//IOREV _notinlined_ purpose not clear to me
	SEQAN_CHECKPOINT;
    (void) file; // When compiled without assertions.
	SEQAN_ASSERT_NOT(_streamEOF(file));
	
	return;
}



//////////////////////////////////////////////////////////////////////////////
// write
//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TStringContainer, typename TSource, typename TSpec>
void _writeImpl(TFile& target, Align<TSource, TSpec>& align, TStringContainer& ids, CGViz) {
//IOREV _batchreading_ _notinlined_ actually writing not reading
	SEQAN_CHECKPOINT

	typedef Align<TSource, TSpec> const TAlign;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Position<typename Rows<TAlign>::Type>::Type TRowsPosition;
	typedef typename Position<TAlign>::Type TPosition;
	TRowsPosition row_count = length(rows(align));
	if (row_count < 2) return;

	unsigned int pair=1;
	unsigned int count=0;
	for(TRowsPosition i=0;i<row_count-1;++i) {
		for(TRowsPosition j=i+1;j<row_count;++j) {
			
			// Print header
			_streamWrite(target, "{DATA dat"); 
			_streamPutInt(target, pair);
			_streamPut(target, '\n');
			_streamWrite(target, "[__GLOBAL__] dimension=2:\n"); 
						
			TPosition begin_ = beginPosition(cols(align));
			TPosition end_ = endPosition(cols(align));
		
			bool match = false;
			while(begin_ < end_) {
				if ((row(align, i)[begin_]==row(align, j)[begin_]) && (row(align, i)[begin_]!='-')) {
					if (!match) {
						match=true;
						_streamPutInt(target, toSourcePosition(row(align,i),begin_+1));
						_streamPut(target, ' ');
						_streamPutInt(target, toSourcePosition(row(align,j),begin_+1));
						_streamPut(target, ' ');
					}
				}
				if ((row(align, i)[begin_]!=row(align, j)[begin_]) || (row(align, i)[begin_]=='-') || (row(align, j)[begin_]=='-')) {
					if (match) {
						_streamPutInt(target, toSourcePosition(row(align,i),begin_));
						_streamPut(target, ' ');
						_streamPutInt(target, toSourcePosition(row(align,j),begin_));
						_streamPut(target, '\n');
						match=false;
					}
				}
				begin_++;
			}
			if (match) {
				_streamPutInt(target, toSourcePosition(row(align,i),begin_));
				_streamPut(target, ' ');
				_streamPutInt(target, toSourcePosition(row(align,j),begin_));
				_streamPut(target, '\n');
				match=false;
			}
			_streamPut(target, '}');
			_streamPut(target, '\n');

			// Write footer
			_streamWrite(target, "{GLYPH Glyph");
			_streamPutInt(target, pair);
			_streamPut(target, '\n');
			_streamWrite(target, "drawerName=Lines\n");
			_streamWrite(target, "lineWidth=3\n");
			_streamPut(target, '}');
			_streamPut(target, '\n');
			_streamWrite(target, "{PANE Pane");
			_streamPutInt(target, pair);
			_streamPut(target, '\n');
			_streamWrite(target, "uLabel=");
			_streamWrite(target, getValue(ids,i));
			_streamPut(target, '\n');
			_streamWrite(target, "uStop=");
			_streamPutInt(target, length(source(row(align,i))));
			_streamPut(target, '\n');
			_streamWrite(target, "vLabel=");
			_streamWrite(target, getValue(ids,j));
			_streamPut(target, '\n');
			_streamWrite(target, "vStop=");
			_streamPutInt(target, length(source(row(align,j))));
			_streamPut(target, '\n');
			_streamPut(target, '}');
			_streamPut(target, '\n');
			_streamWrite(target, "{WINDOW Window");
			_streamPutInt(target, pair);
			_streamPut(target, '\n');
			_streamPut(target, '}');
			_streamPut(target, '\n');
			_streamWrite(target, "{FEEDER Feeder<");
			_streamPutInt(target, pair);
			_streamPut(target, '>');
			_streamPut(target, ' ');
			_streamPutInt(target, count);
			_streamPut(target, ' ');
			_streamPutInt(target, count+1);
			_streamPut(target, '\n');
			_streamPut(target, '}');
			_streamPut(target, '\n');
			++count;
			_streamWrite(target, "{THREADER Threader<");
			_streamPutInt(target, pair);
			_streamPut(target, '>');
			_streamPut(target, ' ');
			_streamPutInt(target, count);
			_streamPut(target, ' ');
			_streamPutInt(target, count+1);
			_streamPut(target, '\n');
			_streamPut(target, '}');
			_streamPut(target, '\n');
			++count;
			_streamWrite(target, "{ANCHOR Anchor<");
			_streamPutInt(target, pair);
			_streamPut(target, '>');
			_streamPut(target, ' ');
			_streamPutInt(target, count);
			_streamPut(target, ' ');
			_streamPutInt(target, count+1);
			_streamPut(target, '\n');
			_streamPut(target, '}');
			_streamPut(target, '\n');
			count+=2;
			++pair;
		}
	}
}


//____________________________________________________________________________

template <typename TFile, typename TSource, typename TSpec>
void write(TFile & file, Align<TSource, TSpec>& align, CGViz) {
//IOREV _notinlined_
	SEQAN_CHECKPOINT
	_writeImpl(file, align, String<String<char> >(), CGViz());
}

//____________________________________________________________________________

template <typename TFile, typename TStringContainer, typename TSource, typename TSpec>
void write(TFile & file, Align<TSource, TSpec> & align, TStringContainer& ids, CGViz) {
//IOREV _notinlined_
	SEQAN_CHECKPOINT
	_writeImpl(file, align, ids, CGViz());
}


//VisualC++ const array bug workaround
template <typename TFile, typename TStringContainer, typename TSource, typename TSpec>
void write(TFile & file, Align<TSource, TSpec>* align, TStringContainer & ids, CGViz) {
//IOREV _notinlined_
	SEQAN_CHECKPOINT
	_writeImpl(file, align, ids, CGViz());
}

//____________________________________________________________________________

template <typename TFile, typename TStringContainer, typename TSource, typename TSpec, typename TMeta>
void write(TFile & file, Align<TSource, TSpec> & align, TStringContainer& ids, TMeta &, CGViz) {
//IOREV _notinlined_
	SEQAN_CHECKPOINT
	_writeImpl(file, align, ids, CGViz());
}



//////////////////////////////////////////////////////////////////////////////
} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
