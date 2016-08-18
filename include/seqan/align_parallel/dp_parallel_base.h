// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_PARALLEL_BASE_H_
#define INCLUDE_SEQAN_ALIGN_PARALLEL_DP_PARALLEL_BASE_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag DPTiled<TBuffer>
// ----------------------------------------------------------------------------

// Tag used to subclass DPScoutState and DPScout.
// T represents the buffer type.
template <typename TBuffer, typename TSimdSpec = void>
struct DPTiled;

// ----------------------------------------------------------------------------
// Class DPTileBuffer
// ----------------------------------------------------------------------------

// The structure owning the horizontal/vertical buffer.
template <typename TDPCellBuff, typename TBuffer = String<TDPCellBuff> >
struct DPTileBuffer
{
    TBuffer horizontalBuffer;
    TBuffer verticalBuffer;
};

// ----------------------------------------------------------------------------
// Class DPScoutState_; DPTiled
// ----------------------------------------------------------------------------

// The overloaded DPScoutState which simply stores the pointers to the corresponding buffer.
template <typename TBuffer>
class DPScoutState_<DPTiled<TBuffer, void> >
{
public:

    TBuffer* ptrHorBuffer;
    TBuffer* ptrVerBuffer;

    DPScoutState_() : ptrHorBuffer(nullptr), ptrVerBuffer(nullptr)
    {}

    DPScoutState_(TBuffer & horBuffer, TBuffer & verBuffer) :
        ptrHorBuffer(&horBuffer),
        ptrVerBuffer(&verBuffer)
    {}
};

// ----------------------------------------------------------------------------
// Class DPScout_; DPTiled
// ----------------------------------------------------------------------------

// Overloaded DPScout to store the corresponding buffer for the current dp tile.
template <typename TDPCell, typename TBuffer>
class DPScout_<TDPCell, DPTiled<TBuffer, void> > :
    public DPScout_<TDPCell, Default>
{
public:
    using TBase = DPScout_<TDPCell, Default>;

    DPScoutState_<DPTiled<TBuffer> > state = DPScoutState_<DPTiled<TBuffer> >{};

    DPScout_() = default;

    DPScout_(DPScoutState_<DPTiled<TBuffer, void> > state) :
        TBase(),
        state(state)
    {}
};

namespace impl
{

namespace debug
{
template <typename TBuffer>
struct DebugSubMatrix
{
    TBuffer hBegin;
    TBuffer vBegin;
    TBuffer hEnd;
    TBuffer vEnd;
    size_t  col;
    size_t  row;
    bool    isSimd;
};

template <typename TBuffer>
struct DebugBuffer
{
    String<String<DebugSubMatrix<TBuffer> > > matrix;


    template <typename TStream>
    void write(TStream & stream)
    {
        auto hBufSize = length(matrix[0][0].hBegin) + 1;
        auto vBufSize = length(matrix[0][0].vBegin) + 1;

//        stream << ',';
//        for (unsigned col = 0; col < length(matrix); ++col)
//        {
//            for (auto c : seqH[col])
//                stream << c << ',';
//            stream << ',';
//        }

        for (unsigned row = 0; row < length(matrix[0]); ++row)
        {
            // Write Block Info.
            for (unsigned col = 0; col < length(matrix); ++col)
            {
                stream << "Block:," << matrix[col][row].col << ',' << matrix[col][row].row << ',' << matrix[col][row].isSimd;
                for (unsigned i = 3; i < hBufSize; ++i)
                    stream << ',';
            }
            stream << '\n';
            // Write hBegin.
            for (unsigned col = 0; col < length(matrix); ++col)
            {
                stream << ',';
                for (unsigned i = 1; i < hBufSize; ++i)
                    stream << matrix[col][row].hBegin[i - 1].i1._score << ',';
            }
            stream << '\n';
            // Write vLines.
            for (unsigned r = 1; r < vBufSize; ++r)
            {
                for (unsigned col = 0; col < length(matrix); ++col)
                {
                    stream << matrix[col][row].vBegin[r - 1].i1._score << ',';
                    for (unsigned c = 1; c < hBufSize - 1; ++c)
                        stream << ',';
                    stream << matrix[col][row].vEnd[r - 1].i1._score << ',';
                }
                stream << '\n';
            }
            // Write hEnd.
            for (unsigned col = 0; col < length(matrix); ++col)
            {
                stream << ',';
                for (unsigned i = 1; i < hBufSize; ++i)
                    stream << matrix[col][row].hEnd[i - 1].i1._score << ',';
            }
            stream << '\n';
        }
    }
};
}  // namespace impl::debug

namespace dp
{

namespace parallel
{

// ----------------------------------------------------------------------------
// Tag LocalTraceStore
// ----------------------------------------------------------------------------

template <typename TSimdVec>
class LocalTraceStore
{
public:

    // Typedefs
    // ----------------------------------------------------------------------------

    using TScalarTraceValue  = typename TraceBitMap_<>::Type;
    using TSimdTraceValue    = typename TraceBitMap_<TSimdVec>::Type;
    using TScalarTraceMatrix = String<TScalarTraceValue>;
    using TSimdTraceMatrix   = String<TSimdTraceValue, Alloc<OverAligned> >;

    // Members
    // ----------------------------------------------------------------------------

    String<TScalarTraceMatrix> mScalarTraceVec;
    String<TSimdTraceMatrix>   mSimdTraceVec;

    // Constructors
    // ----------------------------------------------------------------------------

    // Member functions
    // ----------------------------------------------------------------------------

    inline auto& localScalarTraceMatrix()
    {
        resize(mScalarTraceVec, length(mScalarTraceVec) + 1, Generous());
        return back(mScalarTraceVec);
    }

    inline auto& localSimdTraceMatrix()
    {
        resize(mSimdTraceVec, length(mSimdTraceVec) + 1, Generous());
        return back(mSimdTraceVec);
    }
};

// ----------------------------------------------------------------------------
// Tag TraceProxy
// ----------------------------------------------------------------------------

template <typename TLocalTraceStore_>
class TraceProxy
{
public:
    struct BlockIdentifier
    {
        uint16_t mBlockId  : 1;  // Use first bit for block identifier.
        uint16_t mBlockPos : 15;  // Use remaining bits for pos.
    };

    // Typedefs
    // ----------------------------------------------------------------------------

    using TLocalTraceStore = TLocalTraceStore_;
    using TTraceMatrixIdentifier = std::tuple<TLocalTraceStore*, BlockIdentifier, uint8_t>;

    // Members
    // ----------------------------------------------------------------------------

    size_t                              mHorSize = 0;
    size_t                              mVerSize = 0;
    std::vector<TTraceMatrixIdentifier> mTraceBlockMap;
    
    // Constructors
    // ----------------------------------------------------------------------------
    
    TraceProxy() = default;
    
    template <typename THSize, typename TVSize>
    TraceProxy(THSize const hSize, TVSize const vSize) :
        mHorSize(hSize),
        mVerSize(vSize),
        mTraceBlockMap(mHorSize * mVerSize)
    {}
    
    // Member functions
    // ----------------------------------------------------------------------------
    
    // Thread-safe as long as concurrent threads are not writing into same position.
    template <typename T1, typename T2>
    inline void
    insert(std::pair<T1, T2> const & key,
           TTraceMatrixIdentifier && value)
    {
        SEQAN_ASSERT_LT(key.first * mVerSize + key.second, mTraceBlockMap.size());
        swap(mTraceBlockMap[key.first * mVerSize + key.second], value);
    }
};

}  // namespace impl::dp::parallel
}  // namespace impl::dp
}  // namespace impl


// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ScoutSpecForSimdAlignment_
// ----------------------------------------------------------------------------

template<typename TAlignmentAlgorithm, typename TBuffer>
struct ScoutSpecForAlignmentAlgorithm_<TAlignmentAlgorithm, DPScoutState_<DPTiled<TBuffer, void> > >
{
    using Type = DPTiled<TBuffer, void>;
};

namespace impl
{

namespace dp
{

namespace parallel
{

template <typename TLocalStore, typename TParSpec>
struct ThreadLocalStorage;
}  // namespace parallel
}  // namespace dp
}  // namespace impl

// ============================================================================
// Functions
// ============================================================================

namespace impl
{

namespace debug
{
// Helper function for demangling the symbols.
inline auto demangle(const char * const name)
{
    int status;
    char* result = abi::__cxa_demangle(name, nullptr, nullptr, &status);
    if (status == 0)
    {
        CharString tmp(result);
//        free(result);
        return tmp;
    }
    throw std::runtime_error("Could not demangle name!");
}
}  // namespace debug
}  // namespace impl
}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ALIGN_PARALLEL_DP_PARALLEL_BASE_H_
