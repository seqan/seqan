#ifndef RAZERS_PARALLEL_MISC_H_
#define RAZERS_PARALLEL_MISC_H_

namespace seqan {

template <typename TPos, typename TSize, typename TCount>
void computeSplittersBySlotCount(String<TPos> & splitters, TSize size, TCount count)
{
    resize(splitters, count + 1, Exact());
    splitters[0] = 0;
    TSize blockLength = size / count;
    TSize rest = size % count;
    for (TCount i = 1; i <= count; ++i)
    {
        splitters[i] = splitters[i - 1] + blockLength;
        if (i <= rest)
            splitters[i] += 1;
    }

    SEQAN_ASSERT_EQ(back(splitters), size);
}

template <typename TPos, typename TDataSize, typename TSlotSize, typename TPackageCount>
void computeSplittersBySlotSize(String<TPos> & splitters, TDataSize size, TSlotSize slotSize, TPackageCount maxPackageCount)
{
    // Limit maximal number of verification packages.
    if (maxPackageCount > 0 && (size / slotSize > maxPackageCount))
        slotSize = size / maxPackageCount + (static_cast<TSlotSize>(size % maxPackageCount) > static_cast<TSlotSize>(0));
    // Compute splitters.
    unsigned count = size / slotSize + (static_cast<TSlotSize>(size % slotSize) > static_cast<TSlotSize>(0));
    resize(splitters, count + 1, Exact());
    splitters[0] = 0;
    for (unsigned i = 1; i < count; ++i)
        splitters[i] = splitters[i - 1] + slotSize;
    splitters[count] = size;

    SEQAN_ASSERT_LEQ(splitters[count] - splitters[count - 1], slotSize);
}

template <typename TPos, typename TDataSize, typename TSlotSize>
void computeSplittersBySlotSize(String<TPos> & splitters, TDataSize size, TSlotSize slotSize)
{
    computeSplittersBySlotSize(splitters, size, slotSize, 0u);
}

}  // namespace seqan

#endif  // #ifndef RAZERS_PARALLEL_MISC_H_
