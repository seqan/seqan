// zipstream Library License:
// --------------------------
//
// The zlib/libpng License Copyright (c) 2003 Jonathan de Halleux.
//
// This software is provided 'as-is', without any express or implied warranty. In no event will the authors be held liable for any damages arising from the use of this software.
//
// Permission is granted to anyone to use this software for any purpose, including commercial applications, and to alter it and redistribute it freely, subject to the following restrictions:
//
// 1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
//
// 2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
//
// 3. This notice may not be removed or altered from any source distribution
//
// Altered zipstream library header
// Author: Jonathan de Halleux, dehalleux@pelikhan.com, 2003
// Author: David Weese <david.weese@fu-berlin.de>
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>

#ifndef INCLUDE_SEQAN_STREAM_IOSTREAM_ZIP_IMPL_H_
#define INCLUDE_SEQAN_STREAM_IOSTREAM_ZIP_IMPL_H_

namespace zlib_stream {

// --------------------------------------------------------------------------
// Class basic_zip_streambuf
// --------------------------------------------------------------------------

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::basic_zip_streambuf(
    ostream_reference ostream_,
    size_t level_,
    EStrategy strategy_,
    size_t window_size_,
    size_t memory_level_,
    size_t buffer_size_
    ) :
    m_ostream(ostream_),
    m_output_buffer(buffer_size_, 0),
    m_buffer(buffer_size_, 0)
{
    m_zip_stream.zalloc = (alloc_func)0;
    m_zip_stream.zfree = (free_func)0;

    m_zip_stream.next_in = NULL;
    m_zip_stream.avail_in = 0;
    m_zip_stream.avail_out = 0;
    m_zip_stream.next_out = NULL;

    m_err = deflateInit2(
        &m_zip_stream,
        std::min(9, static_cast<int>(level_)),
        Z_DEFLATED,
        static_cast<int>(window_size_),
        std::min(9, static_cast<int>(memory_level_)),
        static_cast<int>(strategy_)
        );

    this->setp(&(m_buffer[0]), &(m_buffer[m_buffer.size() - 1]));
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::~basic_zip_streambuf()
{
    flush_finalize();
    m_ostream.flush();
    m_err = deflateEnd(&m_zip_stream);
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
int basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::sync()
{
    if (this->pptr() && this->pptr() > this->pbase())
    {
        if (traits_type::eq_int_type(overflow(traits_type::eof()), traits_type::eof()))
            return -1;
    }

    return 0;
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
typename basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::int_type
basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::overflow(
    typename basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::int_type c)
{
    int w = static_cast<int>(this->pptr() - this->pbase());

    if (!traits_type::eq_int_type(c, traits_type::eof()))
    {
        *this->pptr() = c;
        ++w;
    }

    if (zip_to_stream(this->pbase(), w))
    {
        this->setp(this->pbase(), this->epptr() - 1);
        return c;
    }
    else
    {
        return traits_type::eof();
    }
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
bool basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::zip_to_stream(
    typename basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::char_type * buffer_,
    std::streamsize buffer_size_)
{
    std::streamsize written_byte_size = 0, total_written_byte_size = 0;

    m_zip_stream.next_in = (byte_buffer_type)buffer_;
    m_zip_stream.avail_in = static_cast<uInt>(buffer_size_ * sizeof(char_type));
    m_zip_stream.avail_out = static_cast<uInt>(m_output_buffer.size());
    m_zip_stream.next_out = &(m_output_buffer[0]);
    size_t remainder = 0;

    do
    {
        m_err = deflate(&m_zip_stream, 0);

        if (m_err == Z_OK  || m_err == Z_STREAM_END)
        {
            written_byte_size = static_cast<std::streamsize>(m_output_buffer.size()) - m_zip_stream.avail_out;
            total_written_byte_size += written_byte_size;

            // ouput buffer is full, dumping to ostream
            m_ostream.write((const char_type *) &(m_output_buffer[0]),
                            static_cast<std::streamsize>(written_byte_size / sizeof(char_type)));

            // checking if some bytes were not written.
            if ((remainder = written_byte_size % sizeof(char_type)) != 0)
            {
                // copy to the beginning of the stream
                std::memmove(&(m_output_buffer[0]),
                             &(m_output_buffer[written_byte_size - remainder]),
                             remainder);
            }

            m_zip_stream.avail_out = static_cast<uInt>(m_output_buffer.size() - remainder);
            m_zip_stream.next_out = &m_output_buffer[remainder];
        }
    }
    while (m_zip_stream.avail_in != 0 && m_err == Z_OK);

    return m_err == Z_OK;
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
std::streamsize basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::flush(int flush_mode)
{
    int const buffer_size = static_cast<int>(this->pptr() - this->pbase()); // amount of data currently in buffer

    std::streamsize written_byte_size = 0, total_written_byte_size = 0;

    m_zip_stream.next_in = (byte_buffer_type) this->pbase();
    m_zip_stream.avail_in = static_cast<uInt>(buffer_size * sizeof(char_type));
    m_zip_stream.avail_out = static_cast<uInt>(m_output_buffer.size());
    m_zip_stream.next_out = &(m_output_buffer[0]);
    size_t remainder = 0;

    do
    {
        m_err = deflate(&m_zip_stream, flush_mode);
        if (m_err == Z_OK || m_err == Z_STREAM_END)
        {
            written_byte_size = static_cast<std::streamsize>(m_output_buffer.size()) - m_zip_stream.avail_out;
            total_written_byte_size += written_byte_size;

            // ouput buffer is full, dumping to ostream
            m_ostream.write((const char_type *) &(m_output_buffer[0]),
                            static_cast<std::streamsize>(written_byte_size / sizeof(char_type) * sizeof(byte_type)));

            // checking if some bytes were not written.
            if ((remainder = written_byte_size % sizeof(char_type)) != 0)
            {
                // copy to the beginning of the stream
                std::memmove(&(m_output_buffer[0]),
                             &(m_output_buffer[written_byte_size - remainder]),
                             remainder);
            }

            m_zip_stream.avail_out = static_cast<uInt>(m_output_buffer.size() - remainder);
            m_zip_stream.next_out = &m_output_buffer[remainder];
        }
    }
    while (m_err == Z_OK);

    m_ostream.flush();

    return total_written_byte_size;
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
std::streamsize basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::flush()
{
    return flush(Z_SYNC_FLUSH);
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
std::streamsize basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::flush_finalize()
{
    return flush(Z_FINISH);
}

// --------------------------------------------------------------------------
// Class basic_unzip_streambuf
// --------------------------------------------------------------------------

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
basic_unzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::basic_unzip_streambuf(
    istream_reference istream_,
    size_t window_size_,
    size_t read_buffer_size_,
    size_t input_buffer_size_
    ) :
    m_istream(istream_),
    m_input_buffer(input_buffer_size_),
    m_buffer(read_buffer_size_)
{
    // setting zalloc, zfree and opaque
    m_zip_stream.zalloc = (alloc_func)0;
    m_zip_stream.zfree = (free_func)0;

    m_zip_stream.next_in = NULL;
    m_zip_stream.avail_in = 0;
    m_zip_stream.avail_out = 0;
    m_zip_stream.next_out = NULL;

    m_err = inflateInit2(&m_zip_stream, static_cast<int>(window_size_));

    this->setg(&(m_buffer[0]) + 4,  // beginning of putback area
               &(m_buffer[0]) + 4,  // read position
               &(m_buffer[0]) + 4); // end position
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
basic_unzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::~basic_unzip_streambuf()
{
    inflateEnd(&m_zip_stream);
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
typename basic_unzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::int_type
basic_unzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::underflow()
{
    if (this->gptr() && (this->gptr() < this->egptr()))
        return *reinterpret_cast<unsigned char *>(this->gptr());

    int n_putback = static_cast<int>(this->gptr() - this->eback());
    if (n_putback > 4)
        n_putback = 4;

    std::memmove(&(m_buffer[0]) + (4 - n_putback), this->gptr() - n_putback, n_putback * sizeof(char_type));

    int num = unzip_from_stream(&(m_buffer[0]) + 4,
                                static_cast<std::streamsize>((m_buffer.size() - 4) * sizeof(char_type)));

    if (num <= 0)     // ERROR or EOF
        return traits_type::eof();

    // reset buffer pointers
    this->setg(&(m_buffer[0]) + (4 - n_putback),         // beginning of putback area
               &(m_buffer[0]) + 4,                       // read position
               &(m_buffer[0]) + 4 + num);                // end of buffer

    // return next character
    return *reinterpret_cast<unsigned char *>(this->gptr());
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
std::streamsize basic_unzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::unzip_from_stream(
    char_type * buffer_,
    std::streamsize buffer_size_)
{
    m_zip_stream.next_out = (byte_buffer_type)buffer_;
    m_zip_stream.avail_out = static_cast<uInt>(buffer_size_ * sizeof(char_type));
    size_t count = m_zip_stream.avail_in;

    do
    {
        if (m_zip_stream.avail_in == 0)
            count = fill_input_buffer();

        if (m_zip_stream.avail_in)
            m_err = inflate(&m_zip_stream, Z_SYNC_FLUSH);

        if (m_err == Z_STREAM_END)
            inflateReset(&m_zip_stream);
        else if (m_err < 0)
            break;
    }
    while (m_zip_stream.avail_out > 0 && count > 0);

    std::streamsize n_read = buffer_size_ - m_zip_stream.avail_out / sizeof(char_type);

    // check if it is the end
    if (m_zip_stream.avail_out > 0 && m_err == Z_STREAM_END)
        put_back_from_zip_stream();

    return n_read;
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
size_t basic_unzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::fill_input_buffer()
{
    m_zip_stream.next_in = &(m_input_buffer[0]);
    m_istream.read((char_type *)(&(m_input_buffer[0])),
                   static_cast<std::streamsize>(m_input_buffer.size() / sizeof(char_type)));
    return m_zip_stream.avail_in = m_istream.gcount() * sizeof(char_type);
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
void basic_unzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT>::put_back_from_zip_stream()
{
    if (m_zip_stream.avail_in == 0)
        return;

    m_istream.clear(std::ios::goodbit);
    m_istream.seekg(-static_cast<int>(m_zip_stream.avail_in), std::ios_base::cur);

    m_zip_stream.avail_in = 0;
}

} // namespace zlib_stream


#endif // INCLUDE_SEQAN_STREAM_IOSTREAM_ZIP_IMPL_H_
