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

namespace detail {

const int gz_magic[2] = { 0x1f, 0x8b };     // gzip magic header

// gzip flag byte
const int gz_ascii_flag  = 0x01;     // bit 0 set: file probably ascii text
const int gz_head_crc    = 0x02;     // bit 1 set: header CRC present
const int gz_extra_field = 0x04;     // bit 2 set: extra field present
const int gz_orig_name   = 0x08;     // bit 3 set: original file name present
const int gz_comment     = 0x10;     // bit 4 set: file comment present
const int gz_reserved    = 0xE0;     // bits 5..7: reserved

} // namespace detail

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
    m_buffer(buffer_size_, 0),
    m_crc(0)
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
        -static_cast<int>(window_size_),      // <-- changed
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

    // updating crc
    m_crc = crc32(m_crc, m_zip_stream.next_in, m_zip_stream.avail_in);

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

    // updating crc
    m_crc = crc32(m_crc, m_zip_stream.next_in, m_zip_stream.avail_in);

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
    m_buffer(read_buffer_size_),
    m_crc(0)
{
    // setting zalloc, zfree and opaque
    m_zip_stream.zalloc = (alloc_func)0;
    m_zip_stream.zfree = (free_func)0;

    m_zip_stream.next_in = NULL;
    m_zip_stream.avail_in = 0;
    m_zip_stream.avail_out = 0;
    m_zip_stream.next_out = NULL;

    m_err = inflateInit2(&m_zip_stream, -static_cast<int>(window_size_));

    this->setg(&(m_buffer[0]) + 4,  // beginning of putback area
               &(m_buffer[0]) + 4,  // read position
               &(m_buffer[0]) + 4); // end position
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
    }
    while (m_err == Z_OK && m_zip_stream.avail_out != 0 && count != 0);

    // updating crc
    m_crc = crc32(m_crc,
                (byte_buffer_type)buffer_,
                buffer_size_ - m_zip_stream.avail_out / sizeof(char_type));
    std::streamsize n_read = buffer_size_ - m_zip_stream.avail_out / sizeof(char_type);

    // check if it is the end
    if (m_err == Z_STREAM_END)
        put_back_from_zip_stream();

    return n_read;
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

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
int basic_zip_istream<Elem, Tr, ElemA, ByteT, ByteAT>::check_header()
{
    int method;     // method byte
    int flags;      // flags byte
    uInt len;
    int c;
    int err = 0;
    z_stream & zip_stream = this->rdbuf()->get_zip_stream();

    std::basic_istream<Elem, Tr> & istream = this->rdbuf()->get_istream();

    // Check the gzip magic header
    for (len = 0; len < 2; len++)
    {
        c = (int)istream.get();
        if (c != zlib_stream::detail::gz_magic[len])
        {
            if (len != 0)
                istream.unget();
            if (!traits_type::eq_int_type(c, traits_type::eof()))
                istream.unget();

            err = zip_stream.avail_in != 0 ? Z_OK : Z_STREAM_END;
            m_is_gzip = false;
            return err;
        }
    }

    m_is_gzip = true;
    method = (int)istream.get();
    flags = (int)istream.get();
    if (method != Z_DEFLATED || (flags & zlib_stream::detail::gz_reserved) != 0)
    {
        err = Z_DATA_ERROR;
        return err;
    }

    // Discard time, xflags and OS code:
    for (len = 0; len < 6; len++)
        istream.get();

    if ((flags & zlib_stream::detail::gz_extra_field) != 0)
    {
        // skip the extra field
        len  =  (uInt)istream.get();
        len += ((uInt)istream.get()) << 8;
        // len is garbage if EOF but the loop below will quit anyway
        while (len-- != 0 && !traits_type::eq_int_type(istream.get(), traits_type::eof())) ;
    }
    if ((flags & zlib_stream::detail::gz_orig_name) != 0)
    {
        // skip the original file name
        while ((c = istream.get()) != 0 && !traits_type::eq_int_type(c, traits_type::eof())) ;
    }
    if ((flags & zlib_stream::detail::gz_comment) != 0)
    {
        // skip the .gz file comment
        while ((c = istream.get()) != 0 && !traits_type::eq_int_type(c, traits_type::eof())) ;
    }
    if ((flags & zlib_stream::detail::gz_head_crc) != 0)
    {      // skip the header crc
        for (len = 0; len < 2; len++)
            istream.get();
    }
    err = istream.eof() ? Z_DATA_ERROR : Z_OK;

    return err;
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
void basic_zip_istream<Elem, Tr, ElemA, ByteT, ByteAT>::read_footer()
{
    if (m_is_gzip)
    {
        read_uint32(m_gzip_crc);
        read_uint32(m_gzip_data_size);
    }
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
void basic_zip_ostream<Elem, Tr, ElemA, ByteT, ByteAT>::put_uint32(uLong x)
{
    this->write(reinterpret_cast<char *>(&x), 1);
    this->write(reinterpret_cast<char *>(&x + 1), 1);
    this->write(reinterpret_cast<char *>(&x + 2), 1);
    this->write(reinterpret_cast<char *>(&x + 3), 1);
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
void basic_zip_istream<Elem, Tr, ElemA, ByteT, ByteAT>::read_uint32(uLong & x)
{
    this->read(reinterpret_cast<char *>(&x), 1);
    this->read(reinterpret_cast<char *>(&x + 1), 1);
    this->read(reinterpret_cast<char *>(&x + 2), 1);
    this->read(reinterpret_cast<char *>(&x + 3), 1);
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
void basic_zip_ostream<Elem, Tr, ElemA, ByteT, ByteAT>::add_header()
{
          typename Tr::char_type zero = 0;

    this->rdbuf()->get_ostream()
    .put(static_cast<typename Tr::char_type>(detail::gz_magic[0]))
    .put(static_cast<typename Tr::char_type>(detail::gz_magic[1]))
    .put(static_cast<typename Tr::char_type>(Z_DEFLATED))
    .put(zero)         //flags
    .put(zero).put(zero).put(zero).put(zero)         // time
    .put(zero)         //xflags
    .put(static_cast<typename Tr::char_type>(OS_CODE));
}

template <typename Elem,
          typename Tr,
          typename ElemA,
          typename ByteT,
          typename ByteAT>
void basic_zip_ostream<Elem, Tr, ElemA, ByteT, ByteAT>::add_footer()
{
    put_uint32(this->rdbuf()->get_crc());
    put_uint32(this->rdbuf()->get_in_size());
}

} // namespace zlib_stream


#endif // INCLUDE_SEQAN_STREAM_IOSTREAM_ZIP_IMPL_H_
