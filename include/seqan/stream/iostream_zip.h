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

#ifndef INCLUDE_SEQAN_STREAM_IOSTREAM_ZIP_H_
#define INCLUDE_SEQAN_STREAM_IOSTREAM_ZIP_H_

namespace zlib_stream {

// Default gzip buffer size, change this to suite your needs.
const size_t ZIP_DEFAULT_BUFFER_SIZE = 921600;

// --------------------------------------------------------------------------
// Enum EStrategy
// --------------------------------------------------------------------------
// Compression strategy, see zlib doc.

enum EStrategy
{
    StrategyFiltered = 1,
    StrategyHuffmanOnly = 2,
    DefaultStrategy = 0
};

// ===========================================================================
// Classes
// ===========================================================================

// --------------------------------------------------------------------------
// Class basic_zip_streambuf
// --------------------------------------------------------------------------
// A stream decorator that takes raw input and zips it to a ostream.
// The class wraps up the inflate method of the zlib library 1.1.4 http://www.gzip.org/zlib/

template <typename Elem,
          typename Tr = std::char_traits<Elem>,
          typename ElemA = std::allocator<Elem>,
          typename ByteT = unsigned char,
          typename ByteAT = std::allocator<ByteT>
          >
class basic_zip_streambuf :
    public std::basic_streambuf<Elem, Tr>
{
public:
    typedef std::basic_ostream<Elem, Tr> &              ostream_reference;
    typedef ElemA                                       char_allocator_type;
    typedef ByteT                                       byte_type;
    typedef ByteAT                                      byte_allocator_type;
    typedef byte_type *                                 byte_buffer_type;
    typedef Tr                                          traits_type;
    typedef typename Tr::char_type                      char_type;
    typedef typename Tr::int_type                       int_type;
    typedef std::vector<byte_type, byte_allocator_type> byte_vector_type;
    typedef std::vector<char_type, char_allocator_type> char_vector_type;

    // Construct a zip stream
    // More info on the following parameters can be found in the zlib documentation.
    basic_zip_streambuf(ostream_reference ostream_,
                        size_t level_,
                        EStrategy strategy_,
                        size_t window_size_,
                        size_t memory_level_,
                        size_t buffer_size_);

    ~basic_zip_streambuf();

    int sync();
    int_type overflow(int_type c);

    // flushes the zip buffer and output buffer.
    // This method should be called at the end of the compression.
    // Calling flush multiple times, will lower the compression ratio.
    std::streamsize flush();

    // flushes the zip buffer and output buffer and finalize the zip stream
    // This method should be called at the end of the compression.
    std::streamsize flush_finalize();


private:
    bool zip_to_stream(char_type *, std::streamsize);
    size_t fill_input_buffer();
    // flush the zip buffer using a particular mode and flush output buffer
    std::streamsize flush(int flush_mode);

    ostream_reference m_ostream;
    z_stream m_zip_stream;
    int m_err;
    byte_vector_type m_output_buffer;
    char_vector_type m_buffer;
};

// --------------------------------------------------------------------------
// Class basic_unzip_streambuf
// --------------------------------------------------------------------------
// A stream decorator that takes compressed input and unzips it to a istream.
// The class wraps up the deflate method of the zlib library 1.1.4 http://www.gzip.org/zlib/

template <typename Elem,
          typename Tr = std::char_traits<Elem>,
          typename ElemA = std::allocator<Elem>,
          typename ByteT = unsigned char,
          typename ByteAT = std::allocator<ByteT>
          >
class basic_unzip_streambuf :
    public std::basic_streambuf<Elem, Tr>
{
public:
    typedef std::basic_istream<Elem, Tr> &              istream_reference;
    typedef ElemA                                       char_allocator_type;
    typedef ByteT                                       byte_type;
    typedef ByteAT                                      byte_allocator_type;
    typedef byte_type *                                 byte_buffer_type;
    typedef Tr                                          traits_type;
    typedef typename Tr::char_type                      char_type;
    typedef typename Tr::int_type                       int_type;
    typedef std::vector<byte_type, byte_allocator_type> byte_vector_type;
    typedef std::vector<char_type, char_allocator_type> char_vector_type;

    // Construct a unzip stream
    // More info on the following parameters can be found in the zlib documentation.
    basic_unzip_streambuf(istream_reference istream_,
                          size_t window_size_,
                          size_t read_buffer_size_,
                          size_t input_buffer_size_);

    ~basic_unzip_streambuf();

    int_type underflow();

    // returns the compressed input istream
    istream_reference get_istream()  { return m_istream; }
    // returns the zlib stream structure
    z_stream & get_zip_stream()      { return m_zip_stream; }

private:
    void put_back_from_zip_stream();
    std::streamsize unzip_from_stream(char_type *, std::streamsize);
    size_t fill_input_buffer();

    istream_reference m_istream;
    z_stream m_zip_stream;
    int m_err;
    byte_vector_type m_input_buffer;
    char_vector_type m_buffer;
};

// --------------------------------------------------------------------------
// Class basic_zip_ostreambase
// --------------------------------------------------------------------------
// Base class for zip ostreams.
// Contains a basic_zip_streambuf.

template <typename Elem,
          typename Tr = std::char_traits<Elem>,
          typename ElemA = std::allocator<Elem>,
          typename ByteT = unsigned char,
          typename ByteAT = std::allocator<ByteT>
          >
class basic_zip_ostreambase :
    virtual public std::basic_ios<Elem, Tr>
{
public:
    typedef std::basic_ostream<Elem, Tr> &                      ostream_reference;
    typedef basic_zip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT> zip_streambuf_type;

    // Construct a zip stream
    // More info on the following parameters can be found in the zlib documentation.
    basic_zip_ostreambase(ostream_reference ostream_,
                          size_t level_,
                          EStrategy strategy_,
                          size_t window_size_,
                          size_t memory_level_,
                          size_t buffer_size_) :
        m_buf(ostream_, level_, strategy_, window_size_, memory_level_, buffer_size_)
    {
        this->init(&m_buf);
    }

    // returns the underlying zip ostream object
    zip_streambuf_type * rdbuf() { return &m_buf; }

private:
    zip_streambuf_type m_buf;
};

// --------------------------------------------------------------------------
// Class basic_zip_istreambase
// --------------------------------------------------------------------------
// Base class for unzip istreams
// Contains a basic_unzip_streambuf.

template <typename Elem,
          typename Tr = std::char_traits<Elem>,
          typename ElemA = std::allocator<Elem>,
          typename ByteT = unsigned char,
          typename ByteAT = std::allocator<ByteT>
          >
class basic_zip_istreambase :
    virtual public std::basic_ios<Elem, Tr>
{
public:
    typedef std::basic_istream<Elem, Tr> &                        istream_reference;
    typedef basic_unzip_streambuf<Elem, Tr, ElemA, ByteT, ByteAT> unzip_streambuf_type;

    basic_zip_istreambase(istream_reference ostream_,
                          size_t window_size_,
                          size_t read_buffer_size_,
                          size_t input_buffer_size_) :
        m_buf(ostream_, window_size_, read_buffer_size_, input_buffer_size_)
    {
        this->init(&m_buf);
    }

    // returns the underlying unzip istream object
    unzip_streambuf_type * rdbuf() { return &m_buf; }

private:
    unzip_streambuf_type m_buf;
};

// --------------------------------------------------------------------------
// Class basic_zip_ostream
// --------------------------------------------------------------------------
// A zipper ostream
//
// This class is a ostream decorator that behaves 'almost' like any other ostream.
// At construction, it takes any ostream that shall be used to output of the compressed data.
// When finished, you need to call the special method zflush or call the destructor
// to flush all the intermidiate streams.
//
// Example:
//
// // creating the target zip string, could be a fstream
// ostringstream ostringstream_;
// // creating the zip layer
// zip_ostream zipper(ostringstream_);
// // writing data
// zipper<<f<<" "<<d<<" "<<ui<<" "<<ul<<" "<<us<<" "<<c<<" "<<dum;
// // zip ostream needs special flushing...
// zipper.zflush();

template <typename Elem,
          typename Tr = std::char_traits<Elem>,
          typename ElemA = std::allocator<Elem>,
          typename ByteT = unsigned char,
          typename ByteAT = std::allocator<ByteT>
          >
class basic_zip_ostream :
    public basic_zip_ostreambase<Elem, Tr, ElemA, ByteT, ByteAT>,
    public std::basic_ostream<Elem, Tr>
{
public:
    typedef basic_zip_ostreambase<Elem, Tr, ElemA, ByteT, ByteAT> zip_ostreambase_type;
    typedef std::basic_ostream<Elem, Tr>                          ostream_type;
    typedef ostream_type &                                        ostream_reference;

    // Constructs a zipper ostream decorator
    //
    // ostream_ ostream where the compressed output is written
    // is_gzip_ true if gzip header and footer have to be added
    // level_ level of compression 0, bad and fast, 9, good and slower,
    // strategy_ compression strategy
    // window_size_ see zlib doc
    // memory_level_ see zlib doc
    // buffer_size_ the buffer size used to zip data

    basic_zip_ostream(ostream_reference ostream_,
                      size_t level_ = Z_DEFAULT_COMPRESSION,
                      EStrategy strategy_ = DefaultStrategy,
                      size_t window_size_ = 31, // 15 (size) + 16 (gzip header)
                      size_t memory_level_ = 8,
                      size_t buffer_size_ = ZIP_DEFAULT_BUFFER_SIZE) :
        zip_ostreambase_type(ostream_, level_, strategy_, window_size_, memory_level_, buffer_size_),
        ostream_type(this->rdbuf())
    {}

    ~basic_zip_ostream()
    {
        this->flush(); this->rdbuf()->flush_finalize();
    }

    // flush inner buffer and zipper buffer
    basic_zip_ostream<Elem, Tr> & zflush()
    {
        this->flush(); this->rdbuf()->flush(); return *this;
    }

#ifdef _WIN32
private:
    void _Add_vtordisp1() {}  // Required to avoid VC++ warning C4250
    void _Add_vtordisp2() {}  // Required to avoid VC++ warning C4250
#endif
};

// --------------------------------------------------------------------------
// Class basic_zip_istream
// --------------------------------------------------------------------------
// A zipper istream
//
// This class is a istream decorator that behaves 'almost' like any other ostream.
// At construction, it takes any istream that shall be used to input of the compressed data.
//
// Simlpe example:
//
// // create a stream on zip string
// istringstream istringstream_( ostringstream_.str());
// // create unzipper istream
// zip_istream unzipper( istringstream_);
// // read and unzip
// unzipper>>f_r>>d_r>>ui_r>>ul_r>>us_r>>c_r>>dum_r;

template <typename Elem,
          typename Tr = std::char_traits<Elem>,
          typename ElemA = std::allocator<Elem>,
          typename ByteT = unsigned char,
          typename ByteAT = std::allocator<ByteT>
          >
class basic_zip_istream :
    public basic_zip_istreambase<Elem, Tr, ElemA, ByteT, ByteAT>,
    public std::basic_istream<Elem, Tr>
{
public:
    typedef basic_zip_istreambase<Elem, Tr, ElemA, ByteT, ByteAT> zip_istreambase_type;
    typedef std::basic_istream<Elem, Tr>                          istream_type;
    typedef istream_type &                                        istream_reference;
    typedef ByteT                                                 byte_type;
    typedef Tr                                                    traits_type;

    // Construct a unzipper stream
    //
    // istream_ input buffer
    // window_size_
    // read_buffer_size_
    // input_buffer_size_

    basic_zip_istream(istream_reference istream_,
                      size_t window_size_ = 31, // 15 (size) + 16 (gzip header)
                      size_t read_buffer_size_ = ZIP_DEFAULT_BUFFER_SIZE,
                      size_t input_buffer_size_ = ZIP_DEFAULT_BUFFER_SIZE) :
        zip_istreambase_type(istream_, window_size_, read_buffer_size_, input_buffer_size_),
        istream_type(this->rdbuf())
    {}

#ifdef _WIN32
private:
    void _Add_vtordisp1() {}  // Required to avoid VC++ warning C4250
    void _Add_vtordisp2() {}  // Required to avoid VC++ warning C4250
#endif
};

// ===========================================================================
// Typedefs
// ===========================================================================

// A typedef for basic_zip_ostream<char>
typedef basic_zip_ostream<char>     zip_ostream;
// A typedef for basic_zip_ostream<wchar_t>
typedef basic_zip_ostream<wchar_t>  zip_wostream;
// A typedef for basic_zip_istream<char>
typedef basic_zip_istream<char>     zip_istream;
// A typedef for basic_zip_istream<wchart>
typedef basic_zip_istream<wchar_t>  zip_wistream;

} // namespace zlib_stream

#endif // INCLUDE_SEQAN_STREAM_IOSTREAM_ZIP_H_
