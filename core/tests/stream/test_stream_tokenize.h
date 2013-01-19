// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Tests for seqan/stream/tokenize.h
// ==========================================================================

#ifndef TEST_STREAM_TEST_STREAM_TOKENIZING_H_
#define TEST_STREAM_TEST_STREAM_TOKENIZING_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include "test_stream_generic.h"

std::fstream* createFile()
{
    using namespace seqan;

    const char * tempFilename = SEQAN_TEMP_FILENAME();
    char filenameBuffer[1000];
    strcpy(filenameBuffer, tempFilename);

    std::fstream *file = new std::fstream(filenameBuffer, std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    SEQAN_ASSERT(file->is_open());

    char const * STR = "This is a string ....foobar AAAACCCGGGTTT\n\
 TCG > foo.. SJAUDOF78456kapLP345LPL AAAAAAAAAAAAA\n\
\a 123gogogo \r\n\
Tetstetetstststetststetststetsstetetstetststetstdtetstestst    \r\n\
etstetetstststetststetststV           tststetstdtetsteststetstedetstet\r\n\
etstetetstststetststetststetsstetetstetststetstdtetsteststetstedetst|t\r\n\
\n\v\n\r\nACGTACGTACGTACGATACGATCTnn\n\nACGT*";


    file->write(STR, strlen(STR));
    file->seekg(0);
    file->seekp(0);
    return file;
}

// _readHelper function ""readUntil""
SEQAN_DEFINE_TEST(test_stream_tokenizing_readUntil)
{
    using namespace seqan;
    std::fstream * file = createFile();
    RecordReader<std::fstream, SinglePass<void> > reader(*file);

    CharString buf;
    SEQAN_ASSERT_EQ(readUntilWhitespace(buf,reader), 0);
    SEQAN_ASSERT_EQ(buf, "This");
    SEQAN_ASSERT_EQ(value(reader), ' ');

    clear(buf);
    SEQAN_ASSERT_EQ(readUntilChar(buf,reader, 'f'), 0);
    SEQAN_ASSERT_EQ(buf, " is a string ....");
    SEQAN_ASSERT_EQ(value(reader), 'f');

    clear(buf);
    SEQAN_ASSERT_EQ(readUntilBlank(buf,reader), 0);
    SEQAN_ASSERT_EQ(buf, "foobar");
    SEQAN_ASSERT_EQ(value(reader), ' ');

    // check EOF-handling
    clear(buf);
    SEQAN_ASSERT_EQ(readUntilChar(buf,reader, '*'), 0);
    SEQAN_ASSERT_EQ(value(reader), '*');

    clear(buf);
    SEQAN_ASSERT_EQ(readUntilBlank(buf,reader), EOF_BEFORE_SUCCESS);
    SEQAN_ASSERT_EQ(buf, "*");

    file->close();
    delete file;
}

SEQAN_DEFINE_TEST(test_stream_tokenizing_read_until_one_of)
{
    using namespace seqan;
    typedef Stream<CharArray<char const *> > TStream;

    char const * FILE_CONTENTS = "this\tis\ta\ntest\ntest\r";
    TStream stream(FILE_CONTENTS, FILE_CONTENTS + strlen(FILE_CONTENTS));
    RecordReader<TStream, SinglePass<> > reader(stream);
    
    CharString buf;
    
    SEQAN_ASSERT_EQ(readUntilOneOf(buf, reader, '\t'), 0);
    SEQAN_ASSERT_EQ(buf, "this");
    SEQAN_ASSERT_EQ(value(reader), '\t');
    SEQAN_ASSERT_NOT(goNext(reader));
    clear(buf);
    SEQAN_ASSERT_EQ(readUntilOneOf(buf, reader, '\t', '\n'), 0);
    SEQAN_ASSERT_EQ(buf, "is");
    SEQAN_ASSERT_EQ(value(reader), '\t');
    SEQAN_ASSERT_NOT(goNext(reader));
    clear(buf);
    SEQAN_ASSERT_EQ(readUntilOneOf(buf, reader, '\t', '\n', '\v'), 0);
    SEQAN_ASSERT_EQ(buf, "a");
    SEQAN_ASSERT_EQ(value(reader), '\n');
    SEQAN_ASSERT_NOT(goNext(reader));
    clear(buf);
    SEQAN_ASSERT_EQ(readUntilOneOf(buf, reader, '\t', '\n', '\v', ' '), 0);
    SEQAN_ASSERT_EQ(buf, "test");
    SEQAN_ASSERT_EQ(value(reader), '\n');
    SEQAN_ASSERT_NOT(goNext(reader));
    clear(buf);
    SEQAN_ASSERT_EQ(readUntilOneOf(buf, reader, '\t', '\n', '\v', ' ', '\r'), 0);
    SEQAN_ASSERT_EQ(buf, "test");
    SEQAN_ASSERT_EQ(value(reader), '\r');
}

SEQAN_DEFINE_TEST(test_stream_tokenizing_read_digits)
{
    using namespace seqan;
    typedef Stream<CharArray<char const *> > TStream;

    char const * FILE_CONTENTS = "3456.7890";
    TStream stream(FILE_CONTENTS, FILE_CONTENTS + strlen(FILE_CONTENTS));
    RecordReader<TStream, SinglePass<> > reader(stream);
    
    CharString buf;
    
    SEQAN_ASSERT_EQ(readDigits(buf, reader), 0);
    SEQAN_ASSERT_EQ(buf, "3456");

    goNext(reader);
    clear(buf);

    SEQAN_ASSERT_EQ(readDigits(buf, reader), EOF_BEFORE_SUCCESS);
    SEQAN_ASSERT_EQ(buf, "7890");
}

SEQAN_DEFINE_TEST(test_stream_tokenizing_read_alpha_nums)
{
    using namespace seqan;
    typedef Stream<CharArray<char const *> > TStream;

    char const * FILE_CONTENTS = "3a4b5c6.7d8e9f0";
    TStream stream(FILE_CONTENTS, FILE_CONTENTS + strlen(FILE_CONTENTS));
    RecordReader<TStream, SinglePass<> > reader(stream);
    
    CharString buf;
    
    SEQAN_ASSERT_EQ(readAlphaNums(buf, reader), 0);
    SEQAN_ASSERT_EQ(buf, "3a4b5c6");

    goNext(reader);
    clear(buf);

    SEQAN_ASSERT_EQ(readAlphaNums(buf, reader), EOF_BEFORE_SUCCESS);
    SEQAN_ASSERT_EQ(buf, "7d8e9f0");
}

SEQAN_DEFINE_TEST(test_stream_tokenizing_read_float)
{
    using namespace seqan;
    typedef Stream<CharArray<char const *> > TStream;

    char const * FILE_CONTENTS = "3.456 -3.456 +3.456 .1e03 0.1e+03 0.1E-3 0e00";
    TStream stream(FILE_CONTENTS, FILE_CONTENTS + strlen(FILE_CONTENTS));
    RecordReader<TStream, SinglePass<> > reader(stream);
    
    CharString buf;

    // Floating point number without sign.
    SEQAN_ASSERT_EQ(readFloat(buf, reader), 0);
    SEQAN_ASSERT_EQ(buf, "3.456");

    goNext(reader);
    clear(buf);

    // Floating point number with signs.
    SEQAN_ASSERT_EQ(readFloat(buf, reader), 0);
    SEQAN_ASSERT_EQ(buf, "-3.456");
    goNext(reader);
    clear(buf);
    SEQAN_ASSERT_EQ(readFloat(buf, reader), 0);
    SEQAN_ASSERT_EQ(buf, "+3.456");

    // Floating point number with scientific notation.
    goNext(reader);
    clear(buf);
    SEQAN_ASSERT_EQ(readFloat(buf, reader), 0);
    SEQAN_ASSERT_EQ(buf, ".1e03");
    goNext(reader);
    clear(buf);
    SEQAN_ASSERT_EQ(readFloat(buf, reader), 0);
    SEQAN_ASSERT_EQ(buf, "0.1e+03");
    goNext(reader);
    clear(buf);
    SEQAN_ASSERT_EQ(readFloat(buf, reader), 0);
    SEQAN_ASSERT_EQ(buf, "0.1E-3");
    goNext(reader);
    clear(buf);
    SEQAN_ASSERT_EQ(readFloat(buf, reader), EOF_BEFORE_SUCCESS);
    SEQAN_ASSERT_EQ(buf, "0e00");
}

SEQAN_DEFINE_TEST(test_stream_tokenizing_read_until_tab_or_line_break)
{
    using namespace seqan;
    typedef Stream<CharArray<char const *> > TStream;

    char const * FILE_CONTENTS = "this\tis\ta\ntest\n";
    TStream stream(FILE_CONTENTS, FILE_CONTENTS + strlen(FILE_CONTENTS));
    RecordReader<TStream, SinglePass<> > reader(stream);
    
    CharString buf;
    
    SEQAN_ASSERT_EQ(readUntilTabOrLineBreak(buf, reader), 0);
    SEQAN_ASSERT_EQ(buf, "this");
    SEQAN_ASSERT_EQ(value(reader), '\t');
    SEQAN_ASSERT_NOT(goNext(reader));
    clear(buf);
    SEQAN_ASSERT_EQ(readUntilTabOrLineBreak(buf, reader), 0);
    SEQAN_ASSERT_EQ(buf, "is");
    SEQAN_ASSERT_EQ(value(reader), '\t');
    SEQAN_ASSERT_NOT(goNext(reader));
    clear(buf);
    SEQAN_ASSERT_EQ(readUntilTabOrLineBreak(buf, reader), 0);
    SEQAN_ASSERT_EQ(buf, "a");
    SEQAN_ASSERT_EQ(value(reader), '\n');
    SEQAN_ASSERT_NOT(goNext(reader));
    clear(buf);
    SEQAN_ASSERT_EQ(readUntilTabOrLineBreak(buf, reader), 0);
    SEQAN_ASSERT_EQ(buf, "test");
    SEQAN_ASSERT_EQ(value(reader), '\n');
}

// readNChars
SEQAN_DEFINE_TEST(test_stream_tokenizing_readNChars)
{
    using namespace seqan;
    std::fstream * file = createFile();
    RecordReader<std::fstream, SinglePass<void> > reader(*file);

    CharString buf;

    SEQAN_ASSERT_EQ(readNChars(buf,reader, 16), 0);
    SEQAN_ASSERT_EQ(buf, "This is a string");

    // check EOF-handling
    clear(buf);
    SEQAN_ASSERT_EQ(readUntilChar(buf,reader, '*'), 0);
    SEQAN_ASSERT_EQ(value(reader), '*');

    clear(buf);
    SEQAN_ASSERT_EQ(readNChars(buf,reader, 3), EOF_BEFORE_SUCCESS);
    SEQAN_ASSERT_EQ(buf, "*");

    file->close();
    delete file;

}


// _readHelper function "readWhile"
SEQAN_DEFINE_TEST(test_stream_tokenizing_readWhile)
{
    using namespace seqan;
    std::fstream * file = createFile();
    RecordReader<std::fstream, SinglePass<void> > reader(*file);


    CharString buf;
    // skip to where we want to go
    SEQAN_ASSERT_EQ(readUntilChar(buf,reader, 'S'), 0);
    SEQAN_ASSERT_EQ(value(reader), 'S');

    clear(buf);
    SEQAN_ASSERT_EQ(readLetters(buf,reader), 0);
    SEQAN_ASSERT_EQ(buf, "SJAUDOF");

    clear(buf);
    SEQAN_ASSERT_EQ(readAlphaNums(buf,reader), 0);
    SEQAN_ASSERT_EQ(buf, "78456kapLP345LPL");


    // check EOF-handling
    clear(buf);
    SEQAN_ASSERT_EQ(readUntilChar(buf,reader, '*'), 0);
    SEQAN_ASSERT_EQ(value(reader), '*');
    clear(buf);
    SEQAN_ASSERT_EQ(readNChars(buf,reader, 1), 0);
    SEQAN_ASSERT_EQ(buf, "*");

    clear(buf);
    SEQAN_ASSERT_EQ(readLetters(buf,reader), EOF_BEFORE_SUCCESS);
    SEQAN_ASSERT_EQ(buf, "");

    file->close();
    delete file;


}

// _readHelper function with ignored characters
SEQAN_DEFINE_TEST(test_stream_tokenizing_readIgnoring)
{
    using namespace seqan;
    std::fstream * file = createFile();
    RecordReader<std::fstream, SinglePass<void> > reader(*file);

    CharString buf;
    Dna5String buf2;
    // skip to where we want to go
    SEQAN_ASSERT_EQ(readUntilChar(buf,reader, 'A'), 0);
    SEQAN_ASSERT_EQ(value(reader), 'A');

    clear(buf);
    SEQAN_ASSERT_EQ(readDna5IgnoringWhitespaces(buf2, reader), 0);
    SEQAN_ASSERT_EQ(buf2, "AAAACCCGGGTTTTCG");
    SEQAN_ASSERT_EQ(value(reader), '>');

    //double check whitespace implementation
    clear(buf2);
    SEQAN_ASSERT_EQ(readUntilChar(buf,reader, '|'), 0);
    SEQAN_ASSERT_EQ(value(reader), '|');
    clear(buf);
    SEQAN_ASSERT_EQ(readNChars(buf,reader, 1), 0);
    SEQAN_ASSERT_EQ(buf, "|");

    clear(buf);
    SEQAN_ASSERT_EQ(readDna5IgnoringWhitespaces(buf2, reader), 0);
    SEQAN_ASSERT_EQ(buf2, "TACGTACGTACGTACGATACGATCTNNACGT");
    SEQAN_ASSERT_EQ(value(reader), '*');

    // check EOF-handling
    clear(buf2);
    SEQAN_ASSERT_EQ(readNChars(buf,reader, 1), 0);
    SEQAN_ASSERT_EQ(buf, "*");

    clear(buf);
    SEQAN_ASSERT_EQ(readDna5IgnoringWhitespaces(buf, reader), EOF_BEFORE_SUCCESS);
    SEQAN_ASSERT_EQ(buf, "");

    file->close();
    delete file;

}

// readLine
SEQAN_DEFINE_TEST(test_stream_tokenizing_readLine)
{
    using namespace seqan;
    std::fstream * file = createFile();
    RecordReader<std::fstream, SinglePass<void> > reader(*file);


    CharString buf;
    // skip to where we want to go
    SEQAN_ASSERT_EQ(readUntilChar(buf,reader, 'A'), 0);
    SEQAN_ASSERT_EQ(value(reader), 'A');

    clear(buf);
    SEQAN_ASSERT_EQ(readLine(buf, reader), 0);
    SEQAN_ASSERT_EQ(buf, "AAAACCCGGGTTT");
    SEQAN_ASSERT_EQ(value(reader), ' ');

    clear(buf);
    SEQAN_ASSERT_EQ(readLine(buf, reader), 0);
    SEQAN_ASSERT_EQ(value(reader), '\a');

    clear(buf);
    SEQAN_ASSERT_EQ(readUntilChar(buf,reader, '1'), 0);

    clear(buf);
    SEQAN_ASSERT_EQ(readLine(buf, reader), 0);
    SEQAN_ASSERT_EQ(buf, "123gogogo ");
    SEQAN_ASSERT_EQ(value(reader), 'T');

    // check EOF-handling
    clear(buf);
    SEQAN_ASSERT_EQ(readUntilChar(buf,reader, '*'), 0);

    clear(buf);
    SEQAN_ASSERT_EQ(readLine(buf, reader), EOF_BEFORE_SUCCESS);
    SEQAN_ASSERT_EQ(buf, "*");

    file->close();
    delete file;
}

// _skipHelper function "skipUntil"
SEQAN_DEFINE_TEST(test_stream_tokenizing_skipUntil)
{
    using namespace seqan;
    std::fstream * file = createFile();
    RecordReader<std::fstream, SinglePass<void> > reader(*file);

    SEQAN_ASSERT_EQ(skipUntilChar(reader, '.'), 0);
    SEQAN_ASSERT_EQ(value(reader), '.');

    SEQAN_ASSERT_EQ(skipUntilBlank(reader), 0);
    SEQAN_ASSERT_EQ(value(reader), ' ');

    SEQAN_ASSERT_EQ(skipUntilChar(reader, 'A'), 0);
    SEQAN_ASSERT_EQ(value(reader), 'A');

    SEQAN_ASSERT_EQ(skipUntilWhitespace(reader), 0);
    SEQAN_ASSERT_EQ(value(reader), '\n');

    SEQAN_ASSERT_EQ(skipUntilChar(reader, 'A'), 0);
    SEQAN_ASSERT_EQ(value(reader), 'A');

    SEQAN_ASSERT_EQ(skipUntilBlank(reader), 0);
    SEQAN_ASSERT_EQ(value(reader), ' ');

    SEQAN_ASSERT_EQ(skipUntilChar(reader, 'A'), 0);
    SEQAN_ASSERT_EQ(value(reader), 'A');

    SEQAN_ASSERT_EQ(skipUntilWhitespace(reader), 0);
    SEQAN_ASSERT_EQ(value(reader), '\n');

    SEQAN_ASSERT_EQ(skipUntilGraph(reader), 0); // skip over \a
    SEQAN_ASSERT_EQ(value(reader), '1');

    // check EOF-handling
    SEQAN_ASSERT_EQ(skipUntilChar(reader, '*'), 0);
    SEQAN_ASSERT_EQ(value(reader), '*');

    SEQAN_ASSERT_EQ(skipUntilWhitespace(reader), EOF_BEFORE_SUCCESS);

    file->close();
    delete file;
}


// _skipHelper function "skipWhile"
SEQAN_DEFINE_TEST(test_stream_tokenizing_skipWhile)
{
    using namespace seqan;
    std::fstream * file = createFile();
    RecordReader<std::fstream, SinglePass<void> > reader(*file);

    SEQAN_ASSERT_EQ(skipUntilChar(reader, '\n'), 0);
    SEQAN_ASSERT_EQ(value(reader), '\n');

    SEQAN_ASSERT_EQ(skipWhitespaces(reader), 0);
    SEQAN_ASSERT_EQ(value(reader), 'T');

    SEQAN_ASSERT_EQ(skipUntilChar(reader, '\n'), 0);
    SEQAN_ASSERT_EQ(value(reader), '\n');

    SEQAN_ASSERT_EQ(skipWhitespaces(reader), 0); // dont skip over bell char
    SEQAN_ASSERT_EQ(value(reader), '\a');

    SEQAN_ASSERT_EQ(skipUntilChar(reader, 'V'), 0);
    SEQAN_ASSERT_EQ(value(reader), 'V');

    SEQAN_ASSERT_EQ(skipUntilBlank(reader), 0);
    SEQAN_ASSERT_EQ(value(reader), ' ');

    SEQAN_ASSERT_EQ(skipBlanks(reader), 0);
    SEQAN_ASSERT_EQ(value(reader), 't');

    // check EOF-handling
    SEQAN_ASSERT_EQ(skipUntilChar(reader, '*'), 0);
    SEQAN_ASSERT_EQ(value(reader), '*');

    CharString buf;
    SEQAN_ASSERT_EQ(readNChars(buf,reader, 1), 0);
    SEQAN_ASSERT_EQ(buf, "*");

    SEQAN_ASSERT_EQ(skipUntilBlank(reader), EOF_BEFORE_SUCCESS);

    file->close();
    delete file;
}

// skipLine
SEQAN_DEFINE_TEST(test_stream_tokenizing_skipLine)
{
    using namespace seqan;
    std::fstream * file = createFile();
    RecordReader<std::fstream, SinglePass<void> > reader(*file);

    SEQAN_ASSERT_EQ(skipLine(reader), 0);
    SEQAN_ASSERT_EQ(value(reader), ' ');

    SEQAN_ASSERT_EQ(skipLine(reader), 0);
    SEQAN_ASSERT_EQ(value(reader), '\a');

    SEQAN_ASSERT_EQ(skipLine(reader), 0);
    SEQAN_ASSERT_EQ(value(reader), 'T');

    SEQAN_ASSERT_EQ(skipLine(reader), 0);
    SEQAN_ASSERT_EQ(value(reader), 'e');

    SEQAN_ASSERT_EQ(skipLine(reader), 0);
    SEQAN_ASSERT_EQ(value(reader), 'e');

    SEQAN_ASSERT_EQ(skipLine(reader), 0);
    SEQAN_ASSERT_EQ(value(reader), '\n');

    SEQAN_ASSERT_EQ(skipLine(reader), 0);
    SEQAN_ASSERT_EQ(value(reader), '\v');

    SEQAN_ASSERT_EQ(skipLine(reader), 0);
    SEQAN_ASSERT_EQ(value(reader), '\r');

    SEQAN_ASSERT_EQ(skipLine(reader), 0);
    SEQAN_ASSERT_EQ(value(reader), 'A');

    SEQAN_ASSERT_EQ(skipLine(reader), 0);
    SEQAN_ASSERT_EQ(value(reader), '\n');

    SEQAN_ASSERT_EQ(skipLine(reader), 0);
    SEQAN_ASSERT_EQ(value(reader), 'A');

    SEQAN_ASSERT_EQ(skipLine(reader), EOF_BEFORE_SUCCESS);

    file->close();
    delete file;
}

// skipUntilString
SEQAN_DEFINE_TEST(test_stream_tokenizing_skipUntilString)
{
    using namespace seqan;
    std::fstream * file = createFile();
    RecordReader<std::fstream, SinglePass<void> > reader(*file);

    SEQAN_ASSERT_EQ(skipUntilString(reader, "foobar"), 0);
    SEQAN_ASSERT_EQ(value(reader), ' ');

    SEQAN_ASSERT_EQ(skipUntilString(reader, "\r\n"), 0);
    SEQAN_ASSERT_EQ(value(reader), 'T');

    SEQAN_ASSERT_EQ(skipUntilString(reader, "\r\n"), 0);
    SEQAN_ASSERT_EQ(value(reader), 'e');

    SEQAN_ASSERT_EQ(skipUntilString(reader, "stV"), 0);
    SEQAN_ASSERT_EQ(value(reader), ' ');

    SEQAN_ASSERT_EQ(skipUntilString(reader, "\n\nACGT"), 0);
    SEQAN_ASSERT_EQ(value(reader), '*');

    SEQAN_ASSERT_EQ(skipUntilString(reader, "notinhere"), EOF_BEFORE_SUCCESS);

    file->close();
    delete file;
}


// skipUntilLineBeginsWithChar
SEQAN_DEFINE_TEST(test_stream_tokenizing_skipUntilLineBeginsWithChar)
{
    using namespace seqan;
    std::fstream * file = createFile();
    RecordReader<std::fstream, SinglePass<void> > reader(*file);

    SEQAN_ASSERT_EQ(skipUntilLineBeginsWithChar(reader, 'T'), 0);
    SEQAN_ASSERT_EQ(value(reader), 'T');

    SEQAN_ASSERT_EQ(skipUntilLineBeginsWithChar(reader, '1'), 0);
    SEQAN_ASSERT_EQ(value(reader), '1');

    SEQAN_ASSERT_EQ(skipUntilLineBeginsWithChar(reader, 'Z'),
                    EOF_BEFORE_SUCCESS);

    file->close();
    delete file;
}

// skipUntilLineBeginsWithStr
SEQAN_DEFINE_TEST(test_stream_tokenizing_skipUntilLineBeginsWithStr)
{
    using namespace seqan;
    std::fstream * file = createFile();
    RecordReader<std::fstream, SinglePass<void> > reader(*file);

    SEQAN_ASSERT_EQ(skipUntilLineBeginsWithStr(reader, "TCG"), 0);
    SEQAN_ASSERT_EQ(value(reader), ' ');

    SEQAN_ASSERT_EQ(skipUntilLineBeginsWithStr(reader, "123"), 0);
    SEQAN_ASSERT_EQ(value(reader), 'g');

    SEQAN_ASSERT_EQ(skipUntilLineBeginsWithStr(reader, "ACGT"), 0);
    SEQAN_ASSERT_EQ(value(reader), 'A');

    SEQAN_ASSERT_EQ(skipUntilLineBeginsWithStr(reader, "ACGT"), 0);
    SEQAN_ASSERT_EQ(value(reader), '*');

    SEQAN_ASSERT_EQ(skipUntilLineBeginsWithStr(reader, "ACGT"),
                    EOF_BEFORE_SUCCESS);

    file->close();
    delete file;
}

// skipUntilLineBeginsWithOneOfStr
SEQAN_DEFINE_TEST(test_stream_tokenizing_skipUntilLineBeginsWithOneCharOfStr)
{
    using namespace seqan;
    std::fstream * file = createFile();
    RecordReader<std::fstream, SinglePass<void> > reader(*file);

    SEQAN_ASSERT_EQ(skipUntilLineBeginsWithOneCharOfStr(reader, "FGT"), 0);
    SEQAN_ASSERT_EQ(value(reader), 'T');

    SEQAN_ASSERT_EQ(skipUntilLineBeginsWithOneCharOfStr(reader, "987654321"),
                    0);
    SEQAN_ASSERT_EQ(value(reader), '1');

    SEQAN_ASSERT_EQ(skipUntilLineBeginsWithOneCharOfStr(reader, "abcde"), 0);
    SEQAN_ASSERT_EQ(value(reader), 'e');

    SEQAN_ASSERT_EQ(skipUntilLineBeginsWithOneCharOfStr(reader, "UVW"),
                    EOF_BEFORE_SUCCESS);

    file->close();
    delete file;
}

#endif // ndef TEST_STREAM_TEST_STREAM_TOKENIZING_H_
