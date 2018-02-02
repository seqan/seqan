// ==========================================================================
//                   NGS: Regions of Interest Analysis
// ==========================================================================
// Copyright (c) 2012-2018, Bernd Jagla, Institut Pasteur
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Simple and fast creation of PNG files.
//
// This code is based on the Python class PNGCanvas by Rui Carmo.  You can
// find the code with some documentation on his site:
//
//   http://taoofmac.com/space/projects/PNGCanvas
// ==========================================================================

#ifndef SANDBOX_JAGLA_APPS_NGS_ROI_PNG_CANVAS_H_
#define SANDBOX_JAGLA_APPS_NGS_ROI_PNG_CANVAS_H_

#include <seqan/basic.h>  // for uint32_t etc.

#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <zlib.h>

// ============================================================================
// Forwards
// ============================================================================

// ----------------------------------------------------------------------------
// Class PngColor
// ----------------------------------------------------------------------------

// A color for the PngCanvas.

class PngColor
{
public:
    // Red, gree, blue, and alpha.
    unsigned char r, g, b, a;

    // Default constructed colors are transparend black.
    PngColor() : r(0), g(0), b(0), a(0)
    {}

    // Construct directly from values, optionally set alpha to opaque.
    PngColor(unsigned char r, unsigned char g, unsigned char b, unsigned char a = 0xff) :
            r(r), g(g), b(b), a(a)
    {}

    // Alpha-blend the two color using the alpha given by c2.
    PngColor blend(PngColor other) const
    {
        return PngColor((r * (0xff - other.a) + other.r * other.a) >> 8,
                        (g * (0xff - other.a) + other.g * other.a) >> 8,
                        (b * (0xff - other.a) + other.b * other.a) >> 8);
    }

    // Compute new alpha given a 0-0xff intensity.
    PngColor intensity(unsigned char i) const
    {
        return PngColor(r, g, b, (a * i) >> 8);
    }

    // Compute perceptive grayscale value.
    unsigned char grayScale() const
    {
        return static_cast<unsigned char>(r * 0.3 + g * 0.59 + b * 0.11);
    }

    // Return the color white.
    static PngColor WHITE()
    {
        return PngColor(0xff, 0xff, 0xff, 0xff);
    }

    // Return the color black.
    static PngColor BLACK()
    {
        return PngColor(0, 0, 0, 0xff);
    }
};

// ----------------------------------------------------------------------------
// Class PngIhdrChunk
// ----------------------------------------------------------------------------

struct PngIhdrChunk
{
    uint32_t width;
    uint32_t height;
    char bitDepth;
    char colorType;
    char compressionMethod;
    char filterMethod;
    char interlaceMethod;

    PngIhdrChunk() : width(0), height(0), bitDepth(0), colorType(0), compressionMethod(0), filterMethod(0), interlaceMethod(0)
    {}

    // Returns pointer to the payload.
    unsigned char const * payload() const
    {
        return reinterpret_cast<unsigned char const *>(this);
    }

    // Returns number of bytes.
    unsigned size() const
    {
        return 13;
    }

    // Returns the chunk type.
    char const * type() const
    {
        return "IHDR";
    }
};

// ----------------------------------------------------------------------------
// Class PngIdatChunk
// ----------------------------------------------------------------------------

struct PngIdatChunk
{
    std::basic_string<unsigned char> data;

    // Returns pointer to the payload.
    unsigned char const * payload() const
    {
        return reinterpret_cast<unsigned char const *>(&data[0]);
    }

    unsigned size() const
    {
        return data.size();
    }

    char const * type() const
    {
        return "IDAT";
    }
};

// ----------------------------------------------------------------------------
// Class PngIendChunk
// ----------------------------------------------------------------------------

struct PngIendChunk
{
    // Returns pointer to the payload.
    unsigned char const * payload() const
    {
        return 0;
    }

    unsigned size() const
    {
        return 0;
    }

    char const * type() const
    {
        return "IEND";
    }
};

// ----------------------------------------------------------------------------
// Class PngCanvas
// ----------------------------------------------------------------------------

class PngCanvas
{
public:
    // We will operate on a raw vector of PngColor object.
    typedef std::vector<PngColor> TCanvas;

    // The canvas width.
    int width;
    // The canvas height.
    int height;
    // The current background color
    PngColor bgColor;
    // The current foreground color.
    PngColor color;
    // The array of colors that form the canvas and that we will operate upon.
    TCanvas canvas;

    // Default constructor.  The canvas is empty.
    PngCanvas() : width(0), height(0), bgColor(PngColor::WHITE()), color(PngColor::BLACK())
    {}

    // Constructor.
    PngCanvas(int width, int height) :
            width(width), height(height), bgColor(PngColor::WHITE()), color(PngColor::BLACK())
    {
        canvas.resize(width * height, bgColor);
    }

    PngCanvas(int width, int height, PngColor bgColor):
            width(width), height(height), bgColor(bgColor), color(PngColor::BLACK())
    {
        canvas.resize(width * height, bgColor);
    }

    PngCanvas(int width, int height, PngColor bgColor, PngColor color) :
            width(width), height(height), bgColor(bgColor), color(color)
    {
        canvas.resize(width * height, bgColor);
    }

    // Return offset in canvas for the given point.
    inline int _offset(int x, int y) const
    {
        return y * width + x;
    }

    // Draw a point in the current foreground color.
    inline void point(int x, int y, PngColor col)
    {
        if (x < 0 || y < 0 || x > width - 1 || y > height - 1)
            return;  // Make sure that we are on canvas.
        int o = _offset(x, y);
        canvas[o] = canvas[o].blend(col);
    }

    inline void point(int x, int y)
    {
        return point(x, y, color);
    }

    // Ensures that x0 <= x1, y0 <= y1, and swaps values if necessary.
    inline void _rectHelper(int & x0, int & y0, int & x1, int & y1) const
    {
        if (x0 > x1)
            std::swap(x0, x1);
        if (y0 > y1)
            std::swap(y0, y1);

    }

    // Draw vertical gradient from start to end color in rectangle (x0, y0, x1, y1).
    void verticalGradient(int x0, int y0, int x1, int y1, PngColor start, PngColor end)
    {
        _rectHelper(x0, y0, x1, y1);
        std::vector<PngColor> grad;
        gradientList(grad, start, end, y1 - y0);
        for (int x = x0; x <= x1; ++x)
            for (int y = y0; y <= y1; ++y)
                point(x, y, grad[y - y0]);
    }

    // Draw rectangle (x0, y0, x1, y1) in current color.
    void rectangle(int x0, int y0, int x1, int y1)
    {
        _rectHelper(x0, y0, x1, y1);
        std::vector<std::pair<int, int> > points;
        points.push_back(std::make_pair(x0, y0));
        points.push_back(std::make_pair(x1, y0));
        points.push_back(std::make_pair(x1, y1));
        points.push_back(std::make_pair(x0, y1));
        points.push_back(std::make_pair(x0, y0));
        polyline(points);
    }

    // Fill rectangle (x0, y0, x1, y1) in current foreground color.
    void filledRectangle(int x0, int y0, int x1, int y1)
    {
        _rectHelper(x0, y0, x1, y1);
        for (int x = x0; x <= x1; ++x)
            for (int y = y0; y <= y1; ++y)
                point(x, y, color);
    }

    // Copy rectangle (x0, y0, x1, y1) from current canvas to destination canvas with offset dx, dy.
    void copyRect(PngCanvas & other, int x0, int y0, int x1, int y1, int dx, int dy) const
    {
        _rectHelper(x0, y0, x1, y1);
        for (int x = x0; x <= x1; ++x)
            for (int y = y0; y <= y1; ++y)
            {
                int d = other._offset(dx + x - x0, dy + y - y0);
                int o = _offset(x, y);
                other.canvas[d] = canvas[o];
            }
    }

    // Blend rectangle (x0, y0, x1, y1) from current canvas to destination canvas with offset dx, dy and given alpha value.
    void blendRect(PngCanvas & other, int x0, int y0, int x1, int y1, int dx, int dy, unsigned char alpha = 0xff) const
    {
        _rectHelper(x0, y0, x1, y1);
        for (int x = x0; x <= x1; ++x)
            for (int y = y0; y <= y1; ++y)
            {
                int o = _offset(x, y);
                PngColor rgba = canvas[o];
                rgba.a = alpha;
                other.point(dx + x - x0, dy + y - y0, rgba);
            }
    }

    // Draw a line using Xiaolin Wu's antialiasing technique.
    void line(int x0, int y0, int x1, int y1)
    {
        // std::cerr << "line(" << x0 << ", " << y0 << ", " << x1 << ", " << y1 << ")\n";
        // Clean parameters.
        if (y0 > y1)
        {
            std::swap(y0, y1);
            std::swap(x0, x1);
        }
        int dx = x1 - x0;
        int sx = (dx < 0) ? -1 : 1;
        dx *= sx;
        int dy = y1 - y0;

        // Handle the "easy" cases.
        if (dy == 0)
        {
            if (x0 > x1)
            {
                std::swap(x0, x1);
                sx = -sx;
            }
            for (int x = x0; x < x1; x += sx)
                point(x, y0);
            return;
        }
        if (dx == 0)
        {
            for (int y = y0; y < y1; ++y)
                point(x1, y);
            point(x1, y1);
            return;
        }
        if (dx == dy)
        {
            if ((sx < 0) != (x0 > x1))
                std::swap(x0, x1);
            for (int x = x0; x < x1; x += sx, ++y0)
                point(x, y0);
            return;
        }

        // Main loop.
        point(x0, y0);
        int eAcc = 0;
        if (dy > dx)  // vertical displacement, then early exit
        {
            int e = (dx << 16) / dy;
            for (int i = y0; i < y1 - 1; ++i)
            {
                int eAccTemp = eAcc;
                eAcc = (eAcc + e) & 0xffff;
                if (eAcc <= eAccTemp)
                    x0 = x0 + sx;
                int w = 0xff - (eAcc >> 8);
                point(x0, y0, color.intensity(w));
                y0 = y0 + 1;
                point(x0 + sx, y0, color.intensity(0xff - w));
            }
            point(x1, y1);
            return;
        }

        // Handle horizontal displacement.
        int e = (dy << 16) / dx;
        for (int i = x0; i < x1 - sx; i += sx)
        {
            int eAccTemp = eAcc;
            eAcc = (eAcc + e) & 0xffff;
            if (eAcc <= eAccTemp)
                y0 = y0 + 1;
            int w = 0xff - (eAcc >> 8);
            point(x0, y0, color.intensity(w));
            x0 = x0 + sx;
            point(x0, y0 + 1, color.intensity(0xff - w));
        }
        point(x1, y1);
    }

    // Draw polyline.
    void polyline(std::vector<std::pair<int, int> > & points)
    {
        for (unsigned i = 0; i + 1 < points.size(); ++i)
            line(points[i].first, points[i].second,
                 points[i + 1].first, points[i + 1].second);
    }

    // Write to file.
    void write(char const * filename)
    {
        unsigned char const SIGNATURE[8] = {137, 80, 78, 71, 13, 10, 26, 10};

        // Open file.
        std::fstream out(filename, std::ios::binary | std::ios::out);
        // Write header.
        out.write(reinterpret_cast<char const *>(&SIGNATURE[0]), 8);

        // Write IHDR chunk.
        PngIhdrChunk ihdrChunk;
        ihdrChunk.width = swapUInt32(width);
        ihdrChunk.height = swapUInt32(height);
        ihdrChunk.bitDepth = 8;
        ihdrChunk.colorType = 6;
        ihdrChunk.compressionMethod = 0;
        ihdrChunk.filterMethod = 0;
        ihdrChunk.interlaceMethod = 0;
        _writeChunk(out, ihdrChunk);

        // Write IDAT chunk.
        // TODO(holtgrew): Could be improved by not holding everything in memory all the time.
        PngIdatChunk idatChunk;
        idatChunk.data.resize(height + 4 * canvas.size());  // same as input, should be enough

        z_stream zs;
        int status = 0;
        zs.zalloc = NULL;
        zs.zfree = NULL;
        zs.next_in = NULL;
        zs.next_out = &idatChunk.data[0];
        zs.avail_out = idatChunk.data.size();

        // Initialize stream.
        status = deflateInit(&zs, 9);
        if (status != Z_OK)
            throw std::runtime_error("delateInit() failed!");

        // Write data, line-wise.
        for (unsigned y = 0; y < (unsigned)height; ++y)
        {
            // std::cerr << "y = " << y << " height == " << height << "\n";
            // Write '\0'.
            char nofilter[1] = {'\0'};  // filter 0 (no filter)
            zs.next_in = reinterpret_cast<Bytef *>(&nofilter);
            zs.avail_in = 1;
            status = deflate(&zs, 0);
            // std::cerr << "status = " << status << "\n";
            if (status != Z_OK)
                throw std::runtime_error("Writing filter flag failed!");

            // Write sweep line data.
            zs.next_in = reinterpret_cast<Bytef *>(&canvas[_offset(0, y)]);
            zs.avail_in = 4 * width;
            status = deflate(&zs, (y + 1 == (unsigned)height) ? Z_FINISH : 0);
            // std::cerr << "status = " << status << "\n";
            if (y + 1 != (unsigned)height)
            {
                if (status != Z_OK)
                    throw std::runtime_error("Compressing data failed!");
            }
            else
            {
                if (status != Z_STREAM_END)
                    throw std::runtime_error("Last step of compressing data failed!");
            }
        }

        // Finish compression and write chunk.
        status = deflateEnd(&zs);
        if (status != Z_OK)
            throw std::runtime_error("Could not finish compression!");
        idatChunk.data.resize(zs.total_out);
        _writeChunk(out, idatChunk);

        // Write IEND chunk.
        PngIendChunk iendChunk;
        _writeChunk(out, iendChunk);
    }

    uint32_t swapUInt32(uint32_t val)
    {
        val = ((val << 8) & 0xFF00FF00 ) | ((val >> 8) & 0xFF00FF );
        return (val << 16) | (val >> 16);
    }

    int32_t swapInt32(int32_t val)
    {
        val = ((val << 8) & 0xFF00FF00) | ((val >> 8) & 0xFF00FF );
        return (val << 16) | ((val >> 16) & 0xFFFF);
    }

    template <typename TChunk>
    void _writeChunk(std::fstream & out, TChunk const & chunk)
    {
        // Write payload length.
        uint32_t size = swapUInt32(chunk.size());
        out.write(reinterpret_cast<char const *>(&size), 4);

        // Write chunk type.
        out.write(chunk.type(), 4);

        // Write payload.
        out.write(reinterpret_cast<char const *>(chunk.payload()), chunk.size());

        // Compute and write CRC32.
        uint32_t crc = crc32(0, NULL, 0);
        crc = crc32(crc, reinterpret_cast<Bytef const *>(chunk.type()), 4);
        if (chunk.size() > 0u)
            crc = crc32(crc, reinterpret_cast<Bytef const *>(reinterpret_cast<char const *>(chunk.payload())), chunk.size());
        crc = swapUInt32(crc);
        out.write(reinterpret_cast<char const *>(&crc), 4);
    }

    // Load from file.
    void read(char const * /*filename*/)
    {
    }

    // ???
    void defilter()
    {
    }

    // ???
    void chunks()
    {
    }

    // Compute colors for gradients.
    void gradientList(std::vector<PngColor> & result, PngColor start, PngColor end, int steps)
    {
        // PngColor delta(end.r - start.r, end.g - start.g, end.b - start.b, end.a - start.a);
        int deltaR = end.r - start.r;
        int deltaG = end.g - start.g;
        int deltaB = end.b - start.b;
        int deltaA = end.a - start.a;
        for (int i = 0; i <= steps; ++i)
            result.push_back(PngColor(start.r + (deltaR * i) / steps,
                                      start.g + (deltaG * i) / steps,
                                      start.b + (deltaB * i) / steps,
                                      start.a + (deltaA * i) / steps));
    }
};

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_JAGLA_APPS_NGS_ROI_PNG_CANVAS_H_
