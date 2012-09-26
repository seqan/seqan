// FRAGMENT(includes)
#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/graph_types.h>

using namespace seqan;

// FRAGMENT(tags-structs)
struct Newick_;
typedef Tag<Newick_> Newick;

struct NewickBranchLabel
{
    bool isDistanceSet;
    double distance;
    
    NewickBranchLabel() : isDistanceSet(false), distance(0)
    {}
};

// FRAGMENT(read-float)
template <typename TCharSeq, typename TStream, typename TPass>
int _readExponentialPart(TCharSeq & buffer,
                         RecordReader<TStream, TPass> & reader)
{
    // Check preconditions.
    SEQAN_ASSERT_NOT(atEnd(reader));
    SEQAN_ASSERT(value(reader) == 'e' || value(reader) == 'E');
    
    // Read 'e' or 'E';
    appendValue(buffer, value(reader));
    if (goNext(reader))
        return EOF_BEFORE_SUCCESS;
    // Possibly read '+' or '-'.
    if (value(reader) == '+' || value(reader) == '-')
    {
        appendValue(buffer, value(reader));
        if (goNext(reader))
            return EOF_BEFORE_SUCCESS;
    }
    // Read digits.
    if (!isdigit(value(reader)))
        return 1;  // Should have been a digit!
    return readDigits(buffer, reader);
}

template <typename TCharSeq, typename TStream, typename TPass>
int readFloadLiteral(TCharSeq & buffer,
                     RecordReader<TStream, TPass> & reader)
{
    if (atEnd(reader))
        return EOF_BEFORE_SUCCESS;  // Empty field.
        
    // The EBNF for floating point integers is as follows:
    //
    // exponent-indicator     = e | E 
    // exponent-part          = exponent-indicator [+|-]digits
    // floating-point-literal = digits exponent-part
    //                        | digits.[digits][exponent-part]
    //                        | .digits[exponent-part]

    // Read one leading sign if it is there.
    if (value(reader) == '-' || value(reader) == '+')
    {
        appendValue(buffer, value(reader));
        if (goNext(reader))
            return EOF_BEFORE_SUCCESS;  // At end after leading sign.
    }
    
    // Digits or dot?
    if (value(reader) == '.')
    {
        // Dot
        appendValue(buffer, '.');
        if (goNext(reader))
            return EOF_BEFORE_SUCCESS;
        if (!isdigit(value(reader)))
            return 1;  // Invalid format, >= 1 digit have to follow.
        int res = readDigits(buffer, reader);
        if (res != 0)
            return res;  // Error reading digits.
        // Optionally read exponential part.
        if (atEnd(reader))
            return 0;
        if (value(reader) == 'e' || value(reader) == 'E')
            return _readExponentialPart(buffer, reader);
    }
    else
    {
        // Digits
        if (!isdigit(value(reader)))
            return 1;  // >= 1 digit required!
        int res = readDigits(buffer, reader);
        if (res != 0)
            return res;  // Error reading digits.
        if (atEnd(reader))  // Stop if no more data.
            return 0;
        if (value(reader) == '.')
        {
            appendValue(buffer, '.');
            if (goNext(reader))
                return 0;  // End of field.
            if (isdigit(value(reader)))
            {
                res = readDigits(buffer, reader);
                if (res != 0)
                    return res;  // Error reading digits.
            }
            // Optionally read exponential part.
            if (atEnd(reader))
                return 0;
            if (value(reader) == 'e' || value(reader) == 'E')
                return _readExponentialPart(buffer, reader);
        }
        else if (value(reader) == 'e' || value(reader) == 'E')
        {
            return _readExponentialPart(buffer, reader);
        }
    }
    
    return 0;
}

// FRAGMENT(reading)
template <typename TTree, typename TRecordReader, typename TVertexDescriptor>
int _readNewickTree(TTree & tree,
                    String<CharString> & vertexLabels,
                    String<NewickBranchLabel> & branchLabels,
                    TRecordReader & reader,
                    TVertexDescriptor v)
{
    if (atEnd(reader))
        return EOF_BEFORE_SUCCESS;
    CharString buffer;
    int res = 0;
    
#define SKIP_WHITESPACE                            \
    do                                             \
    {                                              \
        int res = skipWhitespaces(reader);         \
        if (res != 0 && res != EOF_BEFORE_SUCCESS) \
            return res;                            \
    } while(false)
    
    if (value(reader) == '(')  // CHILDREN
    {
        if (goNext(reader))
            return EOF_BEFORE_SUCCESS;
        // children
        bool first = true;
        while (true)
        {
            SKIP_WHITESPACE;

            // Skip leading comma.
            if (!first)
            {
                res = skipChar(reader, ',');
                if (res != 0)
                    return res;
            }
            first = false;
            
            SKIP_WHITESPACE;
            
            // Read child.
            TVertexDescriptor x = addChild(tree, v);
            resizeVertexMap(tree, vertexLabels);
            resizeVertexMap(tree, branchLabels);
            res = _readNewickTree(tree, vertexLabels, branchLabels, reader, x);
            if (res != 0)
                return res;
            
            SKIP_WHITESPACE;
            
            // Exit loop.
            if (value(reader) == ')')
                break;
        }
        res = skipChar(reader, ')');
        if (res != 0)
            return res;  // Could not close child list.
        SKIP_WHITESPACE;
    }
    if (value(reader) != ':')  // LABEL
    {
        SKIP_WHITESPACE;
        clear(buffer);
        if (value(reader) == '\'')
        {
            // Read quoted label.
            if (goNext(reader))
                return EOF_BEFORE_SUCCESS;
            while (!atEnd(reader))
            {
                char c = value(reader);
                if (c == '\'')  // Possibly break, if not "''".
                {
                    if (goNext(reader))
                        break;
                    if (value(reader) != '\'')
                        break;
                }
                appendValue(buffer, value(reader));
                if (goNext(reader))
                    return 1;
            }
        }
        else
        {
            // Read unquoted label.
            while (!atEnd(reader))
            {
                char c = value(reader);
                if (isblank(c) || c == '(' || c == ')' || c == '[' ||
                    c == ']' || c == '\'' || c == '.' || c == ';' ||
                    c == ',' || c == ':')
                    break;
                appendValue(buffer, value(reader));
                if (goNext(reader))
                    return 1;
            }
        }
        assignProperty(vertexLabels, v, buffer);
        SKIP_WHITESPACE;
    }
    if (value(reader) == ':')  // DISTANCE
    {
        skipChar(reader, ':');
        SKIP_WHITESPACE;
        clear(buffer);
        res = readFloadLiteral(buffer, reader);
        if (res != 0)
            return res;  // Invalid distance.
        property(branchLabels, v).isDistanceSet = true;
        property(branchLabels, v).distance = lexicalCast<double>(buffer);
        SKIP_WHITESPACE;
    }
    return 0;
}

template <typename TStream, typename TSpec>
int read2(String<Graph<Tree<> > > & forest,
          String<String<CharString> > & vertexLabels,
          String<String<NewickBranchLabel> > & branchLabels,
          RecordReader<TStream, SinglePass<TSpec> > & reader,
          Newick const & /*tag*/)
{
    typedef Graph<Tree<> > TTree;
    typedef typename VertexDescriptor<TTree>::Type TVertexDescriptor;
    int res = 0;

    SKIP_WHITESPACE;
    
    // Read forest.
    while (!atEnd(reader))
    {
        // Allocate graph and maps.
        resize(forest, length(forest) + 1);
        resize(vertexLabels, length(vertexLabels) + 1);
        resize(branchLabels, length(branchLabels) + 1);
        // Allocate root.
        createRoot(back(forest));
        TVertexDescriptor v = root(back(forest));
        resizeVertexMap(back(forest), back(vertexLabels));
        resizeVertexMap(back(forest), back(branchLabels));
        // Read tree.
        res = _readNewickTree(back(forest), back(vertexLabels), back(branchLabels), reader, v);
        if (res != 0)
            return res;
        // Skip trailing semicolon, must be there.
        res = skipChar(reader, ';');
        if (res != 0)
            return res;
        SKIP_WHITESPACE;
    }

#undef SKIP_WHITESPACE

    return 0;
}

// FRAGMENT(writing)
template <typename TStream, typename TTree, typename TVertexDescriptor, typename TVertexLabels,
          typename TBranchLabels>
int _writeNewickRecurse(TStream & stream, TTree & tree, TVertexDescriptor v,
                        TVertexLabels & vertexLabels, TBranchLabels & branchLabels)
{
    if (numChildren(tree, v) > 0u)
    {
        int res = streamPut(stream, '(');
        if (res != 0)
            return res;
        
        typename Iterator<TTree, OutEdgeIterator>::Type it(tree, v);
        bool first = true;
        for (; !atEnd(it); goNext(it))
        {
            if (!first)
            {
                res = streamPut(stream, ',');
                if (res != 0)
                    return res;
            }
            first = false;
            res = _writeNewickRecurse(stream, tree, targetVertex(it), vertexLabels, branchLabels);
            if (res != 0)
                return res;
        }
        
        res = streamPut(stream, ')');
        if (res != 0)
            return res;
    }
    // Write label if any, quoted if required.
    if (length(property(vertexLabels, v)) > 0u)
    {
        bool needsQuoting = false;
        CharString const & label = property(vertexLabels, v);
        typename Iterator<CharString const, Rooted>::Type it = begin(label, Rooted());
        for (; !atEnd(it); ++it)
        {
            if (isblank(*it) || *it == ',' || *it == ';' || *it == '.' ||
                *it == '\'' || *it == '[' || *it == ']' || *it == '(' ||
                *it == ')')
            {
                needsQuoting = true;
                break;
            }
        }
        if (needsQuoting)
        {
            int res = streamPut(stream, '\'');
            if (res != 0)
                return res;
            it = begin(label, Rooted());
            for (; !atEnd(it); ++it)
            {
                if (*it == '\'')
                {
                    res = streamPut(stream, "''");
                    if (res != 0)
                        return res;
                }
                else
                {
                    res = streamPut(stream, *it);
                    if (res != 0)
                        return res;
                }
            }
            res = streamPut(stream, '\'');
            if (res != 0)
                return res;
        }
        else
        {
            int res = streamPut(stream, label);
            if (res != 0)
                return res;
        }
    }
    // Write branch length if any is given.
    if (property(branchLabels, v).isDistanceSet)
    {
        int res = streamPut(stream, ':');
        if (res != 0)
            return res;
        res = streamPut(stream, property(branchLabels, v).distance);
        if (res != 0)
            return res;
    }
    return 0;
}

template <typename TStream>
inline int write2(TStream & stream,
                  Graph<Tree<> > & tree,
                  String<CharString> & vertexLabels,
                  String<NewickBranchLabel> & branchLabels,
                  Newick const & /*tag*/)
{
    // Write <tree>;.
    int res = _writeNewickRecurse(stream, tree, getRoot(tree), vertexLabels, branchLabels);
    if (res != 0)
        return res;
    return streamPut(stream, ';');
}

// FRAGMENT(main)
int main(int argc, char const ** argv)
{
    // Handle arguments, open file.
    if (argc != 2)
    {
        std::cerr << "Incorrect argument count!" << std::endl;
        std::cerr << "USAGE: tutorial_parse_newick INPUT.txt" << std::endl;
        return 1;
    }
    std::fstream stream(argv[1], std::ios::binary | std::ios::in);
    if (!stream.good())
    {
        std::cerr << "Could not open file " << argv[1] << std::endl;
        return 1;
    }
    RecordReader<std::fstream, SinglePass<> > reader(stream);
    
    // Load forest.
    String<Graph<Tree<> > > forest;
    String<String<CharString> > vertexLabels;
    String<String<NewickBranchLabel> > branchLabels;
    int res = read2(forest, vertexLabels, branchLabels, reader, Newick());
    if (res != 0)
    {
        std::cerr << "Could not read Newick file!" << std::endl;
        return res;
    }
    
    // Dump forests.
    for (unsigned i = 0; i < length(forest); ++i)
    {
        res = write2(std::cout, forest[i], vertexLabels[i], branchLabels[i], Newick());
        std::cout << "\n";
        if (res != 0)
        {
            std::cerr << "Error writing to stdout!" << std::endl;
            return 1;
        }
    }
    
    return 0;
}
