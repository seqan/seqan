// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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

#ifndef SEQAN_CORE_INCLUDE_ARG_PARSE_ARGUMENT_PARSER_H_
#define SEQAN_CORE_INCLUDE_ARG_PARSE_ARGUMENT_PARSER_H_

#include <seqan/map.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include <seqan/arg_parse/arg_parse_type_support.h>
#include <seqan/arg_parse/arg_parse_argument.h>
#include <seqan/arg_parse/arg_parse_option.h>

#include <seqan/arg_parse/tool_doc.h>

#include <iostream>
#include <string>
#include <vector>
#include <map>

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================

// friend declaration to make addOption() and hideOption() available
// in ArgumentParser::init()
class ArgumentParser;
class ArgParseOption;
void addOption(ArgumentParser & me, ArgParseOption const & opt);
void hideOption(ArgumentParser & me, std::string const & name, bool hide);
void setValidValues(ArgumentParser & me, std::string const & name, std::string const & values);

// Required in addOption() and addArgument().
inline void hideOption(ArgumentParser & me, std::string const & name, bool hide = true);
inline ArgParseOption & getOption(ArgumentParser & me, std::string const & name);
inline void setValidValues(ArgumentParser & me, std::string const & name, std::vector<std::string> const & values);
inline ArgParseArgument & getArgument(ArgumentParser & me, unsigned position);

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

/*!
 * @class ArgumentParser
 * @headerfile <seqan/arg_parse.h>
 * @brief Parse the command line.
 *
 * @signature class ArgumentParser;
 *
 * Options are stored as @link ArgParseOption @endlink and @link ArgParseArgument @endlink objects.
 *
 * @section Remarks
 *
 * See the documentation of @link ToolDoc @endlink on how to format text.  Wherever possible, formatting is added
 * automatically for you.  You have to use formatting in the following places: (1) usage lines, (2) option help texts,
 * (3) description and additional text sections.
 * 
 * @section Examples
 *
 * The following gives a simple example of how to use the ArgumentParser class.
 *
 * @include demos/arg_parse/argument_parser.cpp
 *
 * @code{.console}
 * $ demo_arg_parse_argument_parser in.fa out.txt --id 0
 * Built target seqan_core
 * Built target demo_arg_parse
 * Verbose:     off
 * Identity:    0
 * Input-File:  in.fa
 * Output-File: out.txt
 * @endcode
 *
 * @see ArgParseArgument
 * @see ArgParseOption
 * @see ToolDoc
 */

/*!
 * @fn ArgumentParser::ArgumentParser
 * @brief Constructor
 *
 * @signature ArgumentParser::ArgumentParser([appName]);
 *
 * @param[in] appName The name of the application (<tt>std::string</tt>), defaults to <tt>argv[0]</tt>.
 */

/*!
 * @enum ArgumentParser::ParseResult
 * @brief Argument parsing result.
 *
 * @signature enum ArgumentParser::ParseResult;
 *
 * @val ArgumentParser::ParseResult ArgumentParser::PARSE_OK;
 * @brief Parsing the program's arguments was successful and no builtin command was triggered.
 *
 * @val ArgumentParser::ParseResult ArgumentParser::PARSE_ERROR;
 * @brief There were errors parsing the arguments.
 *
 * @val ArgumentParser::ParseResult ArgumentParser::PARSE_HELP;
 * @brief Parsing was successful, built-in <tt>--help</tt> option was used.
 *
 * @val ArgumentParser::ParseResult ArgumentParser::PARSE_VERSION;
 * @brief Parsing was successful, built-in <tt>--version</tt> option was used.
 *
 * @val ArgumentParser::ParseResult ArgumentParser::PARSE_WRITE_CTD;
 * @brief Parsing was successful, built-in <tt>--write-ctd</tt> option was used.
 *
 * @val ArgumentParser::ParseResult ArgumentParser::PARSE_EXPORT_HELP;
 * @brief Parsing was successful, built-in <tt>--export-help</tt> option was used.
 */

/**
.Class.ArgumentParser
..cat:Miscellaneous
..summary:Stores multiple @Class.ArgParseOption@ objects and parses the command line arguments for these options.
..signature:ArgumentParser
..include:seqan/arg_parse.h
..remarks:
See the documentation of @Class.ToolDoc@ on how to format text.
Where possible, formatting is added automatically for you.
You have to use formatting in the following places: (1) usage lines, (2) option help texts, (3) description and additional text sections.
..example.text:
The following gives a simple example of how to use the @Class.ArgumentParser@.
..example.code:
ArgumentParser parser("alf");
setShortDescription(parser, "Alignment free sequence comparison");
setVersion(parser, "1.0");
setDate(parser, "Jan 2010");

addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-i\\fP \\fIIN\\fP \\fB-o\\fP \\fIOUT\\fP");

addDescription(parser,
               "ALF can be used to calculate the pairwise similarity of sequences "
               "using alignment-free methods. All methods which are implemented are "
               "based on k-mer counts.");

addOption(parser, ArgParseOption("i", "inputFile", "Name of the multi-FASTA input.",
                                 ArgParseArgument(ArgParseArgument::INPUTFILE, "IN")));
setRequired(parser, "i");

addOption(parser, ArgParseOption("o", "outputFile", "Name of the multi-FASTA input.",
                                 ArgParseArgument(ArgParseArgument::OUTPUTFILE, "OUT")));
setRequired(parser, "o");

addTextSection(parser, "See Also");
addText(parser, "http://www.seqan.de/projects/alf");
..see:Class.ToolDoc
..see:Class.ArgParseArgument
..see:Class.ArgParseOption

.Memfunc.ArgumentParser#ArgumentParser
..class:Class.ArgumentParser
..summary:Constructor
..signature:ArgumentParser ()
..signature:ArgumentParser (applicationName)
..param.applicationName:A std::string containing the name of the application.
..remarks:If the name of the application is not passed to the constructor it will be extracted from the command line.
*/

class ArgumentParser
{
public:

    // ----------------------------------------------------------------------------
    // Enum ParseResult
    // ----------------------------------------------------------------------------

    // will be used as return value of parse(..) to indicate whether parsing worked
    enum ParseResult
    {
        PARSE_OK,
        PARSE_ERROR,
        PARSE_HELP,
        PARSE_VERSION,
        PARSE_WRITE_CTD,
        PARSE_EXPORT_HELP
    };

    // ----------------------------------------------------------------------------
    // Class Typedefs
    // ----------------------------------------------------------------------------

    typedef std::vector<ArgParseOption>   TOptionMap;
    typedef std::vector<ArgParseArgument> TArgumentMap;
    typedef Size<TOptionMap>::Type        TOptionMapSize;
    typedef Size<TArgumentMap>::Type      TArgumentMapSize;

    typedef std::map<std::string, TOptionMapSize> TStringMap;
    typedef std::vector<std::string>              TValueMap;

    // ----------------------------------------------------------------------------
    // Mapping of option names to options
    // ----------------------------------------------------------------------------

    TStringMap   shortNameMap;
    TStringMap   longNameMap;
    TOptionMap   optionMap;
    TArgumentMap argumentList;

    // ----------------------------------------------------------------------------
    // Documentation Members
    // ----------------------------------------------------------------------------

    ToolDoc                  _toolDoc;      // the tool doc for all user specified
                                            // text
    std::vector<std::string> _description;  // the description which we need to
                                            // separate to put it on top of the rest
    std::vector<std::string> _usageText;    // the usage lines as strings, to avoid
                                            // interference with the rest of the doc

    // ----------------------------------------------------------------------------
    // Function init()
    // ----------------------------------------------------------------------------

    void init()
    {
        addOption(*this, ArgParseOption("h", "help", "Displays this help message."));

        // hidden flags used for export of man pages and ctd formats
        addOption(*this, ArgParseOption("",
                                        "write-ctd",
                                        "Exports the app's interface description to a .ctd file.",
                                        ArgParseArgument::OUTPUTFILE));
        hideOption(*this, "write-ctd", true);

        addOption(*this, ArgParseOption("",
                                        "export-help",
                                        "Export help to a format. One of {'html', 'man', 'txt'}.",
                                        ArgParseArgument::STRING,
                                        "FORMAT"));
        hideOption(*this, "export-help", true);
        setValidValues(*this, "export-help", "html man txt");
    }

    // ----------------------------------------------------------------------------
    // Constructors
    // ----------------------------------------------------------------------------

    ArgumentParser()
    {
        init();
    }

    ArgumentParser(std::string const & _appName)
    {
        setName(_toolDoc, _appName);
        init();
    }

};

// ==========================================================================
// Metafunctions
// ==========================================================================

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function hasOption()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#hasOption
 * @headerfile <seqan/arg_parse.h>
 * @brief Query whether a certain option is registered in the parser.
 *
 * @signature bool hasOption(parser, name);
 *
 * @param[in] parser The ArgumentParser to query.
 * @param[in] name   The name to query for (<tt>std::string</tt>).
 *
 * @return bool <tt>true</tt> if there is such an option, <tt>false</tt> otherwise.
 */

/**
.Function.ArgumentParser#hasOption
..class:Class.ArgumentParser
..summary:Returns whether a certain option is registered in the parser.
..cat:Miscellaneous
..signature:hasOption(parser, optionIdentifier)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionIdentifier:A @Shortcut.CharString@ that identifies the option.
..returns:$true$ if the option is registered.
..include:seqan/arg_parse.h
*/

inline bool hasOption(ArgumentParser const & me, std::string const & name)
{
    return hasKey(me.shortNameMap, name) || hasKey(me.longNameMap, name);
}

// ----------------------------------------------------------------------------
// Function addOption()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#addOption
 * @headerfile <seqan/arg_parse.h>
 * @brief Adds an @link ArgParseOption @endlink to an ArgumentParser.
 *
 * @signature void addOption(parser, option);
 *
 * @param[in,out] parser The ArgumentParser to add the option to.
 * @param[in]     option The ArgParseOption to add to <tt>parser</tt>.
 */

/**
.Function.ArgumentParser#addOption
..class:Class.ArgumentParser
..summary:Adds a @Class.ArgParseOption@ object to the @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:addOption(parser, option)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.option:The new @Class.ArgParseOption@ object that should be added.
...type:Class.ArgParseOption
..include:seqan/arg_parse.h
*/

inline void _copyValidValuesToFileExt(ArgumentParser & me, std::string const & name)
{
    // Copy valid values, remove leading dots.
    ArgParseOption & option = getOption(me, name);
    if (isInputFileArgument(option) || isOutputFileArgument(option))
    {
        std::string longName = option.longName.empty() ? option.shortName : option.longName;
        longName += "-file-ext";
        std::vector<std::string> validValues = option.validValues;
        for (unsigned i = 0; i < length(validValues); ++i)
            if (!validValues[i].empty() && validValues[i][0] == '.')
                validValues[i].erase(0, 1);
        setValidValues(me, longName, validValues);
    }
}

inline void addOption(ArgumentParser & me, ArgParseOption const & opt)
{
    // check if an option with the same identifiers was already registered
    SEQAN_CHECK(!hasOption(me, opt.shortName), "There already is an option with the name %s!", toCString(opt.shortName));
    SEQAN_CHECK(!hasOption(me, opt.longName), "There already is an option with the name %s!", toCString(opt.longName));

    // finally append the option
    appendValue(me.optionMap, opt);

    if (!empty(opt.shortName))
        me.shortNameMap.insert(std::make_pair(opt.shortName, length(me.optionMap) - 1));
    if (!empty(opt.longName))
        me.longNameMap.insert(std::make_pair(opt.longName, length(me.optionMap) - 1));

    // handle the case of input and output option: add a string option --${name}-file-ext.
    if (isInputFileArgument(opt) || isOutputFileArgument(opt))
    {
        std::string longName = opt.longName.empty() ? opt.shortName : opt.longName;
        longName += "-file-ext";
        std::string helpText = "Override file extension for --";
        helpText += opt.longName;

        // Add option, copy list argument, number of allowed values.
        addOption(me, ArgParseOption("", longName, helpText, ArgParseOption::STRING, "EXT",
                                     isListArgument(opt), numberOfAllowedValues(opt)));
        getOption(me, longName.c_str()).tags.push_back("file-ext-override");
        getOption(me, longName.c_str()).tags.push_back("gkn-ignore");
        // Hide option.
        hideOption(me, longName);
        // Copy valid values, remove leading dots.
        _copyValidValuesToFileExt(me, opt.longName);
    }
}

// ----------------------------------------------------------------------------
// Function addArgument()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#addArgument
 * @headerfile <seqan/arg_parse.h>
 * @brief Adds an @link ArgParseArgument @endlink to an ArgumentParser.
 *
 * @signature void addArgument(parser, arg);
 *
 * @param[in,out] parser The ArgumentParser to add the argument to.
 * @param[in]     arg    The ArgParseArgument to add to <tt>parser</tt>.
 */

/**
.Function.ArgumentParser#addArgument
..class:Class.ArgumentParser
..summary:Adds a @Class.ArgParseArgument@ object to the @Class.ArgumentParser@.
..cat:Miscellaneous
..signature:addArgument(parser, argument)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.arg:The new @Class.ArgParseArgument@ object that should be added.
...type:Class.ArgParseArgument
..include:seqan/arg_parse.h
*/

inline void _copyValidValuesToFileExt(ArgumentParser & me, unsigned no)
{
    // Copy valid values, remove leading dots.
    ArgParseArgument & arg = getArgument(me, no);
    if (isInputFileArgument(arg) || isOutputFileArgument(arg))
    {
        std::stringstream longNameSS;
        longNameSS << "arg-" << (no + 1) << "-file-ext";
        std::string longName = longNameSS.str();
        std::vector<std::string> validValues = arg.validValues;
        for (unsigned i = 0; i < length(validValues); ++i)
            if (!validValues[i].empty() && validValues[i][0] == '.')
                validValues[i].erase(0, 1);
        setValidValues(me, longName, validValues);
    }
}

inline void addArgument(ArgumentParser & me, ArgParseArgument const & arg)
{
    // check previous arguments
    //  .. lists can only be last argument
    if (!me.argumentList.empty())
    {
        SEQAN_CHECK(!isListArgument(me.argumentList[me.argumentList.size() - 1]),
                    "You cannot add an additional argument after a list argument.");
    }

    // check current argument
    //  .. arguments should not have default values
    SEQAN_CHECK(arg.defaultValue.empty(), "Arguments cannot have default values.");
    SEQAN_CHECK(arg._numberOfValues == 1, "n-Tuple of arguments are not supported.");

    me.argumentList.push_back(arg);

    // handle the case of input and output option: add a string option --${name}-file-ext.
    if (isInputFileArgument(arg) || isOutputFileArgument(arg))
    {
        std::stringstream longNameSS;
        longNameSS << "arg-" << me.argumentList.size() << "-file-ext";
        std::string longName = longNameSS.str();
        std::stringstream helpTextSS;
        helpTextSS << "Override file extension for argument " << me.argumentList.size();
        std::string helpText = helpTextSS.str();

        // Add option, copy list argument, number of allowed values.
        addOption(me, ArgParseOption("", longName, helpText, ArgParseOption::STRING, "EXT",
                                     isListArgument(arg), numberOfAllowedValues(arg)));
        getOption(me, longName.c_str()).tags.push_back("file-ext-override");
        // Hide option.
        hideOption(me, longName);
        // Copy valid values, remove leading dots.
        _copyValidValuesToFileExt(me, me.argumentList.size() - 1);
    }
}

// ----------------------------------------------------------------------------
// Function _getOptionIndex()
// ----------------------------------------------------------------------------
// note that it is assumed that the option exists if this method is called

inline ArgumentParser::TOptionMapSize _getOptionIndex(ArgumentParser const & me,
                                                      std::string const & name)
{
    ArgumentParser::TOptionMapSize option_index;
    if (me.shortNameMap.find(name) != me.shortNameMap.end())
    {
        option_index = me.shortNameMap.find(name)->second;
    }
    else
    {
        option_index = me.longNameMap.find(name)->second;
    }
    return option_index;
}

// ----------------------------------------------------------------------------
// Function getOption()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getOption
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns a reference to the specified option.
 *
 * @signature TOption getOption(parser, name);
 *
 * @param[in] parser The parser to query.
 * @param[in] name   The short or long name of the option (<tt>std::string</tt>).
 *
 * @return TOption Reference to the @link ArgParseOption @endlink with the given short or long name.
 */

/**
.Function.ArgumentParser#getOption
..class:Class.ArgumentParser
..summary:Returns a reference to the specified option.
..cat:Miscellaneous
..signature:getOption(parser, optionName)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionName:The identifier of the command line option.
..returns: a reference to the specified @Class.ArgParseOption@ object.
..include:seqan/arg_parse.h
*/

inline ArgParseOption & getOption(ArgumentParser & me, std::string const & name)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return me.optionMap[_getOptionIndex(me, name)];
}

inline ArgParseOption const & getOption(ArgumentParser const & me, std::string const & name)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return me.optionMap[_getOptionIndex(me, name)];
}

// ----------------------------------------------------------------------------
// Function setRequired()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#setRequired
 * @headerfile <seqan/arg_parse.h>
 * @brief Sets whether or not the option with the givne name is mandatory.
 *
 * @signature void setRequired(parser, name[, required]).
 *
 * @param[in,out] parser   The ArgumentParser to set the flag of.
 * @param[in]     name     The short or long name of the option (<tt>std::string</tt>).
 * @param[in]     required Whether or not the option is required (<tt>bool</tt>, default to <tt>true</tt>).
 */

/**
.Function.ArgumentParser#setRequired
..class:Class.ArgumentParser
..summary:Sets whether or not the option defined by the parameter $name$ (which can be
 either the short or the long name) is mandatory.
..remarks: Note that the empty string is, at least for string options, also a valid string.
Hence setting an option to required does not guarantee that the returned string is not empty.
..cat:Miscellaneous
..signature:setRequired(parser, optionName [, required])
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionName:The identifier of the command line option.
..param.required.remarks: The default value is true.
..include:seqan/arg_parse.h
*/

inline void setRequired(ArgumentParser & me, std::string const & name, bool required = true)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    setRequired(getOption(me, name), required);
}

// ----------------------------------------------------------------------------
// Function hideOption()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#hideOption
 * @headerfile <seqan/arg_parse.h>
 * @brief Hides the ArgParseOption with the given name.
 *
 * @signature void hideOption(parser, name[, hide]).
 *
 * @param[in,out] parser The ArgParseOption to the the hidden flag of.
 * @param[in]     name   The short or long name of the option to modify.
 * @param[in]     hide   Whether or not to hide the flag (<tt>bool</tt>, defaults to <tt>true</tt>).
 */

/**
.Function.ArgumentParser#hideOption
..class:Class.ArgumentParser
..summary:Hides the ArgParseOption defined by the parameter $name$ (which can be
 either the short or the long name) from the help screen.
..cat:Miscellaneous
..signature:hideOption(parser, optionName [, hide])
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionName:The identifier of the command line option.
..param.hide:The new visibility of the option. Default is false.
...type:nolink:bool
..include:seqan/arg_parse.h
*/

inline void hideOption(ArgumentParser & me, std::string const & name, bool hide)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    hideOption(getOption(me, name), hide);
}

// ----------------------------------------------------------------------------
// Function getArgument()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getArgument
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns a reference to the given positional argument.
 *
 * @signature TArgument getArgument(parser, pos);
 *
 * @param[in] parser The ArgumentParser to query.
 * @param[in] pos    The position of the argument to return (<tt>unsigned</tt>, starting at 0).
 *
 * @return TArgument Reference to the @link ArgParseArgument @endlink with the given position.
 */

/**
.Function.ArgumentParser#getArgument
..class:Class.ArgumentParser
..summary:Returns a reference to the specified argument.
..cat:Miscellaneous
..signature:getArgument(parser, argumentPosition)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.argumentPosition:The index of the argument in the argument list.
..returns: a reference to the specified @Class.ArgParseArgument@ object.
..include:seqan/arg_parse.h
*/

inline ArgParseArgument & getArgument(ArgumentParser & me, unsigned position)
{
    SEQAN_CHECK(position < me.argumentList.size(),
                "ArgumentParser: Only %d arguments available", me.argumentList.size());
    return me.argumentList[position];
}

inline ArgParseArgument const & getArgument(ArgumentParser const & me, unsigned position)
{
    SEQAN_CHECK(position < me.argumentList.size(),
                "ArgumentParser: Only %d arguments available", me.argumentList.size());
    return me.argumentList[position];
}

// ----------------------------------------------------------------------------
// Function isSet()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#isSet
 * @headerfile <seqan/arg_parse.h>
 * @brief Query whether an option was set on the command line.
 *
 * @signature bool isSet(parser, name);
 *
 * @param[in] parser The ArgumentParser to query.
 * @param[in] name   The short or long name of the option (<tt>std::string</tt>).
 *
 * @return bool Whether or not the option was set on the command line or not.
 */

/**
.Function.ArgumentParser#isSet
..class:Class.ArgumentParser
..summary:Returns whether an option was set on the parsed command line.
..cat:Miscellaneous
..signature:isSet(parser,optionIdentifier)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionIdentifier:A std::string that identifies the option (either short or long name).
..returns:$true$ if the option was set.
..include:seqan/arg_parse.h
*/

inline bool isSet(ArgumentParser const & me, std::string const & name)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return isSet(getOption(me, name));
}

// ----------------------------------------------------------------------------
// Function hasDefault()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#hasDefault
 * @headerfile <seqan/arg_parse.h>
 * @brief Query whether an option has a default value.
 *
 * @signature bool hasDefault(parser, name);
 *
 * @param[in] parser The ArgumentParser to query.
 * @param[in] name   The short or long name of the option (<tt>std::string</tt>).
 *
 * @return bool Whether or not the option has a default value.
 */

/**
.Function.ArgumentParser#hasDefault
..summary:Returns whether an option has a default value or not.
..cat:Miscellaneous
..signature:hasDefault(parser,optionIdentifier)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionIdentifier:A std::string that identifies the option (either short or long name).
..returns:$true$ if the option has a default value, $false$ otherwise.
..include:seqan/arg_parse.h
*/

inline bool hasDefault(ArgumentParser const & me, std::string const & name)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return hasDefault(getOption(me, name));
}

// ----------------------------------------------------------------------------
// Function _allRequiredSet()
// ----------------------------------------------------------------------------

inline bool _allRequiredSet(ArgumentParser const & me)
{
    for (unsigned o = 0; o < length(me.optionMap); ++o)
        if (!isSet(me.optionMap[o]) && isRequired(me.optionMap[o]))
            return false;

    return true;
}

// ----------------------------------------------------------------------------
// Function _allArgumentsSet()
// ----------------------------------------------------------------------------

inline bool _allArgumentsSet(ArgumentParser const & me)
{
    for (unsigned a = 0; a < me.argumentList.size(); ++a)
        if (!isSet(me.argumentList[a]))
            return false;

    return true;
}

// ----------------------------------------------------------------------------
// Function getOptionValue()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getOptionValue
 * @headerfile <seqan/arg_parse.h>
 * @brief Retrieve the value of an option.
 *
 * @signature bool getOptionValue(dest, parser, name[, pos]);
 *
 * @param[in] dest   The variable to write the result to (the type is a template parameter and the value type of the
 *                    option must be convertible in the type of <tt>dest</tt> for the retrieval to work, also see
 *                    result value).
 * @param[in] parser The ArgumentParser to get the value from.
 * @param[in] name   The short or long name of the option (<tt>std::string</tt>).
 * @param[in] pos    Optional position for multi-value options (<tt>unsigned</tt>, defaults to 0).
 *
 * @return bool <tt>true</tt> if the requested option was given on the command line and could be coverted to the type of
 *              <tt>dest</tt>.
 */

/**
.Function.ArgumentParser#getOptionValue
..class:Class.ArgumentParser
..summary:Retrieves the value of an option given either the short or long name.
..cat:Miscellaneous
..signature:getOptionValue(value, parser, optionIdentifier[, argNo])
..param.value:The variable where the resulting value should be stored.
...remarks:The type of $value$ must be compatible the option type.
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionIdentifier:A std::string that is either the short or long name of the option.
..param.argNo:If the option is list, the $argNo$-th list element is returned.
..returns: $true$ if the requested option is set and has the requested type, $false$ otherwise.
..remarks:The value passed to the method (value) was only updated if the method returns $true$.
..include:seqan/arg_parse.h
*/

template <typename TValue>
inline bool getOptionValue(TValue & val,
                           ArgumentParser const & me,
                           std::string const & name,
                           unsigned argNo)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));

    if (isSet(me, name) || hasDefault(me, name))
        return _convertArgumentValue(val,
                                     getOption(me, name),
                                     getArgumentValue(getOption(me, name), argNo));
    else
        return false;
}

template <typename TValue>
inline bool getOptionValue(TValue & val,
                           ArgumentParser const & me,
                           std::string const & name)
{
    return getOptionValue(val, me, name, 0);
}

// ----------------------------------------------------------------------------
// Function getOptionFileExtension()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getOptionFileExtension
 * @headerfile <seqan/arg_parse.h>
 * @brief Retrieve the file extension of a file option.
 *
 * @signature std::string getOptionFileExtension(parser, name[, pos]);
 *
 * @param[in] parser The ArgumentParser to get the value from.
 * @param[in] name   The short or long name of the option (<tt>std::string</tt>).
 * @param[in] pos    Optional position for multi-value options (<tt>unsigned</tt>, defaults to 0).
 *
 * @return std::string The extension of the option. Empty if not set or no extension.
 *
 * @see ArgumentParser#getArgumentFileExtension
 *
 * @section Overriding File Extension on the Command Line
 *
 * For each option with type <tt>INPUTFILE</tt> and <tt>OUTPUTFILE</tt>, an option with the name
 * <tt>${name}-file-ext</tt> is automatically added to the ArgumentParser (where <tt>${name}</tt> is the name
 * of the original option).  The extension can be overridden by specifying the argument.  Thus, the user of
 * the program could give the value "file.ext" to the parameter "fname" and override the extension on the
 * command line to "ext2" as follows:
 *
 * @code{.console}
 * # program_name --fname file.ext --fname-file-ext ext2
 * @endcode
 */

/**
.Function.ArgumentParser#getOptionFileExtension
..class:Class.ArgumentParser
..summary:Returns the extension of a file option.
..cat:Miscellaneous
..signature:std::string getOptionFileExtension(parser, name[, argNo]);
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParsre
..param.name:The name of the option.
..param.argNo:An optional index for multi-value options.
...type:nolink:$unsigned$
..returns:A $std::string$ with the extension. Empty if no extension or empty value.
*/

inline std::string getOptionFileExtension(ArgumentParser const & me,
                                          std::string const & name,
                                          unsigned argNo = 0)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));

    return getFileExtension(getOption(me, name), argNo);
}

// ----------------------------------------------------------------------------
// Function getOptionValueCount()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getOptionValueCount
 * @headerfile <seqan/arg_parse.h>
 * @brief Query number of values stored for the specified option.
 *
 * @signature unsigned getOptionValueCount(parser, name);
 *
 * @param[in] parser The ArgumentParser to query.
 * @param[in] name   The short or long name of the option (<tt>string</tt>).
 *
 * @return unsigned The number of values for the option with the given name.
 */

/**
.Function.ArgumentParser#getOptionValueCount
..class:Class.ArgumentParser
..summary:Returns the number of values stored in the specified option.
..cat:Miscellaneous
..signature:getOptionValueCount(parser, optionIdentifier)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionIdentifier:A std::string that is either the short or long name of the option.
..returns: The number of values stored for this option.
..include:seqan/arg_parse.h
*/

inline unsigned getOptionValueCount(ArgumentParser const & me, std::string const & name)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return getArgumentValues(getOption(me, name)).size();
}

// ----------------------------------------------------------------------------
// Function getArgumentValueCount()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getArgumentValueCount
 * @headerfile <seqan/arg_parse.h>
 * @brief Query number of values stored for the specified argument.
 *
 * @signature unsigned getArgumentValueCount(parser, pos);
 *
 * @param[in] parser The ArgumentParser to query.
 * @param[in] name   The position of the argument (<tt>unsigned</tt>, 0-based).
 *
 * @return unsigned The number of values for the argument with the given position.
 */

/**
.Function.ArgumentParser#getArgumentValueCount
..class:Class.ArgumentParser
..summary:Returns the number of values stored in the specified option.
..cat:Miscellaneous
..signature:getArgumentValueCount(parser, argumentPosition)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.argumentPosition:The index of the argument in the argument list.
..returns: The number of values stored for the specified argument.
..include:seqan/arg_parse.h
*/

inline unsigned getArgumentValueCount(ArgumentParser const & me, unsigned argumentPosition)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition,
                "Argument Parser has only %d arguments.",
                me.argumentList.size());
    return getArgumentValues(getArgument(me, argumentPosition)).size();
}

// ----------------------------------------------------------------------------
// Function getArgumentValue()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getArgumentValue
 * @headerfile <seqan/arg_parse.h>
 * @brief Retrieves the value of an argument given by its position.
 *
 * @signature bool getArgumentValue(dest, parser, pos[, no]);
 *
 * @param[in] dest   The variable to write the result to (the type is a template parameter and the value type of the
 *                   argument must be convertible in the type of <tt>dest</tt> for the retrieval to work, also see
 *                   result value).
 * @param[in] parser The ArgumentParser to get the value from.
 * @param[in] pos    The position of the argument to get the value of.
 * @param[in] no     Optional position for multi-value arguments (<tt>unsigned</tt>, defaults to 0).
 *
 * @return bool <tt>true</tt> if the retrieval was successful, <tt>false</tt> otherwise.
 */

/**
.Function.ArgumentParser#getArgumentValue
..class:Class.ArgumentParser
..summary:Retrieves the value of an argument given by its position.
..cat:Miscellaneous
..signature:getArgumentValue(value, parser, argumentPosition[, argNo])
..param.value:The variable where the resulting value should be stored.
...remarks:The type of $value$ must be compatible the option type.
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.argumentPosition:The index of the argument in the argument list.
..param.argNo:If the argument is a list, the $argNo$-th list element is returned.
..returns: $true$ if the requested argument is set and has the requested type, $false$ otherwise.
..include:seqan/arg_parse.h
*/

template <typename TValue>
inline bool getArgumentValue(TValue & value,
                             ArgumentParser const & me,
                             unsigned argumentPosition,
                             unsigned argNo)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition,
                "Argument Parser has only %d arguments.",
                me.argumentList.size());
    return _convertArgumentValue(value, getArgument(me, argumentPosition), getArgumentValue(getArgument(me, argumentPosition), argNo));
}

template <typename TValue>
inline bool getArgumentValue(TValue & value,
                             ArgumentParser const & me,
                             unsigned argumentPosition)
{
    return getArgumentValue(value, me, argumentPosition, 0);
}

// ----------------------------------------------------------------------------
// Function getArgumentFileExtension()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getArgumentFileExtension
 * @headerfile <seqan/arg_parse.h>
 * @brief Retrieve the file extension of a file argument.
 *
 * @signature std::string argumentFileExtension(parser, pos[, argNo]);
 *
 * @param[in] parser The ArgumentParser to get the value from.
 * @param[in] pos    The position of the argument to query (<tt>unsigned</tt>).
 * @param[in] argNo  Optional position for multi-value options (<tt>unsigned</tt>, defaults to 0).
 *
 * @return std::string The extension of the argument if any.
 *
 * @see ArgumentParser#getOptionFileExtension
 *
 * @section Overriding File Extensions on the Command Line
 *
 * For each argument with type <tt>INPUTFILE</tt> and <tt>OUTPUTFILE</tt>, an option with the index
 * <tt>arg-${idx}-file-ext</tt> is automatically added to the ArgumentParser (where <tt>${idx}</tt> is the index
 * of the original option).  The extension can be overridden by specifying the argument.  Thus, the user of
 * the program could give the value "file.ext" to the parameter "0" and override the extension on the
 * command line to "ext2" as follows:
 *
 * @code{.console}
 * # program_name file.ext --arg-0-file-ext ext2
 * @endcode
 */

/**
.Function.ArgumentParser#getArgumentFileExtension
..class:Class.ArgumentParser
..summary:Returns the extension of a file argument.
..cat:Miscellaneous
..signature:std::string argumentFileExtension(parser, argPos[, argNo]);
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParsre
..param.argPos:The position of the argument.
...type:nolink:$unsigned$
..param.argNo:An optional index for multi-value arguments.
...type:nolink:$unsigned$
..returns:A $std::string$ with the extension. Empty if no extension or empty value.
*/

inline std::string getArgumentFileExtension(ArgumentParser const & me,
                                            unsigned argumentPosition,
                                            unsigned argNo = 0)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition,
                "Argument Parser has only %d arguments.",
                me.argumentList.size());


    return getFileExtension(getArgument(me, argumentPosition), argNo);
}

// ----------------------------------------------------------------------------
// Function getOptionValues()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getOptionValues
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns all values of an option given on the command line.
 *
 * @signature TVector getOptionValues(parser, name);
 *
 * @param[in] parser The ArgumentParser to query.
 * @param[in] name   The short or long name of the option to get (<tt>std::string</tt>).
 *
 * @return TVector The resulting values (<tt>std::vector&lt;std::string&gt;</tt>).
 */

/**
.Function.ArgumentParser#getOptionValues
..class:Class.ArgumentParser
..summary:Returns all values of an option given on the command line.
..cat:Miscellaneous
..signature:getOptionValues(parser, optionIdentifier)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionIdentifier:A std::string that is either the short or long name of the option.
..returns: A $String<std::string>$ of option values.
..include:seqan/arg_parse.h
*/

inline std::vector<std::string> const & getOptionValues(ArgumentParser const & me,
                                                        std::string const & name)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    return getArgumentValues(getOption(me, name));
}

// ----------------------------------------------------------------------------
// Function getArgumentValues()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getArgumentValues
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns all values of an argument given on the command line.
 *
 * @signature TVector getArgumentValues(parser, pos);
 *
 * @param[in] parser The ArgumentParser to query.
 * @param[in] pos    The position of the argument (<tt>unsigned</tt>, 0-based).
 *
 * @return TVector The resulting values (<tt>std::vector&lt;std::string&gt;</tt>).
 */

/**
.Function.ArgumentParser#getArgumentValues
..class:Class.ArgumentParser
..summary:Returns all values of an option given on the command line.
..cat:Miscellaneous
..signature:getArgumentValues(parser, argumentPosition)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.argumentPosition:The index of the argument in the argument list.
..returns: A $String<std::string>$ of argument values.
..include:seqan/arg_parse.h
*/

inline std::vector<std::string> const & getArgumentValues(ArgumentParser const & me,
                                                          unsigned argumentPosition)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition,
                "Argument Parser has only %d arguments.",
                me.argumentList.size());
    return getArgumentValues(getArgument(me, argumentPosition));
}

// ----------------------------------------------------------------------------
// Function setDefaultValue()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#setDefaultValue
 * @headerfile <seqan/arg_parse.h>
 * @brief Set the default value of an option of an ArgumentParser.
 *
 * @signature void setDefaultValue(parser, name, v);
 *
 * @param[in] parser The ArgumentParser to set the default value to.
 * @param[in] name   The short or long name of the argument (<tt>std::string</tt>).
 * @param[in] v      The value to set (template parameter, must be streamable into a <tt>std::stringstream</tt>).
 */

/**
.Function.ArgumentParser#setDefaultValue
..class:Class.ArgumentParser
..summary:Set default value of an option of an ArgumentParser.
..signature:setDefaultValue(parser, optionName, value)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionName:The identifier of the command line option.
..param.value:The new default value.
..include:seqan/arg_parse.h
*/

template <typename TValue>
inline void setDefaultValue(ArgumentParser & me,
                            std::string const & name,
                            const TValue & value)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    setDefaultValue(getOption(me, name), value);
}

// ----------------------------------------------------------------------------
// Function addDefaultValue()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#addDefaultValue
 * @headerfile <seqan/arg_parse.h>
 * @brief Add/append a value to the default values for an option in an ArgumentParser.
 *
 * @signature void addDefaultValue(parser, name, v);
 *
 * @param[in,out] parser The ArgumentParser to append the default value to.
 * @param[in]     name   The short or long name of the argument (<tt>std::string</tt>).
 * @param[in]     v      The value to append (template parameter, must be streamable into a <tt>std::stringstream</tt>).
 */

/**
.Function.ArgumentParser#addDefaultValue
..class:Class.ArgumentParser
..summary:Add to the default values of an option of an ArgumentParser.
..signature:addDefaultValue(parser, optionName, value)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionName:The identifier of the command line option.
..param.value:The new default value.
..include:seqan/arg_parse.h
*/

template <typename TValue>
inline void addDefaultValue(ArgumentParser & me,
                            std::string const & name,
                            const TValue & value)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    addDefaultValue(getOption(me, name), value);
}

// ----------------------------------------------------------------------------
// Function setMinValue()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#setMinValue
 * @headerfile <seqan/arg_parse.h>
 * @brief Set smallest allowed value for an option or argument of an ArgumentParser.
 *
 * @signature void setMinValue(parser, name, v);
 * @signature void setMinValue(parser, pos, v);
 *
 * @param[in,out] parser The ArgumentParser to set the minimal value for.
 * @param[in]     name   The name of the option to set the minimal value for (<tt>std::string</tt>).
 * @param[in]     pos    The position of the argument to set the minimal value for (<tt>unsigned</tt>, 0-based).
 * @param[in]     v      The minimal value to set (<tt>std::string</tt>).
 *
 * @section Remarks
 *
 * The option/argument must have an integer or double type.
 */

/**
.Function.ArgumentParser#setMinValue
..class:Class.ArgumentParser
..summary:Set smallest allowed value for an option or argument of an ArgumentParser.
..signature:setMinValue(parser,optionName,minValue)
..signature:setMinValue(parser,argumentPosition,minValue)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionName:The identifier of the command line option.
..param.argumentPosition:The index of the argument in the argument list.
..param.minValue:A std::string containing a string representation of the minimum value of the @Class.ArgParseOption@.
..include:seqan/arg_parse.h
*/

inline void setMinValue(ArgumentParser & me,
                        std::string const & name,
                        std::string const & _minValue)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    setMinValue(getOption(me, name), _minValue);
}

inline void setMinValue(ArgumentParser & me,
                        unsigned argumentPosition,
                        std::string const & _minValue)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition,
                "Argument Parser has only %d arguments.",
                me.argumentList.size());
    setMinValue(getArgument(me, argumentPosition), _minValue);
}

// ----------------------------------------------------------------------------
// Function setMaxValue()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#setMaxValue
 * @headerfile <seqan/arg_parse.h>
 * @brief Set largest allowed value for an option or argument of an ArgumentParser.
 *
 * @signature void setMaxValue(parser, name, v);
 * @signature void setMaxValue(parser, pos, v);
 *
 * @param[in,out] parser The ArgumentParser to set the maximal value for.
 * @param[in]     name   The name of the option to set the maximal value for (<tt>std::string</tt>).
 * @param[in]     pos    The position of the argument to set the maximal value for (<tt>unsigned</tt>, 0-based).
 * @param[in]     v      The maximal value to set (<tt>std::string</tt>).
 *
 * @section Remarks
 *
 * The option/argument must have an integer or double type.
 */

/**
.Function.ArgumentParser#setMaxValue
..class:Class.ArgumentParser
..summary:Set largest allowed value for an option or argument of an ArgumentParser.
..signature:setMaxValue(parser,optionName,maxValue)
..signature:setMaxValue(parser,argumentPosition,minValue)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionName:The identifier of the command line option.
..param.argumentPosition:The index of the argument in the argument list.
..param.maxValue:A std::string containing a string representation of the maximum value of the @Class.ArgParseOption@.
..include:seqan/arg_parse.h
*/

inline void setMaxValue(ArgumentParser & me,
                        std::string const & name,
                        std::string const & _maxValue)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    setMaxValue(getOption(me, name), _maxValue);
}

inline void setMaxValue(ArgumentParser & me,
                        unsigned argumentPosition,
                        std::string const & _minValue)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition,
                "Argument Parser has only %d arguments.",
                me.argumentList.size());
    setMaxValue(getArgument(me, argumentPosition), _minValue);
}

// ----------------------------------------------------------------------------
// Function setValidValues()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#setValidValues
 * @headerfile <seqan/arg_parse.h>
 * @brief Set valid values for an argumetn or option of an ArgumentParser.
 *
 * @signature void setValidValues(parser, name, values);
 * @signature void setValidValues(parser, pos, values);
 *
 * @param[in,out] parser The ArgumentParser to set the default values to.
 * @param[in]     name   The name of the option (<tt>std::string</tt>).
 * @param[in]     pos    The position of the argument (<tt>unsigned</tt>, 0-based).
 * @param[in]     values The values to set.  Either a <tt>std::string</tt> with the values as space-separated list
 *                       or a <tt>std::vector&lt;std::string&gt;</tt> with the values.
 */

/**
.Function.ArgumentParser#setValidValues
..class:Class.ArgumentParser
..summary:Set valid values for an argument or option of an ArgumentParser.
..signature:setValidValues(parser,optionName,values)
..signature:setValidValues(parser,argumentPosition,values)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionName:The identifier of the command line option.
..param.argumentPosition:The index of the argument in the argument list.
..param.values:A $std::string$ containing all valid entries for the option.
Alternatively you can pass a string containing all values separated by spaces.
..include:seqan/arg_parse.h
*/

inline void setValidValues(ArgumentParser & me,
                           std::string const & name,
                           std::vector<std::string> const & values)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    setValidValues(getOption(me, name), values);
    _copyValidValuesToFileExt(me, name);
}

inline void setValidValues(ArgumentParser & me,
                           std::string const & name,
                           std::string const & values)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    setValidValues(getOption(me, name), values);
    _copyValidValuesToFileExt(me, name);
}

inline void setValidValues(ArgumentParser & me,
                           unsigned argumentPosition,
                           std::vector<std::string> const & values)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition,
                "Argument Parser has only %d arguments.",
                me.argumentList.size());
    setValidValues(getArgument(me, argumentPosition), values);
    _copyValidValuesToFileExt(me, argumentPosition);
}

inline void setValidValues(ArgumentParser & me,
                           unsigned argumentPosition,
                           std::string const & values)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition,
                "Argument Parser has only %d arguments.",
                me.argumentList.size());
    setValidValues(getArgument(me, argumentPosition), values);
    _copyValidValuesToFileExt(me, argumentPosition);
}

// ----------------------------------------------------------------------------
// Function setHelpText()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#setHelpText
 * @headerfile <seqan/arg_parse.h>
 * @brief Set the help text of an option or argument.
 *
 * @signature void setHelpText(parser, name, text);
 * @signature void setHelpText(parser, pos, text);
 *
 * @param[in,out] parser The ArgumentParser object.
 * @param[in]     name   The name of the option to set the help text for (<tt>std::string</tt>).
 * @param[in]     pos    The position of the argument to set the help text for.
 * @param[in]     text   The string to use for the help text (<tt>std::string</tt>).
 */

/**
.Function.ArgumentParser#setHelpText
..class:Class.ArgumentParser
..summary:Set help text of argument parser.
..signature:setHelpText(parser,optionName,text)
..signature:setHelpText(parser,argumentPosition,text)
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.optionName:The identifier of the command line option.
..param.argumentPosition:The index of the argument in the argument list.
..param.text:A $std::string$ describing the option or argument.
..include:seqan/arg_parse.h
*/

inline void setHelpText(ArgumentParser & me,
                        std::string const & name,
                        std::string const & text)
{
    SEQAN_CHECK(hasOption(me, name), "Unknown option: %s", toCString(name));
    setHelpText(getOption(me, name), text);
}

inline void setHelpText(ArgumentParser & me,
                        unsigned argumentPosition,
                        std::string const & text)
{
    SEQAN_CHECK(me.argumentList.size() > argumentPosition,
                "Argument Parser has only %d arguments.",
                me.argumentList.size());
    setHelpText(getArgument(me, argumentPosition), text);
}

// ----------------------------------------------------------------------------
// Function getFileFormatExtensions()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getFileFormatExtensions
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns file format extension given a format tag.
 *
 * @signature TVector getFormatExtension(tag);
 * @signature TVector getFormatExtension(tagList);
 * @signature TVector getFormatExtension(tagSelector);
 *
 * @param[in] tag         A single file foramt, e.g. <tt>Fastq()</tt>.
 * @param[in] tagList     A list of file format (@link TagList @endlink).
 * @param[in] tagSelector A file format selector (@link TagSelector @endlink).
 *
 * @return TVector A <tt>std::vector&lt;std::string&gt;</tt> with the allowed file format extensions.
 */

/**
.Function.ArgumentParser#getFileFormatExtensions
..class:Class.ArgumentParser
..summary:Returns file format extensions given a format tag.
..signature:getFileFormatExtensions(formatTag)
..signature:getFileFormatExtensions(formatTagList)
..signature:getFileFormatExtensions(formatTagSelector)
..param.format:A single file format, e.g. @Tag.File Format.tag.Fastq@ or @Tag.Sam@.
...type:Tag.Tag
..param.formatTagList:A list of file formats.
...type:Tag.TagList
..param.formatTagSelector:A file format selector.
...type:Class.TagSelector
..returns:A $std::vector<std::string>$ of all extensions supported by a single format or all formats of a list or selector.
..include:seqan/arg_parse.h
*/

template <typename T>
inline std::vector<std::string>
getFileFormatExtensions(T const formatTag)
{
    std::vector<std::string> extensions;
    _getFileFormatExtensions(extensions, formatTag);
    return extensions;
}


}  // namespace seqan

#endif // SEQAN_CORE_INCLUDE_ARG_PARSE_ARGUMENT_PARSER_H_
