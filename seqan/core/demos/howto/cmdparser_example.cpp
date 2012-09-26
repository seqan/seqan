// FRAGMENT(headers)
#include <seqan/basic.h>
#include <seqan/file.h>  // For printing Seqan strings.
#include <seqan/misc/misc_cmdparser.h>

int main(int argc, const char *argv[]) {
    using namespace seqan;

    CommandLineParser parser("<your application name>"); // if you do not pass the name to the c'tor it will guess it from argv

    // FRAGMENT(cmdparser-options)
    addOption(parser, CommandLineOption('d', "double", "a double option", OptionType::Double));
    addOption(parser, CommandLineOption('i', "int", "an integer option", OptionType::Int));
    addOption(parser, CommandLineOption("i2", "int2", "another integer option", OptionType::Int));
    addOption(parser, CommandLineOption('s', "string", "a mandatory string option", (OptionType::String | OptionType::Mandatory)));
    addOption(parser, CommandLineOption('b', "bool", "a boolean option", OptionType::Boolean));
    addOption(parser, CommandLineOption('c', "hidden-bool", "a boolean option (e.g. only for debugging) which is not visible in the help message", (OptionType::Boolean | OptionType::Hidden)));
    addOption(parser, CommandLineOption('l', "list", "integer argument can be given multiple times and is not overwritten", (OptionType::Int | OptionType::List)));

    // FRAGMENT(cmdparser-add-lines)
    addLine(parser, "Here goes your application description!");
    addLine(parser, "And another line of your application description!");

    // FRAGMENT(cmdparser-header)
    addTitleLine(parser, "*************************************");
    addTitleLine(parser, "* My application name               *");
    addTitleLine(parser, "*************************************");

    addUsageLine(parser, "-s <mandatory argument> [Options]");

    // FRAGMENT(cmdparser-arguments)
    // Tell the parser that we expect/require at least one non option argument
    requiredArguments(parser,1);

    // FRAGMENT(cmdparser-parse)
    if(!parse(parser, argc, argv)) {
        return 1;
    } else if(isSetShort(parser,'h')) {
        help(parser,std::cerr);
        return 0;
    }

    // FRAGMENT(cmdparser-get-values)
    // Now, lets get the values from our command line.
    double d = 0.0;
    if (isSetShort (parser, 'd'))
        getOptionValueShort (parser,'d',d);

    // Get a list of options from the '-l' option.
    unsigned lCount = length(getOptionValuesLong(parser, "list"));
    String<int> values;
    resize(values, lCount);
    for (unsigned i = 0; i < lCount; ++i)
        getOptionValueLong(parser, "list", i, values[i]);

    String<char> s = "";

    // W do not have to check if it was set since this option is
    // mandatory if it would not be set by the user then parse would
    // have returned false.
    getOptionValueShort (parser, 's', s);

    // Get the required argument
    String<char> arg = getArgumentValue(parser, 0);

    // Get all arguments.
    String<CharString> args = getArgumentValues(parser);

    // FRAGMENT(main-end)
    return 0;
}
