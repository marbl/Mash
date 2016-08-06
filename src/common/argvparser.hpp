/*
 *   @ingroup utils
 *
 *   C++ command line argument parser
 *
 *   Copyright (C) 2005 by
 *   Michael Hanke        michael.hanke@gmail.com
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 */
#ifndef __ARGVPARSER_H
#define __ARGVPARSER_H

#include <string>
#include <vector>
#include <map>
#include <list>
#include <sstream>

namespace CommandLineProcessing
{

/** Provides parsing and storage of POSIX-like command line arguments (argc, argv).
* To use this class for CLI-option parsing, first define a set of valid options using
* the defineOption() method. An option can have several attributes: it can be required
* itself (its missing is considered as an error) and/or it can require a value.
* Options with optional values can be realized by defining the option to not need a
* value and use the syntax '--option=value' or '-o=value' on the command line.
* Every option can have different alternative labels; see defineOptionAlternative().
* A usage description (see usageDescription()) can be generated from the set of defined options (also see the
* addErrorCode(), setIntroductoryDescription() and setHelpOption() methods).
* \note The implemented parsing algorithm requires that all options have to be given at
* the beginning of the command line. Everything on the commandline after the first
* non-option (or option value) is considered as an argument.
* \attention Short option labels (single letter options) must not be digits
* (the option string itself not a possible value).
* Valid syntaxes are:
* \li program --long-option value -sdfgh -u=5 -i 7 --last=val arg1 arg2 arg3
* Here is a small code example:
* \code
* #include
* ArgvParser cmd;  // the command line parser
*
* // init
* cmd.setIntroductoryDescription("This is foo written by bar.");
*
* //define error codes
* cmd.addErrorCode(0, "Success");
* cmd.addErrorCode(1, "Error");
*
* cmd.setHelpOption("h", "help", "Print this help");
*
* cmd.defineOption("version", ArgvParser::NoOptionAttribute, "Be verbose");
* cmd.defineOptionAlternative("verbose","v");
*
* cmd.defineOption("foo", ArgvParser::OptionRequiresValue, "Fooishness. Default value: 0");
*
* // finally parse and handle return codes (display help etc...)
* int result = cmd.parse(argc, argv);
*
* if (result != ArgvParser::NoParserError)
*   cout << cmd.parseErrorDescription(results);
*   exit(1);
* \endcode
*
* \author Michael Hanke
*/
class ArgvParser
{
public:
    typedef int OptionAttributes;
    typedef int ParserResults;
    typedef std::map<std::string, unsigned int> String2KeyMap;
    typedef std::map<unsigned int, OptionAttributes> Key2AttributeMap;
    typedef std::map<unsigned int, std::string> Key2StringMap;
    typedef std::vector<std::string> ArgumentContainer;

    ArgvParser();
    ~ArgvParser();

    /** Attributes for options. */
    enum
    {
        NoOptionAttribute = 0x00,
        OptionRequiresValue = 0x01,
        OptionRequired = 0x02
    };
    /** Return values of the parser. */
    enum
    {
        NoParserError = 0x00,
        ParserUnknownOption = 0x01,
        ParserMissingValue = 0x02,
        ParserOptionAfterArgument = 0x04,
        ParserMalformedMultipleShortOption = 0x08,
        ParserRequiredOptionMissing = 0x16,
        ParserHelpRequested = 0x32
    };

    /** Defines an option with optional attributes (required, ...) and an
    * additional (also optional) description. The description becomes part of the
    * generated usage help that can be requested by calling the usageDescription()
    * method.
    * \return Returns FALSE if there already is an option with this name
    * OR if a short option string (length == 1) is a digit. In that case no
    * action is peformed.
    */
    bool defineOption(const std::string& _name,
                      const std::string& _description = std::string(),
                      OptionAttributes _attributes = NoOptionAttribute);
    /** Define an alternative name for an option that was previously defined by
    * defineOption().
    * \return Returns FALSE if there already is an option with the alternative
    * name or no option with the original name OR if a short option string
    * (length == 1) is a digit. In that case no action is performed.
    */
    bool defineOptionAlternative(const std::string& _original,
                                 const std::string& _alternative);
    /** Returns whether _name is a defined option. */
    bool isDefinedOption(const std::string& _name) const;
    /** Returns whether _name is an option that was found while parsing
    * the command line arguments with the parse() method. In other word: This
    * method returns true if the string is an option AND it was given on the
    * parsed command line.
    */
    bool foundOption(const std::string& _name) const;
    /** Define a help option. If this option is found a special error code is
    * returned by the parse method.
    * \attention If this method is called twice without an intermediate call
    * to the reset() method the previously set help option will remain a valid
    * option but is not detected as the special help option and will therefore
    * not cause the parse() method to return the special help error code.
    * \return Returns FALSE if there already is an option defined that equals
    * the short or long name.
    */
    bool setHelpOption(const std::string& _longname = "h",
                       const std::string& _shortname = "help",
                       const std::string& _descr = "");
    /** Returns the number of read arguments. Arguments are efined as beeing
    * neither options nor option values and are specified at the end of the
    * command line after all options and their values. */
    unsigned int arguments() const;
    /** Returns the Nth argument. See arguments().
    * \return Argument string or an empty string if there was no argument of
    * that id.
    */
    std::string argument(unsigned int _number) const;
    /** Get the complete argument vector. The order of the arguments in the
    * vector is the same as on the commandline.
    */
    const std::vector<std::string>& allArguments() const;
    /** Add an error code and its description to the command line parser.
    * This will do nothing more than adding an entry to the usage description.
    */
    void addErrorCode(int _code, const std::string& _descr = "");
    /** Set some string as a general description, that will be printed before
    * the list of available options.
    */
    void setIntroductoryDescription(const std::string& _descr);
    /** Parse the command line arguments for all known options and arguments.
    * \return Error code with parsing result.
    * \retval NoParserError Everything went fine.
    * \retval ParserUnknownOption Unknown option was found.
    * \retval ParserMissingValue A value to a given option is missing.
    * \retval ParserOptionAfterArgument Option after an argument detected. All
    * options have to given before the first argument.
    * \retval ParserMalformedMultipleShortOption Malformed short option string.
    * \retval ParserRequiredOptionMissing Required option is missing.
    * \retval ParserHelpRequested Help option detected.
    */
    ParserResults parse(int _argc, char ** _argv);
    /** Return the value of an option.
    * \return Value of a commandline options given by the name of the option or
    * an empty string if there was no such option or the option required no
    * value.
    */
    std::string optionValue(const std::string& _option) const;
    /** Reset the parser. Call this function if you want to parse another set of
    * command line arguments with the same parser object.
    */
    void reset();
    /** Returns the name of the option that was responsible for a parser error.
      * An empty string is returned if no error occured at all.
      */
    const std::string& errorOption() const;
    /** This method can be used to evaluate parser error codes and generate a
    * human-readable description. In case of a help request error code the
    * usage description as returned by usageDescription() is printed.
    */
    std::string parseErrorDescription(ParserResults _error_code) const;
    /** Returns a string with the usage descriptions for all options. The
     * description string is formated to fit into a terminal of width _width.*/
    std::string usageDescription(unsigned int _width = 80) const;

private:
    /** Returns the key of a defined option with name _name or -1 if such option
     * is not defined. */
    int optionKey( const std::string& _name ) const;
    /** Returns a list of option names that are all alternative names associated
     * with a single key value.
     */
    std::list<std::string> getAllOptionAlternatives(unsigned int _key) const;

    /** The current maximum key value for an option. */
    unsigned int max_key;
    /** Map option names to a numeric key. */
    String2KeyMap option2key;

    /** Map option key to option attributes. */
    Key2AttributeMap option2attribute;

    /** Map option key to option description. */
    Key2StringMap option2descr;

    /** Map option key to option value. */
    Key2StringMap option2value;

    /** Map error code to its description. */
    std::map<int, std::string> errorcode2descr;

    /** Vector of command line arguments. */
    ArgumentContainer argument_container;

    /** General description to be returned as first part of the generated help page. */
    std::string intro_description;

    /** Holds the key for the help option. */
    unsigned int help_option;

    /** Holds the name of the option that was responsible for a parser error.
    */
    std::string error_option;
}; // class ArgvParser


// Auxillary functions

/** Returns whether the given string is a valid (correct syntax) option string.
 * It has to fullfill the following criteria:
 *  1. minimum length is 2 characters
 *  2. Start with '-'
 *  3. if if minimal length -> must not be '--'
 *  4. first short option character must not be a digit (to distinguish negative numbers)
 */
bool isValidOptionString(const std::string& _string);

/** Returns whether the given string is a valid (correct syntax) long option string.
 * It has to fullfill the following criteria:
 *  1. minimum length is 4 characters
 *  2. Start with '--'
 */
bool isValidLongOptionString(const std::string& _string);

/** Splits option and value string if they are given in the form 'option=value'.
* \return Returns TRUE if a value was found.
*/
bool splitOptionAndValue(const std::string& _string, std::string& _option,
                         std::string& _value);

/** String tokenizer using standard C++ functions. Taken from here:
 * http://gcc.gnu.org/onlinedocs/libstdc++/21_strings/howto.html#3
 * Splits the string _in by _delimiters and store the tokens in _container.
 */
template <typename Container>
void splitString(Container& _container, const std::string& _in,
                 const char* const _delimiters = " \t\n")
{
    const std::string::size_type len = _in.length();
    std::string::size_type i = 0;

    while ( i < len )
    {
        // eat leading whitespace
        i = _in.find_first_not_of (_delimiters, i);
        if (i == std::string::npos)
            return;   // nothing left but white space

        // find the end of the token
        std::string::size_type j = _in.find_first_of (_delimiters, i);

        // push token
        if (j == std::string::npos)
        {
            _container.push_back (_in.substr(i));
            return;
        }
        else
            _container.push_back (_in.substr(i, j-i));

        // set up for next loop
        i = j + 1;
    }
}

/** Returns true if the character is a digit (what else?). */
bool isDigit(const char& _char);

/** Build a vector of integers from a string of the form:
* '1,3-5,14,25-20'. This string will be expanded to a list of positive
* integers with the following elements: 1,3,4,5,14,25,24,23,22,21,20.
* All of the expanded elements will be added to the provided list.
* \return Returns FALSE if there was any syntax error in the given string
* In that case the function stops at the point where the error occured.
* Only elements processed up to that point will be added to the expanded
* list.
* \attention This function can only handle unsigned integers!
*/
bool expandRangeStringToUInt(const std::string& _string,
                             std::vector<unsigned int>& _expanded);
/** Returns a copy of _str with whitespace removed from front and back. */
std::string trimmedString(const std::string& _str);

/** Formats a string of an arbitrary length to fit a terminal of width
* _width and to be indented by _indent columns.
*/
std::string formatString(const std::string& _string,
                         unsigned int _width,
                         unsigned int _indent = 0);

}
; // namespace CommandLineProcessing

using namespace std;
using namespace CommandLineProcessing;

ArgvParser::ArgvParser()
        : max_key(1),
        help_option(0) // must be smaller than max_key initially

{
    // nothing
}

ArgvParser::~ArgvParser()
{
    // nothing
}

void ArgvParser::reset()
{
    max_key = 1;
    option2key.clear();
    option2attribute.clear();
    option2descr.clear();
    option2value.clear();
    errorcode2descr.clear();
    argument_container.clear();
    intro_description.clear();
    error_option.clear();
    help_option = 0;
}

int ArgvParser::optionKey( const string& _name ) const
{
    String2KeyMap::const_iterator it = option2key.find(_name);

    // if not found
    if (it == option2key.end())
        return(-1);

    return(it->second);
}

bool ArgvParser::isDefinedOption( const string& _name ) const
{
    return(option2key.find(_name) != option2key.end());
}

bool ArgvParser::foundOption( const string & _name ) const
{
    int key = optionKey(_name);

    // not defined -> cannot by found
    if (key == -1)
        return(false);

    // return whether the key of the given option name is in the hash of the
    // parsed options.
    return(option2value.find(key) != option2value.end());
}

string ArgvParser::optionValue(const string& _option) const
{
    int key = optionKey(_option);

    // not defined -> cannot by found
    if (key == -1)
    {
        cerr << "ArgvParser::optionValue(): Requested value of an option the parser did not find or does not know." << endl;
        return("");
    }

    return(option2value.find(key)->second);
}

ArgvParser::ParserResults
ArgvParser::parse(int _argc, char ** _argv)
{
    bool finished_options = false; // flag whether an argument was found (options are passed)

    // loop over all command line arguments
    int i = 1; // argument counter
    while( i< _argc )
    {
        string argument = _argv[i];
        unsigned int key = 0;
        string option; // option name
        string value;  // option value

        // if argument is an option
        if (!isValidOptionString(argument))
        {
            // string is a real argument since values are processed elsewhere
            finished_options=true;
            argument_container.push_back(argument);
        }
        else // can be a long or multiple short options at this point
        {
            // check whether we already found an argument
            if (finished_options)
            {
                error_option = argument;
                return(ParserOptionAfterArgument); // return error code
            }
            // check for long options
            if (isValidLongOptionString(argument))
            {
                // handle long options

                // remove trailing '--'
                argument = argument.substr(2);
                // check for option value assignment 'option=value'
                splitOptionAndValue(argument, option, value);

                if (!isDefinedOption(option)) // is this a known option
                {
                    error_option = option; // store the option that caused the error
                    return(ParserUnknownOption); // return error code if not
                }

                // get the key of this option - now that we know that it is defined
                key = option2key.find(option)->second;
                if (key == help_option) // if help is requested return error code
                    return(ParserHelpRequested);

                // do we need to extract a value
                // AND a value is not already assigned from the previous step
                if ((option2attribute.find(key)->second & OptionRequiresValue) && value.empty())
                {
                    if (i+1 >= _argc) // are there arguments left?
                    {
                        error_option = option; // store the option that caused the error
                        return(ParserMissingValue); // the was no argument left although we need a value
                    }

                    string temp = _argv[i+1]; // get the next element
                    ++i; // increase counter now that we moved forward

                    if (isValidOptionString(temp))
                    {
                        error_option = option; // store the option that caused the error
                        return(ParserMissingValue);  // missing option value
                    }
                    value = temp; // assign value
                }
                // add option-value map entry
                option2value[key] = value;
            }
            else // handle short options
            {
                argument = argument.substr(1);   // remove trailing '-'

                // check for option value assignment 'option=value'
                if (splitOptionAndValue(argument, option, value))
                {
                    // there was an option <- value assignment
                    if (option.length() > 1)
                    {
                        error_option = option; // store the option that caused the error
                        return(ParserMalformedMultipleShortOption); // return error code if option has more than one character
                    }

                    if (!isDefinedOption(option)) // is this a known option
                    {
                        error_option = option; // store the option that caused the error
                        return(ParserUnknownOption); // return error code if not
                    }
                    key = option2key.find(option)->second; // get the key for the extracted option name

                    if (key == help_option) // if help is requested return error code
                        return(ParserHelpRequested);

                    // if value is still empty for some reason: we have an error
                    if ((option2attribute.find(key)->second & OptionRequiresValue) && value.empty())
                    {
                        error_option = option; // store the option that caused the error
                        return(ParserMissingValue);   // missing option value
                    }
                    else
                        // add option-value map entry
                        option2value[key] = value;
                }
                else // no '=' assignment: can be either multiple short options or
                    // something like '-s 4'
                {
                    // handle short options with value like '-s 4'
                    option.clear();
                    value.clear();

                    if (argument.length() == 1) // if a single short option
                    {
                        if (!isDefinedOption(argument)) // is this a known option
                        {
                            error_option = argument; // store the option that caused the error
                            return(ParserUnknownOption); // return error code if not
                        }
                        key = option2key.find(argument)->second; // get the key for the extracted option name

                        if (key == help_option) // if help is requested return error code
                            return(ParserHelpRequested);

                        // check if option needs a value and next arg is not an option
                        if ((option2attribute.find(key)->second & OptionRequiresValue))
                        {
                            if (i+1 >= _argc) // are there arguments left?
                            {
                                error_option = argument; // store the option that caused the error
                                return(ParserMissingValue); // the was no argument left although we need a value
                            }
                            string temp = _argv[i+1]; // get the next element
                            ++i; // increase counter now that we moved forward

                            if (isValidOptionString(temp))
                            {
                                error_option = argument; // store the option that caused the error
                                return(ParserMissingValue);  // missing option value
                            }
                            // add option-value map entry
                            option2value[key] = temp;

                        }
                        else // no value needed
                        {
                            option2value[key] = ""; // assign value
                        }
                    }
                    else // handle multiple short option like '-svxgh'
                    {
                        unsigned int short_option_counter = 0; // position in the multiple short option string
                        while( short_option_counter < argument.length() ) // parse the whole string
                        {
                            option = argument[short_option_counter]; // get the option character

                            if (!isDefinedOption(option)) // is this a known option
                            {
                                error_option = option; // store the option that caused the error
                                return(ParserUnknownOption); // return error code if not
                            }
                            key = option2key.find(option)->second; // get the key for the extracted option name

                            if (key == help_option) // if help is requested return error code
                                return(ParserHelpRequested);

                            option2value[key] = value;

                            ++short_option_counter; // advance one character forward
                        }
                    }
                }
            }
        }
        ++i;
    }

    map<unsigned int, OptionAttributes>::iterator it;
    for( it = option2attribute.begin(); it != option2attribute.end(); it++ )
    {
        // if the current option is required look if we got it
        if (it->second & OptionRequired)
        {
            // is the object missing
            if (option2value.find(it->first) == option2value.end())
            {
                // get the list of alternative names for this option
                list<string> alternatives = getAllOptionAlternatives(it->first);

                unsigned int count = 0;
                for( list<string>::const_iterator alt = alternatives.begin();
                        alt != alternatives.end();
                        ++alt )
                {
                    ++count;
                    // additional '-' for long options
                    if (alt->length() > 1)
                        error_option += "-";

                    error_option += "-" + *alt;

                    // alternatives to come?
                    if (count < alternatives.size())
                        error_option += ", "; // add separator
                }
                return(ParserRequiredOptionMissing);
            }
        }
    }

    return(NoParserError); // everthing went fine -> sucess
}

unsigned int ArgvParser::arguments() const
{
    return(argument_container.size());
}

string ArgvParser::argument(unsigned int _id) const
{
    if (_id >= arguments())
    {
        cerr << "ArgvParser::argument(): Request for non-existing argument." << endl;
        return ("");
    }
    else
        return(argument_container[_id]);
}

const vector<string>& ArgvParser::allArguments() const
{
    return(argument_container);
}

string ArgvParser::usageDescription(unsigned int _width) const
{
    string usage; // the usage description text

    if (intro_description.length())
        usage += formatString(intro_description, _width) + "\n\n";

    if (max_key>1) // if we have some options
        usage += formatString("Available options\n-----------------",_width) + "\n";

    // loop over all option attribute entries (which equals looping over all
    // different options (not option names)
    for (Key2AttributeMap::const_iterator it = option2attribute.begin();
            it != option2attribute.end();
            ++it)
    {
        string os; // temp string for the option

        // get the list of alternative names for this option
        list<string> alternatives = getAllOptionAlternatives(it->first);

        unsigned int count = 0;
        for( list<string>::const_iterator alt = alternatives.begin();
                alt != alternatives.end();
                ++alt )
        {
            ++count;
            // additional '-' for long options
            if (alt->length() > 1)
                os += "-";

            os += "-" + *alt;

            // note if the option requires a value
            if (option2attribute.find(it->first)->second & OptionRequiresValue)
                os += " <value>";

            // alternatives to come?
            if (count < alternatives.size())
                os += ", "; // add separator
        }

        // note if the option is required
        if (option2attribute.find(it->first)->second & OptionRequired)
            os += " [required]";

        usage += formatString(os, _width) + "\n";

        if (option2descr.find(it->first) != option2descr.end())
            usage += formatString(option2descr.find(it->first)->second, _width, 4);
        else
            usage += formatString("(no description)", _width, 4);

        // finally a little gap
        usage += "\n\n";
    }

    if (!errorcode2descr.size()) // if have no errorcodes
        return(usage);

    usage += formatString("Return codes\n-----------------", _width) + "\n";

    //   map<int, string>::const_iterator eit;
    for( std::map<int, std::string>::const_iterator alt = errorcode2descr.begin();
            alt != errorcode2descr.end();
            ++alt )
    {
        ostringstream code;
        code << alt->first;
        string label = formatString(code.str(), _width, 4);
        string descr = formatString(alt->second, _width, 10);
        usage += label + descr.substr(label.length()) + "\n";
    }

    return(usage);
}

const string& ArgvParser::errorOption( ) const
{
    return(error_option);
}

std::string ArgvParser::parseErrorDescription( ParserResults _error_code ) const
{
    string descr;

    switch (_error_code)
    {
    case ArgvParser::NoParserError:
        // no error -> nothing to do
        break;
    case ArgvParser::ParserUnknownOption:
        descr = "Unknown option: '" + errorOption() + "'";
        break;
    case ArgvParser::ParserMissingValue:
        descr = "Missing required value for option: '" + errorOption()+ "'";
        break;
    case ArgvParser::ParserOptionAfterArgument:
        descr = "Misplaced option '" + errorOption() + "' detected. All option have to be BEFORE the first argument";
        break;
    case ArgvParser::ParserMalformedMultipleShortOption:
        descr = "Malformed short-options: '" + errorOption() + "'";
        break;
    case ArgvParser::ArgvParser::ParserRequiredOptionMissing:
        descr = "Required option missing: '" + errorOption() + "'";
        break;
    case ArgvParser::ParserHelpRequested: // help
        descr = usageDescription();
        break;
    default:
        cerr << "ArgvParser::documentParserErrors(): Unknown error code" << endl;
    }

    return(descr);
}

bool ArgvParser::defineOption( const string & _name,
                               const string& _descr,
                               OptionAttributes _attrs)
{
    // do nothing if there already is an option of this name
    if (isDefinedOption(_name))
    {
        cerr << "ArgvParser::defineOption(): The option label equals an already defined option." << endl;
        return(false);
    }

    // no digits as short options allowed
    if (_name.length() == 1 && isDigit(_name[0]))
    {
        cerr << "ArgvParser::defineOption(): Digits as short option labels are not allowd." << endl;
        return(false);
    }

    option2key[_name] = max_key;     // give the option a unique key

    // store the option attributes
    option2attribute[max_key] = _attrs;

    // store the option description if there is one
    if (_descr.length())
        option2descr[max_key] = _descr;

    // inc the key counter
    ++max_key;

    return(true);
}

bool ArgvParser::defineOptionAlternative( const string & _original,
        const string & _alternative )
{
    // do nothing if there already is no option of this name
    if (!isDefinedOption(_original))
    {
        cerr << "ArgvParser::defineOptionAlternative(): Original option label is not a defined option." << endl;
        return(false);
    }

    // AND no digits as short options allowed
    if (_alternative.length() == 1 && isDigit(_alternative[0]))
    {
        cerr << "ArgvParser::defineOptionAlternative(): Digits as short option labels are not allowd." << endl;
        return(false);
    }

    // AND do nothing if there already is an option with the alternativ name
    if (isDefinedOption(_alternative))
    {
        cerr << "ArgvParser::defineOptionAlternative(): The alternative option label equals an already defined option." << endl;
        return(false);
    }

    option2key[_alternative] = optionKey(_original);

    return(true);
}


bool ArgvParser::setHelpOption(const string& _shortname,
                               const string& _longname,
                               const string& _descr)
{
    // do nothing if any name is already in use
    if (isDefinedOption(_shortname) || isDefinedOption(_longname))
    {
        cerr << "ArgvParser::setHelpOption(): Short or long help option label equals an already defined option." << endl;
        return(false);
    }

    // define the help option's short name and the alternative
    // longname
    defineOption(_shortname, _descr, NoOptionAttribute);
    defineOptionAlternative(_shortname, _longname);

    help_option = max_key-1; // store the key in a special member

    return(true);
}

void ArgvParser::addErrorCode(int _code, const string& _descr)
{
    errorcode2descr[_code] = _descr;
}

void ArgvParser::setIntroductoryDescription(const string& _descr)
{
    intro_description = _descr;
}

list<string> ArgvParser::getAllOptionAlternatives( unsigned int _key ) const
{
    // keys go here
    list<string> keys;
    // for all container elements
    for( map<string, unsigned int>::const_iterator it = option2key.begin();
            it != option2key.end();
            ++it )
    {
        if (it->second == _key)
            keys.push_back(it->first);
    }
    return(keys);
}

bool CommandLineProcessing::isDigit(const char& _char)
{
    if (_char == '0' || _char == '1' || _char == '2' || _char == '3'
            || _char == '4' || _char == '5' || _char == '6' || _char == '7'
            || _char == '8' || _char == '9')
        return(true);
    else
        return(false);
}

bool CommandLineProcessing::isValidOptionString(const string& _string)
{
    // minimal short option length is 2
    if (_string.length() < 2)
        return(false);

    // is it an option (check for '-' as first character)
    if (_string.compare(0, 1, "-"))
        return(false);

    // not an option if just '--'
    if (_string.length() == 2 && _string == "--")
        return(false);

    // it might still be a negative number
    // (but not if there is no digit afterwards)
    if (isDigit(_string[1]))
        return(false);

    // let's consider this an option
    return(true);
}

bool CommandLineProcessing::isValidLongOptionString(const string& _string)
{
    if (_string.length() < 4) // must be at least '--??'
        return(false);

    // is it an option (check for '--')
    if (_string.compare(0, 2, "--"))
        return(false);
    else
        return(true);
}

bool CommandLineProcessing::splitOptionAndValue(const string& _string,
        string& _option, string& _value)
{
    // string token container
    std::vector<string> tokens;

    // split string by '=' delimiter
    splitString(tokens, _string, "=");

    // check for option value assignment 'option=value'
    if (tokens.size() < 2)
    {
        _option = _string; // the option is the whole string
        return(false);
    }

    // separate option and value
    _option = tokens[0];

    // concat all remaining tokens to the value string
    for (unsigned int i=1; i<tokens.size(); ++i)
    {
        _value.append(tokens[i]);
    }

    return(true);
}

string CommandLineProcessing::trimmedString( const std::string & _str )
{
    // no string no work
    if(_str.length() == 0)
        return _str;

    string::size_type start_pos = _str.find_first_not_of(" \a\b\f\n\r\t\v");
    string::size_type end_pos = _str.find_last_not_of(" \a\b\f\n\r\t\v");

    // check whether there was any non-whitespace
    if (start_pos == string::npos)
        return("");

    return string(_str, start_pos, end_pos - start_pos + 1);
}

bool CommandLineProcessing::expandRangeStringToUInt( const std::string & _string,
        std::vector< unsigned int > & _expanded )
{
    list<string> tokens;
    // split string by delimiter
    splitString(tokens, _string, ",");

    // loop over all entries
    for(list<string>::const_iterator it = tokens.begin(); it != tokens.end(); it++)
    {
        const string& entry = *it; // convenience reference

#ifdef ARGVPARSER_DEBUG

        cout << "TOKEN: " << entry << endl;
#endif

        // if range was given
        if (entry.find("-") != string::npos)
        {
            // split into upper and lower border
            list<string> range_borders;
            splitString(range_borders, entry, "-");

            // fail if insane range spec
            if (range_borders.size() != 2)
                return(false);

            int first = atoi(range_borders.begin()->c_str());
            int second = atoi((++range_borders.begin())->c_str());

            // write id in increasing order
            if (first <= second)

            {
                for (int j=first; j<=second; ++j)
                {
                    _expanded.push_back(j);
                }
            }
            else // write id in decreasing order
            {
                for (int k=first; k>=second; k--)
                {
                    _expanded.push_back(k);
                }
            }
        }
        else // single number was given
            _expanded.push_back(atoi(entry.c_str())); // store id
    }

    return(true);
}

std::string CommandLineProcessing::formatString(const std::string& _string,
        unsigned int _width,
        unsigned int _indent)
{
    // if insane parameters do nothing
    if (_indent >= _width)
        return(_string);

    // list of lines of the formated string
    list<string> lines;

    // current position in the string
    unsigned int pos = 0;

    // till the end of the string
    while (pos < _string.length())
    {
        // get the next line of the string
        string line = _string.substr(pos, _width - _indent );

#ifdef ARGVPARSER_DEBUG

        cout << "EXTRACT: '" << line << "'" << endl;
#endif

        // check for newlines in the line and break line at first occurence (if any)
        string::size_type first_newline = line.find_first_of("\n");
        if (first_newline != string::npos)
        {
            line = line.substr(0, first_newline);
        }

        // we need to check for possible breaks within words only if the extracted
        // line spans the whole allowed width
        bool check_truncation = true;
        if (line.length() < _width - _indent)
            check_truncation = false;

        // remove unecessary whitespace at front and back
        line = trimmedString(line);

#ifdef ARGVPARSER_DEBUG

        cout << "TRIMMED: '" << line << "'" << endl;
#endif

        // only perform truncation if there was enough data for a full line
        if (!check_truncation)
            pos += line.length() + 1;
        else
        {
            // look for the last whitespace character
            string::size_type last_white_space = line.find_last_of(" \a\b\f\n\r\t\v");

            if (last_white_space != string::npos) // whitespace found!
            {
                // truncated the line at the last whitespace
                line = string(line, 0, last_white_space);
                pos += last_white_space + 1;
            }
            else // no whitespace found
                // rude break! we can leave the line in its current state
                pos += _width - _indent;
        }

        if (!line.empty())
        {
#ifdef ARGVPARSER_DEBUG
            cout << "UNINDEN: '" << line << "'" << endl;
#endif

            if (_indent)
                line.insert(0, _indent, ' ');

#ifdef ARGVPARSER_DEBUG

            cout << "INDENT: '" << line << "'" << endl;
#endif

            lines.push_back(line);
        }
    }

    // concat the formated string
    string formated;
    bool first = true;
    // for all lines
    for (list<string>::iterator it = lines.begin(); it != lines.end(); ++it)
    {
        // prefix with newline if not first
        if (!first)
            formated += "\n";
        else
            first = false;

        formated += *it;
    }
    return(formated);
}



#endif // __CMDLINEPARSER_H
