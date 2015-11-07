/************************************************************************/
/*                                                                      */
/*               Copyright 2014-2015 by Ullrich Koethe                  */
/*                                                                      */
/*    This file is part of the VIGRA2 computer vision library.          */
/*    The VIGRA2 Website is                                             */
/*        http://ukoethe.github.io/vigra2                               */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/

#pragma once

#ifndef VIGRA2_ERROR_HXX
#define VIGRA2_ERROR_HXX

#include <stdexcept>
#include <sstream>
#include <string>
#include "config.hxx"

namespace vigra {

/** \page ErrorReporting Error Reporting
    Exceptions and assertions provided by VIGRA

    <b>\#include</b> \<vigra/error.hxx\>

    VIGRA defines the following exception class:

    \code
    namespace vigra {
        class ContractViolation : public std::runtime_error;
    }
    \endcode

    The following associated macros throw the corresponding exception if
    their PREDICATE evaluates to '<TT>false</TT>':

    \code
    vigra_precondition(PREDICATE, MESSAGE);
    vigra_postcondition(PREDICATE, MESSAGE);
    vigra_invariant(PREDICATE, MESSAGE);
    \endcode

    The MESSAGE is passed to the exception and can be retrieved via
    the overloaded member function '<TT>exception.what()</TT>'. If the compiler
    flag '<TT>NDEBUG</TT>' is <em>not</em> defined, the file name and line number of
    the error are automatically included in the message. The macro

    \code
    vigra_assert(PREDICATE, MESSAGE);
    \endcode

    is identical to <tt>vigra_precondition()</tt> except that it is completely removed
    when '<TT>NDEBUG</TT>' is defined. This is useful for tests that are only needed during
    debugging, such as array index bound checking. The following macro

    \code
    vigra_fail(MESSAGE);
    \endcode

    unconditionally throws a '<TT>std::runtime_error</TT>' constructed from the message
    (along with file name and line number, if NDEBUG is not set).

    <b> Usage:</b>

    Include-File:
    \<vigra/error.hxx\>
    <p>
    Namespace: vigra (except for the macros, of course)

    \code
    int main(int argc, char ** argv)
    {
        try
        {
            const char* input_file_name = argv[1];

            // read input image
            vigra::ImageImportInfo info(input_file_name);

            // fail if input image is not grayscale
            vigra_precondition(info.isGrayscale(), "Input image must be grayscale");

            ...// process image
        }
        catch (std::exception & e)
        {
            std::cerr << e.what() << std::endl; // print message
            return 1;
        }

        return 0;
    }
    \endcode
**/
class ContractViolation
: public std::runtime_error
{
  public:

    static std::string createMessage(char const * prefix, char const * message)
    {
        std::ostringstream s;
        s << "\n" << prefix << "\n" << message << "\n";
        return s.str();
    }

    static std::string createMessage(char const * prefix, char const * message,
                                     char const * file, int line)
    {
        std::ostringstream s;
        s << "\n" << prefix << "\n" << message << "\n("
          << file << ":" << line << ")\n";
        return s.str();
    }

    ContractViolation(char const * prefix, char const * message)
    : std::runtime_error(createMessage(prefix, message))
    {}

    ContractViolation(char const * prefix, char const * message,
                      char const * file, int line)
    : std::runtime_error(createMessage(prefix, message, file, line))
    {}
};


//#ifndef NDEBUG

#if 1

inline
void throw_contract_error(char const * prefix, char const * message,
                          char const * file, int line)
{
    throw vigra::ContractViolation(prefix, message, file, line);
}

inline
void throw_contract_error(char const * prefix, std::string message,
                          char const * file, int line)
{
    throw vigra::ContractViolation(prefix, message.c_str(), file, line);
}

inline
void throw_runtime_error(char const * message, char const * file, int line)
{
    std::ostringstream what;
    what << "\n" << message << "\n(" << file << ":" << line << ")\n";
    throw std::runtime_error(what.str());
}

inline
void throw_runtime_error(std::string message, char const * file, int line)
{
    std::ostringstream what;
    what << "\n" << message << "\n(" << file << ":" << line << ")\n";
    throw std::runtime_error(what.str());
}

// place throw_contract_error() into the else branch of an if
// in order to avoid a potentially expensive message construction
// when there is no error
#define vigra_precondition(PREDICATE, MESSAGE) \
    if((PREDICATE)) {} else \
        vigra::throw_contract_error("Precondition violation!", MESSAGE, __FILE__, __LINE__)

#define vigra_assert(PREDICATE, MESSAGE) \
    vigra_precondition(PREDICATE, MESSAGE)

#define vigra_postcondition(PREDICATE, MESSAGE) \
    if((PREDICATE)) {} else \
        vigra::throw_contract_error("Postcondition violation!", MESSAGE, __FILE__, __LINE__)

#define vigra_invariant(PREDICATE, MESSAGE) \
    if((PREDICATE)) {} else \
        vigra::throw_contract_error("Invariant violation!", MESSAGE, __FILE__, __LINE__)

#define vigra_fail(MESSAGE) \
    vigra::throw_runtime_error(MESSAGE, __FILE__, __LINE__)

#else // NDEBUG

inline
void throw_contract_error(char const * message)
{
    throw vigra::ContractViolation(message);
}

#define vigra_precondition(PREDICATE, MESSAGE) \
    if((PREDICATE)) {} else \
        vigra::throw_contract_error("Precondition violation!", MESSAGE)

#define vigra_assert(PREDICATE, MESSAGE)

#define vigra_postcondition(PREDICATE, MESSAGE) \
    if((PREDICATE)) {} else \
        vigra::throw_contract_error("Postcondition violation!", MESSAGE)

#define vigra_invariant(PREDICATE, MESSAGE) \
    if((PREDICATE)) {} else \
        vigra::throw_contract_error("Invariant violation!", MESSAGE)

#define vigra_fail(MESSAGE) \
    throw std::runtime_error(MESSAGE)

#endif // NDEBUG

} // namespace vigra

#endif // VIGRA2_ERROR_HXX
