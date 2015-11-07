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

#ifndef VIGRA2_TIMING_HXX
#define VIGRA2_TIMING_HXX

#ifndef VIGRA_NO_TIMING

#include <iostream>
#include <sstream>
#include <vector>
#include <chrono>

/** \page TimingMacros  Timing macros for runtime measurements

<b>\#include</b> \<vigra/timing.hxx\>

These macros allow to perform execution speed measurements. Results are reported
in <i>milliseconds</i>. This uses <tt>std::chrono::high_resolution_clock</tt>,
but the actual accuracy is platform dependent.

Basic usage:
\code
   void time_it()
   {
       USETICTOC

       TIC
        ...    code to be timed
       TOC
        ...    untimed code
       TIC
        ...    other code to be timed
       TOC
   }
\endcode

Instead of TOC, which outputs the time difference to std::cerr,
you may use TOCN (the time difference in <i>msec</i> as a double)
or TOCS (the time difference as a std::string).

Alternatively, you can perform nested timing like so:
\code
   void time_it()
   {
       USE_NESTED_TICTOC

       TICPUSH
        ...         code to be timed
           TICPUSH
            ...         nested code to be timed
           TOC          print time for nested code
        ...         more code to be timed
       TOC          print total time
   }
\endcode

*/

/** \file timing.hxx  Timing macros for runtime measurements

  This header defines timing macros for runtime measurements. See \ref TimingMacros for examples.

  \def USETICTOC
  Enable timing using TIC/TOC* pairs. This macro defines temporary storage for the timing data, so it needs to precede the TIC/TOC macros in their context.
  \hideinitializer

  \def USE_NESTED_TICTOC
  Enable timing using TICPUSH/TOC* pairs. This macro defines temporary storage for the timing data, so it needs to precede the TIC/TOC macros in their context.
  \hideinitializer

  \def TIC
  Start timing. Requires USE_TICTOC to be defined in the current context.
  \hideinitializer

  \def TOC
  Stop timing and output result (the time difference w.r.t. the last TIC or TICPUSH
  instance) to std::cerr.
  \hideinitializer

  \def TICPUSH
  Start timing, possibly a nested block of code within some other timed code block.
  Requires USE_NESTED_TICTOC to be defined once in the current context.
  \hideinitializer

  \def TOCN
  Stop timing. This macro evaluates to the time difference (w.r.t. the last TIC
  or TICPUSH) in msec as a double.
  \hideinitializer

  \def TOCS
  Stop timing. This macro evaluates to the time difference (w.r.t. the last TIC
  or TICPUSH) as a std::string (including units).
  \hideinitializer

  \def TICTOCLOOP_BEGIN(inner_repetitions,outer_repetitions)
  Executes the code block up to TICTOCLOOP_END outer_repetitions x
  inner_repetitions times. The measurement is averaged over the
  inner_repetitions, and the best result of the outer_repetitions is
  reported to std::cerr.
  \hideinitializer

  \def TICTOCLOOP_END
  Ends the timing loop started with the TICTOCLOOP_BEGIN macro
  and outputs the result.
  \hideinitializer
*/

namespace timing_detail {

using namespace std::chrono;

typedef high_resolution_clock::time_point time_point;

inline double tic_toc_diff_num(time_point const & tic)
{
    duration<double> span =  duration_cast<duration<double> >(high_resolution_clock::now() - tic);
    return span.count()*1000.0;
}

inline std::string tic_toc_diff_string(time_point const & tic)
{
    double diff = tic_toc_diff_num(tic);
    std::stringstream s;
    s << diff << " msec";
    return s.str();
}

inline void tic_toc_diff(time_point const & tic)
{
    double diff = tic_toc_diff_num(tic);
    std::cerr << diff << " msec" << std::endl;
}

inline double tic_toc_diff_num(std::vector<time_point> & tic)
{
    double res = tic_toc_diff_num(tic.back());
    tic.pop_back();
    return res;
}

inline std::string tic_toc_diff_string(std::vector<time_point> & tic)
{
    std::string res = tic_toc_diff_string(tic.back());
    tic.pop_back();
    return res;
}

inline void tic_toc_diff(std::vector<time_point> & tic)
{
    tic_toc_diff(tic.back());
    tic.pop_back();
}

} // namespace timing_detail

#define USETICTOC timing_detail::time_point tic_timer;
#define TIC  tic_timer = std::chrono::high_resolution_clock::now();
#define TOC  timing_detail::tic_toc_diff       (tic_timer);
#define TOCN timing_detail::tic_toc_diff_num   (tic_timer)
#define TOCS timing_detail::tic_toc_diff_string(tic_timer)
#define USE_NESTED_TICTOC std::vector<timing_detail::time_point> tic_timer;
#define TICPUSH tic_timer.push_back(std::chrono::high_resolution_clock::now());

// TICTOCLOOP runs the body inner_repetitions times, and minimizes the result over a number of outer_repetitions runs,
//  outputting the final minimal average to std::cerr
// We enclose the loop in a dummy do { ... } while(false) in order to make this a true single statement
//  (instead of just a scope).
#define TICTOCLOOP_BEGIN(inner_repetitions,outer_repetitions) \
    do { \
    USETICTOC \
        double tictoc_best_, tictoc_inner_repetitions_=inner_repetitions; size_t tictoc_outer_repetitions_=outer_repetitions; \
        for (size_t tictoccounter_=0; tictoccounter_<tictoc_outer_repetitions_; ++tictoccounter_) { \
        TIC \
        for (size_t tictocinnercounter_=0; tictocinnercounter_<inner_repetitions; ++tictocinnercounter_) { \


#define TICTOCLOOP_END \
                } \
        const double tictoc_cur_ = TOCN; \
                if ((tictoccounter_==0) || (tictoc_cur_ < tictoc_best_)) \
            tictoc_best_ = tictoc_cur_; \
        } \
        std::cerr << tictoc_best_/tictoc_inner_repetitions_ \
             << " msec (best-of-" << tictoc_outer_repetitions_ << ")" << std::endl; \
    } while(false);

#else // VIGRA_NO_TIMING

#define USETICTOC
#define TIC
#define TOC
#define TOCN 0.0
#define TICS ""
#define USE_NESTED_TICTOC
#define TICPUSH
#define TICTOCLOOP_BEGIN(inner_repetitions,outer_repetitions)  do {
#define TICTOCLOOP_END } while(false);

#endif // VIGRA_NO_TIMING

#endif // VIGRA2_TIMING_HXX
