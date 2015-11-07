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

#ifndef VIGRA2_AXISTAGS_HXX_HXX
#define VIGRA2_AXISTAGS_HXX_HXX

#include "config.hxx"
#include "tinyarray.hxx"
#include <string>
#include <iostream>

namespace vigra {

namespace tags {

enum AxisType { Channels = 1,
                Space = 2,
                Angle = 4,
                Time = 8,
                Frequency = 16,
                Edge = 32,
                UnknownAxisType = 64,
                NonChannel = Space | Angle | Time | Frequency | UnknownAxisType,
                AllAxes = 2*UnknownAxisType-1 };

    // tags are arranged in sorting order
enum AxisTag  { axis_c,  // channel axis
                axis_n,  // node map for a graph
                axis_x,  // spatial x-axis
                axis_y,  // spatial y-axis
                axis_z,  // spatial z-axis
                axis_t,  // time axis
                axis_fx, // Fourier transform of x-axis
                axis_fy, // Fourier transform of y-axis
                axis_fz, // Fourier transform of z-axis
                axis_ft, // Fourier transform of t-axis
                axis_e,  // edge map for a graph
                axis_unknown
              };

const char * AxisTagKeys[] = { "c",
                               "n",
                               "x",
                               "y",
                               "z",
                               "t",
                               "x",
                               "y",
                               "z",
                               "t",
                               "e",
                               "?"
                              };

AxisType AxisTagTypes[] = { Channels,                    // c
                            Space,                       // n
                            Space,                       // x
                            Space,                       // y
                            Space,                       // z
                            Time,                        // t
                            AxisType(Space | Frequency), // fx
                            AxisType(Space | Frequency), // fy
                            AxisType(Space | Frequency), // fz
                            AxisType(Time | Frequency),  // ft
                            Edge,                        // e
                            UnknownAxisType              // unknown
                          };
} // namespace tags

using tags::AxisType;
using tags::AxisTag;

class AxisInfo
{
  public:

    AxisInfo(AxisTag tag=tags::axis_unknown,
             double resolution = 0.0, std::string description = "")
    : key_(tags::AxisTagKeys[tag]),
      description_(description),
      resolution_(resolution),
      flags_(tags::AxisTagTypes[tag])
    {}

    AxisInfo(std::string key, AxisType typeFlags = tags::UnknownAxisType,
             double resolution = 0.0, std::string description = "")
    : key_(key),
      description_(description),
      resolution_(resolution),
      flags_(typeFlags)
    {}

    std::string key() const
    {
        return key_;
    }

    std::string description() const
    {
        return description_;
    }

    void setDescription(std::string const & description)
    {
        description_ = description;
    }

    double resolution() const
    {
        return resolution_;
    }

    void setResolution(double resolution)
    {
        resolution_ = resolution;
    }

    AxisType typeFlags() const
    {
        return flags_ == 0
                  ? tags::UnknownAxisType
                  : flags_;
    }

    bool isUnknown() const
    {
        return isType(tags::UnknownAxisType);
    }

    bool isSpatial() const
    {
        return isType(tags::Space);
    }

    bool isTemporal() const
    {
        return isType(tags::Time);
    }

    bool isChannel() const
    {
        return isType(tags::Channels);
    }

    bool isFrequency() const
    {
        return isType(tags::Frequency);
    }

    bool isEdge() const
    {
        return isType(tags::Edge);
    }

    bool isAngular() const
    {
        return isType(tags::Angle);
    }

    bool isType(AxisType type) const
    {
        return (typeFlags() & type) != 0;
    }

    std::string repr() const
    {
        std::string res("AxisInfo: '");
        res += key() + "' (type:";
        if(isUnknown())
        {
            res += " none";
        }
        else
        {
            if(isChannel())
                res += " Channels";
            if(isSpatial())
                res += " Space";
            if(isTemporal())
                res += " Time";
            if(isEdge())
                res += " Edge";
            if(isAngular())
                res += " Angle";
            if(isFrequency())
                res += " Frequency";
        }
        if(resolution_ > 0.0)
        {
            res += ", resolution=";
            res += std::to_string(resolution_);
        }
        res += ")";
        if(description_ != "")
        {
            res += " ";
            res += description_;
        }
        return res;
    }

    AxisInfo toFrequencyDomain(unsigned int size = 0, int sign = 1) const
    {
        AxisType type;
        if(sign == 1)
        {
            vigra_precondition(!isFrequency(),
                "AxisInfo::toFrequencyDomain(): axis is already in the Fourier domain.");
            type = AxisType(tags::Frequency | flags_);
        }
        else
        {
            vigra_precondition(isFrequency(),
                "AxisInfo::fromFrequencyDomain(): axis is not in the Fourier domain.");
            type = AxisType(~tags::Frequency & flags_);
        }
        AxisInfo res(key_, type, 0.0, description_);
        if(resolution_ > 0.0 && size > 0u)
            res.resolution_ = 1.0 / (resolution_ * size);
        return res;
    }

    AxisInfo fromFrequencyDomain(unsigned int size = 0) const
    {
        return toFrequencyDomain(size, -1);
    }

    bool compatible(AxisInfo const & other) const
    {
        return isUnknown() || other.isUnknown() ||
               ((typeFlags() & ~tags::Frequency) == (other.typeFlags() & ~tags::Frequency) &&
                 key() == other.key());
    }

    bool operator==(AxisInfo const & other) const
    {
        return typeFlags() == other.typeFlags() && key() == other.key();
    }

    bool operator!=(AxisInfo const & other) const
    {
        return !operator==(other);
    }

    // the primary ordering is according to axis type:
    //     Channels < Space < Angle < Time < Frequency < Edge < Unknown
    // the secondary ordering is the lexicographic ordering of the keys
    //     "x" < "y" < "z"
    bool operator<(AxisInfo const & other) const
    {
        return (typeFlags() < other.typeFlags()) ||
                (typeFlags() == other.typeFlags() && key() < other.key());
    }

    bool operator<=(AxisInfo const & other) const
    {
        return !(other < *this);
    }

    bool operator>(AxisInfo const & other) const
    {
        return other < *this;
    }

    bool operator>=(AxisInfo const & other) const
    {
        return !(*this < other);
    }

    std::string key_, description_;
    double resolution_;
    AxisType flags_;
};

inline
std::string to_string(AxisInfo const & i)
{
    return i.repr();
}

inline
std::ostream & operator<<(std::ostream & o, AxisInfo const & i)
{
    o << to_string(i);
    return o;
}

// template <int N>
// struct DefaultAxisTags
// {
    // static TinyArray<AxisTag, N> exec()
    // {
        // using namespace tags;
        // static const AxisTag std[] = { t, z, y, x, c };
    // }
// };

} // namespace vigra

#endif // VIGRA2_AXISTAGS_HXX_HXX
