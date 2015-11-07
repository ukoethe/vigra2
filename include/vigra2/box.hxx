/************************************************************************/
/*                                                                      */
/*     Copyright 2009-2016 by Ullrich Koethe and Hans Meine             */
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

#ifndef VIGRA2_BOX_HXX
#define VIGRA2_BOX_HXX

#include "numeric_traits.hxx"
#include "tinyarray.hxx"
#include <type_traits>

namespace vigra {

namespace detail {

// RangePolicy used for floating point coordinate types
struct EndInsidePolicy
{
    template<class VECTOR>
    static inline bool nonEmpty(VECTOR const & b, VECTOR const & e)
    {
        return allLessEqual(b, e);
    }

    template<class VECTOR>
    static inline VECTOR const & pointEnd(VECTOR const & p)
    {
        return p;
    }
};

// RangePolicy used for integer coordinate types
struct EndOutsidePolicy
{
    template<class VECTOR>
    static inline bool nonEmpty(VECTOR const & b, VECTOR const & e)
    {
        return allLess(b, e);
    }

    template<class VECTOR>
    static inline VECTOR pointEnd(VECTOR const & p)
    {
        return p+1;
    }
};

} // namespace vigra::detail

/** \addtogroup RangesAndPoints */
//@{
   /** \brief Represent an n-dimensional box as a (lower, upper) pair.
     * Depending on the value type, upper() is considered to be
     * outside the box (as in the STL, for integer types), or
     * inside (for floating point types).  shape() will always be
     * upper() - lower().
     *
     * <b>\#include</b> \<vigra2/box.hxx\> <br/>
     * Namespace: vigra
     */
template<int DIMENSION=runtime_size, class VALUETYPE=ArrayIndex>
class Box
{
  public:
        /** STL-compatible definition of coordinate valuetype
         */
    typedef VALUETYPE value_type;

        /** Promoted coordinate valuetype, used for volume()
         */
    typedef typename NumericTraits<VALUETYPE>::Promote VolumeType;

        /** Vector type used for lower() and upper()
         */
    typedef TinyArray<VALUETYPE, DIMENSION> Vector;

        /** Ressult type of the erase() operation
         */
    typedef typename std::conditional<DIMENSION==runtime_size,
                                      Box, Box<DIMENSION-1, VALUETYPE> >::type
            EraseResult;

        /** Ressult type of the insert() operation
         */
    typedef typename std::conditional<DIMENSION==runtime_size,
                                      Box, Box<DIMENSION+1, VALUETYPE> >::type
            InsertResult;

    enum { Dimension = DIMENSION };

  protected:
    Vector lower_, upper_;

        /** Range policy (EndInsidePolicy/EndOutsidePolicy, depending on valuetype)
         */
    typedef typename std::conditional<std::is_integral<VALUETYPE>::value,
                        detail::EndOutsidePolicy,
                        detail::EndInsidePolicy>::type RangePolicy;

  public:
        /** Construct an empty box (empty() will return true).
         * (Internally, this will initialize all dimensions with the
         * empty range [1..0].)
         */
    Box()
    {
        if(DIMENSION != runtime_size)
            lower_ += 1;
    }

        /** Construct a box representing the given range.  Depending
         * on the value type, upper() is considered to be outside the
         * box (as in the STL, for integer types), or inside (for
         * floating point types).
         */
    Box(Vector const &lower, Vector const &upper)
    : lower_(lower)
    , upper_(upper)
    {
        VIGRA_ASSERT_RUNTIME_SIZE(DIMENSION, lower.size() == upper.size(),
            "Box(): size mismatch between lower and upper bound.");
    }

        /** Construct a box of given shape at the origin (i.e. upper() ==
         * shape()).
         */
    explicit Box(Vector const &shape)
    : lower_(shape, 0)
    , upper_(shape)
    {}

        /** Get the dimension of the box. If the box is statically
         *  sized, this equals the <tt>DIMENSION</tt> template paramter.
         *  For a dynamically sized box (<tt>DIMENSION=runtime_dimension</tt>)
         *  it equals the size of the bound vectors.
         */
    int ndim() const
    {
        return lower_.size();
    }

        /** Get lower bound (i.e. smallest coordinates for each
         * dimension).  This is the first point (scan-order wise)
         * which is considered to be "in" the box.
         */
    Vector const & lower() const
    {
        return lower_;
    }

        /** Get upper bound (i.e. coordinates higher than lower() in
         * each dimension for non-empty boxes).  This is lower() +
         * shape(), and depending on the valuetype (float/int), this is
         * the last point within or the first point outside the box,
         * respectively.
         */
    Vector const & upper() const
    {
        return upper_;
    }


        /** Change bounds. If the box has <tt>DIMENSION=runtime_size</tt>,
            the function can also change the dimension of the box.
         */
    void set(Vector const &lower, Vector const &upper)
    {
        VIGRA_ASSERT_RUNTIME_SIZE(DIMENSION, lower.size() == upper.size(),
            "Box::set(): size mismatch between lower and upper bound.");
        lower_ = lower;
        upper_ = upper;
    }

        /** Change lower() without changing upper(), changing shape()
         * accordingly.
         */
    void setLower(Vector const &lower)
    {
        VIGRA_ASSERT_RUNTIME_SIZE(DIMENSION, lower.size() == ndim(),
            "Box::setLower(): size mismatch between ndim() and new lower bound.");
        lower_ = lower;
    }

        /** Change upper() without changing lower(), which will change
         * the shape() most probably.
         */
    void setUpper(Vector const &upper)
    {
        VIGRA_ASSERT_RUNTIME_SIZE(DIMENSION, upper.size() == ndim(),
            "Box::setUpper(): size mismatch between ndim() and new upper bound.");
        upper_ = upper;
    }

        /** Move the whole box so that the given point will be
         * lower() afterwards.
         */
    void moveTo(Vector const &newBegin)
    {
        upper_ += newBegin - lower_;
        lower_  = newBegin;
    }

        /** Move the whole box by the given offset.
         * (Equivalent to operator+=)
         */
    void moveBy(Vector const &offset)
    {
        lower_ += offset;
        upper_ += offset;
    }

        /** Create a Box with the given dimension erased.
         */
    EraseResult erase(int dim) const
    {
        return EraseResult(lower_.erase(dim), upper_.erase(dim));
    }

        /** Create a Box with a new dimension added at the given index.
         *  The dimensions at and above that index are shifted up by one.
         */
    InsertResult insert(int index, VALUETYPE lower, VALUETYPE upper) const
    {
        return InsertResult(lower_.insert(index, lower), upper_.insert(index, upper));
    }

        /** Determine and return the area of this box. That is,
         * if this rect empty(), returns zero, otherwise returns the
         * product of the extents in each dimension.
         */
    VolumeType volume() const
    {
        if(empty())
            return 0;
        return prod(shape());
    }

        /** Determine and return the shape of this box. The shape
         * might be zero or even negative in one or more dimensions,
         * and if so, empty() will return true.
         */
    Vector shape() const
    {
        return upper_ - lower_;
    }

        /** Reshape this box to the given extents. This will
         * change upper() only.
         */
    void setShape(Vector const &shape)
    {
        upper_ = lower_ + shape;
    }

        /** Increase the shape of the box by the given
         * offset. This will move upper() only. (If any of offset's
         * components is negative, the box will get smaller
         * accordingly.)
         */
    void addShape(Vector const &offset)
    {
        upper_ += offset;
    }

        /** Adds a border of the given width around the box. That
         * means, lower()'s components are moved by -borderWidth
         * and upper()'s by borderWidth. (If borderWidth is
         * negative, the box will get smaller accordingly.)
         */
    void addBorder(VALUETYPE borderWidth)
    {
        lower_ -= borderWidth;
        upper_ += borderWidth;
    }

        /** Adds a border of the given width around the box. That
         * means, lower()'s components are moved by -borderWidth
         * and upper()'s by borderWidth. (If borderWidth is
         * negative, the box will get smaller accordingly.)
         */
    void addBorder(const Vector & borderWidth)
    {
        lower_ -= borderWidth;
        upper_ += borderWidth;
    }

        /// equality check
    bool operator==(Box const &r) const
    {
        return (lower_ == r.lower_) && (upper_ == r.upper_);
    }

        /// inequality check
    bool operator!=(Box const &r) const
    {
        return (lower_ != r.lower_) || (upper_ != r.upper_);
    }

        /** Return whether this box is considered empty. It is
         * non-empty if all upper() coordinates are greater than (or
         * equal, for floating point valuetypes) the corresponding
         * lower() coordinates. Uniting an empty box with something
         * will return the bounding box of the 'something', and
         * intersecting any box with an empty box will again yield an
         * empty box.
         */
    bool empty() const
    {
        return ndim() == 0 || !RangePolicy::nonEmpty(lower_, upper_);
    }

        /** Return whether this box contains the given point.
         * That is, if the point lies within the range [lower, upper] in
         * each dimension (excluding upper() itself for integer valuetypes).
         */
    bool contains(Vector const &p) const
    {
        return ndim() > 0 && allLessEqual(lower_, p) && RangePolicy::nonEmpty(p, upper_);
    }

        /** Return whether this box contains the given
         * one. <tt>r1.contains(r2)</tt> returns the same as
         * <tt>r1 == (r1|r2)</tt> (but is of course more
         * efficient). That also means, a box (even an empty one!)
         * contains() any empty box.
         */
    bool contains(Box const &r) const
    {
        return r.empty() || (contains(r.lower()) && allLessEqual(r.upper(), upper()));
    }

        /** Return whether this box overlaps with the given
         * one. <tt>r1.intersects(r2)</tt> returns the same as
         * <tt>!(r1&r2).empty()</tt> (but is of course much more
         * efficient).
         */
    bool intersects(Box const &r) const
    {
        return !empty() && !r.empty() && (ndim() == r.ndim()) &&
               (RangePolicy::nonEmpty(r.lower(), upper()) ||
                RangePolicy::nonEmpty(lower(), r.upper()));
    }

        /** Modifies this box by including the given point.
         * The result will be the bounding box of the box and the
         * point.  If empty() returns true on the original box, the
         * union will be a box containing only the given point.
         */
    Box &operator|=(Vector const &p)
    {
        if(empty())
        {
            lower_ = p;
            upper_ = RangePolicy::pointEnd(p);
        }
        else
        {
            lower_ = min(lower_, p);
            upper_ = max(upper_, RangePolicy::pointEnd(p));
        }
        return *this;
    }

        /** Returns the union of this box and the given point.
         * The result will be the bounding box of the box and the
         * point.  If empty() returns true on the original box, the
         * union will be a box containing only the given point.
         */
    Box operator|(Vector const &p) const
    {
        Box result(*this);
        result |= p;
        return result;
    }

        /** Modifies this box by uniting it with the given one.
         * The result will be the bounding box of both boxes. If one of
         * the boxes is empty(), the union will be the other one.
         */
    Box &operator|=(Box const &r)
    {
        if(r.empty())
            return *this;
        if(empty())
            return this->operator=(r);

        lower_ = min(lower_, r.lower());
        upper_ = max(upper_, r.upper());
        return *this;
    }

        /** Returns the union of this box and the given one.
         * The result will be the bounding box of both boxs. If one of
         * the boxes empty(), the union will be the other one.
         */
    Box operator|(Box const &r) const
    {
        Box result(*this);
        result |= r;
        return result;
    }

        /** Modifies this box by intersecting it with the given one.
         * The result will be the maximal box contained in both
         * original ones. Intersecting with an empty box will yield
         * again an empty box.
         */
    Box &operator&=(Box const &r)
    {
        if(empty())
            return *this;
        if(r.empty())
            return this->operator=(r);

        lower_ = max(lower_, r.lower());
        upper_ = min(upper_, r.upper());
        return *this;
    }

        /** Intersects this box with the given one.
         * The result will be the maximal box contained in both
         * original ones.  Intersecting with an empty box will yield
         * again an empty box.
         */
    Box operator&(Box const &r) const
    {
        Box result(*this);
        result &= r;
        return result;
    }

        /**
         * Scale box by scalar multiply-assignment.  The same scalar
         * multiply-assignment operation will be performed on both
         * lower() and upper().
         */
    Box &operator*=(double scale)
    {
        lower_ *= scale;
        upper_ *= scale;
        return *this;
    }

        /**
         * Return box scaled by given factor.  The same scalar
         * multiplication will be performed on both lower() and upper().
         */
    Box operator*(double scale)const
    {
        Box result(*this);
        result *= scale;
        return result;
    }

        /**
         * Scale box by scalar divide-assignment.  The same scalar
         * divide-assignment operation will be performed on both
         * lower() and upper().
         */
    Box &operator/=(double scale)
    {
        lower_ /= scale;
        upper_ /= scale;
        return *this;
    }

        /**
         * Return box scaled by inverse of given factor.  The same scalar
         * division will be performed on both lower() and upper().
         */
    Box operator/(double scale)const
    {
        Box result(*this);
        result /= scale;
        return result;
    }

        /**
         * Translate box by vector addition-assignment.  The same vector
         * addition-assignment operation will be performed on both
         * lower() and upper().
         */
    Box &operator+=(const Vector &offset)
    {
        VIGRA_ASSERT_RUNTIME_SIZE(DIMENSION, offset.size() == ndim(),
            "Box::operator+=(): size mismatch between ndim() and offset.");
        lower_ += offset;
        upper_ += offset;
        return *this;
    }

        /**
         * Translate box by vector addition.  The same vector addition
         * operation will be performed on both lower() and upper().
         */
    Box operator+(const Vector &offset)const
    {
        Box result(*this);
        result += offset;
        return result;
    }

        /**
         * Translate box by vector subtract-assignment.  The same vector
         * subtract-assignment operation will be performed on both
         * lower() and upper().
         */
    Box &operator-=(const Vector &offset)
    {
        VIGRA_ASSERT_RUNTIME_SIZE(DIMENSION, offset.size() == ndim(),
            "Box::operator-=(): size mismatch between ndim() and offset.");
        lower_ -= offset;
        upper_ -= offset;
        return *this;
    }

        /**
         * Translate box by vector subtract.  The same vector subtract
         * operation will be performed on both lower() and upper().
         */
    Box operator-(const Vector &offset) const
    {
        Box result(*this);
        result -= offset;
        return result;
    }
};

template<int DIMENSION, class VALUETYPE>
std::ostream& operator<< (std::ostream& stream, const Box<DIMENSION, VALUETYPE> & box) {
    stream<<"["<<box.lower()<<", "<<box.upper()<<" ]";
    return stream;
}

//@}

} // namespace vigra

#endif // VIGRA2_BOX_HXX
