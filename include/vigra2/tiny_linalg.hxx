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

#ifndef VIGRA_TINY_LINALG_HXX
#define VIGRA_TINY_LINALG_HXX

#include "config.hxx"
#include "numeric_traits.hxx"
#include "tinyarray.hxx"

namespace vigra {

    // vector-matrix product
template <class V1, class D1, class V2, class D2, int N1, int N2>
inline
TinyArray<Promote<V1, V2>, N2>
dot(TinyArrayBase<V1, D1, N1> const & l,
    TinyArrayBase<V2, D2, N1, N2> const & r)
{
    TinyArray<Promote<V1, V2>, N2> res;
    for(int j=0; j < N2; ++j)
    {
        res[j] = l[0] * r(0,j);
        for(int i=1; i < N1; ++i)
            res[j] += l[i] * r(i,j);
    }
    return res;
}

    // vector - symmetric matrix product
template <class V1, class D1, class V2, int N>
inline
TinyArray<Promote<V1, V2>, N>
dot(TinyArrayBase<V1, D1, N> const & l,
    TinySymmetricView<V2, N> const & r)
{
    TinyArray<Promote<V1, V2>, N> res;
    for(int j=0; j < N; ++j)
    {
        res[j] = l[0] * r(0,j);
        for(int i=1; i < N; ++i)
            res[j] += l[i] * r(i,j);
    }
    return res;
}

    // matrix-vector product
template <class V1, class D1, class V2, class D2, int N1, int N2>
inline
TinyArray<Promote<V1, V2>, N1>
dot(TinyArrayBase<V1, D1, N1, N2> const & l,
    TinyArrayBase<V2, D2, N2> const & r)
{
    TinyArray<Promote<V1, V2>, N1> res;
    for(int i=0; i < N1; ++i)
    {
        res[i] = l(i,0) * r[0];
        for(int j=1; j < N2; ++j)
            res[i] += l(i,j) * r[j];
    }
    return res;
}

    // symmetric matrix - vector product
template <class V1, class V2, class D2, int N>
inline
TinyArray<Promote<V1, V2>, N>
dot(TinySymmetricView<V1, N> const & l,
    TinyArrayBase<V2, D2, N> const & r)
{
    TinyArray<Promote<V1, V2>, N> res;
    for(int i=0; i < N; ++i)
    {
        res[i] = l(i,0) * r[0];
        for(int j=1; j < N; ++j)
            res[i] += l(i,j) * r[j];
    }
    return res;
}

    // matrix-matrix product
template <class V1, class D1, class V2, class D2, 
          int N1, int N2, int N3>
inline
TinyArray<Promote<V1, V2>, N1, N3>
dot(TinyArrayBase<V1, D1, N1, N2> const & l,
    TinyArrayBase<V2, D2, N2, N3> const & r)
{
    TinyArray<Promote<V1, V2>, N1, N3> res;
    for(int i=0; i < N1; ++i)
    {
        for(int j=0; j < N3; ++j)
        {
            res(i,j) = l(i,0) * r(0,j);
            for(int k=1; k < N2; ++k)
                res(i,j) += l(i,k) * r(k,j);
        }
    }
    return res;
}

    // matrix - symmetric matrix product
template <class V1, class D1, class V2,
          int N1, int N2>
inline
TinyArray<Promote<V1, V2>, N1, N2>
dot(TinyArrayBase<V1, D1, N1, N2> const & l,
    TinySymmetricView<V2, N2> const & r)
{
    TinyArray<Promote<V1, V2>, N1, N2> res;
    for(int i=0; i < N1; ++i)
    {
        for(int j=0; j < N2; ++j)
        {
            res(i,j) = l(i,0) * r(0,j);
            for(int k=1; k < N2; ++k)
                res(i,j) += l(i,k) * r(k,j);
        }
    }
    return res;
}

    // matrix - symmetric matrix product
template <class V1, class D1, class V2,
          int N1, int N2>
inline
TinyArray<Promote<V1, V2>, N1, N2>
dot(TinySymmetricView<V2, N1> const & l,
    TinyArrayBase<V1, D1, N1, N2> const & r)
{
    TinyArray<Promote<V1, V2>, N1, N2> res;
    for(int i=0; i < N1; ++i)
    {
        for(int j=0; j < N2; ++j)
        {
            res(i,j) = l(i,0) * r(0,j);
            for(int k=1; k < N1; ++k)
                res(i,j) += l(i,k) * r(k,j);
        }
    }
    return res;
}

    // symmetric matrix - symmetric matrix product
template <class V1, class V2, int N>
inline
TinyArray<Promote<V1, V2>, N*(N+1)/2>
dot(TinySymmetricView<V1, N> const & l,
    TinySymmetricView<V2, N> const & r)
{
    TinyArray<Promote<V1, V2>, N*(N+1)/2> res;
    TinySymmetricView<Promote<V1, V2>, N> view(res.data());
    for(int i=0; i < N; ++i)
    {
        for(int j=i; j < N; ++j)
        {
            view(i,j) = l(i,0) * r(0,j);
            for(int k=1; k < N; ++k)
                view(i,j) += l(i,k) * r(k,j);
        }
    }
    return res;
}


} // namespace vigra

#endif // VIGRA_TINY_LINALG_HXX
