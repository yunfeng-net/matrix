/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

//////////////////////////////////////////////////////////////////////
//
// $Id: vector.h,v 1.1 2003/09/21 07:54:34 lyf1998 Exp $
// content: interface for the CVector class.
// Writed by Yunfeng Liu, 
//
// The file contains the class template for the vector
// regardless of the base type
//
//////////////////////////////////////////////////////////////////////

#if !defined(_VECTOR_H_)
#define _VECTOR_H_

#include "basic.h"
#include <vector>

_NAMESPACE_FCA_BEGIN_

/// Vector in linear algebra, which integrates the basic algorithm
template <typename TYPE>
class CVector : public vector<TYPE>
{
public:
	// value type declaratioin
	typedef TYPE value_type;
	typedef value_type* pointer;
	typedef const value_type* const_pointer;
	typedef value_type& reference;
	typedef const value_type& const_reference;
	// iterator declaration
	typedef typename vector<TYPE>::iterator iterator;
	typedef typename vector<TYPE>::const_iterator const_iterator;
	typedef typename vector<TYPE>::reverse_iterator reverse_iterator;
	typedef typename vector<TYPE>::const_reverse_iterator const_reverse_iterator;
	typedef typename vector<TYPE>::size_type size_type;
	typedef typename vector<TYPE>::difference_type difference_type;
	
	explicit CVector(int nLen=0) : vector<TYPE>(nLen) { }
	CVector(const CVector& A) { *this = A; }
	~CVector() {}

	/** @name Self-modification by a vector: =, +=, -=
	Two vectors should have the same dimension with this object
	all the A item will be applied on the corresponding item of this object
	*/
	//@{
	/// Copy from another vector
	template<class Vector> CVector& operator =(const Vector& A)
	{
		typename Vector::const_iterator q = A.begin();
		for(iterator p = begin(); p!=end(); p++) *p = *q++;
		return *this;
	}
	/// Add up another vector 
	template<class Vector> CVector& operator +=(const Vector& A)
	{
		typename Vector::const_iterator q = A.begin();
		for(iterator p = begin(); p!=end(); p++) *p += *q++;
		return *this;
	}
	/// Substract another vector 
	template<class Vector> CVector& operator -=(const Vector& A)
	{
		typename Vector::const_iterator q = A.begin();
		for(iterator p = begin(); p!=end(); p++) *p -= *q++;
		return *this;
	}
	//@}
	/// @name Self-modification by a value: fill, scale
	//@{
	/// Fill in vector with a value
	CVector& fill(const TYPE A) { _STD_::fill(begin(), end(), A); return *this; }
	/// Scale each component with a value
	CVector& scale(const TYPE A) { _FCA_::scale(begin(), end(), A); return *this; }
	//@}

	/// 1, 2 & infinite norm
	double norm(int nNorm = -1) // the norm of vector
	{
		double fNorm = 0;
		if(nNorm<1) 	// the infinite norm : maximum item
        {
		  for(iterator p = begin(); p!=end(); p++)
             if(_FCA_::norm(*p)>fNorm) fNorm = _FCA_::norm(*p);
        }
		else if(nNorm==1)  // the sum of the absolute value of the items
		  for(iterator p = begin(); p!=end(); p++) fNorm += _FCA_::norm(*p);
		else
		{
		  for(iterator p = begin(); p!=end(); p++) fNorm += sqr(*p);
		  fNorm = sqrt(fNorm);
		}
		return fNorm;
	}
};

// the usual vector derived from the template
/// the integer vector
typedef CVector<int> CIntVector;
/// the real vector in double precision
typedef CVector<double> CDblVector;
/// the complex vector in double precision
typedef CVector<Complex> CCmplVector;

/// Map a vector into another one in given order
class CMapping : public CIntVector
{
public:
	CMapping(int nLen = 0) : CIntVector(nLen) { reset(); }
	~CMapping() {}

	/// Combine two mappings
	void combine(const CMapping& A, CMapping &B) const { A.translate(*this, B); }
	/// Inverse-translate a vector into another one
	template<class Vector>
	void inverse(const Vector& A, Vector& B) const // by index
	{
		typename Vector::const_iterator q = A.begin();
		for(const_iterator p = begin(); p!=end(); p++)  B[*p] = *q++;
	}
	/// Translate a vector into another one
	template<class Vector>
	void translate(const Vector& A, Vector& B) const // by index
	{
		typename Vector::iterator q = B.begin();
		for(const_iterator p = begin(); p!=end(); p++)  *q++ = A[*p];
	}
	/// Restore the identity mapping
	void reset()
	{
		int i = 0;
		for(iterator p = begin(); p!=end(); *p++ = i++);
	}
	/// Swap two items
	void swap(int i, int j) { if(i!=j) _STD_::swap((*this)[i], (*this)[j]); }
	/// Change the size
	void resize(int size) { CIntVector::resize(size); reset(); }	
};

/**
	Output vector to stream
	Format: N '\n' { ' ' A[i] } '\n'
*/
template <class TYPE>
_STD_::ostream& operator<<(_STD_::ostream &s, const CVector<TYPE> &A)
{
    int N = A.size();
    s <<  N << endl;
    for(typename CVector<TYPE>::const_iterator p = A.begin(); p!=A.end(); p++) s << " " << *p;
    s << endl;
    return s;
}

/**
	Input vector from stream
	Format: N '\n' { ' ' A[i] } '\n'
 */
template <class TYPE>
_STD_::istream & operator>>(_STD_::istream &s, CVector<TYPE> &A)
{
    int N;
    s >> N;
    if(N!=(int)A.size()) A.resize(N);
    for(int i=0; i<N; i++)  s >> A[i];
    return s;
}

_NAMESPACE_FCA_END_

#endif // _VECTOR_H_
