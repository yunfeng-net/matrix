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
// $Id: basic.h,v 1.4 2003/10/03 05:37:35 lyf1998 Exp $
// content: Numerical Analysis Template Object
// Writed by Yunfeng Liu, 2003/04
//
//////////////////////////////////////////////////////////////////////

#if !defined(_BASIC_H_)
#define _BASIC_H_

// limits
#include <limits>
#include <complex>
#include <functional>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <ctime>
#include <complex>
#include <iostream>
#include <cassert>

//////////////////////////////////////////////////////////////////////
// Physical/Electrical Constant
//////////////////////////////////////////////////////////////////////

#define BOLTZMAN	(1.3806226e-23)
#define CHARGE		(1.6021918e-19)
#define EPSILON		(10e-9/36/PI) /* Permitivity */
#define MU			(4*PI*10e-7) /* Permeability */
#define PI			(3.14159265358979323846)
#define REF_TEMP	(300.15) /* 27 degrees C */
#define VC			(3e+8)   /* Light Velocity */
#define ZERO_TEMP	(273.15) /* 0 degrees C */

#define _NAMESPACE_FCA_BEGIN_ namespace fca {
#define _NAMESPACE_FCA_END_   }
#define _STD_   std
#define _FCA_   fca

using namespace std;

_NAMESPACE_FCA_BEGIN_

/// Timer for performance evaluation
class CTimer
{
public:
	CTimer() { restart(); }

	void restart() { tick = clock(); }
	double time() { return (double)(clock()-tick)/CLOCKS_PER_SEC; }
protected:
	clock_t tick;
};

/// Complex with double precision component
typedef _STD_::complex<double> Complex;

// Common operation between real and complex
inline double abs(double X) { return fabs(X); }
inline double norm(double X) { return fabs(X); }
inline double conj(double X) { return X; }
inline double sqr(double X) { return X*X; }

inline double abs(Complex X) { return _STD_::abs(X); }
inline double norm(Complex X) { return fabs(real(X))+fabs(imag(X)); }
inline double sqr(Complex X) { return real(X)*real(X)+imag(X)*imag(X); }

/// Random value generator
class CRandomizer
{
public:
	CRandomizer(double fB = 0, double fE = 1.0, bool bNorm = false)
		{ m_fB = fB; m_fE = fE; m_bNorm = bNorm; }
	virtual ~CRandomizer() { }

	int Get(int& nA)		  const	{ double fA = nA; Get(fA); return (nA = (int)fA);	}
	double Get(double& fA)	  const { return (fA = rand()*(m_fE-m_fB)/5000); }
	Complex& Get(Complex& cA) const	{ double fA; return (cA = Complex(Get(fA), Get(fA))); }

protected:
	double m_fB, m_fE;
	bool m_bNorm;
};

/// Range for the double precision comparison
class CRange
{
public:
	explicit CRange(double fMin = 1e-15, double fRelMin = 1e-15)
		{ m_fMin = fMin, m_fRelMin = fRelMin; }
	CRange(const CRange& A) { *this = A; }
	/// Judge if the value is zero
	template<typename TYPE>	bool IsZero(TYPE c) const { return norm(c)<m_fMin; }
	/// Judge if two value is equal
	template<typename TYPE>
	bool IsEqual(TYPE c1, TYPE c2) const{ return norm(c1-c2)<m_fRelMin*norm(c1+c2); }
	CRange& operator = (const CRange& A) // assignment
		{ m_fMin = A.m_fMin, m_fRelMin = A.m_fRelMin; return *this; }
	double LowerBound() const { return m_fMin; }

protected:
	double m_fMin, m_fRelMin;
};

/// Boundary condition for solver iteration
class CIteRange : public CRange// Addition of the iteration limit
{
public:
	explicit CIteRange(int nMax = -1, double fMin = 1e-15, double fRelMin = 1e-15)
		: CRange(fMin, fRelMin)
			{ m_nMax = (nMax<0)? 0x8fffffff : nMax; }
	CIteRange(const CIteRange& A) { *this = A; }
	~CIteRange() { }

	int count()  const { return m_nMax; }
	CIteRange& operator = (const CIteRange& A) // assignment
	{ 	
		CRange::operator=(A);	m_nMax = A.m_nMax;
		return *this;
	}

protected:
	int m_nMax; // max iteration number
};

/// Rotation for Givens transform
template<typename TYPE>
class CRotation // Iteration limit
{
public:
	CRotation(TYPE a, TYPE b) { set(a, b); }
	/// Compute the c & s of the transform
	void set(TYPE a, TYPE b);
	/// Rotate a & b by the pre-computed c & s 
	void operator()(TYPE &a, TYPE& b) const
	{
		TYPE tmp = m_C*a + m_S*b;
		b = m_C*b - conj(m_S)*a;	a = tmp;
	}
protected:
	double m_C;	TYPE m_S;
};

template<> void CRotation<Complex>::set(Complex a, Complex b)
{
	if(!norm(b)) // simple case
	{
		m_C = 1;	m_S = 0.0;
	}
	else // 
	{
		double fScale = abs(a) + abs(b), fA = abs(a)/fScale, fB = abs(b)/fScale;
		double fNorm = fScale*sqrt(fA*fA + fB*fB);
		m_C = abs(a)/fNorm;	m_S = a/abs(a)*conj(b)/fNorm;
	}
}
template<> void CRotation<double>::set(double a, double b)
{
	if(!b) // simple case
	{
		m_C = 1;	m_S = 0;
	}
	else if(fabs(a)>fabs(b)) 
	{
		double tmp = b/a;	m_C = 1/sqrt(1+tmp*tmp);	m_S = m_C*tmp;
	}
	else // if(fabs(a)<=fabs(b))
	{
		double tmp = a/b;	m_S = 1/sqrt(1+tmp*tmp);	m_C = m_S*tmp;
	}
}

/// Erase the pointers and their objects
template<class vector>
inline void erase_ptr_vector(vector& ptr_vector)
{
	for(typename vector::iterator p = ptr_vector.begin(); p!=ptr_vector.end(); p++)
		delete *p;
	ptr_vector.resize(0);
}

/// Scale the number container with a given value
template<class iterator>
inline void scale(iterator begin, iterator end, const typename iterator::value_type V)
{
	while(begin!=end) *begin++ *= V; // scale the value
}

/// Rotate-shift the value in the vector to right
template<class iterator>
inline void shift_right(iterator begin, iterator end)
{
	for(iterator q = end-1; q!=begin; q--) _STD_::swap(*q, *(q-1));
}

/// Interproduct of two vectors
template<class Vector, class Vector2>
typename Vector::value_type dot(const Vector& A, const Vector2& B) // Inner product
{
	assert((unsigned)A.size()==(unsigned)B.size());
	typename Vector::value_type fSum = 0.0;
	typename Vector::const_iterator p = A.begin();
	typename Vector2::const_iterator q = B.begin();
	for(; p!=A.end(); p++) fSum += *p * conj(*q++);
	return fSum;
}

_NAMESPACE_FCA_END_

#endif // _BASIC_H_
