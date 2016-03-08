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
// $Id: algorithm.h,v 1.2 2003/10/02 08:20:30 lyf1998 Exp $
// content: Interface for the CVector class.
// Writed by Yunfeng Liu, 2003/9
//
// The file contains the algorithm template for the vector
// regardless of the base type
//
//////////////////////////////////////////////////////////////////////

#if !defined(_ALGORITHM_H_)
#define _ALGORITHM_H_

#include "matrix.h"
#include "spmatrix.h"
#include <iterator>

_NAMESPACE_FCA_BEGIN_

/// Delayed iterator for a unary function in the delayed vector
template<class Iterator, class Functor>
struct unary_iterator : public iterator<random_access_iterator_tag, typename Functor::result_type>
{
	typedef iterator<random_access_iterator_tag, typename Functor::result_type> base_type;
	typedef typename base_type::value_type              value_type;
	
	unary_iterator(Iterator p, Functor func) : m_pointer(p), m_func(func) { }

	value_type operator *() const { return m_func(*m_pointer); } // operation
	unary_iterator& operator++() { m_pointer++; return *this; }
	unary_iterator operator++(int) { return unary_iterator(m_pointer++, m_func); }
	unary_iterator operator+(int N) const { return unary_iterator(m_pointer+N, m_func); }
	bool operator !=(unary_iterator& A) const { return m_pointer!=A.m_pointer; }
protected:
	Iterator m_pointer;
	Functor m_func;
};

/// Delayed vector for a unary function
template<class Vector, class Functor>
struct unary_vector
{
	typedef typename Functor::result_type value_type;
	typedef unary_iterator<typename Vector::const_iterator, Functor> const_iterator;

	unary_vector(const Vector &A, Functor func) : m_vA(A), m_func(func) { }
	
	// iterator production
	const_iterator begin() const { return const_iterator(m_vA.begin(), m_func); }
	const_iterator end() const { return const_iterator(m_vA.end(), m_func); }
	int size() const { return m_vA.size(); }
protected:
	const Vector &m_vA;
	Functor m_func;
};

/// Delayed iterator for a binary function in the delayed vector
template<class Iterator, class Iterator2, class Functor>
struct binary_iterator : public iterator<random_access_iterator_tag, typename Functor::result_type>
{
	typedef iterator<random_access_iterator_tag, typename Functor::result_type> base_type;
	typedef typename base_type::value_type              value_type;
	
	binary_iterator(Iterator p, Iterator2 p2, Functor func)// constructor
		: m_pointer(p), m_pointer2(p2), m_func(func) { }

	value_type operator *() const { return m_func(*m_pointer, *m_pointer2); } // operation
	binary_iterator& operator++() { m_pointer++; m_pointer2++; return *this; }
	binary_iterator operator++(int)
		{ return binary_iterator(m_pointer++, m_pointer2++, m_func); }
	binary_iterator operator+(int N) const
		{ return binary_iterator(m_pointer+N, m_pointer2+N, m_func); }
	bool operator !=(binary_iterator& A) const
		{ return pointer1!=A.m_pointer && m_pointer2!=A.m_pointer2; }
protected:
	Iterator m_pointer;	Iterator2 m_pointer2;
	Functor m_func;
};

/// Delayed vector for a binary function
template<class Vector, class Vector2, class Functor>
struct binary_vector
{
	typedef typename Functor::result_type value_type;
	typedef binary_iterator<typename Vector::const_iterator,
		typename Vector2::const_iterator, Functor> const_iterator;

	binary_vector(const Vector &A, const Vector2 &B, Functor func)  // constructor
		: m_vA(A), m_vB(B), m_func(func) { }

	// iterator production
	const_iterator begin() const { return const_iterator(m_vA.begin(), m_vB.begin(), m_func); }
	const_iterator end() const { return const_iterator(m_vA.end(), m_vB.end(), m_func); }
	int size() const { return m_vA.size(); }

protected:
	const Vector &m_vA;	const Vector2 &m_vB;
	Functor m_func;
};
/// Delayed iterator for a matrix times a vector
template<class Matrix, class Iterator>
struct matrix_iterator : public iterator<random_access_iterator_tag,
	typename Matrix::value_type>
{
	typedef typename Matrix::value_type 					value_type;
	typedef typename Matrix::vector_type::const_iterator	const_iterator;
	
	matrix_iterator(const Matrix& Mtx, const_iterator p, Iterator p2)// constructor
		: m_Mtx(Mtx), m_pointer(p), m_pointer2(p2), m_bTrans(Mtx.m_bTrans)	{ }

	value_type operator *() const  // operation
	{
		value_type sum(0);
		int i;	const_iterator p = m_pointer;
		if(m_bTrans)
			for(i = 0; i<m_Mtx.m_nRow; i++, p+=m_Mtx.m_nCol) sum += *p*(*(m_pointer2+i));
		else
			for(i = 0; i<m_Mtx.m_nCol; i++, p++) sum += *p*(*(m_pointer2+i));
		return sum; 
	}
	matrix_iterator& operator++()
		{ if(m_bTrans) m_pointer++; else m_pointer+=m_Mtx.m_nCol; return *this; }
	matrix_iterator operator++(int)
	{ 
		const_iterator p = m_pointer;
		if(m_bTrans) m_pointer++; else m_pointer+=m_Mtx.m_nCol;
		return matrix_iterator(m_Mtx, p, m_pointer2);
	}
	matrix_iterator operator+(int N) const
	{ 
		return matrix_iterator(m_Mtx,
			(m_bTrans) ? m_pointer+N : m_pointer+m_Mtx.m_nCol*N, m_pointer2);
	}
	bool operator !=(matrix_iterator& A) const { return m_pointer!=A.m_pointer; }

protected:
	const Matrix& m_Mtx;	const_iterator m_pointer;	const Iterator m_pointer2;
	bool m_bTrans;
};

/// Delayed vector of a matrix times a vector
template<class Matrix, class Vector>
struct matrix_vector
{
	typedef typename Matrix::value_type value_type;
	typedef matrix_iterator<Matrix, typename Vector::const_iterator> const_iterator;

	matrix_vector(const Matrix &Mtx, const Vector &V) : m_Mtx(Mtx), m_vA(V)
	{ 
		int M, N; m_Mtx.size(M, N);
		assert(!m_Mtx.m_bTrans && N==(int)m_vA.size() || m_Mtx.m_bTrans && M==(int)m_vA.size());
	}

	// iterator production
	const_iterator begin() const
		{ return const_iterator(m_Mtx, m_Mtx.data().begin(), m_vA.begin()); }
	const_iterator end() const
	{
		return const_iterator(m_Mtx, (!m_Mtx.m_bTrans)?m_Mtx.data().end():
			m_Mtx.data().begin()+m_Mtx.m_nCol, m_vA.begin());
	}
	int size() const { return m_vA.size(); }

protected:
	const Matrix &m_Mtx;	const Vector &m_vA;
};

/// Delayed iterator for a sparse matrix times a vector
template<class ItemPtrIterator, class Iterator>
struct spmatrix_iterator : public iterator<random_access_iterator_tag,
	typename Iterator::value_type>
{
	typedef typename Iterator::value_type value_type;
	typedef typename ItemPtrIterator::value_type ptr_type;
	
	spmatrix_iterator(ItemPtrIterator p, Iterator p2, bool trans)// constructor
		: m_pointer(p), m_pointer2(p2), m_bTrans(trans)	{ }

	value_type operator *() const  // operation
	{
		value_type sum(0);
		if(m_bTrans)
			for(ptr_type p = *m_pointer; p; p = p->down)
				sum += *(m_pointer2+p->m_nRow) * p->value;
		else
			for(ptr_type p = *m_pointer; p; p = p->next)
				sum += *(m_pointer2+p->m_nCol) * p->value;
		return sum; 
	}
	spmatrix_iterator& operator++() { m_pointer++; return *this; }
	spmatrix_iterator operator++(int)
		{ return spmatrix_iterator(m_pointer++, m_pointer2, m_bTrans); }
	spmatrix_iterator operator+(int N) const
		{ return spmatrix_iterator(m_pointer+N, m_pointer2, m_bTrans); }
	bool operator !=(spmatrix_iterator& A) const { return m_pointer!=A.m_pointer; }

protected:
	ItemPtrIterator m_pointer;	const Iterator m_pointer2;
	bool m_bTrans;
};

/// Delayed vector of a sparse matrix times a vector
template<class Matrix, class Vector>
struct spmatrix_vector
{
	typedef typename Matrix::value_type value_type;
	typedef spmatrix_iterator<typename Matrix::CItemVector::const_iterator,
		typename Vector::const_iterator> const_iterator;

	spmatrix_vector(const Matrix &Mtx, const Vector &V) : m_Mtx(Mtx), m_vA(V)
	{ 
		int M, N; m_Mtx.size(M, N);
		assert(!m_Mtx.m_bTrans && N==(int)m_vA.size() || m_Mtx.m_bTrans && M==(int)m_vA.size() );
	}

	// iterator production
	const_iterator begin() const
	{
		return const_iterator((!m_Mtx.m_bTrans)?m_Mtx.m_aRow.begin():m_Mtx.m_aCol.begin(),
			m_vA.begin(), m_Mtx.m_bTrans);
	}
	const_iterator end() const
	{
		return const_iterator((!m_Mtx.m_bTrans)?m_Mtx.m_aRow.end():m_Mtx.m_aCol.end(),
			m_vA.begin(), m_Mtx.m_bTrans);
	}
	int size() const { return m_vA.size(); }

protected:
	const Matrix &m_Mtx;	const Vector &m_vA;
};

/// Scale a vector for the delayed computation
template<class Vector>
inline unary_vector<Vector, // times
	_STD_::binder2nd< _STD_::multiplies< typename Vector::value_type> > >
scaled(const Vector &A, typename Vector::value_type alpha)
{
	typedef _STD_::multiplies<typename Vector::value_type> multiplies;
	typedef _STD_::binder2nd<multiplies> times;
	return unary_vector<Vector, times>(A, times(multiplies(), alpha));
}

/// Add vectors for the delayed computation
template<class Vector, class Vector2>
inline binary_vector<Vector, Vector2, _STD_::plus<typename Vector::value_type> >
add(const Vector &A, const Vector2& B)
{ 
	typedef _STD_::plus<typename Vector::value_type> plus;
	return binary_vector<Vector, Vector2, plus>(A, B, plus());
}

/// Substract vectors for the delayed computation
template<class Vector, class Vector2>
inline binary_vector<Vector, Vector2, _STD_::minus<typename Vector::value_type> >
sub(const Vector &A, const Vector2 &B)
{ 
	typedef typename Vector::value_type value_type;
	typedef typename _STD_::minus<value_type> minus;
	return binary_vector<Vector, Vector2, minus>(A, B, minus());
}

/// Multiply a matrix with vectors for the delayed computation
template<typename TYPE, class Vector>
inline matrix_vector<CMatrix<TYPE>, Vector>
mult(const CMatrix<TYPE> &A, const Vector &B)
{ 
	return matrix_vector<CMatrix<TYPE>, Vector>(A, B);
}

/// Multiply a sparse matrix with vectors for the delayed computation
template<typename TYPE, class Vector>
inline spmatrix_vector<CSpMatrix<TYPE>, Vector>
mult(const CSpMatrix<TYPE> &A, const Vector &B)
{ 
	return spmatrix_vector<CSpMatrix<TYPE>, Vector>(A, B);
}

_NAMESPACE_FCA_END_

#endif // _ALGORITHM_H_
