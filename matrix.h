/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

//////////////////////////////////////////////////////////////////////
// $Id: matrix.h,v 1.1 2003/09/21 07:54:33 lyf1998 Exp $
// content: interface for the CMatrix class.
// Writed by Yunfeng Liu, 2002/01
//
// The file contains the class template for the matrix
// regardless of the base type
//
//////////////////////////////////////////////////////////////////////

#if !defined(_MATRIX_H_)
#define _MATRIX_H_

#include "vector.h"

_NAMESPACE_FCA_BEGIN_

template<class Matrix, class Vector> struct matrix_vector;	// for M*V algorithm
template<class Matrix, class Vector> struct matrix_iterator;

/** Dense matrix
	that integrates the dimension-sizeing, item-accessing and matrix-solving
*/
template <class TYPE>
class CMatrix
{
public:
	typedef TYPE value_type;
	typedef CVector<TYPE>							vector_type;
	typedef typename CVector<TYPE>::iterator		iterator;
	typedef typename CVector<TYPE>::const_iterator	const_iterator;

	explicit CMatrix(int nRow=0, int nCol=0, CRange range = CRange())
		: m_nRow(0), m_nCol(0), m_bTrans(false), m_Zero(0), m_Range(range)
			{ resize(nRow, nCol); }
	CMatrix(const CMatrix& A) { *this = A; }
	virtual ~CMatrix() {}

	/// @name Item-accessing
	//@{
	/// Position an item to read
	TYPE operator() (int nRow, int nCol) const
	{
		assert(nRow<m_nRow && nCol<m_nCol);
		if(m_bTrans) _STD_::swap(nRow, nCol); 
		return (nRow<0 || nCol<0) ? m_Zero : (*this)[nRow][nCol];
	}
	/// Position an item to write
	TYPE& operator() (int nRow, int nCol)
	{
		assert(nRow<m_nRow && nCol<m_nCol);
		if(m_bTrans) _STD_::swap(nRow, nCol); 
		return (nRow<0 || nCol<0) ? m_Zero : (*this)[nRow][nCol];
	}
	/// Clear all items
	void clear() { m_aItem.fill(0); }
	/// copy the sub-matrix
	template<class Matrix>
	void copy(int left, int right, Matrix& Mtx, int org) const;
	//@}
	/// @name Dimension-sizeing
	//@{
	/// Set the size of the matrix
	void resize(int M, int N)
		{ m_aItem.resize(M*N); m_nRow = M; m_nCol = N; m_bTrans = false; }
	/// Get the size of the matrix
	void size(int &M, int &N) const	{ M = m_nRow; N = m_nCol; if(m_bTrans) _STD_::swap(M,N); }
	/// Get the size of the symmetric matrix
	int  size() const { assert(m_nRow==m_nCol); return m_nRow; }
	//@}
	/// @name matrix-solving
	//@{
	/// Get the condition number after LU factorization
	double condition(const CMapping &aIndex, double fNorm) const;
	/// Factor to Unit Lower * Upper and store the changed order by pivoting
	bool LU(CMapping &aIndex);
	/// Hermite transform
	void H();
	/// Get the norm of the matrix
	double norm(int nNorm = -1) const;
	/// Set the range of items
	void SetRange(const CRange& range) { m_Range = range; }
	/// Solve linear equations by the matrix factored in LU way
	void SolveLU(const CMapping &aIndex, const CVector<TYPE>& b, CVector<TYPE>& x) const;
	/// tlu solver for GMRES : T(L) * x = z
	void SolveTL(const vector_type& b, int N, vector_type& x) const;
	/// Transpose transform
	void T() { m_bTrans = !m_bTrans; }
	//@}

	/// @name matrix-transforming
	//@{
	/// Apply the Givens transform on two columns
	void Givens_col(const CRotation<TYPE> &rot, int i, int j, unsigned N);
	/// Apply the Givens transform on two rows
	void Givens_row(const CRotation<TYPE> &rot, int i, int j, unsigned N);
	//@}
	
	/// @name items-accessing by vector
	//@{
	/// fetch the item for write
	iterator operator [] (int nRow)
		{ return m_aItem.begin()+nRow*m_nCol; }
	/// fetch the item for read
	const_iterator operator [] (int nRow) const
		{ return m_aItem.begin()+nRow*m_nCol; }
	/// fecth the vector to write
	vector_type& data() { return m_aItem; }
	/// fecth the vector to read
	const vector_type& data() const { return m_aItem; }
	//@}

protected:
	template<class Matrix, class Vector>
	friend struct matrix_vector;	// for M*V algorithm
	template<class Matrix, class Vector>
	friend struct matrix_iterator;
	
	CVector<TYPE> m_aItem;
	int m_nRow, m_nCol;		bool m_bTrans;
	TYPE m_Zero;
	CRange m_Range;		// value range
	
	/// @name protected member functions
	//@{
	/// get the pivot, the largest item in a column
	int pivot(int nIdx);
	/// lu solver : A * x = b
	void solve_lu(const CMapping &aIndex, const vector_type& b, vector_type& x) const;
	/// hlu solver : H(L) * x = z
	void solve_hl(const CMapping &aIndex, vector_type& z, vector_type& x) const;
	/// swap two rows
	void swap_row(int nIdx, int nPivot);
	/// swap two columns
	void swap_col(int nIdx, int nPivot);
	//@}
};

// the common class derived from the template
/// the real matrix in double precision
typedef CMatrix<double> CDblMatrix;
/// the complex matrix in double precision
typedef CMatrix<Complex> CCmplMatrix;

//////////////////////////////////////////////////////////////////////
// CMatrix Class Template
//////////////////////////////////////////////////////////////////////

template <class TYPE> template<class Matrix>
void CMatrix<TYPE>::copy(int left, int right, Matrix& Mtx, int org) const
{
	Mtx.clear();
	const_iterator p = (*this)[left];
	for(int i = left; i<right; i++, p+=m_nCol) // for each row
		for(int j = left; j<right; j++) 
				Mtx(org+i-left, org+j-left) = p[j];
}

// partial instantiation for double precision 
void CMatrix<double>::H() { T(); }

// partial instantiation for double precision 
void CMatrix<Complex>::H()
{
	T();
	for(iterator p = m_aItem.begin(); p<m_aItem.end(); p++) *p = conj(*p);
}

// Select the column pivot and swap the correponding rows
template <class TYPE> inline
int CMatrix<TYPE>::pivot(int nIdx)
{
	assert(nIdx<m_nCol);
	iterator p = m_aItem.begin() + nIdx*m_nCol + nIdx;
	
	double Max = _FCA_::norm(*p); // previously set n to pivot
	int i, nPivot = nIdx;
	// find the row with max value at the column
	for(i = nIdx+1, p += m_nCol; i<m_nRow; i++, p += m_nCol)
		if(_FCA_::norm(*p)>Max) { nPivot = i; Max = _FCA_::norm(*p); }
	if(m_Range.IsZero(Max)) return -1;
	if(nIdx!=nPivot) swap_row(nIdx, nPivot); // swap the row and column
	return nPivot;
}

template <class TYPE>
void CMatrix<TYPE>::swap_row(int nIdx, int nPivot)
{
	if(nIdx==nPivot) return;
	iterator p = (*this)[nIdx], q = (*this)[nPivot];
	for(int i = 0; i<m_nCol; i++) _STD_::swap(*p++, *q++);
}

template <class TYPE>
void CMatrix<TYPE>::swap_col(int nIdx, int nPivot)
{
	if(nIdx==nPivot) return;
	iterator p = m_aItem.begin()+nIdx;
	iterator q = m_aItem.begin()+nPivot;
	for(int i = 0; i<m_nRow; i++, p+=m_nCol, q+=m_nCol) _STD_::swap(*p, *q);
}

/*  
	LU decomposition in the Crout way with the pivot choice

	P * A = L * U
	P is the permutation matrix for the pivot choice and is returned by pIndex

	L11 * U11 = a11, L11 = 1, U11 = a11
	L11 * U12 = a12, U12 = a12
	L21 * U11 = a21, L21 = a21/a11
	L21 * U12 + L22 * U22 = a22, L22*U22 = a22 - L21*U12
*/
template <class TYPE>
bool CMatrix<TYPE>::LU(CMapping &aIndex)
{
	int N = size(); // must be a square
	aIndex.resize(N); // Set the dimension
	iterator q = m_aItem.begin();
	for(int i = 0; i<N; i++, q+=N) // for each row
	{
		int j, iPivot = pivot(i);  // select the pivot
		if(iPivot<0) return false;
		aIndex.swap(i, iPivot);
		// L21 // Unit lower matrix
		iterator p = q + N + i;
		for(j = i+1; j<N; j++, p+=N) *p /= q[i];
		
		// Generate L22*U22 = A22 - A21*A12
		// A21 = pA->down->down, ...; A12 = pB->next->next,...
		p = q + N;
		for(j = i+1; j<N; j++, p+=N) // for each row
			for(int k = i+1; k<N; k++) p[k] -= q[k]*p[i];
	}
	return true;
}

template <class TYPE>
void CMatrix<TYPE>::SolveLU(
	const CMapping &aIndex,	const CVector<TYPE> &b, CVector<TYPE> &x) const
{
	if(!m_bTrans) { solve_lu(aIndex, b, x); return; }
	
	// H(A) * x = b => H(U) * H(L) * y = b where x = H(P) * y, A = H(P) * L * U
	vector_type z = b;
	// Forward substitude : H(U) * z = b
	const_iterator p = m_aItem.begin();
	for(int j = 0; j<m_nRow; j++, p+=m_nCol) // for each row
	{
		z[j] /= p[j];
		for(int i = j+1; i<m_nCol; i++)	z[i] -= p[i]*z[j];
	}
	solve_hl(aIndex, z, x);
}

template <class TYPE>
void CMatrix<TYPE>::solve_lu(
	const CMapping &aIndex,	const CVector<TYPE> &b, CVector<TYPE> &x) const
{
	int i, j, M = size(); // Unit lower matrix
	aIndex.translate(b, x); // L * U * x = P * b

	// Forward substitude : L * y = b
	const_iterator p = m_aItem.begin()+M;
	for(j = 0; j<M-1; j++, p+=M+1) // for each column
	{
		const_iterator q = p;
		for(i = j+1; i<M; i++, q+=M) x[i] -= *q*x[j];
	}
	// Backward substitude : U * x = y
	p = m_aItem.end()-1;
	for(j = M-1; j>=0; j--) // for each column
	{
		const_iterator q = m_aItem.begin()+j;
		x[j] /= *p;
		for(i = 0; i<j; i++, q+=M) x[i] -= *q*x[j];
		p -= M+1;  // diagonal item
	}
}

template <class TYPE>
void CMatrix<TYPE>::solve_hl(const CMapping &aIndex, vector_type& z, vector_type& x) const
{
	// Backward substitude : H(L) * y = z
	const_iterator p = m_aItem.end()-m_nCol; // the last row
	for(int j = size()-1; j>0; j--, p-=m_nCol) // for each row
		for(int i = 0; i<j; i++) z[i] -= p[i]*z[j];
	aIndex.inverse(z, x);
}

template <class TYPE>
void CMatrix<TYPE>::SolveTL(const vector_type& b, int N, vector_type& x) const
{
	// Backward substitude : T(L) * x = b
	x = b;
	const_iterator p = m_aItem.begin()+m_nCol*(N-1); // the last row
	for(int j = N-1; j>=0; j--, p-=m_nCol) // for each row
	{
		x[j] /= p[j]; // non-unit lower triagonal
		for(int i = 0; i<j; i++) x[i] -= p[i]*x[j];
	}
}

template <class TYPE>
double CMatrix<TYPE>::norm(int nNorm) const
{
	double fNorm = 0.0;
	const_iterator p = m_aItem.begin();
	if(m_bTrans && abs(nNorm)==1) nNorm = -nNorm;
	if(nNorm==1) // 1st norm
		for(int i = 0; i<m_nCol; i++) // maximum the sum of a column value
		{
			double fSum = 0.0;
			p = m_aItem.begin()+i; // begin with the ith column vector 
			for(int j = 0; j<m_nRow; j++, p+=m_nCol) fSum += _FCA_::norm(*p);
			if(fSum>fNorm) fNorm = fSum;
		}
	else if(nNorm==2)
	{
		while(p!=m_aItem.end()) fNorm += sqr(*p++);
		fNorm = sqrt(fNorm);
	}
	else // infinite norm
		while(p!=m_aItem.end()) // maximum the sum of a row value
		{
			double fSum = 0.0;
			for(int j = 0; j<m_nCol; j++) fSum += _FCA_::norm(*p++);
			if(fSum>fNorm) fNorm = fSum;
		}
	return fNorm;
}

/*
	return the condition number || A || * || inv(A) || for A * x = b
	|| inv(A) || = 1/sigma ~= max(|| y || / || d ||)
	where H(A) * d = b, A * y = d
*/
template <class TYPE>
double CMatrix<TYPE>::condition(const CMapping &aIndex, double fNorm) const
{
	int i, j, N = size();
	vector_type y(N), d(N), b(N), b2(N), b3(N);
	b.fill(0);

	// H(A) * x = b => H(U) * H(L) * y = b
	// where x = H(P) * y, P * A = L * U
	// Forward substitude : H(U) * z = b
	const_iterator q = m_aItem.begin();
	for(j = 0; j<N; j++, q+=m_nCol) // for each row
	{
		TYPE xplus = (1.0- b[j])/q[j];
		TYPE xminus = (-1.0- b[j])/q[j];
		double splus = _FCA_::norm(xplus), sminus = _FCA_::norm(xminus);
		for(i = j+1; i<N; i++)  // compute the influence of choosing x
		{
			b2[i] = b[i] + q[i]*xplus;
			b3[i] = b[i] + q[i]*xminus;
			splus += _FCA_::norm(b2[i]);
			sminus += _FCA_::norm(b3[i]);
		}
		if(splus>sminus)
			{ b[j] = xplus; _STD_::copy(b2.begin()+j+1, b2.end(), b.begin()+j+1); }
		else
			{ b[j] = xminus; _STD_::copy(b3.begin()+j+1, b3.end(), b.begin()+j+1); }
	}
	solve_hl(aIndex, b, d);
	solve_lu(aIndex, d, y);
	
	i = (m_bTrans) ? 1 : -1; // keep same order as fNorm
	return fNorm * y.norm(i) / d.norm(i);
}

template <class TYPE>
void CMatrix<TYPE>::Givens_col(const CRotation<TYPE> &rot, int i, int j, unsigned N)
{
	iterator p = (*this)[i]+i; // begin with the left-upper corner
	iterator q = (*this)[i]+j;
	for(unsigned k = i; k<N; k++, p+=m_nCol, q+=m_nCol) rot(*p, *q);
}

template <class TYPE>
void CMatrix<TYPE>::Givens_row(const CRotation<TYPE> &rot, int i, int j, unsigned N)
{
	const_iterator p = (*this)[i]+i; // begin with the left-upper corner
	const_iterator q = (*this)[j]+i;
	for(unsigned k = i; k<N; k++) rot(*p++, *q++);
}
/**
	Output the matrix to a stream	Syntax: M N { { A[i][j] }j '\n' }i
*/
template <class TYPE>
std::ostream& operator<<(std::ostream &s, const CMatrix<TYPE>&A)
{
    int M, N;
	A.size(M, N);
    s << M << " " << N << "\n";
    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N; j++) s << A[i][j] << " ";
        s << endl;
    }
    return s;
}

/**
	Iutput a matrix from a stream	Syntax: M N { { A[i][j] }j '\n' }i
*/
template <class TYPE>
std::istream& operator>>(std::istream &s, CMatrix<TYPE>&A)
{
    int M, N, AM, AN;
	A.size(AM, AN);
    s >> M >> N;
    if ( M!=AM || N!=AN) A.resize(M,N);
    for (int i=0; i<M; i++)
        for (int j=0; j<N; j++)  s >> A[i][j];
    return s;
}

_NAMESPACE_FCA_END_

#endif // _MATRIX_H_
