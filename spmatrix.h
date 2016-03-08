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
// $Id: spmatrix.h,v 1.1 2003/09/21 07:54:34 lyf1998 Exp $
// content: interface and implementation of the CSpMatrix class.
// Writed by Yunfeng Liu, 2002/05
//
// The file contains the class template for the sparse matrix
// regardless of the base type
//
//////////////////////////////////////////////////////////////////////

#ifndef _SPARSE_MATRIX_H_
#define _SPARSE_MATRIX_H_

#include "vector.h"

_NAMESPACE_FCA_BEGIN_

template<class Matrix, class Vector> struct spmatrix_vector;
template<class ItemPtrIterator, class Iterator> struct spmatrix_iterator;

/** Sparse matrix
	that integrates the dimension-sizeing, item-accessing and matrix-solving
*/
template <typename TYPE>
class CSpMatrix
{
public:
	typedef TYPE	 				value_type;
	typedef CVector<TYPE>	 		vector_type;
	
	explicit CSpMatrix(int nRow = 0, int nCol = 0, CRange range = CRange())
		: m_nRow(0), m_nCol(0), m_bTrans(false), m_Zero(0), m_Range(range)
			{ resize(nRow, nCol); }
	~CSpMatrix() { destroy(); }

	/// @name Item-accessing
	//@{
	/// Position an item to read
	TYPE operator() (int nRow, int nCol) const;
	/// Position an item to write
	TYPE& operator() (int nRow, int nCol) { return set(nRow, nCol); }
	/// Clear all items
	void clear();
	/// copy the sub-matrix
	template<class Matrix>
	void copy(int left, int right, Matrix& Mtx, int org) const;
	/// Position an item to write
	TYPE& set(int nRow, int nCol);
	//@}
	/// @name Dimension-sizeing
	//@{
	/// Set the size of the matrix
	void resize(int M, int N);
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
	/// Get the norm of the matrix
	double norm(int nNorm = -1) const;
	/// Hermite transform
	void H();
	/// Set the range of items
	void SetRange(const CRange& range) { m_Range = range; }
	/// Solve linear equations by the matrix factored in LU way
	void SolveLU(const CMapping &aIndex, const vector_type& b, vector_type& x) const;
	/// Transpose transform
	void T() { m_bTrans = !m_bTrans; }
	//@}

	/// @name i/o function
	//@{
	/// Print to an output stream 
	void print(ostream &output) const
	{
		output << m_nRow << ' ' << m_nCol << endl;
		for(int i = 0; i<m_nRow; i++, output << endl)
			for(CItem *pT = m_aRow[i]; pT; pT = pT->next)
				if(!m_bTrans)
					output << " " << pT->m_nRow << " " << pT->m_nCol << " " << pT->value;
				else
					output << " " << pT->m_nCol << " " << pT->m_nRow << " " << pT->value;
		output << m_nRow << endl; // end flag
	}
	/// Scan from an output stream 
	void scan(istream &input)
	{
		int M, N;
		input >> M >> N;
		resize(M, N);
		while(!input.eof())
		{
			int nRow, nCol;	TYPE value;
			input >> nRow;
			if(nRow>=m_nRow) return; // end of spmatrix
			input >> nCol >> value;
			(*this)(nRow, nCol) = value;
		}
	}
	//@}
protected:
	template<class Matrix, class Vector>
	friend struct spmatrix_vector;	// for M*V algorithm
	template<class ItemPtrIterator, class Iterator>
	friend struct spmatrix_iterator;
	
	struct CItem
	{
		TYPE value;
		CItem(int nRow, int nCol) { m_nRow = nRow; m_nCol = nCol; }
		int m_nRow, m_nCol;
		CItem *next, *down;
	};
	typedef vector<CItem*> CItemVector;

	int m_nRow, m_nCol;		bool m_bTrans;
	CItemVector m_aRow, m_aCol, m_aDiag; // index for items
	TYPE m_Zero;		// for null items
	CRange m_Range;		// value range
	
	/// @name protected member functions
	//@{
	/// free all items
	void destroy();
	/// link the diagonal item to the index
	void link_diag(int nIdx)
	{
		CItem *pT;
		for(pT = m_aCol[nIdx]; pT && pT->m_nRow<nIdx; pT = pT->down); // for each row
		m_aDiag[nIdx] = (pT && pT->m_nRow==nIdx) ? pT : NULL; 
	}
	/// get the pivot, the largest item in a column
	int pivot(int nIdx);
	/// lu solver : A * x = b
	void solve_lu(const CMapping &aIndex, const vector_type& b, vector_type& x) const;
	/// hlu solver : H(L) * x = z
	void solve_hl(const CMapping &aIndex, vector_type& z, vector_type& x) const;
	/// swap two rows
	void swap_row(int nIdx, int nPivot);
	/// insert an item in a row
	inline void insert_item_in_row(CItem *pA, int nIdx);
	/// move an item in a column
	inline void move_item_in_column(CItem *pA, int nIdx, int nPivot);
	/// move a row to another row's tail
	inline void move_row(CItem **ppA, CItem **ppB, int nIdx, int nPivot);
	/// integrity assertion
	bool integrity() const
	{
		assert((int)m_aDiag.size()==_STD_::min(m_nRow, m_nCol));
		unsigned i;	CItem *pT;
		for(i = 0; i<m_aRow.size(); i++)
			for(pT = m_aRow[i]; pT; pT = pT->next) assert(pT->m_nRow==(int)i); // for each row
		for(i = 0; i<m_aCol.size(); i++)
			for(pT = m_aCol[i]; pT; pT = pT->down) assert(pT->m_nCol==(int)i); // for each row
		for(i = 0; i<m_aDiag.size(); i++) // for each diagonal item
		{ 
			pT = m_aDiag[i]; assert(!pT || pT->m_nCol==(int)i && pT->m_nRow==(int)i); 
		}
		return m_nRow==(int)m_aRow.size() && m_nCol==(int)m_aCol.size();
	}
	//@}
};

// the usual sparse matrix derived from the template
/// the real sparse matrix in double precision
typedef CSpMatrix<double> CDblSpMatrix;
/// the complex sparse matrix in double precision
typedef CSpMatrix<Complex> CCmplSpMatrix;

//////////////////////////////////////////////////////////////////////
// CSpMatrix Class Template
//////////////////////////////////////////////////////////////////////

template <class TYPE>
void CSpMatrix<TYPE>::resize(int nRow, int nCol)
{
	assert(nRow>=0 && nCol>=0);
	if(nRow==m_nRow && m_nCol==nCol) return;
	destroy(); // delete the previous

	m_aRow.resize(nRow);	m_aCol.resize(nCol); // Length
	m_aDiag.resize(min(nRow,nCol));

	int i;
	for(unsigned j = 0; j<m_aDiag.size(); j++) m_aDiag[j] = NULL;
	for(i = 0; i<nRow; i++) m_aRow[i] = NULL;
	for(i = 0; i<nCol; i++) m_aCol[i] = NULL; // Set NULL

	m_nRow = nRow; m_nCol = nCol; // member initiating
	m_Zero = (TYPE)0.0;
	m_bTrans = false;
	assert(integrity());
}

template <class TYPE>
void CSpMatrix<TYPE>::destroy()
{
	for(unsigned i = 0; i<m_aRow.size(); i++) // delete all rows
		for(CItem *pB, *pA = m_aRow[i]; pA; pA = pB)	// delete all items
	{
		pB = pA->next;	delete pA;
	}
	m_aRow.resize(0);	m_aCol.resize(0);	m_aDiag.resize(0);
}

template <class TYPE>
TYPE CSpMatrix<TYPE>::operator () (int nRow, int nCol) const
{
	CItem *pA;
	if(m_bTrans) _STD_::swap(nRow, nCol);
	for(pA = m_aCol[nCol]; pA && pA->m_nRow<nRow; pA = pA->down);
;
	if(pA && pA->m_nRow==nRow) return pA->value; else return (TYPE)0.0;
}

template <class TYPE>
void CSpMatrix<TYPE>::clear()
{
	for(unsigned i = 0; i<m_aRow.size(); i++) // for all rows
		for(CItem *pA = m_aRow[i]; pA; pA = pA->next)	// for all items
			pA->value = (TYPE)0.0;
}

template <class TYPE> template<class Matrix>
void CSpMatrix<TYPE>::copy(int left, int right, Matrix& Mtx, int org) const
{
	Mtx.clear();
	for(int i = left; i<right; i++) //
		for(CItem *pA = m_aRow[i]; pA && pA->m_nCol<right; pA = pA->next)
			if(pA->m_nCol>=left)
				Mtx(org+pA->m_nRow-left, org+pA->m_nCol-left) = pA->value;
}

template <class TYPE>
TYPE& CSpMatrix<TYPE>::set(int nRow, int nCol)
{
	assert(nRow<m_nRow && nCol<m_nCol);
	if(nRow<0 || nCol<0) return m_Zero;
	if(nRow==nCol)
		{ if(m_aDiag[nRow]) return m_aDiag[nRow]->value; } // diagonoal item
	else
		if(m_bTrans) _STD_::swap(nRow, nCol); // hermite-transformed

	CItem **ppA = &m_aRow[nRow];	// find if it exists
	while(ppA[0] && ppA[0]->m_nCol<nCol) ppA = &(ppA[0]->next);
	if(ppA[0] && ppA[0]->m_nCol==nCol) return ppA[0]->value;

	// new an item and insert into the row, column and diagonal indices
	CItem *pB = new CItem(nRow,nCol);
	pB->next = ppA[0]; ppA[0] = pB;
	
	// insert into the column
	ppA = &m_aCol[nCol];
	while(ppA[0] && ppA[0]->m_nRow<nRow) ppA = &(ppA[0]->down);
	pB->down = ppA[0];	ppA[0] = pB;

	// insert into the diagonal line
	if(nRow==nCol) m_aDiag[nRow] = pB;

	assert(integrity());
	return pB->value;
}

// partial instantiation for double precision
template<> void CSpMatrix<double>::H() { T(); }

// partial instantiation for double precision 
template<> void CSpMatrix<Complex>::H()
{
	T();
	for(unsigned i = 0; i<m_aRow.size(); i++) // for all rows
		for(CItem *pA = m_aRow[i]; pA; pA = pA->next)	// for all items
			pA->value = conj(pA->value);
}

template <class TYPE>
inline void CSpMatrix<TYPE>::insert_item_in_row(CItem *pA, int nIdx)
{
	assert(pA->m_nRow==nIdx);
	CItem **ppA = &m_aRow[nIdx];
	for(; ppA[0] && ppA[0]->m_nCol<pA->m_nCol; ppA = &(ppA[0]->next));
	pA->next = ppA[0]; ppA[0] = pA;
}

template <class TYPE>
inline void CSpMatrix<TYPE>::move_item_in_column(CItem *pA, int nIdx, int nPivot)
{
	pA->m_nRow = nPivot;
	if(nIdx<nPivot && (!pA->down || pA->down->m_nRow<nPivot)) return;
	CItem **ppA = &m_aCol[pA->m_nCol];
	for(; ppA[0]->m_nRow<_STD_::min(nIdx, nPivot); ppA = &(ppA[0]->down)); // pointer to pA
	if(nIdx<nPivot)
	{
		ppA[0] = pA->down; // delete
		for(; ppA[0] && ppA[0]->m_nRow<nPivot; ppA = &(ppA[0]->down));
		pA->down = ppA[0];	ppA[0] = pA; // insert;
	}
	else // nIdx>nPivot
	{
		CItem *pB = pA->down;
		pA->down = ppA[0];	ppA[0] = pA; // insert
		for(ppA = &(ppA[0]->down); ppA[0]!=pA; ppA = &(ppA[0]->down));
		ppA[0] = pB; // link the previous pA's tail
	}
}

template <class TYPE>
inline void CSpMatrix<TYPE>::move_row(CItem **ppA, CItem **ppB, int nIdx, int nPivot)
{
	ppB[0] = ppA[0]; ppA[0] = NULL; // move the ppA up to ppB
	for( ; ppB[0]; ppB = &(ppB[0]->next)) // for each item in the row
		move_item_in_column(ppB[0], nIdx, nPivot);
	assert(integrity());
}

template <class TYPE>
void CSpMatrix<TYPE>::swap_row(int nIdx, int nPivot)
{
	if(nIdx==nPivot) return; else if(nIdx>nPivot) _STD_::swap(nIdx, nPivot);
	CItem **ppA = &m_aRow[nIdx], **ppB = &m_aRow[nPivot]; 

	while(ppA[0] || ppB[0])
	{
		if(!ppA[0])
		{
			 move_row(ppB, ppA, nPivot, nIdx); // exchange items in each column of ppB
			 break;
		}
		else if(!ppB[0])
		{
			 move_row(ppA, ppB, nIdx, nPivot); // exchange items in each column of ppA
			 break;
		}
		else if(ppA[0]->m_nCol==ppB[0]->m_nCol)
		{
			_STD_::swap(ppA[0]->value, ppB[0]->value);
			ppA = &(ppA[0]->next);	ppB = &(ppB[0]->next);
		}
		else if(ppA[0]->m_nCol<ppB[0]->m_nCol)
		{
			CItem *pA = ppA[0];
			ppA[0] = ppA[0]->next;	// delete from row nIdx
			move_item_in_column(pA, nIdx, nPivot);
			insert_item_in_row(pA, nPivot); // insert into row nPivot
			assert(ppB[0]==pA);	ppB = &(ppB[0]->next);
		}
		else // ppA[0]->m_nCol>ppB[0]->m_nCol
		{
			CItem *pB = ppB[0];
			ppB[0] = ppB[0]->next;	// delete from row nPivot
			move_item_in_column(pB, nPivot, nIdx);
			insert_item_in_row(pB, nIdx); // insert into row nIdx
			assert(ppA[0]==pB);	ppA = &(ppA[0]->next);
		}
	}

	link_diag(nIdx); link_diag(nPivot); // link the diagonal item to the index
	assert(integrity());
}


template <class TYPE>
int CSpMatrix<TYPE>::pivot(int nIdx)
{
	double Max = _FCA_::norm((!m_aDiag[nIdx]) ? 0.0 : m_aDiag[nIdx]->value);
	const CItem *pB, *pA = m_aDiag[nIdx];
	if(pA) { if(!pA->down) return nIdx; }
	else
	{
 		for(pA = m_aCol[nIdx]; pA && pA->m_nRow<nIdx; pA = pA->down);
		if(!pA) return -1; // error for null column
		Max = _FCA_::norm(pA->value);
	}
	
	for(pB = pA->down; pB; pB = pB->down)	// find the row with max value at the column
		if(_FCA_::norm(pB->value)>Max) { pA = pB; Max = _FCA_::norm(pB->value); }
	if(m_Range.IsZero(Max)) return -1;
	if(nIdx!=pA->m_nRow) swap_row(nIdx, pA->m_nRow); // swap the row and column
	return pA->m_nRow;
}

template <class TYPE>
bool CSpMatrix<TYPE>::LU(CMapping &aIndex)
{
	int N = size(); // must be a square
	aIndex.resize(N); // Set the dimension
	for(int i = 0; i<N; i++)
	{
		int iPivot = pivot(i);  // select the pivot
		if(iPivot<0) return false;
		aIndex.swap(i, iPivot);
		CItem *pA = m_aDiag[i]->down; // Select the column
		// L21 // Unit lower matrix
		for(; pA; pA = pA->down) pA->value /= m_aDiag[i]->value;
		
		// Generate L22*U22 = A22 - A21*A12
		// A21 = pA->down->down, ...; A12 = pB->next->next,...
		for(pA = m_aDiag[i]->down; pA; pA = pA->down)
			for(CItem *pC = pA, *pB = m_aDiag[i]->next; pB; pB = pB->next)
			{
				while(pC && pC->m_nCol<pB->m_nCol) pC = pC->next;
				TYPE &Value = (pC && pC->m_nCol==pB->m_nCol) ? // if it exists
					pC->value : set(pA->m_nRow, pB->m_nCol);
				Value -= pA->value * pB->value;
			}
	}
	return true;
}

template <class TYPE>
void CSpMatrix<TYPE>::SolveLU(
	const CMapping &aIndex,	const CVector<TYPE> &b, CVector<TYPE> &x) const
{
	if(!m_bTrans) { solve_lu(aIndex, b, x); return; }
	
	CVector<TYPE> z = b;
	// H(A) * x = b => H(U) * H(L) * y = b where x = H(P) * y, A = H(P) * L * U
	// Forward substitude : H(U) * z = b
	for(int j = 0; j<size(); j++)
	{
		z[j] /= m_aDiag[j]->value;
		for(CItem *pA = m_aDiag[j]->next; pA; pA = pA->next)
			z[pA->m_nCol] -= pA->value*z[j];
	}
	solve_hl(aIndex, z, x);
}

template <class TYPE>
void CSpMatrix<TYPE>::solve_lu(
	const CMapping &aIndex,	const CVector<TYPE> &b, CVector<TYPE> &x) const
{
	int j, M = size(); // Unit lower matrix
	aIndex.translate(b, x); // L * U * x = P * b
	// Forward substitude : L * y = b
	for(j = 0; j<M-1; j++)
		for(CItem *pA = m_aDiag[j]->down; pA; pA = pA->down)
			x[pA->m_nRow] -= pA->value*x[j];
	// Backward substitude : U * x = y
	for(j = M-1; j>=0; j--)
	{
		x[j] /= m_aDiag[j]->value;
		for(CItem *pA = m_aCol[j]; pA && pA->m_nRow<j; pA = pA->down)
			x[pA->m_nRow] -= pA->value*x[j];
	}
}

template <class TYPE>
void CSpMatrix<TYPE>::solve_hl(
	const CMapping &aIndex,	CVector<TYPE> &z, CVector<TYPE> &x) const
{
	// Backward substitude : H(L) * y = z
	for(int j = size()-1; j>0; j--)
		for(CItem *pA = m_aRow[j]; pA && pA->m_nCol<j; pA = pA->next)
			z[pA->m_nCol] -= pA->value*z[j];
	aIndex.inverse(z, x);
}

template <class TYPE>
double CSpMatrix<TYPE>::norm(int nNorm) const
{
	double fNorm = 0.0;
	int i;	CItem *pA;
	if(m_bTrans && abs(nNorm)==1) nNorm = -nNorm;
	if(nNorm==1)// 1st norm
		for(i = 0; i<m_nCol; i++) // maximum the sum of a column value
		{
			double fSum = 0.0;
			for(pA = m_aCol[i]; pA; pA = pA->down) fSum += _FCA_::norm(pA->value);
			if(fSum>fNorm) fNorm = fSum;
		}
	else if(nNorm==2)
	{
		for(i = 0; i<m_nCol; i++) // maximum the sum of a column value
			for(pA = m_aCol[i]; pA; pA = pA->down) fNorm += sqr(pA->value);
		fNorm = sqrt(fNorm);
	}
	else // infinite norm
		for(i = 0; i<m_nRow; i++) // maximum the sum of a row value
		{
			double fSum = 0.0;
			for(pA = m_aRow[i]; pA; pA = pA->next) fSum += _FCA_::norm(pA->value);
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
double CSpMatrix<TYPE>::condition(const CMapping &aIndex, double fNorm) const
{
	int i, j, N = size();
	CVector<TYPE> y(N), d(N), b(N), b2(N), b3(N);
	b.fill(0);

	// H(A) * x = b => H(U) * H(L) * y = b
	// where x = H(P) * y, P * A = L * U
	// Forward substitude : H(U) * z = b
	for(j = 0; j<N; j++)
	{
		TYPE xplus = (1.0- b[j])/m_aDiag[j]->value;
		TYPE xminus = (-1.0- b[j])/m_aDiag[j]->value;
		double splus = _FCA_::norm(xplus), sminus = _FCA_::norm(xminus);
		for(CItem *pA = m_aDiag[j]->next; pA; pA = pA->next)
		{
			b2[pA->m_nCol] = b[pA->m_nCol] + pA->value*xplus;
			b3[pA->m_nCol] = b[pA->m_nCol] + pA->value*xminus;
			splus += _FCA_::norm(b2[pA->m_nCol]);
			sminus += _FCA_::norm(b3[pA->m_nCol]);
		} // compute the influence of choosing x
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

/// Output spmatrix to stream
template <class TYPE>
_STD_::ostream& operator<<(_STD_::ostream &s, const CSpMatrix<TYPE> &A)
{
	A.print(s);		return s;
}

///	Input spmatrix from stream
template <class TYPE>
_STD_::istream & operator>>(_STD_::istream &s, CSpMatrix<TYPE> &A)
{
	A.scan(s);		return s;
}

_NAMESPACE_FCA_END_

#endif // _SPARSE_MATRIX_H_
