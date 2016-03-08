//////////////////////////////////////////////////////////////////////
// $Id: polynomial.h,v 1.1 2003/09/21 07:54:33 lyf1998 Exp $
// content: interface for the CPolynomial class.
// Writed by Yunfeng Liu, 2002/05
//
// CTable, CPolynomial, CRational
// 
//////////////////////////////////////////////////////////////////////

#if !defined(_POLYNOMIAL_H_)
#define _POLYNOMIAL_H_

#include "matrix.h"
_NAMESPACE_FCA_BEGIN_

/// Tabletemplate<class IDX, class ITEM>
class CTable
{
public:
	typedef enum { eStep, eLinear } EInterpol;

	// Constructor and destructor
	CTable(int nLen=0, EInterpol eInterpol=eStep) { SetDim(nLen,eInterpol); }
	CTable(const CTable& A) { 	m_aIndex = A.m_aIndex; m_aItem = A.m_aItem;	}
	virtual ~CTable() {}

    // Dimension operation
	int GetDim(void) const { return m_aIndex.GetDim(); }
    void SetDim(int nDim, EInterpol eInterpol=eStep)
	{ 
		m_eInterpol = eInterpol;
		m_aIndex.SetDim(nDim); m_aItem.SetDim(nDim); // automatic growth
	}  
	CTable& Add(CVector<IDX>& A, CVector<ITEM>& B) // append the items
	{
		ASSERT(A.GetDim()==B.GetDim());
		for(int i = 0; i<A.GetDim(); i++) { m_aIndex.Add(A[i]); m_aItem.Add(B[i]); }
		return *this;
	}

	// Data access
	CVector<IDX>& GetIndex() { return m_aIndex; }
	CVector<ITEM>& GetItem() { return m_aItem; }

	IDX& Index(int nIndex) { return m_aIndex[nIndex]; }
	ITEM& Item(int nIndex) { return m_aItem[nIndex]; }
	ITEM& operator ()(const IDX& X); 

protected:
	EInterpol m_eInterpol;
	CVector<IDX> m_aIndex; // 1 dimension
	CVector<ITEM> m_aItem; // different data type

	int FindIndex(const IDX x) // By bisection
	{
		int a = 0, b = GetDim()-1, c = (a+b)/2;
		while(m_aIndex[c]!=x)
		{
			if(x>m_aIndex[c]) a = c+1; else b = c-1; // Switch to another half region
			if(a>b) return b;
			c = (a+b)/2; // Divide into 2 pieces
		}
		return c;
	}
};

// the common classes derived from CTable
typedef CTable<double,double> CDblTable;
typedef CTable<double,Complex> CCmplTable;
typedef CTable<double,CDblMatrix> CDblMatTable;

class CCmplMatTable : public CTable<double,CCmplMatrix>
{
public:
	CCmplMatTable() {}
	virtual ~CCmplMatTable() {}

	CCmplMatrix& operator ()(const double& X)
	{ 
		ASSERT(GetDim()>1);
		int i = FindIndex(X);
		if(i<0) // out of the left bound
			i = 0;
		else if(i>=m_aItem.GetDim()-1) // out of the right bound
			i = (m_eInterpol==eStep) ? m_aItem.GetDim()-1 : m_aItem.GetDim()-2;

		if(m_eInterpol==eStep) return m_aItem[i]; // Step
		else
		{
			CCmplMatrix item = m_aItem[i+1]; 
			item -= m_aItem[i];
			item *= (Complex)(X-m_aIndex[i]);
			return (item +=	m_aItem[i]); // Linear
		}
	}

};

/// Polynomial operation
template <class TYPE>
class CPolynomial
{
public: // constructor and destructor
	CPolynomial(int nLen=1, bool init=false) : m_aItem(nLen, init) {}
	CPolynomial(const CPolynomial& A) : m_aItem(A.m_aItem) {}
	virtual ~CPolynomial() {}

    // Operation about dimension
	CPolynomial& Add(const TYPE& A) { m_aItem.Add(A); return *this; }
	int GetDim(void) const { return m_aItem.GetDim(); }
    void SetDim(int nDim, bool init=false) { m_aItem.SetDim(nDim, init); }
	CPolynomial& Reduce(CLimit& lmt = CLimit()) // Reduce the zero item of the high order
	{ 
		for(int N = m_aItem.GetDim(); N>1 && lmt.IsZero(m_aItem[GetDim()-1]); N--);
		SetDim(N);
		return *this; 
	}
	
	// Access the item
    TYPE& operator [](int nIndex) { return m_aItem[nIndex]; }
    const TYPE& operator [](int nIndex) const { return m_aItem[nIndex]; }
	CVector<TYPE>& GetItem() { return m_aItem; }

	TYPE operator()(const TYPE& x) const // return the computation of X
	{ 
		TYPE sum=0; 
		for(int i=GetDim()-1; i>=0; i--) sum = x*sum + m_aItem[i];
		return sum;
	}
	CPolynomial Transform(const CPolynomial& X) // Replace the basic variable
	{
		CPolynomial<TYPE> Sum(1, true); Sum[0] = m_aItem[GetDim()-1];
		for(int i = GetDim()-2; i>=0; i--) // Recursive computation
			(Sum *= X)[0] += m_aItem[i]; // sum = X*sum + coeff[i]
		return Sum.Reduce();
	}
	CPolynomial& Derive()
	{
		for(int i = 0; i<GetDim()-1; i++) m_aItem[i] = m_aItem[i+1]*(i+1);
	    m_aItem.RemoveAt(GetDim()-1,1);
		return *this;
	}
		
	// Self-modification with a polynomial
	CPolynomial& operator =(const CPolynomial<TYPE>& A) // use the vector's operation
		{ m_aItem = A.m_aItem; return *this; }
	CPolynomial& operator +=(const CPolynomial<TYPE>& A);
	CPolynomial& operator -=(const CPolynomial<TYPE>& A);
	CPolynomial& operator *=(const CPolynomial<TYPE>& A);
	CPolynomial& operator /=(const CPolynomial<TYPE>& A);
	CPolynomial& operator %=(const CPolynomial<TYPE>& A);

	// Self-modification with an item
	CPolynomial& operator *=(const TYPE& A) { m_aItem *= A; return *this; }

protected:
	 CVector<TYPE> m_aItem;

};

// the common classes derived from CPolynomial
typedef CPolynomial<int> CIntPolynomial;
typedef CPolynomial<double> CDblPolynomial;
typedef CPolynomial<Complex> CCmplPolynomial;

/// Rational operation
template <class TYPE>
class CRational
{
public: // constructors
	CRational() {}
	CRational(const int A) { }
	CRational(const CRational& A) { (*this) = A; }
	CRational(const CPolynomial<TYPE>& A, const CPolynomial<TYPE>& B)
	{ 
		m_aNum = A;  m_aDenom = B; 
	}
	virtual ~CRational() {}

	// For the polynomial operation
	CPolynomial<TYPE>& GetNum() { return m_aNum; }
	CPolynomial<TYPE>& GetDenom() { return m_aDenom; }

	TYPE operator()(const TYPE& x) const // return the computation of X
	{	return m_aNum(x)/m_aDenom(x);	}
	CRational& Transform(const CRational& X); // Replace the basic variable
	CRational& Derive(); // Change to the derivative

	// Self-modification with rational
	CRational& operator = (const CRational& A)
		{ m_aNum = A.m_aNum; m_aDenom = A.m_aDenom; return *this; }
	CRational& operator +=(const CRational<TYPE>& A);
	CRational& operator -=(const CRational<TYPE>& A);
	CRational& operator *=(const CRational<TYPE>& A)
		{ m_aNum *= A.m_aNum; m_aDenom *= A.m_aDenom; return *this; }
	CRational& operator /=(const CRational<TYPE>& A)
		{ m_aNum *= A.m_aDenom; m_aDenom *= A.m_aNum; return *this; }

	// Self-modification with an item
	CRational& operator *=(const TYPE& A) { m_aNum *= A; return *this; }

protected:
	CPolynomial<TYPE> m_aNum, m_aDenom;

};

// the common classes derived from CRational
typedef CRational<int>CIntRational;
typedef CRational<double>CDblRational;
typedef CRational<Complex>CCmplRational;
inline CCmplRational& Random(CCmplRational& aCR) { return aCR; }

//--CTable-------------------------------------------------------------------

// Return the value by the interpolation
template<class IDX, class ITEM>
ITEM& CTable<IDX,ITEM>::operator ()(const IDX& X) 
{ 
	ASSERT(GetDim()>1);
	int i = FindIndex(X);
	if(i<0) // out of the left bound
		i = 0;
	else if(i>=m_aItem.GetDim()-1) // out of the right bound
		i = (m_eInterpol==eStep) ? m_aItem.GetDim()-1 : m_aItem.GetDim()-2;

	if(m_eInterpol==eStep) return m_aItem[i]; // Step
	else
	{
		ITEM item = m_aItem[i+1]; 
		item -= m_aItem[i];
		item *= (IDX)(X-m_aIndex[i]);
		return (item +=	m_aItem[i]); // Linear
	}
}

// Output to a stream
template<class IDX, class ITEM>
std::ostream& operator<<(std::ostream &s, CTable<IDX, ITEM> &A)
{ 	
	for(int i = 0; i<A.GetDim(); i++)	s << A.Index(i) << A.Item(i) << endl;  return s;
}

//--CPolynomial-------------------------------------------------------------------

/*
	C[i] = C[i] + A[i] where i<C.GetDim() && i<A.GetDim()	<1>
		 = C[i] where i<C.GetDim() && i>=A.GetDim()			<2> No need
		 = A[i] where i<A.GetDim() && i>=C.GetDim()			<3>
*/
template <class TYPE>
CPolynomial<TYPE>& CPolynomial<TYPE>::operator +=(const CPolynomial<TYPE>& A)
{
	int M = GetDim(); // Keep the original dimension
	int i, N = A.GetDim();
	SetDim(max(M,N)); // Change to the maximum dimension
	for(i = 0; i<min(M,N); i++) m_aItem[i] += A[i]; // region <1>
	for(i = min(M,N); i<max(M,N); i++) 	if(N>M)m_aItem[i] = A[i]; // region <3>
	return *this;
}

/*
	C[i] = C[i] - A[i] where i<C.GetDim() && i<A.GetDim()	<1>
		 = C[i] where i<C.GetDim() && i>=A.GetDim()			<2> No need
		 = -A[i] where i<A.GetDim() && i>=C.GetDim()		<3>
*/
template <class TYPE>
CPolynomial<TYPE>& CPolynomial<TYPE>::operator -=(const CPolynomial<TYPE>& A)
{
	int M = GetDim(); // Keep the original dimension
	int i, N = A.GetDim();
	SetDim(max(M,N)); // Change to the maximum dimension
	for(i = 0; i<min(M,N); i++) m_aItem[i] -= A[i]; // region <1>
	for(i = min(M,N); i<max(M,N); i++) 	if(N>M)m_aItem[i] = -A[i]; // region <3>
	return *this;
}

//	this *= A
template <class TYPE>
CPolynomial<TYPE>& CPolynomial<TYPE>::operator *=(const CPolynomial<TYPE>& A)
{
	int N = GetDim(); // Keep the original dimension
	CPolynomial<TYPE> C(N+A.GetDim()-1, true);
	for(int i = C.GetDim()-1; i>=0; i--) // travel all the items
		for(int j = min(i,N-1); i-j<A.GetDim() && j>=0; j--) // travel the pairs
			C[i] += (*this)[j] * A[i-j];
	return (*this = C);
}

/*
	C = 0 if A.GetDim()<B.GetDim()
	C = sigma(A[i]) * sigma(B[j]) = sigma( sigma(A[i])*B[j] )
*/
template <class TYPE>
CPolynomial<TYPE>& CPolynomial<TYPE>::operator /=(const CPolynomial<TYPE>& A)
{
	int M = GetDim(), N = A.GetDim();
	CPolynomial<TYPE> C(max(1, M-N+1), true);
	if(M<N) return (*this = C); // return 0;
	for(int i = M-1; i>=N-1; i--)
	{
		C[i+1-N] = m_aItem[i] / A[N-1];
		for(int j = N-2; j>=0; j--) // Travel A[j]
			m_aItem[i+j-(N-1)] -= C[i+1-N] * A[j]; // Eliminate
	}
	return (*this = C);
}

// modulator
template <class TYPE>
CPolynomial<TYPE>& CPolynomial<TYPE>::operator %=(const CPolynomial<TYPE>& A)
{
	if(GetDim()<A.GetDim()) return *this; // No need to compute

	int M = GetDim(), N = A.GetDim();
	for(int i = M-1; i>=N-1; i--)
		for(int j = N-1; j>=0; j--) // Travel A[j]
			m_aItem[i+j-(N-1)] -= m_aItem[j] / A[N-1] * m_aItem[i]; // Eliminate
	SetDim(N-1);
	return *this;
}

// Output to a stream
template <class TYPE>
std::ostream& operator<<(std::ostream &s, CPolynomial<TYPE> &A)
{ 	
	s << A.GetItem();  return s;
}

// Input from a stream
template <class TYPE>
std::istream & operator>>(std::istream &s, CPolynomial<TYPE> &A)
{	
	s >> A.m_aItem;	return *this;
}

//--CRational-------------------------------------------------------------------
template <class TYPE>
CRational<TYPE>& CRational<TYPE>::operator +=(const CRational<TYPE>& A)
{ 
	CPolynomial<TYPE> D = m_aDenom; D *= A.m_aNum;
	(m_aNum *= A.m_aDenom) += D; // Num*A.Denum + Denom*A.Num
	m_aDenom *= A.m_aDenom; // Denom*A.Denom
	return *this;
}

template <class TYPE>
CRational<TYPE>& CRational<TYPE>::operator -=(const CRational<TYPE>& A)
{ 
	CPolynomial<TYPE> D = m_aDenom; D *= A.m_aNum;
	(m_aNum *= A.m_aDenom) -= D; // Num*A.Denum - Denom*A.Num
	m_aDenom *= A.m_aDenom; // Denom*A.Denom
	return *this;
}

template <class TYPE>
CRational<TYPE>& CRational<TYPE>::Transform(const CRational& X) // Replace the basic variable
{
	int i;
	CPolynomial<TYPE> SumN(1,true), SumD(1,true);
	SumN[0] = m_aNum[m_aNum.GetDim()-1];
	SumD[0] = 1;
	CRational<TYPE> SumNum(SumN, SumD); // Compute the numerator at first
	for(i = m_aNum.GetDim()-2; i>=0; i--) // Recursive computation
	{
		SumN[0] = m_aNum[i];
		(SumNum *= X) += CRational<TYPE>(SumN,SumD); // sum = X*sum + coeff[i]
	}
	SumN[0] = m_aDenom[m_aDenom.GetDim()-1];
	CRational<TYPE> SumDenom(SumN, SumD); // Compute the numerator at first
	for(i = m_aDenom.GetDim()-2; i>=0; i--) // Recursive computation
	{
		SumN[0] = m_aNum[i];
		(SumDenom *= X) += CRational<TYPE>(SumN,SumD); // sum = X*sum + coeff[i]
	}
	return (*this = SumNum) /= SumDenom;
}

// Change to the derivative
template <class TYPE>
CRational<TYPE>& CRational<TYPE>::Derive()
{
	return *this;
}

// Output to stream
template <class TYPE>
std::ostream& operator<<(std::ostream &s, CRational<TYPE> &A)
{
	s << A.GetNum() << '/' << A.GetDenom();	return s;
}

// Input from stream
template <class TYPE>
std::istream & operator>>(std::istream &s, CRational<TYPE> &A)
{
	s >> A.GetNum() >> '/' >> A.GetDenom();	return s;
}

_NAMESPACE_FCA_END_

#endif // _POLYNOMIAL_H_
