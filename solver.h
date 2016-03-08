/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

//////////////////////////////////////////////////////////////////////
// $Id: solver.h,v 1.5 2003/10/03 05:37:35 lyf1998 Exp $
// content: Solvers in Flexible Computing Architecture
// Writed by Yunfeng Liu, 2003/04
//
//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------
// 
// Solvers solves the functor by the aid of the other solvers
//
// TimeFunctor -> CNonlinearFunctor -> Matrix
//				   |                    |
// COdeSolver   <- CNonlinearSolver  <- CLinearSolver
//
//-----------------------------------------------------------

#if !defined(_SOLVER_H_)
#define _SOLVER_H_

#include "algorithm.h"
#include "functor.h"

_NAMESPACE_FCA_BEGIN_

//-----------------------------------------------------------
// Linear Solver : CLuSolver, CGmrSolver
//-----------------------------------------------------------
/// Interface for the linear solver
template<class Matrix>
class CLinearSolver
{
public: // easy to use and combine
	typedef typename Matrix::value_type value_type;
	typedef CVector<value_type> 		vector_type;

	CLinearSolver(Matrix *pA) { m_pA = pA; }
	virtual ~CLinearSolver() {}

	bool valid() { return m_bValid; } // validation
	bool set(Matrix *pA) { m_pA = pA; return (m_bValid = m_pA!=NULL && prepare()); };
	int size() const { return m_pA->size(); }
	
	virtual bool solve(vector_type& aB, vector_type& aX)=0;

protected:
	Matrix *m_pA; // coefficient
	bool m_bValid;
	virtual bool prepare() = 0; // prepare to solve
};

/// Solver of the linear equation by the LU decomposition
template<class Matrix>
class CLuSolver : public CLinearSolver<Matrix>
{
public: // by the LU decomposition
	typedef typename Matrix::value_type value_type;
	typedef CVector<value_type> 		vector_type;

	/// Constructor
	CLuSolver(Matrix *pA = NULL) : CLinearSolver<Matrix>(pA) { prepare(); }
	virtual ~CLuSolver() {}
	
	/// Solving function
	virtual bool solve(CVector<value_type>& aB, CVector<value_type>& aX)
		{ m_pA->SolveLU(m_aIndex, aB, aX); return true;}
protected:
	CMapping m_aIndex;
	virtual bool prepare()
		{ 	return m_pA && m_pA->LU(m_aIndex);	}
};

/// Solver of the linear equation by the GMRES
template<class Matrix>
class CGmrSolver : public CLinearSolver<Matrix>
{
public: // by the GMRES iteration
	typedef typename Matrix::value_type value_type;
	typedef CVector<value_type> 		vector_type;
	
	/// Constructor
	CGmrSolver(Matrix *pA, CIteRange EB) : CLinearSolver<Matrix>(pA)
		{ m_EB = EB; prepare(); }
	virtual ~CGmrSolver() {}

	/// Solving function
	virtual bool solve(vector_type& aB, vector_type& aX)
	{
		const int M = m_EB.count();
		m_vR = sub(aB, mult(*m_pA, aX));
		double w = m_vR.norm(2);
		if(m_EB.IsZero(w)) return true;
		m_Q[0] = scaled(m_vR, 1.0/w);

		for(m_nIte = 0; m_nIte<M; m_nIte++)
		{
			int N = GetHQ(); // construct H and the orthogonal base
//			assert(TestHQ(N)); // disable for schmit orthogonalizeation
			GetLS(N, w, aX); // least square problem : min|| Hy-||r0||e1 ||2
			m_vR = sub(aB, mult(*m_pA, aX));
			w = m_vR.norm(2);
			if(m_EB.IsZero(w)) return true;
			m_Q[0] = scaled(m_vR, 1.0/w);
		}
		return false;
	}
	int count() const{ return m_nIte; }
protected:
	int m_nIte;
	CIteRange m_EB;
	vector_type m_vR, m_vY;
	vector<vector_type> m_Q;
	CMatrix<value_type> m_H;

	/// prepare the iteration matrix and vector
	virtual bool prepare()
	{
		if(!m_pA) return false;
		const int N = m_pA->size(), M = m_EB.count();
		m_vR.resize(N);	m_vY.resize(M);
		m_Q.resize(M+1);	m_H.resize(M,M+1);
		for(int i = 0; i<M+1; i++) m_Q[i].resize(N);
		return true;
	}
	/// Construct the Hessenberg matrix and the orthogonal bases
	int GetHQ()
	{
		for(int i = 0; i<m_EB.count(); i++) // construct H & Q
		{
			m_Q[i+1] = mult(*m_pA, m_Q[i]);
			for(int k = 0; k<=i; k++) m_H[i][k] = dot(m_Q[i+1], m_Q[k]);
			mult_eq(_STD_::minus<value_type>(), m_H[i], i+1, m_Q[i+1]);
			m_H[i][i+1] = m_Q[i+1].norm(2);
			if(m_EB.IsZero(m_H[i][i+1]) || i==m_EB.count()-1) return i+1;
			m_Q[i+1] = scaled(m_Q[i+1], (value_type)1.0/m_H[i][i+1]);
		}
		return m_EB.count();
	}
	/// Test the Q matrix
	bool TestHQ(int N)
	{
		cout << " H = " << m_H << endl;
		CMatrix<value_type> mtA(N+1,N+1);
		mtA.clear();
		for(int k = 0; k<N; k++)
			for(int r = 0; r<N+1; r++) cout << " Q*Q'= " << dot(m_Q[k], m_Q[r]) << endl;
		for(int i = 0; i<N; i++) for(int j = 0; j<i+2; j++)
			mtA[i][j] = m_H[i][j] - dot(m_Q[i], mult(*m_pA, m_Q[j]));
		cout << " E = " << mtA << endl;
		return true;
		return m_EB.IsZero(mtA.norm());
	}
	/// Get the solution of the LS problem by the LQ method
	void GetLS(int N, double w, vector_type& aX)
	{
		m_vR.fill(0);	m_vR[0] = w;
		for(int k = 0; k<N-1; k++) // work out R
		{
			CRotation<value_type> rot(m_H[k][k], m_H[k][k+1]);
			m_H.Givens_col(rot, k, k+1, N);
			rot(m_vR[k], m_vR[k+1]);
		}
		m_H.SolveTL(m_vR, N, m_vY); // solve R * Y = ||r0||*e1
		mult_eq(_STD_::plus<value_type>(), m_vY.begin(), N, aX);
	}
	/// Multiply a vector array with a sized vector
	template<class iterator, class functor>
	void mult_eq(functor func, iterator F, int N, vector_type& aX)
	{
		assert(m_Q.size()>0);
		for(unsigned i = 0; i<m_Q[0].size(); i++)
			for(int j = 0; j<N; j++) // work out the new aX
				aX[i] = func(aX[i], *(F+j)*m_Q[j][i]);
	}

};

//-----------------------------------------------------------
// Non-Linear Solver : CNrSolver, CHmSolver
//-----------------------------------------------------------
/// Interface for the linear solver
template <class Functor>
class CNonLinearSolver
{
public: // Newton-Ralphson
	typedef typename Functor::value_type	value_type;
	typedef typename Functor::vector_type	vector_type;
	typedef typename Functor::matrix_type	matrix_type;
	typedef CLinearSolver<matrix_type>		LinearSolver;

	CNonLinearSolver(LinearSolver *pSolver, Functor *pFunc, CIteRange EB)
		{ m_pSolver = pSolver; m_EB = EB; set(pFunc); }
	virtual ~CNonLinearSolver() { delete m_pSolver; }

	virtual bool solve(CVector<value_type>& aX) = 0;
	
	CVector<value_type>& result() { return m_pFunc->result(); }

	bool valid() const { return m_bValid; }
	bool set(Functor *pFunc) // validation
	{ 
		m_pFunc = pFunc;
		m_bValid = (m_pFunc!=NULL) && m_pSolver->set(&m_pFunc->derivative());
		return m_bValid;
	}
	int size() const { return m_pFunc->size(); }
	int count() const { return m_nIte; }

protected:
	Functor* m_pFunc; // function to solve
	LinearSolver *m_pSolver;
	CIteRange m_EB;
	bool m_bValid;
	int m_nIte;
};

/// Solver of the nonlinear equations in a damping-NR way
template <class Functor>
class CNrSolver : public CNonLinearSolver<Functor>
{
public: // Newton-Ralphson
	typedef typename CNonLinearSolver<Functor>::value_type value_type;
	typedef typename CNonLinearSolver<Functor>::vector_type vector_type;
	typedef typename CNonLinearSolver<Functor>::matrix_type matrix_type;
	typedef typename CNonLinearSolver<Functor>::LinearSolver LinearSolver;
	/// constructor
	CNrSolver(LinearSolver *pSolver, Functor *pFunc, CIteRange EB)
		: CNonLinearSolver<Functor>(pSolver, pFunc, EB)
			, m_aF(pFunc->size()) {}
	virtual ~CNrSolver() {}
	/// solve from an initial point
	virtual bool solve(vector_type& aX) // initial value
	{
		double fFMax = 1000;
		
		for(m_nIte = 0; m_nIte<m_EB.count(); m_nIte++)
		{
			vector_type &aF = m_pFunc->result(aX); // f0 = f(x)
			if(m_EB.IsZero(aF.norm())) return true;
			if(!set(m_pFunc) || !m_pSolver->solve(aF, m_aF)) return false; // non-convergence
			double fNorm = m_aF.norm();
			if(fNorm>fFMax)
				aX -= scaled(m_aF, 0.5 * fFMax / fNorm); // not converge
			else
				aX -= m_aF, fFMax = fNorm; // converge
		}
		return false;
	}
protected:
	vector_type m_aF;
};

/// Solver of the nonliear equations in a homotopy way
template <class Functor>
class CHmSolver : public CNonLinearSolver<Functor>
{
public: // Homotopy solver
	typedef typename CNonLinearSolver<Functor>::value_type value_type;
	typedef typename CNonLinearSolver<Functor>::vector_type vector_type;
	typedef typename CNonLinearSolver<Functor>::matrix_type matrix_type;
	typedef typename CNonLinearSolver<Functor>::LinearSolver LinearSolver;
	/// constructor
	CHmSolver(LinearSolver *pSolver, Functor *pFunc, CIteRange& EB)
		: CNonLinearSolver<Functor>(Func, Solver, EB) {}
	virtual ~CHmSolver() {}
	
	/// solve from a initial point
	virtual bool solve(vector_type& aX) // initial value
	{
		return true;
	}
};

/// Solver of the nonliear equations mixed with ODE equations
template<class TimeFunctor>
class COdeSolver 
{
public:
	typedef typename TimeFunctor::value_type	value_type;
	typedef typename TimeFunctor::vector_type	vector_type;
	typedef typename TimeFunctor::matrix_type	matrix_type;
	typedef CNonLinearSolver<TimeFunctor>		NLnrSolver;

	COdeSolver(NLnrSolver *pSolver, TimeFunctor *pFunc)
		{ m_pSolver = pSolver; m_pFunc = pFunc; m_pSolver->set(m_pFunc);}
	virtual ~COdeSolver() { delete m_pSolver; }

	// fStep is the estimated step
	/// solve in the given range
	bool solve(double fStep, double fStop, double fStart = 0)
	{
		double fDelta = fStep, fGamma; // the real step & the Predicted step
		m_pFunc->SetTime(fStart, fStep, fStop);

		m_nIte = m_nPoint = 0;
		for(m_fTime = fStart; m_fTime<=fStop; m_fTime+=fGamma) // iterate
			for(;;) // loop for non-convergence
			{
				if(!m_pSolver->solve(m_pFunc->Predict(fDelta)))
					{ if(!m_pFunc->Shrink(fDelta)) return false; } // error
				else // solved
				{
					fGamma = fDelta;	// keep the used delta
					if(m_pFunc->Accept(fDelta))
					{
						m_nPoint++;
						m_nIte += m_pSolver->count();
						break; // converge
					}
				}
			}
		return true;
	}
	
	vector_type& result() { return m_pFunc->result(); }
	int count(int nTotal = 0) const { return (!nTotal)?m_nPoint:m_nIte; }
	bool valid() const { return m_bValid; }
	/// set the functor to solve
	bool set(TimeFunctor *pFunc) // validation
	{ 
		m_pFunc = pFunc;
		m_bValid = (m_pFunc!=NULL) && m_pSolver->set(pFunc);
		return m_bValid;
	}
protected:
	TimeFunctor *m_pFunc; // function to solve
	NLnrSolver *m_pSolver;
	bool m_bValid;
	int m_nPoint, m_nIte;
	double m_fTime;
};

_NAMESPACE_FCA_END_

#endif // _SOLVER_H_
