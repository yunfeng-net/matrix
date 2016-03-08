//////////////////////////////////////////////////////////////////////
// $Id: functor.h,v 1.4 2003/10/03 05:37:35 lyf1998 Exp $
// content: interface for the VecFunctor class of CVector<TYPE> 
// Writed by Yunfeng Liu, 2003/01
//
//////////////////////////////////////////////////////////////////////
#if !defined(_FUNCTOR_H_)
#define _FUNCTOR_H_

_NAMESPACE_FCA_BEGIN_

/// Nonlinear VecFunctor
template<class Matrix>
class CVecFunctor
{
public:
	typedef typename Matrix::value_type value_type;
	typedef CVector<value_type> 		vector_type;
	typedef Matrix						matrix_type;
	/// constructor
	CVecFunctor(int N) : m_vResult(N), m_Deriv(N, N) { }
	virtual ~CVecFunctor() { }
	/// dimension of vector	int size() const { return m_vResult.size(); }
	/// the cached result
	vector_type& result() { return m_vResult; }
	/// the cached derivative
	Matrix& derivative() { return m_Deriv; }	// f'(x)
	/// compute the derivative and result by the para
	virtual  vector_type& result(const CDblVector& para)=0;	// f(x)

protected:
	vector_type m_vResult;
	Matrix m_Deriv;
};

/// Nonlinear VecFunctor with the time parameter
template<class Matrix>
class COdeVecFunctor : public CVecFunctor<Matrix>
{
public:
	typedef typename CVecFunctor<Matrix>::value_type	value_type;
	typedef typename CVecFunctor<Matrix>::vector_type	vector_type;
	typedef typename CVecFunctor<Matrix>::matrix_type	matrix_type;
	/// constructor
	COdeVecFunctor(int N, int iOrd = 2, double fEPS = 1e-6)
		: CVecFunctor<Matrix>(N), m_aHistory(4), m_EB(fEPS, fEPS)
	{
		assert(iOrd<=2);	m_iMaxOrd = iOrd;	m_iOrd = iOrd;	// default order
		for(int i = 0; i<4; i++ ) m_aHistory[i] = new vector_type(N);
	}
	/// destructor
	virtual ~COdeVecFunctor()
	{
		erase_ptr_vector(m_aHistory);
	}
	/// accept the new step and data, Predict data at the next step
	virtual bool Accept(double& fNewStep)
	{
		if(!CanAccept(fNewStep)) return false;
		m_fTime += m_aDelta[0];
		shift_right(m_aDelta, m_aDelta+3);
		shift_right(m_aHistory.begin(), m_aHistory.end());
		return true;	// always accept the step under the default condition
	}
	/// compute the data at the given time
	virtual vector_type& Predict(double fStep)
	{
		m_aDelta[0] = fStep;
		double alpha = fStep/m_aDelta[1];		// linear Prediction
		*m_aHistory[0] = add(scaled(*m_aHistory[1], -alpha), scaled(*m_aHistory[2], 1+alpha));
		return *m_aHistory[0];
	}
	/// set the time range of solution
	virtual void SetTime(double fStart, double fStep, double fStop)
	{
		m_fStart = fStart; m_fStep = fStep; m_fStop = fStop;  m_fTime = m_fStart;
		_STD_::fill(m_aDelta, m_aDelta+3, m_fStep);
	}
	/// set the initial vector
	void SetZeroPoint(const vector_type &aX)
	{
		for(int i = 0; i<4; i++) (*m_aHistory[i]) = aX;
	}
	/// shrink the time step when non-convergence
	virtual bool Shrink(double& fStep) { fStep /= 4; return fStep>m_fStep*1e-6; }

protected:
	_STD_::vector<vector_type*> m_aHistory;
	double m_fStart, m_fStep, m_fStop, m_aDelta[3], m_fTime;
	int m_iOrd, m_iMaxOrd;
	CRange m_EB;
	
	/// Judge LPE and predict the new step
	bool CanAccept(double &fNewStep) const
	{
		double fError = GetLPE();
		double fRatio = (m_iOrd==1) ? exp(log(m_EB.LowerBound()/fError)/3)
			: sqrt(m_EB.LowerBound()/fError);
		if(fRatio<0.9) // shrink the step
		{
			fNewStep = _STD_::max(fRatio, 0.1)*m_aDelta[0];
		}
		else // relax the step
		{
			fNewStep = _STD_::min(fRatio*m_aDelta[0], m_aDelta[0]+m_aDelta[1]);
		}
		return fRatio>=0.9;
	}
	/// compute LPE
	double GetLPE() const
	{
		double fError = 1e-15;
		double aDD[3][3];
		
		if(m_iOrd>1) // high order
			for(int i = 0; i<size(); i++)
			{
				aDD[0][1] = ((*m_aHistory[0])[i]-(*m_aHistory[1])[i])/m_aDelta[0];
				aDD[1][2] = ((*m_aHistory[1])[i]-(*m_aHistory[2])[i])/m_aDelta[1];
				aDD[2][3] = ((*m_aHistory[2])[i]-(*m_aHistory[3])[i])/m_aDelta[2];
				aDD[0][2] = (aDD[0][1]-aDD[1][2])/(m_aDelta[0]+m_aDelta[1]);
				aDD[1][3] = (aDD[1][2]-aDD[2][3])/(m_aDelta[1]+m_aDelta[2]);
				aDD[0][3] = (aDD[0][2]-aDD[1][3])/(m_aDelta[0]+m_aDelta[1]+m_aDelta[2]);
				double fErr = (aDD[0][2]/6+aDD[0][3]*m_aDelta[0]/24);
				fErr = m_aDelta[0]*m_aDelta[0]*m_aDelta[0];
				if(fabs(fErr)>fError) fError = fabs(fErr);
			}
		else // low order
			for(int i = 0; i<size(); i++)
			{
				aDD[0][1] = ((*m_aHistory[0])[i]-(*m_aHistory[1])[i])/m_aDelta[0];
				aDD[1][2] = ((*m_aHistory[1])[i]-(*m_aHistory[2])[i])/m_aDelta[1];
				aDD[0][2] = (aDD[0][1]-aDD[1][2])/(m_aDelta[0]+m_aDelta[1]);
				double fErr = aDD[0][2]/6*m_aDelta[0]*m_aDelta[0]*m_aDelta[0];
				if(fabs(fErr)>fError) fError = fabs(fErr);
			}
		return fError;
	}
};

_NAMESPACE_FCA_END_

#endif // _FUNCTOR_H_
