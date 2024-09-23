//PJRN Scalatron class
//for use with EMTG and CSALT
//Jacob Englander 8/16/2018

#include <cmath>
#include "ScalatronPJRN.h"

namespace Scalatron
{	
	//constructor with arguments
	ScalatronPJRN::ScalatronPJRN(const std::vector<double>& X0,
		const std::vector<double>& Xlowerbounds,
		const std::vector<double>& Xupperbounds,
		const std::vector<double>& F0,
		const std::vector<double>& Flowerbounds,
		const std::vector<double>& Fupperbounds,
		const std::vector<size_t>& iGfun,
		const std::vector<size_t>& jGvar,
		const std::vector<double>& G0,
		size_t ObjectiveIndex)
		: ScalatronBase::ScalatronBase(X0, Xlowerbounds, Xupperbounds, F0, Flowerbounds, Fupperbounds, iGfun, jGvar, G0, ObjectiveIndex)
	{}

	//compute scale factors
	void ScalatronPJRN::compute_scale_factors()
	{
		//scale objective
		this->compute_objective_scaling();

		//scale equality constraints
		this->compute_equality_constraint_scaling();

		//scale inequality constraints
		this->compute_inequality_constraint_scaling();
	}//end compute_scale_factors()

	void ScalatronPJRN::compute_objective_scaling()
	{
		size_t Findex = this->ObjectiveIndex;

		double Kf_thisConstraint = 0.0;

		for (size_t Gindex = 0; Gindex < this->nG; ++Gindex)
		{
			if (this->iGfun[Gindex] == Findex)
			{
				size_t Xindex = this->jGvar[Gindex];

				double candidate_Kf = std::fabs(this->G0[Gindex] / this->Kx[Xindex]);

				Kf_thisConstraint = (candidate_Kf > Kf_thisConstraint) ? candidate_Kf : Kf_thisConstraint;
			}
		}//end loop over Jacobian entries

		this->Kf[Findex] = Kf_thisConstraint;
	}//end compute_objective_scaling

	void ScalatronPJRN::compute_inequality_constraint_scaling()
	{
		//PJRN method
		//from Sagliano, "Performance analysis of linear and nonlinear techniques for automatic scaling of discretized control problems"
		//Operations Research Letters 42 (2014) 213-216
		//https://core.ac.uk/download/pdf/31014768.pdf
		//Equation 7

		for (size_t Findex : this->InequalityConstraintIndices)
		{
			double Kf_thisConstraint = 0.0;

			for (size_t Gindex = 0; Gindex < this->nG; ++Gindex)
			{
				if (this->iGfun[Gindex] == Findex)
				{
					size_t Xindex = this->jGvar[Gindex];

					double candidate_Kf = std::fabs(this->G0[Gindex] / this->Kx[Xindex]);

					Kf_thisConstraint = (candidate_Kf > Kf_thisConstraint) ? candidate_Kf : Kf_thisConstraint;
				}
			}//end loop over Jacobian entries

			this->Kf[Findex] = Kf_thisConstraint;
		}//end loop over constraints
	}//end compute_inequality_constraint_scaling()
}//end namespace Scalatron