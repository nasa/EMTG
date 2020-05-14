//PJRN Scalatron class
//for use with EMTG and CSALT
//Jacob Englander 8/16/2018

#pragma once
#include "ScalatronBase.h"

namespace Scalatron
{
	class ScalatronPJRN : public ScalatronBase
	{
	public:
		//default constructor
		ScalatronPJRN() {};

		//constructor with arguments
		ScalatronPJRN(const std::vector<double>& X0,
			const std::vector<double>& Xlowerbounds,
			const std::vector<double>& Xupperbounds,
			const std::vector<double>& F0,
			const std::vector<double>& Flowerbounds,
			const std::vector<double>& Fupperbounds,
			const std::vector<size_t>& iGfun,
			const std::vector<size_t>& jGvar,
			const std::vector<double>& G0,
			size_t ObjectiveIndex);

		//compute scale factors
		virtual void compute_scale_factors();

	protected:
		virtual void compute_objective_scaling();

		virtual void compute_inequality_constraint_scaling();
	};
}//end namespace Scalatron