// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2013 - 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Other Rights Reserved.

// Licensed under the NASA Open Source License (the "License"); 
// You may not use this file except in compliance with the License. 
// You may obtain a copy of the License at:
// https://opensource.org/licenses/NASA-1.3
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
// express or implied.   See the License for the specific language
// governing permissions and limitations under the License.

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