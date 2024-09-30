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

//Scalatron factory
//Jacob Englander 8/16/2018

#include "ScalatronFactory.h"
#include <iostream>
#include <exception>
#include <string>

namespace Scalatron
{
	ScalatronBase* Create_Scalatron(const ScalatronType& myScalatronType)
	{
		switch (myScalatronType)
		{
		case ScalatronType::PJRN:
			return new ScalatronPJRN();

		default:
            throw std::invalid_argument("Your chosen Scalatron type does not exist. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
		}
	}//end Create_Scalatron()

	ScalatronBase* Create_Scalatron(const std::vector<double>& X0,
		const std::vector<double>& Xlowerbounds,
		const std::vector<double>& Xupperbounds,
		const std::vector<double>& F0,
		const std::vector<double>& Flowerbounds,
		const std::vector<double>& Fupperbounds,
		const std::vector<std::size_t>& iGfun,
		const std::vector<std::size_t>& jGvar,
		const std::vector<double>& G0,
		std::size_t ObjectiveIndex,
		const ScalatronType& myScalatronType)
	{
		switch (myScalatronType)
		{
		case ScalatronType::PJRN:
			return new ScalatronPJRN(X0, Xlowerbounds, Xupperbounds, F0, Flowerbounds, Fupperbounds, iGfun, jGvar, G0, ObjectiveIndex);

		default:
            throw std::invalid_argument("Your chosen Scalatron type does not exist. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));;
		}
	}//end Create_Scalatron()
}//close namespace Scalatron
