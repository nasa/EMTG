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
