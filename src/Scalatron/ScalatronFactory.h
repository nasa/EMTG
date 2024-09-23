//Scalatron factory
//Jacob Englander 8/16/2018

#pragma once

#include "ScalatronPJRN.h"

namespace Scalatron
{
	enum ScalatronType { PJRN };

	ScalatronBase* Create_Scalatron(const ScalatronType& myScalatronType);

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
		const ScalatronType& myScalatronType);
}//close namespace Scalatron
