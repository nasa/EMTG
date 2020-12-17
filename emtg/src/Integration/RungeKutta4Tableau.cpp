#include "RungeKutta4Tableau.h"

namespace EMTG {
	namespace Integration
	{
		// standard constructor
		RungeKutta4Tableau::RungeKutta4Tableau()
		{
			this->setHasVariableStepCoefficients(false); // fixed-step method
			this->num_stages = 4;
			this->resizeArrays(); // set sizes of coefficient arrays

			// populate A
			this->A.assign_zeros(); // mostly zeros
			this->A(1, 0) = 1.0 / 2.0;
			this->A(2, 1) = 1.0 / 2.0;
			this->A(3, 2) = 1.0;

			// populate bUpper
			this->bUpper(0) = 1.0 / 6.0;
			this->bUpper(1) = 1.0 / 3.0;
			this->bUpper(2) = 1.0 / 3.0;
			this->bUpper(3) = 1.0 / 6.0;

			// populate c
			this->c(0, 0) = 0.0;
			this->c(1, 0) = 1.0 / 2.0;
			this->c(2, 0) = 1.0 / 2.0;
			this->c(3, 0) = 1.0;
			
		}
	} // end integration namespace
} // end EMTG namespace