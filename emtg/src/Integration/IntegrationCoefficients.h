#ifndef INTEGRATIONCOEFFICIENTS_H
#define INTEGRATIONCOEFFICIENTS_H

#include "EMTG_Matrix.h"

namespace EMTG {
	namespace Integration
	{
		class IntegrationCoefficients
		{

		public:
			// constructors
			IntegrationCoefficients();

			// destructor
			~IntegrationCoefficients();

            // getters for private variables
            inline bool getHasVariableStepCoefficients() const { return this->has_variable_step_coeffs; }

            // setters for private variables
            inline void setHasVariableStepCoefficients(const bool & has_variable_step_coeffs) { this->has_variable_step_coeffs = has_variable_step_coeffs; }

		private:
			bool has_variable_step_coeffs;

		}; // end class definition

	} // end integration namespace
} // end EMTG namespace


#endif
