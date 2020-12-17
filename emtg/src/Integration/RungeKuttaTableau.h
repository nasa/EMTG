#ifndef RUNGEKUTTATABLEAU_H
#define RUNGEKUTTATABLEAU_H

#include <iostream>
#include <fstream>
#include <iomanip>

#include "IntegrationCoefficients.h"

namespace EMTG {
	namespace Integration
	{
		class RungeKuttaTableau : public IntegrationCoefficients
		{

		public:
			// constructors
			RungeKuttaTableau();

			// destructor
			~RungeKuttaTableau();

			// methods
            void resizeArrays();

            // getters for private variables
            inline size_t getNumStages() const { return this->num_stages; }
            inline math::Matrix<double> getA() const { return this->A; }
            inline math::Matrix<double> getBUpper() const { return this->bUpper; }
            inline math::Matrix<double> getBLower() const { return this->bLower; }
            inline math::Matrix<double> getC() const { return this->c; }

            // setters
            //inline void setNumStages(const size_t & num_stages) { this->num_stages = num_stages; }
            //inline void setA(const math::Matrix<double> & A) { this->A = A; }
            //inline void setBUpper(const math::Matrix<double> & bUpper) { this->bUpper = bUpper; }
            //inline void setBLower(const math::Matrix<double> & bLower) { this->bLower = bLower; }
            //inline void setC(const math::Matrix<double> & c) { this->c = c; }

			// write to file (for debugging)
			void writeTableauToFile(const std::string & fileName);

			// fields

		protected:
			// fields
			size_t num_stages; // number of stages in the method
			math::Matrix<double> A; // A matrix
			math::Matrix<double> bUpper; // b vector for higher-order method if variable-step. b vector for only method if fixed step
			math::Matrix<double> bLower; // b vector for lower-order method if variable-step. not used if fixed step
			math::Matrix<double> c; // c vector

		}; // end class definition

	} // end integration namespace
} // end EMTG namespace


#endif
