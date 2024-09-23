#include "RungeKuttaTableau.h"

namespace EMTG {
	namespace Integration
	{
		// standard constructor
		RungeKuttaTableau::RungeKuttaTableau()
		{
		}

		// destructor
		RungeKuttaTableau::~RungeKuttaTableau()
		{
		}

		// methods
		void RungeKuttaTableau::resizeArrays()
		{
			this->A.resize(this->getNumStages(), this->getNumStages());
			this->bUpper.resize(this->getNumStages(), 1);
			this->bLower.resize(this->getNumStages(), 1);
			this->c.resize(this->getNumStages(), 1);
		}

		void RungeKuttaTableau::writeTableauToFile(const std::string & fileName)
		{
			std::ofstream f;
			f.open(fileName);

            size_t num_stages = this->getNumStages();

			// A 
			for (size_t i = 0; i < num_stages; ++i)
			{
				for (size_t j = 0; j < num_stages; ++j)
				{
					f << std::fixed << std::setprecision(16) << "A(" << i << ", " << j << ") = " << this->getA()(i, j) << "\n";
				}
			}

			// bUpper
			for (size_t i = 0; i < num_stages; ++i)
			{
				f << "bUpper(" << i << ") = " << this->getBUpper()(i, 0) << "\n";
			}

			// bLower
			if (this->getHasVariableStepCoefficients())
			{
				for (size_t i = 0; i < num_stages; ++i)
				{
					f << "bLower(" << i << ") = " << this->getBLower()(i, 0) << "\n";
				}
			}

			// c 
			for (size_t i = 0; i < num_stages; ++i)
			{
				f << "c(" << i << ") = " << this->getC()(i, 0) << "\n";
			}

			f.close();
		}
	} // end integration namespace
} // end EMTG namespace