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

//Propulator
//contains a vector of PropulatorJourney objects

#include "missionoptions.h"
#include "universe.h"
#include "Spacecraft.h"
#include "EMTG_Matrix.h"

#ifdef SPLINE_EPHEM
#include "SplineEphem_universe.h"
#endif


#include "PropulatorJourney.h"

#include <vector>
#include <string>

#include "boost/filesystem.hpp"
#include "boost/ptr_container/ptr_vector.hpp"

#ifdef PROPULATOR_PYTHON_INTERFACE
#include "boost/python.hpp"
#include "boost/python/list.hpp"
#include "boost/python/extract.hpp"
#endif

namespace EMTG
{
    namespace Propulator
    {
        class Propulator
        {
        public:
            Propulator();
            Propulator(const std::string& options_file_name);

            ~Propulator();

            void initialize(const std::string& options_file_name);

            void setJourneyIndex(const size_t& journeyIndex) { this->journeyIndex = journeyIndex; }

            void setStepSize(const double& stepSize) { this->myPropulatorJourneys[this->journeyIndex].setStepSize(stepSize); }

            math::Matrix<double> getSTM() { return this->myPropulatorJourneys[this->journeyIndex].getSTM(); }

            void propulate(const math::Matrix<double>& InputState,
                const math::Matrix<double>& ControlVector,
                const double& DutyCycle,
                const double& propagationTime,
                math::Matrix<double>& OutputState,
                const bool& needSTM);

            void propulate(const math::Matrix<double>& InputState,
                const double& propagationTime,
                math::Matrix<double>& OutputState,
                const bool& needSTM);

#ifdef PROPULATOR_PYTHON_INTERFACE
			boost::python::list propulateThrust(const boost::python::list& InputState,
				const boost::python::list& ControlVector,
				const double& DutyCycle,
				const double& propagationTime,
				const bool& needSTM);

			boost::python::list propulateCoast(const boost::python::list& InputState,
				const double& propagationTime,
				const bool& needSTM);

			boost::python::list getpythonSTM();
#endif

        protected:
            size_t journeyIndex;
            boost::ptr_vector<PropulatorJourney> myPropulatorJourneys;

            missionoptions myOptions;
            std::vector<Astrodynamics::universe > TheUniverse;
            HardwareModels::Spacecraft mySpacecraft;
            std::vector<::boost::filesystem::path> SPICE_files_required;

#ifdef SPLINE_EPHEM
            SplineEphem::universe* mySplineUniverse;
#endif

        };


#ifdef PROPULATOR_PYTHON_INTERFACE
		BOOST_PYTHON_MODULE(propulator)
		{
			boost::python::class_<Propulator>("propulator", boost::python::init<std::string>())
											  .def("propulateThrust", &Propulator::propulateThrust)
											  .def("propulateCoast", &Propulator::propulateCoast)
											  .def("getSTM", &Propulator::getpythonSTM)
											  .def("setJourneyIndex", &Propulator::setJourneyIndex)
                                              .def("setStepSize", &Propulator::setStepSize);

		}
#endif
    }//close namespace Propulator
}//close namespace EMTG
