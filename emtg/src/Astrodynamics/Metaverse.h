// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2014 - 2017 United States Government as represented by the
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

//EMTGv9 Metaverse class
//this is just a container for a vector of universes and a SplineEphemUniverse object
//Jacob Englander 1-5-2018

#include "missionoptions.h"
#include "universe.h"

#ifdef SPLINE_EPHEM
#include "SplineEphem_universe.h"
#endif

namespace EMTG
{
    namespace Astrodynamics
    {
        class Metaverse
        {
        public:
            Metaverse() {};
            Metaverse(missionoptions& options);

            ~Metaverse();

            void initialize(missionoptions& options);

            universe* getUniversePointer(const size_t& j) { return &this->TheUniverse[j]; }

        protected:
            SplineEphem::universe SplineUniverse;
            std::vector<universe> TheUniverse;
            missionoptions* myOptions;
        };
    }//close namespace Astrodynamics
}//close namespace EMTG