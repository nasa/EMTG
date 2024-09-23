// EMTG: Evolutionary Mission Trajectory Generator
// An open-source global optimization tool for preliminary mission design
// Provided by NASA Goddard Space Flight Center
//
// Copyright (c) 2013 - 2020 United States Government as represented by the
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

#include "CartesianStateRepresentation.h"
#include "SphericalRADECStateRepresentation.h"
#include "SphericalAZFPAStateRepresentation.h"
#include "COEStateRepresentation.h"
#include "MEEStateRepresentation.h"
#include "IncomingBplaneStateRepresentation.h"
#include "OutgoingBplaneStateRepresentation.h"

//#include "IncomingAsymptoteStateRepresentation.h"
#include "EMTG_enums.h"

namespace EMTG
{
    namespace Astrodynamics
    {

        StateRepresentationBase* CreateStateRepresentation(StateRepresentation state_representation, const double& mu)
		{
			if (state_representation == StateRepresentation::Cartesian)
			{
				return new CartesianStateRepresentation(mu);
			}
			else if (state_representation == StateRepresentation::SphericalAZFPA)
			{
				return new SphericalAZFPAStateRepresentation(mu);
			}
            else if (state_representation == StateRepresentation::SphericalRADEC)
            {
                return new SphericalRADECStateRepresentation(mu);
            }
            else if (state_representation == StateRepresentation::COE)
            {
                return new COEStateRepresentation(mu);
            }
            else if (state_representation == StateRepresentation::MEE)
            {
                return new MEEStateRepresentation(mu);
            }
            else if (state_representation == StateRepresentation::IncomingBplane)
            {
                return new IncomingBplaneStateRepresentation(mu);
            }
            else if (state_representation == StateRepresentation::OutgoingBplane)
            {
                return new OutgoingBplaneStateRepresentation(mu);
            }
            else
            {
                throw std::runtime_error("State representation '" + StateRepresentationStrings[state_representation] + "' has not been implemented.");
            }
		}

    }//end namespace HardwareModels
}//end namespace EMTG