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

#ifndef DETIC_ELEVATION_FROM_GROUND_STATION_CONSTRAINT_H
#define DETIC_ELEVATION_FROM_GROUND_STATION_CONSTRAINT_H

#include "SpecializedBoundaryConstraintBase.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        class BoundaryEventBase;

        namespace SpecializedConstraints
        {

            class BoundaryDeticElevationFromGroundStationConstraint : virtual public SpecializedBoundaryConstraintBase
            {
            public:
                //constructors
                BoundaryDeticElevationFromGroundStationConstraint() {};

                BoundaryDeticElevationFromGroundStationConstraint(const std::string& name,
                    const size_t& journeyIndex,
                    const size_t& phaseIndex,
                    const size_t& stageIndex,
                    Astrodynamics::universe* Universe,
                    HardwareModels::Spacecraft* mySpacecraft,
                    missionoptions* myOptions,
                    BoundaryEventBase* myBoundaryEvent,
                    const std::string& constraintDefinition);

                //destructor
                virtual ~BoundaryDeticElevationFromGroundStationConstraint() {};


                //public methods
                virtual void output(std::ofstream& outputfile);

                virtual void calcbounds();

                //process
                virtual void process_constraint(const std::vector<doubleType>& X,
                    size_t& Xindex,
                    std::vector<doubleType>& F,
                    size_t& Findex,
                    std::vector<double>& G,
                    const bool& needG);

            protected:
                //fields
                doubleType detic_elevation;
                doubleType sin_detic_elevation;
                double ground_station_latitude;
                double ground_station_longitude;
                double ground_station_altitude;
                Astrodynamics::body * ground_station_body;
                std::string ground_station_body_name;
            
                // position of spacecraft w.r.t. central body
                math::Matrix<doubleType> r_cb2sc_ICRF;

                // position of ground station body w.r.t. central body
                math::Matrix<doubleType> r_cb2gsbody_ICRF;
                math::Matrix<double> d_r_cb2gsbody_dt_ICRF;

                // position of ground station w.r.t. central body
                math::Matrix<doubleType> r_cb2gs_ICRF;

                // position of ground station w.r.t. ground station body
                math::Matrix<doubleType> r_gsbody2gs_ICRF;
                math::Matrix<doubleType> r_gsbody2gs_BCF;

                // position of spacecraft w.r.t. ground station
                math::Matrix<doubleType> r_gs2sc_ICRF;
                math::Matrix<doubleType> r_gs2sc_BCF;
                math::Matrix<doubleType> d_r_gs2sc_dt_BCF;

                // zenith vector at ground station
                math::Matrix<doubleType> surface_normal_vec;


                std::vector<size_t> Gindex_constraint_wrt_time_variables;
                std::vector< std::vector<size_t> > Gindex_constraint_wrt_PositionBeforeEvent_variables;//stateIndex, Gindex
                std::vector< std::vector<size_t> > Gindex_constraint_wrt_PositionBeforeEvent_time_variables;//stateIndex, Gindex
                                                   
                std::vector< std::vector<size_t> > dIndex_with_respect_to_PositionBeforeEvent;//stateIndex, dIndex
                std::vector< std::vector<size_t> > dIndex_with_respect_to_PositionBeforeEvent_wrt_Time;//stateIndex, dIndex
            };
        } // end namespace SpecializedConstraints
    } // end namespace BoundaryEvents
} // end namespace EMTG

#endif // DETIC_ELEVATION_FROM_GROUND_STATION_CONSTRAINT_H