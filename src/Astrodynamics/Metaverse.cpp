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

#include "Metaverse.h"

#include <sstream>

namespace EMTG
{
    namespace Astrodynamics
    {
        Metaverse::Metaverse(missionoptions& options)
        {
            this->initialize(options);
        }//end constructor

        void Metaverse::initialize(missionoptions& options)
        {
            this->myOptions = &options;

            //if SplineEphem is enabled, create an empty SplineEphem universe
#ifdef SPLINE_EPHEM
            std::vector< std::tuple<int, int, int, double> > SplineUniverse_keyList;

            this->SplineUniverse = SplineEphem::universe(SplineUniverse_keyList);
#endif

            //create a vector of universes for each journey
            this->TheUniverse.clear();
            this->myOptions->TU = 0;
            for (int j = 0; j < this->myOptions->number_of_journeys; ++j)
            {
#ifdef SPLINE_EPHEM
                TheUniverse.push_back(EMTG::Astrodynamics::universe(j, this->myOptions->universe_folder + "//" + this->myOptions->Journeys[j].journey_central_body + ".emtg_universe", options, &this->SplineUniverse));
#else
                TheUniverse.push_back(EMTG::Astrodynamics::universe(j, this->myOptions->universe_folder + "//" + this->myOptions->Journeys[j].journey_central_body + ".emtg_universe", options));
#endif
                std::stringstream universenamestream;

                universenamestream << this->myOptions->Journeys[j].journey_central_body + "_Journey_" << j << ".universe_output";

                if (TheUniverse[j].TU > this->myOptions->TU)
                    this->myOptions->TU = TheUniverse[j].TU;
            }

            for (int j = 0; j < this->myOptions->number_of_journeys; ++j)
            {
                if (j > 0)
                {
                    TheUniverse[j - 1].set_nextUniverse(TheUniverse[j]);
                }
            }

            //now that we have a Universe vector, we can use it to populate the SplineEphem::universe
            //add every body that will we used in the mission to the SplineUniverse
#ifdef SPLINE_EPHEM
            SplineUniverse_keyList.clear();
            //if (this->myOptions->ephemeris_source == 2)
            {
                for (size_t j = 0; j < this->myOptions->number_of_journeys; ++j)
                {
                    std::vector<int> body_index_array;
                    if (this->myOptions->Journeys[j].destination_list[0] > 0)
                        body_index_array.push_back(this->myOptions->Journeys[j].destination_list[0] - 1);
                    if (this->myOptions->Journeys[j].destination_list[1] > 0)
                        body_index_array.push_back(this->myOptions->Journeys[j].destination_list[1] - 1);
                    for (size_t b = 0; b < this->myOptions->Journeys[j].sequence_input.size(); ++b)
                        if (this->myOptions->Journeys[j].sequence_input[b] > 0)
                            body_index_array.push_back(this->myOptions->Journeys[j].sequence_input[b] - 1);
                    for (size_t b = 0; b < TheUniverse[j].perturbation_menu.size(); ++b)
                        body_index_array.push_back(TheUniverse[j].perturbation_menu[b]);

                    for (size_t b = 0; b < body_index_array.size(); ++b)
                    {
                        //do we already have this body?
                        bool body_in_keylist = false;
                        for (size_t k = 0; k < SplineUniverse_keyList.size(); ++k)
                        {
                            if (std::get<0>(SplineUniverse_keyList[k]) == TheUniverse[j].bodies[body_index_array[b]].spice_ID
                                && std::get<1>(SplineUniverse_keyList[k]) == TheUniverse[j].central_body_SPICE_ID)
                            {
                                body_in_keylist = true;
                                break;
                            }
                        }

                        if (!body_in_keylist && body_index_array[b] >= 0)
                        {
                            SplineUniverse_keyList.push_back(std::make_tuple(
                                TheUniverse[j].bodies[body_index_array[b]].spice_ID,
                                TheUniverse[j].central_body_SPICE_ID,
                                360.0,
                                TheUniverse[j].mu));
                        }
                    }//end loop over bodies in the universe

                     //is this universe's central body the sun? If not, let's add this body with respect to the sun. Let's add extra ephemeris points, too.
                    if (!(TheUniverse[j].central_body_SPICE_ID == 10))
                    {
                        SplineUniverse_keyList.push_back(std::make_tuple(
                            TheUniverse[j].central_body_SPICE_ID,
                            10,
                            1000,
                            1.32712440018e+11));
                    }
                }
                this->SplineUniverse.reinitialize(SplineUniverse_keyList,
                    this->myOptions->launch_window_open_date + this->myOptions->Journeys.front().journey_wait_time_bounds[0] - 10.0 * 86400.0,
                    100000.0 * 86400.0);
            }
#endif
        }//end initialize
    }//close namespace Astrodynamics
}//close namespace EMTG