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

//Edelbaum spiral
//Jacob Englander 10-24-2017

#include "EdelbaumSpiral.h"

namespace EMTG
{
    namespace BoundaryEvents
    {
        EdelbaumSpiral::EdelbaumSpiral(const std::string& name,
            const size_t& journeyIndex,
            const size_t& phaseIndex,
            size_t& stageIndex,
            Astrodynamics::universe* Universe,
            Astrodynamics::body* myBody,
            HardwareModels::Spacecraft* mySpacecraft,
            missionoptions* myOptions)
        {
            this->name = name;
            this->journeyIndex = journeyIndex;
            this->phaseIndex = phaseIndex;
            this->stageIndex = stageIndex;
            this->myUniverse = Universe;
            this->myBody = myBody;
            this->mySpacecraft = mySpacecraft;
            this->myOptions = myOptions;
            this->myJourneyOptions = &this->myOptions->Journeys[this->journeyIndex];

            //size the storage vectors and matrices
            this->spiralRadius = std::vector<double>(this->myOptions->spiral_segments + 1, 0.0);
            this->state_after_spiral = math::Matrix<doubleType>(8, 1, 0.0);

            //create the spiral segments
            for (size_t segmentIndex = 0; segmentIndex < this->myOptions->spiral_segments; ++segmentIndex)
            {
                std::string segmentName = this->name + "Segment" + std::to_string(segmentIndex);


                this->SpiralSegments.push_back(EdelbaumSpiralSegment(segmentName,
                    this->journeyIndex,
                    this->phaseIndex,
                    this->stageIndex,
                    segmentIndex,
                    this->myUniverse,
                    this->myBody,
                    this->mySpacecraft,
                    this->myOptions,
                    this));

                EdelbaumSpiralSegment* previousSegment;
                if (segmentIndex == 0)
                    previousSegment = NULL;
                else
                    previousSegment = &this->SpiralSegments.back();

                this->SpiralSegments.back().set_previous_segment(previousSegment);

                this->SpiralSegments.back().set_radius_before_segment(this->spiralRadius[segmentIndex]);

                this->SpiralSegments.back().set_radius_after_segment(this->spiralRadius[segmentIndex + 1]);
            }
        }//end constructor

        void EdelbaumSpiral::initialize_spiral_segments()
        {
            this->TotalDeltav = fabs(sqrt(this->myBody->mu / this->SpiralEndRadius) - sqrt(this->myBody->mu / this->SpiralStartRadius));
            for (size_t segmentIndex = 0; segmentIndex < this->myOptions->spiral_segments; ++segmentIndex)
            {
                this->SpiralSegments[segmentIndex].set_deltav(this->TotalDeltav / this->myOptions->spiral_segments);
            }
        }//end initialize_spiral_segments()

        void EdelbaumSpiral::setup_calcbounds(
            std::vector<double>* Xupperbounds,
            std::vector<double>* Xlowerbounds,
            std::vector<double>* X_scale_factors,
            std::vector<double>* Fupperbounds,
            std::vector<double>* Flowerbounds,
            std::vector<double>* F_scale_factors,
            std::vector<std::string>* Xdescriptions,
            std::vector<std::string>* Fdescriptions,
            std::vector<size_t>* iGfun,
            std::vector<size_t>* jGvar,
            std::vector<std::string>* Gdescriptions,
            std::vector<size_t>* iAfun,
            std::vector<size_t>* jAvar,
            std::vector<std::string>* Adescriptions,
            std::vector<double>* A)
        {
            this->prefix = this->name + ": ";
            this->Xupperbounds = Xupperbounds;
            this->Xlowerbounds = Xlowerbounds;
            this->X_scale_factors = X_scale_factors;
            this->Xdescriptions = Xdescriptions;
            this->Fupperbounds = Fupperbounds;
            this->Flowerbounds = Flowerbounds;
            this->F_scale_factors = F_scale_factors;
            this->Fdescriptions = Fdescriptions;
            this->iGfun = iGfun;
            this->jGvar = jGvar;
            this->Gdescriptions = Gdescriptions;
            this->iAfun = iAfun;
            this->jAvar = jAvar;
            this->Adescriptions = Adescriptions;
            this->A = A;

            //spiral segments
            for (size_t segmentIndex = 0; segmentIndex < this->myOptions->spiral_segments; ++segmentIndex)
            {
                this->SpiralSegments[segmentIndex].setup_calcbounds(
                    Xupperbounds,
                    Xlowerbounds,
                    X_scale_factors,
                    Fupperbounds,
                    Flowerbounds,
                    F_scale_factors,
                    Xdescriptions,
                    Fdescriptions,
                    iGfun,
                    jGvar,
                    Gdescriptions,
                    iAfun,
                    jAvar,
                    Adescriptions,
                    A);
            }
        }//end setup_calcbounds()


        void EdelbaumSpiral::calcbounds(std::vector<size_t> Xindices_leftEpoch)
        {
            this->Xindices_SpiralEndEpoch = Xindices_leftEpoch;
            for (size_t segmentIndex = 0; segmentIndex < this->myOptions->spiral_segments; ++segmentIndex)
            {
                this->SpiralSegments[segmentIndex].calcbounds();
                this->Xindices_SpiralEndEpoch.push_back(this->SpiralSegments[segmentIndex].get_Xindex_spiralSegmentTime());
            }

            //set up the vectors of derivatives
            //mass has a derivative with respect to final segment's mass after segment
            this->Derivatives_of_StateAfterSpiral.push_back({ this->SpiralSegments.back().get_Xindex_mass_after_segment(), 6, 1.0 });

            //everything else has a derivative with respect to time
            for (size_t Xindex_SpiralEndEpoch : this->Xindices_SpiralEndEpoch)
            {
                //initialize state derivatives to zero
                for (size_t stateIndex : {0, 1, 2, 3, 4, 5})
                {
                    this->Derivatives_of_StateAfterSpiral_wrt_Time.push_back({ Xindex_SpiralEndEpoch, stateIndex, 0.0 });
                }

                //initialize epoch derivative to 1.0
                this->Derivatives_of_StateAfterSpiral_wrt_Time.push_back({ Xindex_SpiralEndEpoch, 7, 1.0 });
            }
        }//end EdelbaumSpiral()

        void EdelbaumSpiral::process_spiral(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            //Step 1: set the initial radius
            this->spiralRadius[0] = this->SpiralStartRadius;

            //Step 2: compute the total delta-v
            double v0 = sqrt(this->myBody->mu / this->SpiralStartRadius);
            double vf = sqrt(this->myBody->mu / this->SpiralEndRadius);
            this->TotalDeltav = fabs(vf - v0);

            //Step 3: compute the spiral segments
            for (EdelbaumSpiralSegment& mySegment : this->SpiralSegments)
            {
                mySegment.set_deltav(this->TotalDeltav / this->myOptions->spiral_segments);

                mySegment.process(X, Xindex, F, Findex, G, needG);
            }

            //Step 4: compute state_after_spiral
            //Step 4.1: when are we?
            this->process_spiral_end_epoch(X, Xindex, F, Findex, G, needG);

            //Step 4.2: where are we?
            doubleType temp_body_state[12];

            this->myBody->locate_body(this->state_after_spiral(7),
                temp_body_state,
                needG,
                *this->myOptions);

            for (size_t stateIndex = 0; stateIndex < 6; ++stateIndex)
                this->state_after_spiral(stateIndex) = temp_body_state[stateIndex];

            //Step 4.3: how heavy are we?
            this->state_after_spiral(6) = this->SpiralSegments.back().get_mass_after_segment();

            //Step 5: derivatives of state_after_spiral
            if (needG)
            {
                //mass with respect to mass (should be) already encoded

                //position and velocity with respect to time
                for (std::tuple<size_t, size_t, double>& derivativeEntry : this->Derivatives_of_StateAfterSpiral_wrt_Time)
                {
                    size_t stateIndex = std::get<1>(derivativeEntry);

                    if (stateIndex < 6)
                        std::get<2>(derivativeEntry) = temp_body_state[6 + stateIndex]_GETVALUE;
                }

                //epoch with respect to time is already encoded
            }
        }//end process_spiral()

        void EdelbaumSpiral::process_spiral_end_epoch(const std::vector<doubleType>& X,
            size_t& Xindex,
            std::vector<doubleType>& F,
            size_t& Findex,
            std::vector<double>& G,
            const bool& needG)
        {
            this->LaunchDate = X[this->Xindices_SpiralEndEpoch[0]];

            this->state_after_spiral(7) = 0.0;
            for (size_t Xindex_SpiralEndEpoch : this->Xindices_SpiralEndEpoch)
            {
                this->state_after_spiral(7) += X[Xindex_SpiralEndEpoch];
            }
        }//end process_spiral_end_epoch

        void EdelbaumSpiral::output(std::ofstream& outputfile,
            const double& launchdate,
            size_t& eventcount)
        {
            for (EdelbaumSpiralSegment& mySegment : this->SpiralSegments)
                mySegment.output(outputfile, eventcount);
        }//end output()
    }//close namespace BoundaryEvents
}//close namespace EMTG