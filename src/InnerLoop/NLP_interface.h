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

//header file for NLP interface abstract base class
#pragma once

#include "problem.h"
#include "EMTG_enums.h"
#include "NLPoptions.h"

namespace EMTG
{
    namespace Solvers
    {
        class NLP_interface
        {

        public:
            //constructor
            NLP_interface();
            NLP_interface(problem* myProblem, 
                const NLPoptions& myOptions);

            //get/set
            inline void setX0_scaled(const std::vector<doubleType>& X0_scaled_in) {this->X0_scaled = X0_scaled_in; }
            inline void setX0_unscaled(const std::vector<doubleType>& X0_unscaled_in) { this->X0_unscaled = X0_unscaled_in; }
            

            inline std::vector<doubleType> getX_scaled() const { return this->X_scaled; }
            inline std::vector<doubleType> getX_unscaled() const { return this->X_unscaled; }
            inline std::vector<doubleType> getF() const { return this->F; }
            inline std::vector<double> getFlowerbounds() const { return this->Flowerbounds; }
            inline std::vector<double> getFupperbounds() const { return this->Flowerbounds; }
            inline std::vector<size_t> getiGfun() const { return this->iGfun; }
            inline std::vector<size_t> getjGvar() const { return this->jGvar; }
            inline std::vector<double> getG() const { return this->G; }
            inline size_t getnX() const { return this->nX; }
            inline size_t getnF() const { return this->nF; }
            inline size_t getnG() const { return this->nG; }
            inline std::vector<doubleType> getX_NLP_incumbent_unscaled() const { return this->X_NLP_incumbent_unscaled; }
            inline std::vector<doubleType> getX_NLP_incumbent_scaled() const { return this->X_NLP_incumbent_scaled; }
            inline std::vector<doubleType> getF_NLP_incumbent() const { return this->F_NLP_incumbent; }
            inline std::vector<double> getG_NLP_incumbent() const { return this->G_NLP_incumbent; }
            inline doubleType getfeasibility_metric() const { return this->feasibility_metric; }
            inline doubleType getfeasibility_metric_NLP_incumbent() const { return this->feasibility_metric_NLP_incumbent; }


            //run NLP
            virtual void run_NLP(const bool& X0_is_scaled = true) = 0;

        protected:
            //scaling functions
            inline void scaleX0()
            {
                for (size_t Xindex = 0; Xindex < this->nX; ++Xindex)
                    X0_scaled[Xindex] = (X0_unscaled[Xindex] - this->myProblem->Xlowerbounds[Xindex]) / this->myProblem->X_scale_factors[Xindex];
            }

            inline void scaleX()
            {
                for (size_t Xindex = 0; Xindex < this->nX; ++Xindex)
                    X_scaled[Xindex] = (X_unscaled[Xindex] - this->myProblem->Xlowerbounds[Xindex]) / this->myProblem->X_scale_factors[Xindex];
            }

            inline void scaleX_NLP_incumbent()
            {
                for (size_t Xindex = 0; Xindex < this->nX; ++Xindex)
                    X_NLP_incumbent_scaled[Xindex] = (X_NLP_incumbent_unscaled[Xindex] - this->myProblem->Xlowerbounds[Xindex]) / this->myProblem->X_scale_factors[Xindex];
            }

            inline void unscaleX0()
            {
                for (size_t Xindex = 0; Xindex < this->nX; ++Xindex)
                    X0_unscaled[Xindex] = X0_scaled[Xindex] * this->myProblem->X_scale_factors[Xindex] + this->myProblem->Xlowerbounds[Xindex];
            }

            inline void unscaleX()
            {
                for (size_t Xindex = 0; Xindex < this->nX; ++Xindex)
                    X_unscaled[Xindex] = X_scaled[Xindex] * this->myProblem->X_scale_factors[Xindex] + this->myProblem->Xlowerbounds[Xindex];
            }

            inline void unscaleX_NLP_incumbent()
            {
                for (size_t Xindex = 0; Xindex < this->nX; ++Xindex)
                    X_NLP_incumbent_unscaled[Xindex] = X_NLP_incumbent_scaled[Xindex] * this->myProblem->X_scale_factors[Xindex] + this->myProblem->Xlowerbounds[Xindex];
            }


            //fields
            problem* myProblem;
            NLPoptions myOptions;
            time_t NLP_start_time;
            size_t movie_frame_count = 0;

            std::vector<doubleType> X0_scaled;
            std::vector<doubleType> X0_unscaled;
            std::vector<doubleType> X_scaled;
            std::vector<doubleType> X_unscaled;
            std::vector<double> Xupperbounds;
            std::vector<double> Xlowerbounds;
            std::vector<double> G;
            std::vector<double> A;
            std::vector<doubleType> F;
            std::vector<double> Fupperbounds;
            std::vector<double> Flowerbounds;
            size_t nX;
            size_t nF;
            size_t nG;
            size_t nA;

            std::vector<size_t> iGfun;
            std::vector<size_t> jGvar;
            std::vector<size_t> iAfun;
            std::vector<size_t> jAvar;

            std::vector<size_t> original_G_indices_of_filament_critical_inequality_constraints;

            //chaperone
            doubleType J_NLP_incumbent;
            doubleType feasibility_metric_NLP_incumbent;
            doubleType feasibility_metric;
            doubleType normalized_feasibility_metric;
            size_t worst_constraint;
            size_t worst_decision_variable;
            doubleType decision_vector_feasibility_metric;
            doubleType distance_from_equality_filament;
            std::vector<doubleType> X_NLP_incumbent_scaled;
            std::vector<doubleType> X_NLP_incumbent_unscaled;
            std::vector<doubleType> F_NLP_incumbent;
            std::vector<double> G_NLP_incumbent;
			
			bool first_feasibility;
        };//end class NLP_interface
}//end namespace Solvers
}//end namespace EMTG