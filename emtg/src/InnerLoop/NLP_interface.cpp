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

//NLP interface abstract base class

#include "NLP_interface.h"

namespace EMTG
{
    namespace Solvers
    {
        //constructor
        EMTG::Solvers::NLP_interface::NLP_interface() :
            nX(1),
            nF(1)
        {}

        EMTG::Solvers::NLP_interface::NLP_interface(problem* myProblem_in,
            const NLPoptions& myOptions) :
            myProblem(myProblem_in),
            myOptions(myOptions),
            nX(myProblem->total_number_of_NLP_parameters),
            nF(myProblem->total_number_of_constraints),
            nG(myProblem->Gdescriptions.size()),
            nA(myProblem->Adescriptions.size()),
            J_NLP_incumbent(math::LARGE),
            feasibility_metric_NLP_incumbent(math::LARGE)
        {

            this->X0_scaled = std::vector<doubleType>(nX, 0.0);
            this->X0_unscaled = std::vector<doubleType>(nX, 0.0);
            this->X_scaled = std::vector<doubleType>(nX, 0.0);
            this->X_unscaled = std::vector<doubleType>(nX, 0.0);
            this->F = std::vector<doubleType>(nF, 0.0);
            this->G = std::vector<double>(nG, 0.0);
            this->X_NLP_incumbent_scaled = std::vector<doubleType>(nX, 0.0);
            this->X_NLP_incumbent_unscaled = std::vector<doubleType>(nX, 0.0);
            this->F_NLP_incumbent = std::vector<doubleType>(nF, 0.0);
            this->G_NLP_incumbent = std::vector<double>(nG, 0.0);
            this->Xupperbounds = std::vector<double>(this->nX);
            this->Xlowerbounds = std::vector<double>(this->nX);

            for (size_t Xindex = 0; Xindex < this->nX; ++Xindex)
            {
                this->Xupperbounds[Xindex] = this->myProblem->Xupperbounds[Xindex] / this->myProblem->X_scale_factors[Xindex];
                this->Xlowerbounds[Xindex] = this->myProblem->Xlowerbounds[Xindex] / this->myProblem->X_scale_factors[Xindex];
            }
            this->Fupperbounds = myProblem->Fupperbounds;
            this->Flowerbounds = myProblem->Flowerbounds;
            //if we are in filament finder mode, then we need to construct our own version of F, G, iGfun/jGvar that are not the same as the Problem's
            if (this->myOptions.get_SolverMode() == NLPMode::FilamentFinder)
            {
                //for now, let's use the same linear constraints
                this->A = myProblem->A;
                this->iAfun = myProblem->iAfun;
                this->jAvar = myProblem->jAvar;

                //the filament finding problem is effectively unconstrained
                this->nF = 1;
                this->nG = this->nX;
                this->iGfun.resize(this->nG, 0);
                this->jGvar.resize(this->nG, 1); //we have to make everybody influence the filament otherwise SNOPT gets confused, even if not really true
                this->Flowerbounds.resize(1, 0.0);
                this->Fupperbounds.resize(1, math::LARGE);

                for (size_t Gindex = 0; Gindex < this->nG; ++Gindex)
                {
                    this->iGfun[Gindex] = 0;
                    this->jGvar[Gindex] = Gindex;
                }

                //add critical inequality constraints
                this->myProblem->locate_filament_critical_inequality_constraints();
                for (size_t Findex = 0; Findex < this->myProblem->F_indices_of_filament_critical_inequality_constraints.size(); ++Findex)
                {
                    ++this->nF;
                    this->Flowerbounds.push_back(this->myProblem->Flowerbounds[this->myProblem->F_indices_of_filament_critical_inequality_constraints[Findex]]);
                    this->Fupperbounds.push_back(this->myProblem->Fupperbounds[this->myProblem->F_indices_of_filament_critical_inequality_constraints[Findex]]);

                    for (size_t Gindex = 0; Gindex < this->myProblem->iGfun.size(); ++Gindex)
                    {
                        if (this->myProblem->iGfun[Gindex] == this->myProblem->F_indices_of_filament_critical_inequality_constraints[Findex])
                        {
                            this->original_G_indices_of_filament_critical_inequality_constraints.push_back(Gindex);
                            ++this->nG;
                            this->iGfun.push_back(this->nF - 1);
                            this->jGvar.push_back(this->myProblem->jGvar[Gindex]);
                        }
                    }
                }

                //size F and G appropriately
                this->F.resize(this->nF, 1.0e+100);
                this->G.resize(this->nG, 0.0);
            }
            else
            {
                this->A = myProblem->A;
                this->iAfun = myProblem->iAfun;
                this->jAvar = myProblem->jAvar;
                this->iGfun = myProblem->iGfun;
                this->jGvar = myProblem->jGvar;
            }
        }
    }//end namespace Solvers
}//end namespace EMTG