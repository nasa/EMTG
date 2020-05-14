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

//EMTGv9 Problem class
//Jacob Englander 6-22-2017


#pragma once

#include <string>
#include <vector>

#include "missionoptions.h"
#include "SpacecraftOptions.h"
#include "LaunchVehicleOptions.h"

#include "doubleType.h"

#include "ScalatronBase.h"

namespace EMTG
{
    //enum for type of solution to write out
    enum SolutionOutputType { FAILURE, HOP, SUCCESS };

    class problem
    {
    public:
        //constructor
        problem();

        //destructor
        virtual ~problem();

        //methods

        //optimize function
        //return 0 for success, 1 for failure
        bool optimize();

        //bounds calculation function
        virtual void calcbounds() = 0;

        //function to output X and F bounds, descriptions
        void output_problem_bounds_and_descriptions();
        void output_problem_bounds_and_descriptions(std::string filestring);

        //function to determine which constraints are equalities vs inequalities
        void locate_equality_constraints();   

        //function to find constraints that are considered critical to the filament walker
        virtual void locate_filament_critical_inequality_constraints() {};
        
        //function to rescale problem to [0, 1]
        virtual void scale_to_unit_hypercube();

        //function to output the Jacobian sparsity information
        virtual void output_Jacobian_sparsity_information(std::string filestring);

        //function to check feasibility via SNOPT's metric
        void check_feasibility(const std::vector<doubleType>& X,
            const std::vector<doubleType>& F,
            size_t& worst_decision_variable,
            size_t& worst_constraint,
            double& max_constraint_violation,
            double& normalized_max_constraint_violation, 
            double& distance_from_equality_filament,
            double& decision_variable_infeasibility,
            bool calledFromNLP = false);

        //function to construct an initial guess from inputs
        virtual std::vector<double> construct_initial_guess() = 0;
    
        //functions to check inputs (X) and outputs (F) for NaN, inf
        bool check_inputs_for_invalid_entries(const std::vector<doubleType>& X, const bool& PrintAnyway);
        bool check_outputs_for_invalid_entries(const std::vector<doubleType>& X, const std::vector<doubleType>& F, const bool& PrintAnyway);
        bool check_derivatives_for_invalid_entries(const std::vector<doubleType>& X, const std::vector<double>& G, const bool& PrintAnyway);

        //virtual function templates
        virtual void evaluate(const std::vector<doubleType>& X, 
            std::vector<doubleType>& F, 
            std::vector<double>& G, 
            const bool& needG) = 0;

        void evaluate(const std::vector<doubleType>& X);

        virtual void output(const std::string& outputfilename) = 0;

        //method to determine the output file name
        void what_the_heck_am_I_called(const SolutionOutputType& mySolutionType, const size_t number_of_hops = 0);
		
		virtual doubleType getUnscaledObjective() = 0;

        //performance characteristics function
        //used to extract various pieces of mission data for a multi-objective GA
        virtual void extract_outerloop_objective_values(std::vector<double>& objective_functions) {};

        //fields

        //encapsulated options structure
        missionoptions options;

		//encapsulated Scalatron
		Scalatron::ScalatronBase* myScalatron;

        //solver parameters
        std::vector<doubleType> X0; //initial guess, may be supplied or generated at random
        std::vector<doubleType> X; //current decision vector
        std::vector<doubleType> Xopt; //optimal decision vector
        std::vector<doubleType> F; //constraint vector
        std::vector<double> G; //nonlinear Jacobian vector
        std::vector<double> A; //linear Jacobian vector
        double best_cost;
        size_t total_number_of_constraints; //total number of nonlinear constraints
        size_t total_number_of_NLP_parameters; //total number of NLP parameters
        std::vector<double> Xupperbounds, Xlowerbounds, Fupperbounds, Flowerbounds;
        std::vector<double> X_scale_factors;
		std::vector<double> F_scale_factors;
        std::vector<std::string> Xdescriptions, Fdescriptions;
        std::vector<size_t> iAfun, jAvar, iGfun, jGvar;
        std::vector<std::string> Adescriptions, Gdescriptions;
        size_t number_of_solutions;

        std::vector<bool> F_equality_or_inequality; //true for equality, false for inequality, 1 entry for each entry in F[1:]
        std::vector<size_t> F_indices_of_filament_critical_inequality_constraints;

        //vector of synodic periods for any periodic variables
        std::vector<double> synodic_periods;
    };

} /* namespace EMTG */
