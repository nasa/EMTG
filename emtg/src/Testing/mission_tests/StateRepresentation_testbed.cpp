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

//frame testbed class

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <boost/ptr_container/ptr_vector.hpp>

#include "frame.h"
#include "EMTG_math.h"
#include "EMTG_Matrix.h"
#include "EMTG_Tensor.h"

#include "StateRepresentation_testbed.h"

#include "SphericalAZFPAStateRepresentation.h"
#include "SphericalRADECStateRepresentation.h"
#include "COEStateRepresentation.h"
#include "CartesianStateRepresentation.h"
#include "MEEStateRepresentation.h"
#include "IncomingBplaneStateRepresentation.h"
#include "OutgoingBplaneStateRepresentation.h"
#include "StateRepresentationFactory.h"

void StateRepresentation_testbed(EMTG::missionoptions& options,
    std::vector< EMTG::Astrodynamics::universe > TheUniverse,
    std::mt19937 RNG,
    std::uniform_real_distribution<> UniformDouble)
{
    size_t GSADindex = 0;

    //make six random entries
    EMTG::math::Matrix<doubleType> inputState(6, 1, 0.0);
    for (size_t i = 0; i < 6; ++i)
    {
        inputState(i).setValue(UniformDouble(RNG));
        inputState(i).setDerivative(GSADindex++, 1.0);
    }

    double mu = 1.0;


    std::ofstream outputfile("tests/StateRepresentation_output.csv", std::ios::trunc);

    //test all of the representations
    EMTG::Astrodynamics::StateRepresentationBase* cartesianStateRepresentation = EMTG::Astrodynamics::CreateStateRepresentation(EMTG::StateRepresentation::Cartesian, mu);
    testRepresentation(cartesianStateRepresentation, inputState, outputfile);
    delete cartesianStateRepresentation;

    EMTG::Astrodynamics::StateRepresentationBase* SphericalRADECStateRepresentation = EMTG::Astrodynamics::CreateStateRepresentation(EMTG::StateRepresentation::SphericalRADEC, mu);
    testRepresentation(SphericalRADECStateRepresentation, inputState, outputfile);
    delete SphericalRADECStateRepresentation;

    EMTG::Astrodynamics::StateRepresentationBase* SphericalAZFPAStateRepresentation = EMTG::Astrodynamics::CreateStateRepresentation(EMTG::StateRepresentation::SphericalAZFPA, mu);
    testRepresentation(SphericalAZFPAStateRepresentation, inputState, outputfile);
    delete SphericalAZFPAStateRepresentation;

    EMTG::Astrodynamics::StateRepresentationBase* COEStateRepresentation = EMTG::Astrodynamics::CreateStateRepresentation(EMTG::StateRepresentation::COE, mu);
    testRepresentation(COEStateRepresentation, inputState, outputfile);
    delete COEStateRepresentation;

    EMTG::Astrodynamics::StateRepresentationBase* MEEStateRepresentation = EMTG::Astrodynamics::CreateStateRepresentation(EMTG::StateRepresentation::MEE, mu);
    testRepresentation(MEEStateRepresentation, inputState, outputfile);
    delete MEEStateRepresentation;

    EMTG::math::Matrix<doubleType> noble_test_values(6, 1, { 7000.0, 1000.0, -5000.0 , -1.0, 22.0, 7.0 });
    for (size_t i = 0; i < 6; ++i)
    {
        noble_test_values(i).setDerivative(i, 1.0);
    }

    EMTG::Astrodynamics::StateRepresentationBase* IncomingBplaneStateRepresentation = EMTG::Astrodynamics::CreateStateRepresentation(EMTG::StateRepresentation::IncomingBplane, 3.986e+5);
    testRepresentation(IncomingBplaneStateRepresentation, noble_test_values, outputfile);
    delete IncomingBplaneStateRepresentation;

    noble_test_values.assign_all({ -1389.182993860459, -6830.153393624416, -1183.630020141228, 10.080073981938, -3.269814704898709, 1.400774599839738 });
    for (size_t i = 0; i < 6; ++i)
    {
        noble_test_values(i).setDerivative(i, 1.0);
    }

    EMTG::Astrodynamics::StateRepresentationBase* OutgoingBplaneStateRepresentation = EMTG::Astrodynamics::CreateStateRepresentation(EMTG::StateRepresentation::OutgoingBplane, 3.986e+5);
    testRepresentation(OutgoingBplaneStateRepresentation, noble_test_values, outputfile);
    delete OutgoingBplaneStateRepresentation;

    outputfile.close();
}


void testRepresentation(EMTG::Astrodynamics::StateRepresentationBase* myStateRepresentation,
    EMTG::math::Matrix<doubleType>& inputState,
    std::ofstream& outputfile)
{
    std::vector<std::string> CartesianStateNames({ "x", "y", "z", "vx", "vy", "vz" });
    std::vector<std::string> RepresentationStateNames = myStateRepresentation->getStateNames();
    

    try
    {

        //test the transformation from cartesian to the state representation

        EMTG::math::Matrix<doubleType> outputState2 = myStateRepresentation->convertFromCartesianToRepresentation(inputState, true);
        EMTG::math::Matrix<doubleType> TransformationMatrix = myStateRepresentation->getCartesianToRepresentationTransitionMatrix();



        outputfile << "testing derivatives of " << "Cartesian" << " to " << myStateRepresentation->getName() << " transformation." << std::endl;
        outputfile << "name, analytical, algorithmic, absolute error, relative error, fabs(relerror), an/al, al/an" << std::endl;
        for (size_t i : {0, 1, 2, 3, 4, 5})
        {
            for (size_t j : {0, 1, 2, 3, 4, 5})
            {
                double analytical = TransformationMatrix(i, j).getValue();
                double algorithmic = outputState2(i).getDerivative(j);
                outputfile << "d_" << RepresentationStateNames[i] << "_" << myStateRepresentation->getName() << "_d_" << "cartesian" << "_" << CartesianStateNames[j] << ","
                    << analytical << ","
                    << algorithmic << ","
                    << analytical - algorithmic << ","
                    << (analytical - algorithmic) / algorithmic << ","
                    << fabs((analytical - algorithmic) / algorithmic) << ","
                    << analytical / algorithmic << ","
                    << algorithmic / analytical << std::endl;
            }
        }

        //test the transformation from the state representation to cartesian
        EMTG::math::Matrix<doubleType> input2(6, 1, 0.0);
        for (size_t i : {0, 1, 2, 3, 4, 5})
        {
            input2(i).setValue(outputState2(i).getValue());
            input2(i).setDerivative(i, 1.0);
        }
        EMTG::math::Matrix<doubleType> outputState = myStateRepresentation->convertFromRepresentationToCartesian(input2, true);
        TransformationMatrix = myStateRepresentation->getRepresentationToCartesianTransitionMatrix();



        outputfile << "testing derivatives of " << myStateRepresentation->getName() << " to Cartesian transformation." << std::endl;
        outputfile << "name, analytical, algorithmic, absolute error, relative error, fabs(relerror), an/al, al/an" << std::endl;
        for (size_t i : {0, 1, 2, 3, 4, 5})
        {
            for (size_t j : {0, 1, 2, 3, 4, 5})
            {
                double analytical = TransformationMatrix(i, j).getValue();
                double algorithmic = outputState(i).getDerivative(j);
                outputfile << "d_" << CartesianStateNames[i] << "_" << "cartesian" << "_d_" << RepresentationStateNames[j] << "_" << myStateRepresentation->getName() << ","
                    << analytical << ","
                    << algorithmic << ","
                    << analytical - algorithmic << ","
                    << (analytical - algorithmic) / algorithmic << ","
                    << fabs((analytical - algorithmic) / algorithmic) << ","
                    << analytical / algorithmic << ","
                    << algorithmic / analytical << std::endl;
            }
        }

        //make sure that the output state is reversible
        std::cout << myStateRepresentation->getName() << " -> " << "cartesian" << " -> " << myStateRepresentation->getName() << std::endl;
        //EMTG::math::Matrix<doubleType> checkState = myStateRepresentation->convertFromCartesianToRepresentation(outputState);
        (inputState - outputState).element_divide(inputState).print_to_screen();

        outputfile << std::endl;
    }
    catch (std::exception& error)
    {
        std::cout << error.what() << std::endl;
    }
}
//end frameTestbed()