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

//event testbed class

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>

#include "mission.h"

//needed for GSAD A5
//GSAD::adouble GSAD::adouble::temp;
//std::vector<size_t>::size_type GSAD::adouble::point = 0;

void missionTestbed(EMTG::missionoptions& options,
                  std::vector< EMTG::Astrodynamics::universe > TheUniverse,
                  EMTG::HardwareModels::Spacecraft& mySpacecraft,
                  EMTG::HardwareModels::LaunchVehicle& myLaunchVehicle,
                  std::mt19937 RNG,
                  std::uniform_real_distribution<> UniformDouble)
{

    //********************************************************************************mission setup

    //create the mission
    EMTG::Mission myMission(options, 
        TheUniverse,
        myLaunchVehicle, 
        mySpacecraft);
    //myMission.scale_to_unit_hypercube();

    //********************************************************************************create decision and constraint vectors

    std::vector<double> newX = myMission.construct_initial_guess();
    for (size_t Xindex = 0; Xindex < myMission.Xdescriptions.size(); ++Xindex)
    {
        myMission.X[Xindex] = newX[Xindex];
        myMission.X[Xindex].setDerivative(Xindex, myMission.X_scale_factors[Xindex]);
    }

    //********************************************************************************evaluate!
    size_t Xindex = 0;
    size_t Findex = 0; 

    for (size_t ntest = 0; ntest < 1; ++ntest)
        myMission.evaluate(myMission.X, myMission.F, myMission.G, true);

    //********************************************************************************print stuff
    //print XF
    myMission.output_problem_bounds_and_descriptions("tests/Mission_XFout.csv");

    //print (sparse) nonlinear Jacobian
    std::ofstream Gout("tests/Mission_Gout.csv", std::ios::trunc);
    Gout << "index, iGfun, jGvar, Name, Analytical, Algorithmic, Absolute error, Relative error, abs(Relative error), an/al, al/an, log10(analytical)" << std::endl;
    Gout.precision(15);
    for (size_t Gindex = 0; Gindex < myMission.Gdescriptions.size(); ++Gindex)
    {
        size_t Findex = myMission.iGfun[Gindex];
        size_t Xindex = myMission.jGvar[Gindex];
        Gout << Gindex << "," << Findex << "," << Xindex << "," << myMission.Gdescriptions[Gindex];
        double Galgorithmic = myMission.F[Findex].getDerivative(Xindex);
        double Ganalytical = myMission.G[Gindex];
        double abserror = Ganalytical - Galgorithmic;
        double relerror = abserror / Galgorithmic;
        double an_al = Ganalytical / Galgorithmic;
        double al_an = Galgorithmic / Ganalytical;
        double log10_analytical = log10(fabs(Ganalytical));
        Gout << "," << Ganalytical << "," << Galgorithmic << "," << abserror << "," << relerror << "," << fabs(relerror) << "," << an_al << "," << al_an << "," << log10_analytical << std::endl;
    }
    Gout.close();

    //print (sparse) linear Jacobian
    std::ofstream Aout("tests/Mission_Aout.csv", std::ios::trunc);
    Aout << "index, iAfun, jAvar, Name, Analytical, Algorithmic, Absolute error, Relative error, abs(Relative error), an/al, al/an" << std::endl;
    for (size_t Aindex = 0; Aindex < myMission.Adescriptions.size(); ++Aindex)
    {
        size_t Findex = myMission.iAfun[Aindex];
        size_t Xindex = myMission.jAvar[Aindex];
        Aout << Aindex << "," << Findex << "," << Xindex << "," << myMission.Adescriptions[Aindex];
        double Aalgorithmic = myMission.F[Findex].getDerivative(Xindex);
        double Aanalytical = myMission.A[Aindex];
        double abserror = Aanalytical - Aalgorithmic;
        double relerror = abserror / Aalgorithmic;
        double an_al = Aanalytical / Aalgorithmic;
        double al_an = Aalgorithmic / Aanalytical;
        Aout << "," << Aanalytical << "," << Aalgorithmic << "," << abserror << "," << relerror << "," << fabs(relerror) << "," << an_al << "," << al_an << std::endl;
    }
    Aout.close();

    //search for missing Jacobian entries
    std::vector< std::tuple<size_t, size_t, double> > missingDerivatives; //Findex, Xindex, magnitude
    for (size_t Findex = 0; Findex < myMission.Fdescriptions.size(); ++Findex)
    {
        std::vector<size_t> derivativeIndices = myMission.F[Findex].getDerivativeIndicies();

        for (size_t Xindex : derivativeIndices)
        {
            bool found = false;
            for (size_t Gindex = 0; Gindex < myMission.Gdescriptions.size(); ++Gindex)
            {
                if (myMission.jGvar[Gindex] == Xindex && myMission.iGfun[Gindex] == Findex)
                {
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                for (size_t Aindex = 0; Aindex < myMission.Adescriptions.size(); ++Aindex)
                {
                    if (myMission.jAvar[Aindex] == Xindex && myMission.iAfun[Aindex] == Findex)
                    {
                        found = true;
                        break;
                    }
                }
            }

            if (!found)
            {
                double magnitude = myMission.F[Findex].getDerivative(Xindex);
                missingDerivatives.push_back({ Findex, Xindex, magnitude });
            }
        }
    }

    std::ofstream Mout("tests/Mission_MissingEntries.csv", std::ios::trunc);
    Mout << "Findex, Xindex, description, magnitude" << std::endl;

    for (std::tuple<size_t, size_t, double> missingDerivative : missingDerivatives)
    {
        size_t Findex = std::get<0>(missingDerivative);
        size_t Xindex = std::get<1>(missingDerivative);
        Mout << Findex << "," << Xindex << ",";
        Mout << "Derivative of " << myMission.Fdescriptions[Findex] << " F[" << Findex << "] with respect to X[" << Xindex << "]: " << myMission.Xdescriptions[Xindex];
        Mout << "," << std::get<2>(missingDerivative);
        Mout << std::endl;
    }
    Mout.close();


    //print
    myMission.Xopt = myMission.X;
    myMission.output("tests/outputcheck.emtg");
    mySpacecraft.getSpacecraftOptions().write_output_file("tests/Spacecraft.emtg_spacecraftopt");
}//end main