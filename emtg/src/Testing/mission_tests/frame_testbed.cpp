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
#include <vector>

#include "frame.h"
#include "EMTG_math.h"
#include "EMTG_Matrix.h"
#include "EMTG_Tensor.h"

#include "frame_testbed.h"

void frameTestbed(EMTG::missionoptions& options,
    std::vector< EMTG::Astrodynamics::universe > TheUniverse,
    std::mt19937 RNG,
    std::uniform_real_distribution<> UniformDouble)
{
    //mimic Noble Hatten's Python testbed
    size_t GSADindex = 0;
    std::ofstream framefile("tests/frame_output.csv", std::ios::trunc);
    std::vector<std::string> statenames({ "x", "y", "z", "vx", "vy","vz","t" });

    // set parameters of ellipsoid
    double a = 10.0;
    double b = 8.0;
    double c = 6.0;

    // set Euler angles relating BCF frame to PA frame(3 - 1 - 3)
    // Note on comparisons with STK :
    // The Body_Axes_Body_Fixed system was created as a "Fixed in Axes" type set of axes.
    // The reference axes is principal axis frame.
    // Theta1 corresponds to Euler Angle A
    // Theta2 corresponds to Euler Angle B
    // Theta3 corresponds to Euler Angle C
    // However, all the signs of the angles need to be reversed because the rotations are going the other direction.
    double theta1_0 = -21.0 * EMTG::math::deg2rad;
    double theta2_0 = -9.0 * EMTG::math::deg2rad; 
    double theta3_0 = -14.0 * EMTG::math::deg2rad;

    // assume the theta's can vary linearly in time with constant rates
    double theta1dot_0 = 0.2 * EMTG::math::deg2rad; //deg/day
    double theta2dot_0 = 0.1 * EMTG::math::deg2rad; //deg/day 
    double theta3dot_0 = 0.4  * EMTG::math::deg2rad; //deg/day

    // calculate a point on that ellipsoid
    EMTG::math::Matrix<doubleType> r_pa(3, 1, 0.0);
    r_pa(0) = 0.2*a;
    r_pa(1) = -0.4*b;
    //r_pa[2] = 0.5*c
    //r_pa[0] = np.sqrt(a**2 * (1.0 - r_pa[1] * *2 / b * *2 - r_pa[2] * *2 / c * *2))
    //r_pa[1] = np.sqrt(b**2 * (1.0 - r_pa[0] * *2 / a * *2 - r_pa[2] * *2 / c * *2))
    r_pa(2) = sqrt(c*c * (1.0 - r_pa(0)*r_pa(0) / (a*a) - r_pa(1)*r_pa(1) / (b * b)));

    for (size_t i : {0, 1, 2})
        r_pa(i).setDerivative(GSADindex++, 1.0);

    // set velocity in BCF frame
    EMTG::math::Matrix<doubleType> vBCF(3, 1, std::vector<doubleType>({ 0.2, -0.1, -0.4 }));
    for (size_t i : {0, 1, 2})
        vBCF(i).setDerivative(GSADindex++, 1.0);

    // set rotation of BCF w.r.t.BCI
    double W0 = 18.0 * EMTG::math::deg2rad;// deg
    double Wdot = 1.5 * EMTG::math::deg2rad; // deg / s

    // set the time
    doubleType t = 0.4; // s
    t.setDerivative(GSADindex++, 1.0);

    //extract other angles from the Universe
    double alpha0 = TheUniverse[0].LocalFrame.get_alpha0();
    double alphadot = TheUniverse[0].LocalFrame.get_alphadot();
    double delta0 = TheUniverse[0].LocalFrame.get_delta0();
    double deltadot = TheUniverse[0].LocalFrame.get_deltadot();

    //define some tests
    std::vector< std::tuple<EMTG::ReferenceFrame, EMTG::ReferenceFrame> > testFramePairs;
    testFramePairs.push_back({ EMTG::ReferenceFrame::TrueOfDate_BCF, EMTG::ReferenceFrame::ICRF });
    testFramePairs.push_back({ EMTG::ReferenceFrame::PrincipleAxes, EMTG::ReferenceFrame::TrueOfDate_BCF });
    testFramePairs.push_back({ EMTG::ReferenceFrame::PrincipleAxes, EMTG::ReferenceFrame::ICRF });
    testFramePairs.push_back({ EMTG::ReferenceFrame::Topocentric, EMTG::ReferenceFrame::TrueOfDate_BCF });
    testFramePairs.push_back({ EMTG::ReferenceFrame::Topocentric, EMTG::ReferenceFrame::ICRF });
    testFramePairs.push_back({ EMTG::ReferenceFrame::PrincipleAxes, EMTG::ReferenceFrame::Topocentric });
    testFramePairs.push_back({ EMTG::ReferenceFrame::TrueOfDate_BCF, EMTG::ReferenceFrame::Topocentric });

    EMTG::math::Matrix<doubleType> inputVector = r_pa;
    EMTG::math::Matrix<doubleType> referenceVector = r_pa;

    //run the tests
    for (std::tuple<EMTG::ReferenceFrame, EMTG::ReferenceFrame> FramePair : testFramePairs)
    {
        testRotation(std::get<0>(FramePair),
            std::get<1>(FramePair),
            framefile,
            alpha0,
            alphadot,
            delta0,
            deltadot,
            W0,
            Wdot,
            theta1_0,
            theta1dot_0,
            theta2_0,
            theta2dot_0,
            theta3_0,
            theta3dot_0,
            a,
            b,
            c,
            t,
            inputVector,
            referenceVector);
    }

    framefile.close();
}


void testRotation(EMTG::ReferenceFrame InputFrame,
                  EMTG::ReferenceFrame OutputFrame,
                  std::ofstream& outputfile,
                  const double& alpha0,
                  const double& alphadot,
                  const double& delta0,
                  const double& deltadot,
                  const double& W0,
                  const double& Wdot,
                  const double& theta1_0,
                  const double& theta1dot_0,
                  const double& theta2_0,
                  const double& theta2dot_0,
                  const double& theta3_0,
                  const double& theta3dot_0,
                  const double& semi_axis_a,
                  const double& semi_axis_b,
                  const double& semi_axis_c,
                  const doubleType& ETepoch,
                  EMTG::math::Matrix<doubleType> inputVector,
                  EMTG::math::Matrix<doubleType> referenceVector)
{
    //vector of frame names
    std::vector<std::string> ReferenceFrameNames({ "ICRF", "J2000_BCI", "J2000_BCF", "TrueOfDate_BCI", "TrueOfDate_BCF", "PrincipleAxes", "Topocentric", "Polar", "SAM", "ObjectRefereneced" });

    //vector of state names
    std::vector<std::string> statenames({ "x", "y", "z", "vx", "vy","vz","t" });

    //make a frame
    EMTG::Astrodynamics::frame myFrame(alpha0,
        alphadot,
        delta0,
        deltadot,
        W0,
        Wdot,
        theta1_0,
        theta1dot_0,
        theta2_0,
        theta2dot_0,
        theta3_0,
        theta3dot_0);

    myFrame.set_semi_axis_a(semi_axis_a);
    myFrame.set_semi_axis_b(semi_axis_b);
    myFrame.set_semi_axis_c(semi_axis_c);

    //construct the various matrices
    myFrame.construct_rotation_matrices_J2000();
    myFrame.construct_rotation_matrices(ETepoch, true);
    myFrame.construct_PrincipleAxes_to_BCF_rotation(ETepoch, true);
    myFrame.construct_Topocentric_to_BCF_rotation(ETepoch, referenceVector, true);

    //perform the rotation
    EMTG::math::Matrix<doubleType> dstate_dt(3, 1, 0.0);
    EMTG::math::Matrix<doubleType> dreferenceVector_dt(3, 1, 0.0);
    EMTG::math::Matrix<doubleType> Rotated3Vector;
    EMTG::math::Matrix<doubleType> dRotated3Vector_dt;
    EMTG::math::Matrix<doubleType> dRotated3Vector_dOriginal3Vector;
    EMTG::math::Matrix<doubleType> dRotated3Vector_dreferenceVector;

    myFrame.rotate_frame_to_frame(InputFrame,//OriginalReferenceFrame
        inputVector,//state
        dstate_dt,//dstate_dt - this state is ENCODED, so this is always zero
        referenceVector,//referenceVector,
        dreferenceVector_dt,//dreferenceVector_dt
        OutputFrame,//RotatedReferenceFrame
        Rotated3Vector,//Rstate
        dRotated3Vector_dOriginal3Vector,//dRstate_dstate
        dRotated3Vector_dt,//dRstate_dt
        dRotated3Vector_dreferenceVector,//dRstate_dreferenceVector
        ETepoch,//ReferenceEpochJD
        true);//GenerateDerivatives

    //test the rotation
    outputfile << "testing derivatives of " << ReferenceFrameNames[InputFrame] << " to " << ReferenceFrameNames[OutputFrame] << " rotation." << std::endl;
    outputfile << "name, analytical, algorithmic, absolute error, relative error, fabs(relerror), an/al, al/an" << std::endl;
    for (size_t i : {0, 1, 2})
    {
        for (size_t j : {0, 1, 2})
        {
            double analytical = dRotated3Vector_dOriginal3Vector(i, j).getValue();
            double algorithmic = Rotated3Vector(i).getDerivative(j);
            outputfile << "d_" << statenames[i] << "_" << ReferenceFrameNames[OutputFrame] << "_d_" << statenames[j] << "_" << ReferenceFrameNames[InputFrame] << ","
                << analytical << ","
                << algorithmic << ","
                << analytical - algorithmic << ","
                << (analytical - algorithmic) / algorithmic << ","
                << fabs((analytical - algorithmic) / algorithmic) << ","
                << analytical / algorithmic << ","
                << algorithmic / analytical << std::endl;
        }
        //wrt time
        {
            size_t j = 6;
            double analytical = dRotated3Vector_dt(i).getValue();
            double algorithmic = Rotated3Vector(i).getDerivative(j);
            outputfile << "d_" << statenames[i] << "_" << ReferenceFrameNames[OutputFrame] << "_d_" << statenames[j] << ","
                << analytical << ","
                << algorithmic << ","
                << analytical - algorithmic << ","
                << (analytical - algorithmic) / algorithmic << ","
                << fabs((analytical - algorithmic) / algorithmic) << ","
                << analytical / algorithmic << ","
                << algorithmic / analytical << std::endl;
        }
    }
    outputfile << std::endl;
}
//end frameTestbed()