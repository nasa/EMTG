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

// Acceleration model accuracy verification against MIRAGE data

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include "boost/algorithm/string/split.hpp"                                    
#include "boost/algorithm/string.hpp"
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"

#include "doubleType.h"
#include "EMTG_Matrix.h"
#include "MIRAGE_verification.h"
#include "mission.h"

#include "SpacecraftAccelerationModel.h"

#include "SpiceUsr.h"

void MIRAGE_verification(EMTG::missionoptions& options,
	                     std::vector< EMTG::Astrodynamics::universe > TheUniverse,
	                     EMTG::HardwareModels::Spacecraft& mySpacecraft,
	                     EMTG::HardwareModels::LaunchVehicle& myLaunchVehicle)
{

	//********************************************************************************mission setup
	options.outputfile = "tests/MIRAGE_error_comparison.emtg";

	//create the mission
	EMTG::Mission myMission(options,
		                    TheUniverse,
		                    myLaunchVehicle,
		                    mySpacecraft);


	std::string MIRAGE_input_file = "accel_all.txt";

	// parse the MIRAGE data file 
	std::ifstream inputfile(MIRAGE_input_file.c_str());

	if (!inputfile.is_open())
	{
		throw std::runtime_error("Cannot open MIRAGE input file. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
		return;
	}


	std::vector<std::string> line_data;
	std::vector<std::string> MIRAGE_headers;
	std::string linebuffer;
	/*std::copy(std::istream_iterator<std::string>(inputfile),
		      std::istream_iterator<std::string>(),
		      std::back_inserter(MIRAGE_data));
	*/


	std::vector<MIRAGE_data_row> MIRAGE_data;
	size_t data_row_index = 0;
	while (inputfile.is_open() && std::getline(inputfile, linebuffer))
	{
		if (linebuffer.length() != 0)
		{
			boost::algorithm::split(line_data, linebuffer, boost::is_any_of("\t"), boost::token_compress_on);

			size_t data_col_index = 0;
			size_t header_index = 0;
			// we are reading headers
			if (line_data[0] == "TIME")
			{
				MIRAGE_headers = line_data;
				continue;
			}
			else // we are reading data
			{
				// create a new empty row of data
				MIRAGE_data.push_back(MIRAGE_data_row ());

				// all other MIRAGE data is a 3-vector
				EMTG::math::Matrix<double> temp_3_vec(3, 1, 0.0);

				// sweep across MIRAGE data columns
				while (data_col_index < line_data.size())
				{

                    // take care of the singletons TIME and MASS
                    if (MIRAGE_headers[header_index] == "TIME")
                    {
                        MIRAGE_data[data_row_index].epoch = line_data[data_col_index++];
                    }
                    else if (MIRAGE_headers[header_index] == "MASS")
                    {
                        MIRAGE_data[data_row_index].mass = atof(line_data[data_col_index++].c_str());
                    }
                    else
                    {
                        temp_3_vec(0) = atof(line_data[data_col_index++].c_str());
                        temp_3_vec(1) = atof(line_data[data_col_index++].c_str());
                        temp_3_vec(2) = atof(line_data[data_col_index++].c_str());
                        MIRAGE_data[data_row_index].acceleration_components.emplace(MIRAGE_headers[header_index], temp_3_vec);
                    }
					++header_index;
				}
			}
		}
		else
		{
			continue;
		}

        /*
		std::cout << MIRAGE_data[data_row_index].time << " " << MIRAGE_data[data_row_index].acceleration_components["RCP"](0) << " " 
                                                             << MIRAGE_data[data_row_index].acceleration_components["RCP"](1) << " " 
                                                             << MIRAGE_data[data_row_index].acceleration_components["RCP"](2) << " " << std::endl << std::endl;
        */

		// on to the next row of MIRAGE data
		++data_row_index;

	}


	size_t state_vector_size = 10;
	size_t STM_start_index = 9;
	size_t num_STM_rows = 14;
	size_t num_STM_columns = 14;
    std::vector<std::string> stringthing;
	EMTG::Astrodynamics::SpacecraftAccelerationModel test_acceleration_model(&myMission.options,
		                                                                     &myMission.options.Journeys[0],
		                                                                     &TheUniverse[0],
                                                                             &stringthing,
		                                                                     &mySpacecraft,
		                                                                     num_STM_rows);

    test_acceleration_model.setDutyCycle(myMission.options.engine_duty_cycle);


    //test_acceleration_model.setLaunchEpoch(2.0);
    
    
    
    test_acceleration_model.setEpoch(2.0);

    // state is {x, y, z, vx, vy, vz, m, chem fuel, elec. fuel/oxidizer}

    // for each MIRAGE state, compute the EMTG acceleration
    //create the working directory
    try
    {
        boost::filesystem::path p("MIRAGE_verification");
        boost::filesystem::create_directories(p);
    }
    catch (std::exception &e)
    {
        std::cerr << "Error " << e.what() << ": Directory creation failed" << std::endl;
    }
    std::ofstream MIRAGE_stream("MIRAGE_verification/spacecraft_acceleration_model.csv", std::ios::trunc);
    for (MIRAGE_data_row data_row : MIRAGE_data)
    {
        EMTG::math::Matrix<doubleType> state(9, 1, 0.0);
        state(0) = data_row.acceleration_components["RCP"](0);
        state(1) = data_row.acceleration_components["RCP"](1);
        state(2) = data_row.acceleration_components["RCP"](2);
        state(3) = data_row.acceleration_components["VCP"](0);
        state(4) = data_row.acceleration_components["VCP"](1);
        state(5) = data_row.acceleration_components["VCP"](2);
        state(6) = data_row.mass;
        state(7) = 0.0;
        state(8) = 0.0;

        double epoch_seconds_past_J2000;
        str2et_c(data_row.epoch.c_str(), &epoch_seconds_past_J2000);

        test_acceleration_model.setEpoch(epoch_seconds_past_J2000 + 51544.5 * 86400.0);


        test_acceleration_model.computeAcceleration(state, false);

        std::map<std::string, std::tuple<double, EMTG::math::Matrix<doubleType>>> EMTG_gravity_accel_sources = test_acceleration_model.getGravityAccelSources();

        MIRAGE_stream << data_row.epoch << std::endl;
        MIRAGE_stream << "Acceleration source, EMTG value, MIRAGE value, absolute error, relative error" << std::endl;

        MIRAGE_stream << std::setprecision(20);

        // compare 
        if (EMTG_gravity_accel_sources.count("Mercury") && data_row.acceleration_components["MEP"].norm() != 0.0)
        {
            MIRAGE_stream << "Mercury x,"
                << std::get<1>(EMTG_gravity_accel_sources["Mercury"])(0) << "," << data_row.acceleration_components["MEP"](0) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Mercury"])(0) - data_row.acceleration_components["MEP"](0)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Mercury"])(0) - data_row.acceleration_components["MEP"](0)) / data_row.acceleration_components["MEP"](0)) << std::endl;

            MIRAGE_stream << "Mercury y,"
                << std::get<1>(EMTG_gravity_accel_sources["Mercury"])(1) << "," << data_row.acceleration_components["MEP"](1) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Mercury"])(1) - data_row.acceleration_components["MEP"](1)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Mercury"])(1) - data_row.acceleration_components["MEP"](1)) / data_row.acceleration_components["MEP"](1)) << std::endl;

            MIRAGE_stream << "Mercury z,"
                << std::get<1>(EMTG_gravity_accel_sources["Mercury"])(2) << "," << data_row.acceleration_components["MEP"](2) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Mercury"])(2) - data_row.acceleration_components["MEP"](2)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Mercury"])(2) - data_row.acceleration_components["MEP"](2)) / data_row.acceleration_components["MEP"](2)) << std::endl;
        }

        if (EMTG_gravity_accel_sources.count("Venus") && data_row.acceleration_components["VEP"].norm() != 0.0)
        {
            MIRAGE_stream << "Venus x,"
                << std::get<1>(EMTG_gravity_accel_sources["Venus"])(0) << "," << data_row.acceleration_components["VEP"](0) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Venus"])(0) - data_row.acceleration_components["VEP"](0)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Venus"])(0) - data_row.acceleration_components["VEP"](0)) / data_row.acceleration_components["VEP"](0)) << std::endl;

            MIRAGE_stream << "Venus y,"
                << std::get<1>(EMTG_gravity_accel_sources["Venus"])(1) << "," << data_row.acceleration_components["VEP"](1) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Venus"])(1) - data_row.acceleration_components["VEP"](1)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Venus"])(1) - data_row.acceleration_components["VEP"](1)) / data_row.acceleration_components["VEP"](1)) << std::endl;

            MIRAGE_stream << "Venus z,"
                << std::get<1>(EMTG_gravity_accel_sources["Venus"])(2) << "," << data_row.acceleration_components["VEP"](2) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Venus"])(2) - data_row.acceleration_components["VEP"](2)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Venus"])(2) - data_row.acceleration_components["VEP"](2)) / data_row.acceleration_components["VEP"](2)) << std::endl;
        }

        if (EMTG_gravity_accel_sources.count("Earth") && data_row.acceleration_components["EAP"].norm() != 0.0)
        {
            MIRAGE_stream << "Earth x,"
                << std::get<1>(EMTG_gravity_accel_sources["Earth"])(0) << "," << data_row.acceleration_components["EAP"](0) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Earth"])(0) - data_row.acceleration_components["EAP"](0)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Earth"])(0) - data_row.acceleration_components["EAP"](0)) / data_row.acceleration_components["EAP"](0)) << std::endl;

            MIRAGE_stream << "Earth y,"
                << std::get<1>(EMTG_gravity_accel_sources["Earth"])(1) << "," << data_row.acceleration_components["EAP"](1) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Earth"])(1) - data_row.acceleration_components["EAP"](1)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Earth"])(1) - data_row.acceleration_components["EAP"](1)) / data_row.acceleration_components["EAP"](1)) << std::endl;

            MIRAGE_stream << "Earth z,"
                << std::get<1>(EMTG_gravity_accel_sources["Earth"])(2) << "," << data_row.acceleration_components["EAP"](2) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Earth"])(2) - data_row.acceleration_components["EAP"](2)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Earth"])(2) - data_row.acceleration_components["EAP"](2)) / data_row.acceleration_components["EAP"](2)) << std::endl;
        }

        if (EMTG_gravity_accel_sources.count("Mars") && data_row.acceleration_components["MAP"].norm() != 0.0)
        {
            MIRAGE_stream << "Mars x,"
                << std::get<1>(EMTG_gravity_accel_sources["Mars"])(0) << "," << data_row.acceleration_components["MAP"](0) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Mars"])(0) - data_row.acceleration_components["MAP"](0)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Mars"])(0) - data_row.acceleration_components["MAP"](0)) / data_row.acceleration_components["MAP"](0)) << std::endl;

            MIRAGE_stream << "Mars y,"
                << std::get<1>(EMTG_gravity_accel_sources["Mars"])(1) << "," << data_row.acceleration_components["MAP"](1) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Mars"])(1) - data_row.acceleration_components["MAP"](1)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Mars"])(1) - data_row.acceleration_components["MAP"](1)) / data_row.acceleration_components["MAP"](1)) << std::endl;

            MIRAGE_stream << "Mars z,"
                << std::get<1>(EMTG_gravity_accel_sources["Mars"])(2) << "," << data_row.acceleration_components["MAP"](2) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Mars"])(2) - data_row.acceleration_components["MAP"](2)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Mars"])(2) - data_row.acceleration_components["MAP"](2)) / data_row.acceleration_components["MAP"](2)) << std::endl;
        }

        if (EMTG_gravity_accel_sources.count("Jupiter") && data_row.acceleration_components["JUP"].norm() != 0.0)
        {
            MIRAGE_stream << "Jupiter x,"
                << std::get<1>(EMTG_gravity_accel_sources["Jupiter"])(0) << "," << data_row.acceleration_components["JUP"](0) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Jupiter"])(0) - data_row.acceleration_components["JUP"](0)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Jupiter"])(0) - data_row.acceleration_components["JUP"](0)) / data_row.acceleration_components["JUP"](0)) << std::endl;

            MIRAGE_stream << "Jupiter y,"
                << std::get<1>(EMTG_gravity_accel_sources["Jupiter"])(1) << "," << data_row.acceleration_components["JUP"](1) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Jupiter"])(1) - data_row.acceleration_components["JUP"](1)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Jupiter"])(1) - data_row.acceleration_components["JUP"](1)) / data_row.acceleration_components["JUP"](1)) << std::endl;

            MIRAGE_stream << "Jupiter z,"
                << std::get<1>(EMTG_gravity_accel_sources["Jupiter"])(2) << "," << data_row.acceleration_components["JUP"](2) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Jupiter"])(2) - data_row.acceleration_components["JUP"](2)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Jupiter"])(2) - data_row.acceleration_components["JUP"](2)) / data_row.acceleration_components["JUP"](2)) << std::endl;
        }

        if (EMTG_gravity_accel_sources.count("Saturn") && data_row.acceleration_components["SAP"].norm() != 0.0)
        {
            MIRAGE_stream << "Saturn x,"
                << std::get<1>(EMTG_gravity_accel_sources["Saturn"])(0) << "," << data_row.acceleration_components["SAP"](0) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Saturn"])(0) - data_row.acceleration_components["SAP"](0)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Saturn"])(0) - data_row.acceleration_components["SAP"](0)) / data_row.acceleration_components["SAP"](0)) << std::endl;

            MIRAGE_stream << "Saturn y,"
                << std::get<1>(EMTG_gravity_accel_sources["Saturn"])(1) << "," << data_row.acceleration_components["SAP"](1) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Saturn"])(1) - data_row.acceleration_components["SAP"](1)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Saturn"])(1) - data_row.acceleration_components["SAP"](1)) / data_row.acceleration_components["SAP"](1)) << std::endl;

            MIRAGE_stream << "Saturn z,"
                << std::get<1>(EMTG_gravity_accel_sources["Saturn"])(2) << "," << data_row.acceleration_components["SAP"](2) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Saturn"])(2) - data_row.acceleration_components["SAP"](2)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Saturn"])(2) - data_row.acceleration_components["SAP"](2)) / data_row.acceleration_components["SAP"](2)) << std::endl;
        }

        if (EMTG_gravity_accel_sources.count("Uranus") && data_row.acceleration_components["URP"].norm() != 0.0)
        {
            MIRAGE_stream << "Uranus x,"
                << std::get<1>(EMTG_gravity_accel_sources["Uranus"])(0) << "," << data_row.acceleration_components["URP"](0) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Uranus"])(0) - data_row.acceleration_components["URP"](0)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Uranus"])(0) - data_row.acceleration_components["URP"](0)) / data_row.acceleration_components["URP"](0)) << std::endl;

            MIRAGE_stream << "Uranus y,"
                << std::get<1>(EMTG_gravity_accel_sources["Uranus"])(1) << "," << data_row.acceleration_components["URP"](1) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Uranus"])(1) - data_row.acceleration_components["URP"](1)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Uranus"])(1) - data_row.acceleration_components["URP"](1)) / data_row.acceleration_components["URP"](1)) << std::endl;

            MIRAGE_stream << "Uranus z,"
                << std::get<1>(EMTG_gravity_accel_sources["Uranus"])(2) << "," << data_row.acceleration_components["URP"](2) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Uranus"])(2) - data_row.acceleration_components["URP"](2)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Uranus"])(2) - data_row.acceleration_components["URP"](2)) / data_row.acceleration_components["URP"](2)) << std::endl;
        }

        if (EMTG_gravity_accel_sources.count("Neptune") && data_row.acceleration_components["NEP"].norm() != 0.0)
        {
            MIRAGE_stream << "Neptune x,"
                << std::get<1>(EMTG_gravity_accel_sources["Neptune"])(0) << "," << data_row.acceleration_components["NEP"](0) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Neptune"])(0) - data_row.acceleration_components["NEP"](0)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Neptune"])(0) - data_row.acceleration_components["NEP"](0)) / data_row.acceleration_components["NEP"](0)) << std::endl;

            MIRAGE_stream << "Neptune y,"
                << std::get<1>(EMTG_gravity_accel_sources["Neptune"])(1) << "," << data_row.acceleration_components["NEP"](1) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Neptune"])(1) - data_row.acceleration_components["NEP"](1)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Neptune"])(1) - data_row.acceleration_components["NEP"](1)) / data_row.acceleration_components["NEP"](1)) << std::endl;

            MIRAGE_stream << "Neptune z,"
                << std::get<1>(EMTG_gravity_accel_sources["Neptune"])(2) << "," << data_row.acceleration_components["NEP"](2) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Neptune"])(2) - data_row.acceleration_components["NEP"](2)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Neptune"])(2) - data_row.acceleration_components["NEP"](2)) / data_row.acceleration_components["NEP"](2)) << std::endl;
        }

        if (EMTG_gravity_accel_sources.count("Pluto") && data_row.acceleration_components["PLP"].norm() != 0.0)
        {
            MIRAGE_stream << "Pluto x,"
                << std::get<1>(EMTG_gravity_accel_sources["Pluto"])(0) << "," << data_row.acceleration_components["PLP"](0) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Pluto"])(0) - data_row.acceleration_components["PLP"](0)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Pluto"])(0) - data_row.acceleration_components["PLP"](0)) / data_row.acceleration_components["PLP"](0)) << std::endl;

            MIRAGE_stream << "Pluto y,"
                << std::get<1>(EMTG_gravity_accel_sources["Pluto"])(1) << "," << data_row.acceleration_components["PLP"](1) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Pluto"])(1) - data_row.acceleration_components["PLP"](1)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Pluto"])(1) - data_row.acceleration_components["PLP"](1)) / data_row.acceleration_components["PLP"](1)) << std::endl;

            MIRAGE_stream << "Pluto z,"
                << std::get<1>(EMTG_gravity_accel_sources["Pluto"])(2) << "," << data_row.acceleration_components["PLP"](2) << ","
                << fabs(std::get<1>(EMTG_gravity_accel_sources["Pluto"])(2) - data_row.acceleration_components["PLP"](2)) << ","
                << fabs((std::get<1>(EMTG_gravity_accel_sources["Pluto"])(2) - data_row.acceleration_components["PLP"](2)) / data_row.acceleration_components["PLP"](2)) << std::endl;
        }

        if (data_row.acceleration_components["SRP"].norm() != 0.0)
        {

            EMTG::math::Matrix<doubleType> EMTG_SRP_vec = test_acceleration_model.getSRPAccelerationVec();

            MIRAGE_stream << "SRP x,"
                << EMTG_SRP_vec(0) << "," << data_row.acceleration_components["SRP"](0) << ","
                << fabs(EMTG_SRP_vec(0) - data_row.acceleration_components["SRP"](0)) << ","
                << fabs((EMTG_SRP_vec(0) - data_row.acceleration_components["SRP"](0)) / data_row.acceleration_components["SRP"](0)) << std::endl;

            MIRAGE_stream << "SRP y,"
                << EMTG_SRP_vec(1) << "," << data_row.acceleration_components["SRP"](1) << ","
                << fabs(EMTG_SRP_vec(1) - data_row.acceleration_components["SRP"](1)) << ","
                << fabs((EMTG_SRP_vec(1) - data_row.acceleration_components["SRP"](1)) / data_row.acceleration_components["SRP"](1)) << std::endl;

            MIRAGE_stream << "SRP z,"
                << EMTG_SRP_vec(2) << "," << data_row.acceleration_components["SRP"](2) << ","
                << fabs(EMTG_SRP_vec(2) - data_row.acceleration_components["SRP"](2)) << ","
                << fabs((EMTG_SRP_vec(2) - data_row.acceleration_components["SRP"](2)) / data_row.acceleration_components["SRP"](2)) << std::endl;
        }


        MIRAGE_stream << std::endl << std::endl;
    }
    MIRAGE_stream.close();
    
}//end main