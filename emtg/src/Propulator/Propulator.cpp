//Propulator
//contains a vector of PropulatorJourney objects

#include "Propulator.h"

#include "file_utilities.h"
#include "SpacecraftOptionsFactory.h"

#include "SpiceUsr.h"

namespace EMTG
{
    namespace Propulator
    {
        Propulator::Propulator() :
            journeyIndex(0)
        {}//empty default constructor except for initialization list

        Propulator::Propulator(const std::string& options_file_name)
            : Propulator()
        {
            this->initialize(options_file_name);
        }//end constructor

        void Propulator::initialize(const std::string& options_file_name)
        {
            std::cout << "Loading " << options_file_name << std::endl;
            this->myOptions = missionoptions(options_file_name);

            //set up Universe stuff
            {

                //load all ephemeris data if using SPICE
                std::vector<::boost::filesystem::path> SPICE_files_initial;
                std::vector<::boost::filesystem::path> SPICE_files_not_required;
                this->SPICE_files_required.clear();
                std::vector<int> SPICE_bodies_required;
                std::string filestring;
                if (this->myOptions.ephemeris_source >= 1)
                {
                    //load all BSP files
                    EMTG::file_utilities::get_all_files_with_extension(::boost::filesystem::path(this->myOptions.universe_folder + "/ephemeris_files/"), ".bsp", SPICE_files_initial);

                    for (size_t k = 0; k < SPICE_files_initial.size(); ++k)
                    {
                        filestring = this->myOptions.universe_folder + "/ephemeris_files/" + SPICE_files_initial[k].string();
                        furnsh_c(filestring.c_str());
                        std::cout << filestring << std::endl;
                    }

                    //disable quit-on-SPICE-error so that we can see what happens if the leap second and/or frame kernels don't load properly
                    erract_c((SpiceChar*)"SET", 100, (SpiceChar*)"RETURN");

                    //SPICE reference frame kernel
                    std::string leapsecondstring = this->myOptions.universe_folder + "/ephemeris_files/" + this->myOptions.SPICE_leap_seconds_kernel;
                    std::string referenceframestring = this->myOptions.universe_folder + "/ephemeris_files/" + this->myOptions.SPICE_reference_frame_kernel;
                    furnsh_c(leapsecondstring.c_str());
                    furnsh_c(referenceframestring.c_str());

                    //disable SPICE error printing. This is because we can, and will often, go off the edge of an ephemeris file.
                    errprt_c((SpiceChar*)"SET", 100, (SpiceChar*)"NONE");

                    SPICE_files_required = SPICE_files_initial;

                    std::cout << "Completed loading SPICE kernels." << std::endl;
                }
#ifdef SPLINE_EPHEM
                std::vector< std::tuple<int, int, int, double> > SplineUniverse_keyList;

                this->mySplineUniverse = new SplineEphem::universe(SplineUniverse_keyList);
#endif



                for (int j = 0; j < this->myOptions.number_of_journeys; ++j)
                {
#ifdef SPLINE_EPHEM
                    TheUniverse.push_back(EMTG::Astrodynamics::universe(j, this->myOptions.universe_folder + "//" + this->myOptions.Journeys[j].journey_central_body + ".emtg_universe", this->myOptions, this->mySplineUniverse));
#else
                    TheUniverse.push_back(EMTG::Astrodynamics::universe(j, options.universe_folder + "//" + options.Journeys[j].journey_central_body + ".emtg_universe", options));
#endif

                    if (TheUniverse[j].TU > this->myOptions.TU)
                        this->myOptions.TU = TheUniverse[j].TU;
                }

                for (int j = 0; j < this->myOptions.number_of_journeys; ++j)
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
                try
                {
                    //double earliest_possible_epoch = options.launch_window_open_date + options.Journeys.front().wait_time_bounds[0];
                    //double latest_possible_epoch = options.latestPossibleEpoch * 86400.0;

                    size_t number_of_journeys_to_spline = std::min(this->myOptions.number_of_journeys, this->myOptions.stop_after_journey + 1);
                    for (size_t j = 0; j < number_of_journeys_to_spline; ++j)
                    {
                        std::vector<int> body_index_array;

                        //first boundary point
                        if (this->myOptions.Journeys[j].destination_list[0] > 0)
                            body_index_array.push_back(this->myOptions.Journeys[j].destination_list[0] - 1);

                        //last boundary point
                        if (this->myOptions.Journeys[j].destination_list[1] > 0)
                            body_index_array.push_back(this->myOptions.Journeys[j].destination_list[1] - 1);

                        //sequence
                        for (int body : this->myOptions.Journeys[j].sequence)
                            if (body > 0)
                                body_index_array.push_back(body - 1);

                        //perturbation list
                        if (this->myOptions.perturb_thirdbody)
                        {
                            for (size_t b = 0; b < TheUniverse[j].perturbation_menu.size(); ++b)
                                body_index_array.push_back(TheUniverse[j].perturbation_menu[b]);
                        }

                        //distance constraint list
                        for (std::string& constraint : this->myOptions.Journeys[j].PhaseDistanceConstraintDefinitions)
                        {
                            std::vector<std::string> ConstraintDefinitionCell;
                            boost::split(ConstraintDefinitionCell,
                                constraint,
                                boost::is_any_of("_"),
                                boost::token_compress_on);

                            if (boost::to_lower_copy(ConstraintDefinitionCell[1]) != "cb")
                            {
                                int bodyIndex = std::stoi(ConstraintDefinitionCell[1]) - 1;

                                body_index_array.push_back(bodyIndex);
                            }
                        }

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
                                    this->myOptions.SplineEphem_points_per_period,
                                    TheUniverse[j].mu));
                            }
                        }//end loop over bodies in the universe

                        //is this universe's central body the sun? If not, let's add this body with respect to the sun. Let's add extra ephemeris points, too.
                        if (!(TheUniverse[j].central_body_SPICE_ID == 10))
                        {
                            bool body_in_keylist = false;
                            for (size_t k = 0; k < SplineUniverse_keyList.size(); ++k)
                            {
                                if (std::get<0>(SplineUniverse_keyList[k]) == TheUniverse[j].central_body_SPICE_ID
                                    && std::get<1>(SplineUniverse_keyList[k]) == 10)
                                {
                                    body_in_keylist = true;
                                    break;
                                }
                            }

                            if (!body_in_keylist)
                            {
                                SplineUniverse_keyList.push_back(std::make_tuple(
                                    TheUniverse[j].central_body_SPICE_ID,
                                    10,
                                    this->myOptions.SplineEphem_non_central_body_sun_points_per_period,
                                    1.32712440018e+11));
                            }
                        }

                        ////do we need to update the earliest or latest possible epoch?
                        //if (options.Journeys[j].arrival_class == EMTG::BoundaryClass::FreePoint)
                        //{
                        //    earliest_possible_epoch = options.Journeys[j].arrival_elements_reference_epoch < earliest_possible_epoch ? options.Journeys[j].arrival_elements_reference_epoch : earliest_possible_epoch;
                        //    latest_possible_epoch = options.Journeys[j].arrival_elements_reference_epoch > latest_possible_epoch ? options.Journeys[j].arrival_elements_reference_epoch : latest_possible_epoch;
                        //}
                        //if (options.Journeys[j].departure_class == EMTG::BoundaryClass::FreePoint && j == 0)
                        //{
                        //    earliest_possible_epoch = options.Journeys[j].departure_elements_reference_epoch < earliest_possible_epoch ? options.Journeys[j].departure_elements_reference_epoch : earliest_possible_epoch;
                        //    latest_possible_epoch = options.Journeys[j].departure_elements_reference_epoch > latest_possible_epoch ? options.Journeys[j].departure_elements_reference_epoch : latest_possible_epoch;
                        //}
                    }

                    double earliestPossibleEpoch = this->myOptions.earliestPossibleEpoch * 86400.0;
                    double latestPossibleEpoch = this->myOptions.latestPossibleEpoch * 86400.0;

                    if (this->myOptions.SplineEphem_truncate_ephemeris_at_maximum_mission_epoch
                        && latestPossibleEpoch < (this->myOptions.launch_window_open_date + this->myOptions.Journeys.front().wait_time_bounds[1] + this->myOptions.total_flight_time_bounds[1]))
                        latestPossibleEpoch = this->myOptions.launch_window_open_date + this->myOptions.Journeys.front().wait_time_bounds[1] + this->myOptions.total_flight_time_bounds[1] * 86400.0;
                    /*if (earliest_possible_epoch > this->myOptions.earliestPossibleEpoch * 86400.0)
                        earliest_possible_epoch = this->myOptions.earliestPossibleEpoch * 86400.0;*/
                    this->mySplineUniverse->reinitialize(SplineUniverse_keyList,
                        earliestPossibleEpoch - 10.0 * 86400.0,
                        latestPossibleEpoch + 10.0 * 86400.0);
                }
                catch (std::exception &myError)
                {
                    std::cout << "Failure while configuring SplineEphem." << std::endl;
                    std::cout << myError.what() << std::endl;
                    std::cout << "Submit this error message to the EMTG development team, along with your .emtgopt, .emtg_universe file(s), your hardware model files, any relevant ephemeris files, and which branch you are using. This information will allow us to properly help you." << std::endl;
#ifndef BACKGROUND_MODE //macro overrides if statement
                    std::cout << "Press enter to close window." << std::endl;
                    std::cin.ignore();
#endif
                    throw;
                }
#endif
            }//end setting up universe stuff

            //spacecraft
            EMTG::HardwareModels::SpacecraftOptions mySpacecraftOptions = EMTG::HardwareModels::CreateSpacecraftOptions(this->myOptions);
            this->mySpacecraft = EMTG::HardwareModels::Spacecraft(mySpacecraftOptions);


            //set up PropulatorJourneys
            size_t stageIndex = 0;
            this->myPropulatorJourneys.clear();

            for (size_t journeyIndex = 0; journeyIndex < this->myOptions.number_of_journeys; ++journeyIndex)
            {
                this->myPropulatorJourneys.push_back(new PropulatorJourney(journeyIndex,
                    stageIndex,
                    &this->myOptions,
                    &this->TheUniverse[journeyIndex],
                    &this->mySpacecraft));
            }
        }//end initialize()

        Propulator::~Propulator()
        {
            if (this->myOptions.ephemeris_source >= 1)
            {
                for (size_t k = 0; k < this->SPICE_files_required.size(); ++k)
                {
                    std::string filestring = this->myOptions.universe_folder + "ephemeris_files/" + this->SPICE_files_required[k].string();
                    unload_c(filestring.c_str());
                }

                unload_c((this->myOptions.universe_folder + "ephemeris_files/" + this->myOptions.SPICE_leap_seconds_kernel).c_str());
                unload_c((this->myOptions.universe_folder + "ephemeris_files/" + this->myOptions.SPICE_reference_frame_kernel).c_str());
            }
#ifdef SPLINE_EPHEM
            delete this->mySplineUniverse;
#endif
        }//end destructor

        void Propulator::propulate(const math::Matrix<double>& InputState,
            const math::Matrix<double>& ControlVector,
            const double& DutyCycle,
            const double& propagationTime,
            math::Matrix<double>& OutputState,
            const bool& needSTM)
        {
            this->myPropulatorJourneys[this->journeyIndex].propulate(InputState, 
                ControlVector,
                DutyCycle, 
                propagationTime, 
                OutputState,
                needSTM);
        }//end propulate()

        void Propulator::propulate(const math::Matrix<double>& InputState,
            const double& propagationTime,
            math::Matrix<double>& OutputState,
            const bool& needSTM)
        {
            this->myPropulatorJourneys[this->journeyIndex].propulate(InputState,
                propagationTime,
                OutputState,
                needSTM);
        }//end propulate()



#ifdef PROPULATOR_PYTHON_INTERFACE
        boost::python::list Propulator::propulateThrust(const boost::python::list& InputState,
            const boost::python::list& ControlVector,
            const double& DutyCycle,
            const double& propagationTime,
            const bool& needSTM)
        {
            //Step 1: wrap inputs
            size_t nstates = boost::python::len(InputState);
            size_t ncontrols = boost::python::len(ControlVector);
            EMTG::math::Matrix<double> inputStateLine(10, 1, 0.0);
            EMTG::math::Matrix<double> inputControlLine(3, 1, 0.0);

            //position, velocity, and mass
            for (size_t stateIndex = 0; stateIndex < 7; ++stateIndex)
                inputStateLine(stateIndex) = boost::python::extract<double>(InputState[stateIndex]);

            //epoch
            double JDepoch = boost::python::extract<double>(InputState[7]);
            inputStateLine(7) = (JDepoch - 2400000.5) * 86400.0;

            //leave the tanks zeroed

            //control
            for (size_t controlIndex = 0; controlIndex < ncontrols; ++controlIndex)
                inputControlLine(controlIndex) = boost::python::extract<double>(ControlVector[controlIndex]);

            //Step 2: call propulator
            EMTG::math::Matrix<double> outputStateLine(10, 1, 0.0);
            this->propulate(inputStateLine, 
                inputControlLine,
                DutyCycle, 
                propagationTime,
                outputStateLine,
                needSTM);

            //Step 3: wrap outputs
            boost::python::list outputState;
            for (size_t stateIndex = 0; stateIndex < nstates; ++stateIndex)
                outputState.append(outputStateLine(stateIndex));

            //Step 4: return the output state
            return outputState;
        }//end python propulateThrust()

        boost::python::list Propulator::propulateCoast(const boost::python::list& InputState,
            const double& propagationTime,
            const bool& needSTM)
        {
            //Step 1: wrap inputs
            size_t nstates = boost::python::len(InputState);
            EMTG::math::Matrix<double> inputStateLine(10, 1, 0.0);

            //position, velocity, and mass
            for (size_t stateIndex = 0; stateIndex < 7; ++stateIndex)
                inputStateLine(stateIndex) = boost::python::extract<double>(InputState[stateIndex]);

            //epoch
            double JDepoch = boost::python::extract<double>(InputState[7]);
            inputStateLine(7) = (JDepoch - 2400000.5) * 86400.0;

            //leave the tanks zeroed

            //Step 2: call propulator
            EMTG::math::Matrix<double> outputStateLine(10, 1, 0.0);
            this->propulate(inputStateLine,
                propagationTime,
                outputStateLine,
                needSTM);

            //Step 3: wrap outputs
            boost::python::list outputState;
            for (size_t stateIndex = 0; stateIndex < nstates; ++stateIndex)
                outputState.append(outputStateLine(stateIndex));

            //Step 4: return the output state
            return outputState;
        }//end python propulateCoast()
        
        boost::python::list Propulator::getpythonSTM()
        {
            boost::python::list STMlist;

            math::Matrix<double> STM = this->getSTM();
            size_t n = STM.get_n();

            for (size_t i = 0; i < n; ++i)
            {
                boost::python::list STMrow;
                for (size_t j = 0; j < n; ++j)
                {
                    STMrow.append(STM(i, j));
                }

                STMlist.append(STMrow);
            }

            return STMlist;
        }//end getpythonSTM
#endif
    }//close namespace Propulator
}//close namespace EMTG