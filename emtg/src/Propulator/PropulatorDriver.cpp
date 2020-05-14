#include "Propulator.h"

#include "EMTG_Matrix.h"

#include <iostream>

int main(int argc, char* argv[])
{
    try
    {
        std::cout << "Propulation starting" << std::endl;

        //parse the options file
        std::string options_file_name;
        std::string input_states_file_name;
        std::string output_states_file_name;
        std::string output_STM_file_name;
        bool printSTM = false;
        if (argc >= 4)
        {
            options_file_name.assign(argv[1]);
            input_states_file_name.assign(argv[2]);
            output_states_file_name.assign(argv[3]);

            if (argc == 5)
            {
                output_STM_file_name.assign(argv[4]);
                printSTM = true;
            }
        }
        else
        {
            std::cout << "Syntax is PropagatorDriver OptionsFilePath InputStatesFilePath OutputStatesFilePath [OutputSTMFilePath]" << std::endl;
            std::cout << "input columns are:" << std::endl;
            std::cout << "#JD, journey, duty cycle, step size [sec], x [km/s], y [km/s], z [km/s], xdot [km/s], ydot [km/s], zdot [km/s], mass [kg], propagation time [sec]" << std::endl;
            return 0;
        }

        //parse the input file
        std::vector< EMTG::math::Matrix<double> > InputStateLines;
        std::vector< EMTG::math::Matrix<double> > ControlVectorLines;
        std::vector<double> PropagationTimes;
        std::string inputLine;
        std::ifstream inputfile(input_states_file_name);
        std::vector<double> dutyCycle;
        std::vector<double> stepSize;
        std::vector<int> journeyIndex;

        while (EMTG::file_utilities::safeGetline(inputfile, inputLine))
        {
            if (inputLine.size() > 0)
            {
                if (!(inputLine.front() == *"#"))
                {
                    std::vector<std::string> linecell;
                    boost::split(linecell, inputLine, boost::is_any_of(" ,"), boost::token_compress_on);

                    if (linecell.size() >= 12)
                    {
                        //epoch, x, y, z, xdot, ydot, zdot, mass, u_x, u_y, u_z, propagationTime
                        size_t cellIndex = 0;
                        EMTG::math::Matrix<double> inputStateLine(10, 1, 0.0);
                        EMTG::math::Matrix<double> inputControlLine(3, 1, 0.0);

                        //epoch
                        double JDepoch = std::stod(linecell[cellIndex++]);
                        inputStateLine(7) = (JDepoch - 2400000.5) * 86400.0;

                        //journey, duty cycle, and step size
                        journeyIndex.push_back(std::stoi(linecell[cellIndex++]));
                        dutyCycle.push_back(std::stod(linecell[cellIndex++]));
                        stepSize.push_back(std::stod(linecell[cellIndex++]));

                        //state
                        for (size_t stateIndex = 0; stateIndex < 7; ++stateIndex)
                            inputStateLine(stateIndex) = std::stod(linecell[cellIndex++]);

                        InputStateLines.push_back(inputStateLine);


                        if (linecell.size() == 15)
                        {
                            for (size_t controlIndex = 0; controlIndex < 3; ++controlIndex)
                                inputControlLine(controlIndex) = std::stod(linecell[cellIndex++]);

                            ControlVectorLines.push_back(inputControlLine);
                        }
                        else
                            ControlVectorLines.push_back(EMTG::math::Matrix<double>(1, 1, 0.0));

                        PropagationTimes.push_back(std::stod(linecell[cellIndex++]));

                    }
                    else
                    {
                        std::cout << "invalid line, '" << inputLine << "'" << std::endl;
                    }
                }
            }
        }//end loop over lines in file

        //instantiate a propulator
        EMTG::Propulator::Propulator myPropulator(options_file_name);

        //run the propulator
        EMTG::math::Matrix<double> OutputState(10, 1, 0.0);
        std::ofstream outputfile(output_states_file_name, std::ios::trunc);
        outputfile << "#States output from propulator" << std::endl;

        outputfile << std::endl;
        outputfile << "JD, x [km/s], y [km/s], z [km/s], xdot [km/s], ydot [km/s], zdot [km/s], mass [kg]" << std::endl;

        std::ofstream STMfile;
        if (printSTM)
        {
            STMfile = std::ofstream(output_STM_file_name, std::ios::trunc);
            STMfile << "#STM output from propulator" << std::endl;
            STMfile << "row ordered, i.e. each line goes [first row, second row, ...]" << std::endl;
            STMfile << std::endl;
        }

        for (size_t caseIndex = 0; caseIndex < InputStateLines.size(); ++caseIndex)
        {

            myPropulator.setJourneyIndex(journeyIndex[caseIndex]);
            if (stepSize[caseIndex] > 0.0)
                myPropulator.setStepSize(stepSize[caseIndex]);

            if (ControlVectorLines[caseIndex].get_n() > 1) //propagate with control
            {
                if (fabs(dutyCycle[caseIndex] < 1.0e-10))
                {
                    myPropulator.propulate(InputStateLines[caseIndex],
                        PropagationTimes[caseIndex],
                        OutputState,
                        printSTM);
                }
                else
                {
                    myPropulator.propulate(InputStateLines[caseIndex],
                        ControlVectorLines[caseIndex],
                        dutyCycle[caseIndex],
                        PropagationTimes[caseIndex],
                        OutputState,
                        printSTM);
                }
            }
            else //propagate without control
            {
                myPropulator.propulate(InputStateLines[caseIndex],
                    PropagationTimes[caseIndex],
                    OutputState,
                    printSTM);
            }
            outputfile.precision(14);
            outputfile << OutputState(7) / 86400.0 + 2400000.5;
            for (size_t stateIndex = 0; stateIndex < 7; ++stateIndex)
                outputfile << ", " << OutputState(stateIndex);
            outputfile << std::endl;

            if (printSTM)
            {
                EMTG::math::Matrix<double> STM = myPropulator.getSTM();
                for (size_t i = 0; i < STM.get_n(); ++i)
                {
                    for (size_t j = 0; j < STM.get_m(); ++j)
                    {
                        STMfile << STM(i, j);

                        if (!(i == STM.get_n() - 1 && j == STM.get_m() - 1)) //if we're not the last entry, write a comma
                            STMfile << ", ";
                    }
                }
                STMfile << std::endl;
            }
        }//end loop over cases


        std::cout << "Propulation complete." << std::endl;
    }
    catch (std::exception exception)
    {
        std::cout << "Propulator failed with error:" << std::endl;
        std::cout << exception.what() << std::endl;
    }
    return 0;
}