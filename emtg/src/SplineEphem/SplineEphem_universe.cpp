//SplineEphem Universe container
//all this does is hold all of your bodies
//Jacob Englander 11-8-2016

#include "SplineEphem_universe.h"

#include <iostream>

namespace SplineEphem
{
    //constructors
    universe::universe()
    {
        this->number_of_bodies = 0;
    }

    universe::universe(const std::vector< std::tuple<int, int, int, double> >& keyList, 
        const double& ephemeris_window_open,
        const double& ephemeris_window_close)
    {
        this->initialize(keyList,
            ephemeris_window_open,
            ephemeris_window_close);
    }

    //destructor
    universe::~universe()
    {
        this->clear();
    }

    //reinitialize function
    void universe::reinitialize(const std::vector< std::tuple<int, int, int, double> >& keyList,
        const double& ephemeris_window_open,
        const double& ephemeris_window_close)
    {
        this->clear();
        this->initialize(keyList,
            ephemeris_window_open,
            ephemeris_window_close);
    }

    //get functions
    void universe::getBodyPosition(const int& SPICE_ID, const int& ReferenceBody_SPICE_ID, const double& epoch, double* PositionArray)
    {
        int body_index = this->getBodyIndex(SPICE_ID, ReferenceBody_SPICE_ID);

        if (body_index < 0)
        {
            throw std::invalid_argument("User has requested a (body, reference) pair from SplineEphem that SplineEphem did not load. The requested pair was ("
                + std::to_string(SPICE_ID) + ", " +  std::to_string(ReferenceBody_SPICE_ID) + ")."
                + "This could be user error or a bug in SplineEphem. Please report this message to the developers. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
        }

        this->bodies[body_index]->getPosition(epoch, PositionArray);
    }
    void universe::getBodyVelocity(const int& SPICE_ID, const int& ReferenceBody_SPICE_ID, const double& epoch, double* VelocityArray)
    {
        int body_index = this->getBodyIndex(SPICE_ID, ReferenceBody_SPICE_ID);

        if (body_index < 0)
        {
            throw std::invalid_argument("User has requested a (body, reference) pair from SplineEphem that SplineEphem did not load. The requested pair was ("
                + std::to_string(SPICE_ID) + ", " + std::to_string(ReferenceBody_SPICE_ID) + ")."
                + "This could be user error or a bug in SplineEphem. Please report this message to the developers. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
        }

        this->bodies[body_index]->getVelocity(epoch, VelocityArray);
    }
    void universe::getBody6State(const int& SPICE_ID, const int& ReferenceBody_SPICE_ID, const double& epoch, double* StateArray)
    {
        int body_index = this->getBodyIndex(SPICE_ID, ReferenceBody_SPICE_ID);

        if (body_index < 0)
        {
            throw std::invalid_argument("User has requested a (body, reference) pair from SplineEphem that SplineEphem did not load. The requested pair was ("
                + std::to_string(SPICE_ID) + ", " + std::to_string(ReferenceBody_SPICE_ID) + ")."
                + "This could be user error or a bug in SplineEphem. Please report this message to the developers. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
        }

        this->bodies[body_index]->get6State(epoch, StateArray);
    }
    void universe::getBodyPositionDerivative(const int& SPICE_ID, const int& ReferenceBody_SPICE_ID, const double& epoch, double* PositionDerivativeArray)
    {
        int body_index = this->getBodyIndex(SPICE_ID, ReferenceBody_SPICE_ID);

        if (body_index < 0)
        {
            throw std::invalid_argument("User has requested a (body, reference) pair from SplineEphem that SplineEphem did not load. The requested pair was ("
                + std::to_string(SPICE_ID) + ", " + std::to_string(ReferenceBody_SPICE_ID) + ")."
                + "This could be user error or a bug in SplineEphem. Please report this message to the developers. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
        }

        this->bodies[body_index]->getPositionDerivative(epoch, PositionDerivativeArray);
    }
    void universe::getBodyVelocityDerivative(const int& SPICE_ID, const int& ReferenceBody_SPICE_ID, const double& epoch, double* VelocityDerivativeArray)
    {
        int body_index = this->getBodyIndex(SPICE_ID, ReferenceBody_SPICE_ID);

        if (body_index < 0)
        {
            throw std::invalid_argument("User has requested a (body, reference) pair from SplineEphem that SplineEphem did not load. The requested pair was ("
                + std::to_string(SPICE_ID) + ", " + std::to_string(ReferenceBody_SPICE_ID) + ")."
                + "This could be user error or a bug in SplineEphem. Please report this message to the developers. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
        }

        this->bodies[body_index]->getVelocityDerivative(epoch, VelocityDerivativeArray);
    }
    void universe::getBody6StateDerivative(const int& SPICE_ID, const int& ReferenceBody_SPICE_ID, const double& epoch, double* StateDerivativeArray)
    {
        int body_index = this->getBodyIndex(SPICE_ID, ReferenceBody_SPICE_ID);

        if (body_index < 0)
        {
            throw std::invalid_argument("User has requested a (body, reference) pair from SplineEphem that SplineEphem did not load. The requested pair was ("
                + std::to_string(SPICE_ID) + ", " + std::to_string(ReferenceBody_SPICE_ID) + ")."
                + "This could be user error or a bug in SplineEphem. Please report this message to the developers. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
        }

        this->bodies[body_index]->get6StateDerivative(epoch, StateDerivativeArray);
    }
    void universe::getBody6StateAndDerivative(const int& SPICE_ID, const int& ReferenceBody_SPICE_ID, const double& epoch, double* StateAndDerivativeArray)
    {
        int body_index = this->getBodyIndex(SPICE_ID, ReferenceBody_SPICE_ID);

        if (body_index < 0)
        {
            throw std::invalid_argument("User has requested a (body, reference) pair from SplineEphem that SplineEphem did not load. The requested pair was ("
                + std::to_string(SPICE_ID) + ", " + std::to_string(ReferenceBody_SPICE_ID) + ")."
                + "This could be user error or a bug in SplineEphem. Please report this message to the developers. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
        }

        this->bodies[body_index]->get6StateAndDerivative(epoch, StateAndDerivativeArray);
    }


    double universe::getEphemerisWindowOpen(const int& SPICE_ID, const int& ReferenceBody_SPICE_ID)
    {
        int body_index = this->getBodyIndex(SPICE_ID, ReferenceBody_SPICE_ID);
        if (body_index >= 0)
            return this->bodies[body_index]->getEphemerisWindowOpen();
		else
			return 0;
    }

    double universe::getEphemerisWindowClose(const int& SPICE_ID, const int& ReferenceBody_SPICE_ID)
    {
        int body_index = this->getBodyIndex(SPICE_ID, ReferenceBody_SPICE_ID);
        if (body_index >= 0)
            return this->bodies[body_index]->getEphemerisWindowClose();
		else
			return 0;
    }

    //initialize method
    void universe::initialize(const std::vector< std::tuple<int, int, int, double> >& keyList,
        const double& ephemeris_window_open,
        const double& ephemeris_window_close)
    {
        for (size_t b = 0; b < keyList.size(); ++b)
        {
            try
            {
                this->bodies.push_back(std::unique_ptr<SplineEphem::body>(new SplineEphem::body(std::get<0>(keyList[b]),
                    std::get<1>(keyList[b]),
                    std::get<3>(keyList[b]), 
                    std::get<2>(keyList[b]),
                    ephemeris_window_open,
                    ephemeris_window_close)));
            }
            catch (std::exception &error)
            {
                std::cout << error.what() << std::endl;
                std::cout << "error in configuring body " << b << ", " << std::get<0>(keyList[b]) << " with respect to " << std::get<1>(keyList[b]) << std::endl;
                throw;
            }
        }

        this->number_of_bodies = this->bodies.size();
    }

    //clear method
    void universe::clear()
    {
        this->bodies.clear();
        this->number_of_bodies = 0;
    }

    int universe::getBodyIndex(const int& SPICE_ID, const int& ReferenceBody_SPICE_ID)
    {
        int BodyIndex = -1;//this will return -1 if the body is not found

        for (size_t b = 0; b < this->number_of_bodies; ++b)
        {
            if (SPICE_ID == this->bodies[b]->getSPICE_ID() && ReferenceBody_SPICE_ID == this->bodies[b]->getReferenceBody_SPICE_ID())
            {
                BodyIndex = b;
                break;
            }
        }    

        return BodyIndex;
    }
}//close namespace SplineEphem