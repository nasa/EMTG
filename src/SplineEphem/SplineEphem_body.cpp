//SplineEphem body
//this is where the real work gets done
//Jacob Englander 11-8-2016

#include "SplineEphem_body.h"

#include <iostream>
#include <cmath>
#include <string>

namespace SplineEphem
{
    //default constructor - we actually don't want to use this so it will throw an error
    body::body()
    {
        throw std::invalid_argument("Don't call the SplineEphem::body default constructor. It doesn't do anything. Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
    }

    //actually useful constructor
    body::body(const int& SPICE_ID,
        const int& reference_body_SPICE_ID,
        const double& reference_body_mu,
        const size_t& number_of_steps_per_period,
        const double& tLowerBound,
        const double& tUpperBound)
    {
        this->initialize(SPICE_ID,
            reference_body_SPICE_ID,
            reference_body_mu,
            number_of_steps_per_period,
            tLowerBound,
            tUpperBound);
    }

    //destructor
    body::~body()
    {
        this->free();
    }

    //reinitialize function
    void body::reinitialize(const int& SPICE_ID,
        const int& reference_body_SPICE_ID,
        const double& reference_body_mu,
        const size_t& number_of_steps_per_period,
        const double& tLowerBound,
        const double& tUpperBound)
    {
        this->free();

        this->initialize(SPICE_ID,
            reference_body_SPICE_ID,
            reference_body_mu,
            number_of_steps_per_period,
            tLowerBound,
            tUpperBound);
    }

    //initialize function
    void body::initialize(const int& SPICE_ID,
        const int& reference_body_SPICE_ID,
        const double& reference_body_mu,
        const size_t& number_of_steps_per_period,
        const double& tLowerBound,
        const double& tUpperBound)
    {
        //ingest inputs
        this->SPICE_ID = SPICE_ID;
        this->reference_body_SPICE_ID = reference_body_SPICE_ID;
        this->reference_body_mu = reference_body_mu;
        this->ephemeris_window_open = tLowerBound;
        this->ephemeris_window_close = tUpperBound;

        //find my coverage window
        this->getCoverageWindow();
        

        this->ephemeris_window_close -= 86400.0;
        this->ephemeris_window_open += 86400.0;
        this->ephemeris_window_width = this->ephemeris_window_close - this->ephemeris_window_open;
        this->number_of_steps_per_period = number_of_steps_per_period;
		this->getPeriodFromSPICE();
        //different rules for bounded and unbounded orbits
        double SPICE_state[6], elements[8];
        double LT_dump;
        spkez_c(this->SPICE_ID, this->ephemeris_window_open - (51544.5 * 86400.0), "J2000", "NONE", this->reference_body_SPICE_ID, SPICE_state, &LT_dump);
        oscelt_c(SPICE_state,
            this->ephemeris_window_open - (51544.5 * 86400.0),
            this->reference_body_mu,
            elements);
        if (elements[0] / (1.0 - elements[1]) > 0) //SMA check for bounded orbits
        {
            this->number_of_periods = std::ceil(this->ephemeris_window_width / (this->Period * 86400.0));
            this->number_of_steps = this->number_of_periods * this->number_of_steps_per_period;
            this->time_step_width = this->ephemeris_window_width / this->number_of_steps;
        }
        else //unbounded orbits
        {
            this->number_of_periods = 1;
            this->number_of_steps = this->number_of_periods * this->number_of_steps_per_period;
            this->time_step_width = this->ephemeris_window_width / this->number_of_steps;
        }

        //size vectors
        this->t.resize(this->number_of_steps + 1);
        this->x.resize(this->number_of_steps + 1);
        this->y.resize(this->number_of_steps + 1);
        this->z.resize(this->number_of_steps + 1);
        this->xdot.resize(this->number_of_steps + 1);
        this->ydot.resize(this->number_of_steps + 1);
        this->zdot.resize(this->number_of_steps + 1);

        //read in SPICE data
        for (size_t tindex = 0; tindex <= this->number_of_steps; ++tindex)
        {
            this->t[tindex] = this->ephemeris_window_open + this->time_step_width * tindex;
            double SPICE_state[6];

            // if (this->t[tindex] > this->ephemeris_window_close)

            double LT_dump;
            spkez_c(this->SPICE_ID, this->t[tindex] - (51544.5 * 86400.0), "J2000", "NONE", this->reference_body_SPICE_ID, SPICE_state, &LT_dump);

            if (failed_c())
            {
                throw std::overflow_error("Error producing spline! Cannot find body " + std::to_string(this->SPICE_ID) + " with respect to "
                    + std::to_string(this->reference_body_SPICE_ID) + " on epoch " + std::to_string(this->t[tindex] / 86400.0) + ".\n"
                    + "Ephemeris window bounds for this body, from your .bsp, are [" + std::to_string(this->ephemeris_window_open / 86400.0) + ", " + std::to_string(this->ephemeris_window_close / 86400.0) + "]."
                    + "Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }    
            
            this->x[tindex] = SPICE_state[0];
            this->y[tindex] = SPICE_state[1];
            this->z[tindex] = SPICE_state[2];
            this->xdot[tindex] = SPICE_state[3];
            this->ydot[tindex] = SPICE_state[4];
            this->zdot[tindex] = SPICE_state[5];
        }

        //set up GSL accelerators and splines
        this->SplineAccelerator = gsl_interp_accel_alloc();
        this->Spline_x = gsl_spline_alloc(gsl_interp_cspline, this->number_of_steps + 1);
        this->Spline_y = gsl_spline_alloc(gsl_interp_cspline, this->number_of_steps + 1);
        this->Spline_z = gsl_spline_alloc(gsl_interp_cspline, this->number_of_steps + 1);
        this->Spline_xdot = gsl_spline_alloc(gsl_interp_cspline, this->number_of_steps + 1);
        this->Spline_ydot = gsl_spline_alloc(gsl_interp_cspline, this->number_of_steps + 1);
        this->Spline_zdot = gsl_spline_alloc(gsl_interp_cspline, this->number_of_steps + 1);

        //initialize GSL splines
        gsl_spline_init(this->Spline_x, t.data(), x.data(), this->number_of_steps + 1);
        gsl_spline_init(this->Spline_y, t.data(), y.data(), this->number_of_steps + 1);
        gsl_spline_init(this->Spline_z, t.data(), z.data(), this->number_of_steps + 1);
        gsl_spline_init(this->Spline_xdot, t.data(), xdot.data(), this->number_of_steps + 1);
        gsl_spline_init(this->Spline_ydot, t.data(), ydot.data(), this->number_of_steps + 1);
        gsl_spline_init(this->Spline_zdot, t.data(), zdot.data(), this->number_of_steps + 1);

        //clear the data arrays, because GSL splines make copies of their own data arrays for some reason and so we no longer need the originals
        t.clear();
        x.clear();
        y.clear();
        z.clear();
        xdot.clear();
        ydot.clear();
        zdot.clear();
    }

    //get coverage window
    void body::getCoverageWindow()
    {
        //Step 1: get coverage window for this body
        {
            //use spkcov to get coverage window
            const size_t  FILSIZ = 256;
            const size_t  LNSIZE = 81;
            const size_t  MAXCOV = 100000;
            const size_t  WINSIZ = (2 * MAXCOV);
            const size_t  TIMLEN = 51;
            const size_t MAXOBJ = 10000;

            SPICEDOUBLE_CELL(cover, WINSIZ);
            SPICEINT_CELL(ids, MAXOBJ);

            SpiceBoolean found;
            bool FoundBody = false;

            SpiceChar file[FILSIZ];
            SpiceChar source[FILSIZ];
            SpiceChar type[LNSIZE];

            SpiceDouble b = 0;
            SpiceDouble e = 0;

            SpiceInt count;
            SpiceInt handle;
            SpiceInt i;
            SpiceInt idcode = this->SPICE_ID;
            SpiceInt niv;

            /*
                Find out how many kernels are loaded.Loop over the
                kernels : for each loaded SPK file, add its coverage
                for `idcode', if any, to the coverage window.
            */
            ktotal_c("SPK", &count);

            double EarliestEpoch = this->ephemeris_window_open;
            double LatestEpoch = this->ephemeris_window_open;
            double temp_latest = this->ephemeris_window_open;

            for (i = 0; i < count; ++i)
            {
                kdata_c(i, "SPK", FILSIZ, LNSIZE, FILSIZ,
                    file, type, source, &handle, &found);


                //let's see if the object we want is in this file
                scard_c(0, &ids);
                spkobj_c(file, &ids);

                for (size_t j = 0; j < card_c(&ids); ++j)
                {
                    if (SPICE_CELL_ELEM_I(&ids, j) == this->SPICE_ID)
                    {
                        FoundBody = true;
                        scard_c(0, &cover);
                        spkcov_c(file, idcode, &cover);

                        niv = wncard_c(&cover);
                        FoundBody = true;
                        for (size_t k = 0; k < niv; ++k)
                        {
                            //Get the endpoints of the ith interval.
                            wnfetd_c(&cover, k, &b, &e);

                            if (b + (51544.5 * 86400.0) < EarliestEpoch)
                                EarliestEpoch = b + (51544.5 * 86400.0);
                            if (e + (51544.5 * 86400.0) > LatestEpoch)
                                LatestEpoch = e + (51544.5 * 86400.0);
                        }
                    }
                }

                this->ephemeris_window_open = fmax(this->ephemeris_window_open, EarliestEpoch);
                temp_latest = fmax(temp_latest, LatestEpoch);
            }

            this->ephemeris_window_close = fmin(this->ephemeris_window_close, temp_latest);

            if (!FoundBody)
            {
                throw std::invalid_argument("The SPICE kernel pool does not contain body " + std::to_string(this->SPICE_ID) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
        }//end coverage for the body of interest
        

        //Step 2: coverage for the reference body
        {
            //use spkcov to get coverage window
            const size_t  FILSIZ = 256;
            const size_t  LNSIZE = 81;
            const size_t  MAXCOV = 100000;
            const size_t  WINSIZ = (2 * MAXCOV);
            const size_t  TIMLEN = 51;
            const size_t MAXOBJ = 10000;

            SPICEDOUBLE_CELL(cover, WINSIZ);
            SPICEINT_CELL(ids, MAXOBJ);

            SpiceBoolean found;
            bool FoundBody = false;

            SpiceChar file[FILSIZ];
            SpiceChar source[FILSIZ];
            SpiceChar type[LNSIZE];

            SpiceDouble b = 0;
            SpiceDouble e = 0;

            SpiceInt count;
            SpiceInt handle;
            SpiceInt i;
            SpiceInt idcode = this->reference_body_SPICE_ID;
            SpiceInt niv;

            /*
                Find out how many kernels are loaded.Loop over the
                kernels : for each loaded SPK file, add its coverage
                for `idcode', if any, to the coverage window.
            */
            ktotal_c("SPK", &count);

            double EarliestEpoch = this->ephemeris_window_open;
            double LatestEpoch = this->ephemeris_window_open;
            double temp_latest = this->ephemeris_window_open;

            for (i = 0; i < count; ++i)
            {
                kdata_c(i, "SPK", FILSIZ, LNSIZE, FILSIZ,
                    file, type, source, &handle, &found);


                //let's see if the object we want is in this file
                scard_c(0, &ids);
                spkobj_c(file, &ids);

                for (size_t j = 0; j < card_c(&ids); ++j)
                {
                    if (SPICE_CELL_ELEM_I(&ids, j) == this->reference_body_SPICE_ID)
                    {
                        FoundBody = true;
                        scard_c(0, &cover);
                        spkcov_c(file, idcode, &cover);

                        niv = wncard_c(&cover);
                        FoundBody = true;
                        for (size_t k = 0; k < niv; ++k)
                        {
                            //Get the endpoints of the ith interval.
                            wnfetd_c(&cover, k, &b, &e);

                            if (b + (51544.5 * 86400.0) < EarliestEpoch)
                                EarliestEpoch = b + (51544.5 * 86400.0);
                            if (e + (51544.5 * 86400.0) > LatestEpoch)
                                LatestEpoch = e + (51544.5 * 86400.0);
                        }
                    }
                }

                this->ephemeris_window_open = fmax(this->ephemeris_window_open, EarliestEpoch);
                temp_latest = fmax(temp_latest, LatestEpoch);
            }

            this->ephemeris_window_close = fmin(this->ephemeris_window_close, temp_latest);

            if (!FoundBody)
            {
                throw std::invalid_argument("The SPICE kernel pool does not contain body " + std::to_string(this->reference_body_SPICE_ID) + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }
        }//end coverage of the reference body
    }//end getCoverageWindow()

    //function to get the period of the body with respect to the reference
    void body::getPeriodFromSPICE()
    {
        //we'll make the assumption that the period at the beginning of the coverage window is good enough to drive our sampling of ephemeris points

        //first we need the state vector of the body with respect to the reference body, 1 second after the coverage window opens
        double SPICE_state[6];

        if (failed_c())
            reset_c();

        double LT_dump;
        spkez_c(this->SPICE_ID, this->ephemeris_window_open - (51544.5 * 86400.0) + 1.0, "J2000", "NONE", this->reference_body_SPICE_ID, SPICE_state, &LT_dump);

        if (failed_c())
        {
            throw std::overflow_error("Error producing spline! Could not find a SPICE state for body " + std::to_string(this->SPICE_ID) + " with respect to "
                + std::to_string(this->reference_body_SPICE_ID) + " on epoch " + std::to_string((this->ephemeris_window_open + 1.0) / 86400.0) + ", and therefore could not compute the period."
                + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
        }     

        //now we can calculate orbit elements
        double SPICE_elements[8];
        oscelt_c(SPICE_state, this->ephemeris_window_open - (51544.5 * 86400.0) + 1.0, this->reference_body_mu, SPICE_elements);

        //and finally we can compute the period
        double SMA = SPICE_elements[0] / (1 - SPICE_elements[1]);
        this->Period = 2.0 * 3.1415926535897932 * sqrt(SMA * SMA * SMA / this->reference_body_mu) / 86400.0;
    }

    //free function
    void body::free()
    {
        //free GSL accelerators and splines
        gsl_spline_free(this->Spline_x);
        gsl_spline_free(this->Spline_y);
        gsl_spline_free(this->Spline_z);
        gsl_spline_free(this->Spline_xdot);
        gsl_spline_free(this->Spline_ydot);
        gsl_spline_free(this->Spline_zdot);
        gsl_interp_accel_free(this->SplineAccelerator);
    }

    //evaluate functions
    void body::get6State(const double& epoch, double* StateArray)
    {
        this->getPosition(epoch, StateArray);
        this->getVelocity(epoch, StateArray + 3);
    }

    void body::get6StateDerivative(const double& epoch, double* StateDerivativeArray)
    {
        this->getPositionDerivative(epoch, StateDerivativeArray);
        this->getVelocityDerivative(epoch, StateDerivativeArray + 3);
    }

    void body::get6StateAndDerivative(const double& epoch, double* StateAndDerivativeArray)
    {
        this->get6State(epoch, StateAndDerivativeArray);
        this->get6StateDerivative(epoch, StateAndDerivativeArray + 6);
    }
    
}//end namespace SplineEphem