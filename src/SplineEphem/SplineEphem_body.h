//header file for SplineEphem body
//this is where the real work gets done
//Jacob Englander 11-8-2016

//relies on the GSL, which is NOT compatible with the NOSA
//therefore do not EVER distribute GSL with this, users must get their own

#ifdef SPLINE_EPHEM
#ifndef SPLINEPHEM_BODY
#define SPLINEPHEM_BODY

#include "SpiceUsr.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_errno.h" //I'm not really sure what this is but the GSL examples all have it, probably error handling

#include <vector>
#include <iostream>
#include <string>

namespace SplineEphem
{
    class body
    {
    public:
        //constructors
        body();
        body(const int& SPICE_ID,
            const int& reference_body_SPICE_ID,
            const double& reference_body_mu,
            const size_t& number_of_steps_per_period = 100,
            const double& tLowerBound = 30000.0*86400.0, 
            const double& tUpperBound = 100000.0*86400.0);
        //destructor
        ~body();

        //methods
        void reinitialize(const int& SPICE_ID,
            const int& reference_body_SPICE_ID,
            const double& reference_body_mu,
            const size_t& number_of_steps_per_period = 100,
            const double& tLowerBound = 30000.0*86400.0,
            const double& tUpperBound = 100000.0*86400.0);

        inline void getPosition(const double& epoch, double* PositionArray)
        {
            if (epoch > this->ephemeris_window_close || epoch < this->ephemeris_window_open)
            {
                throw std::runtime_error("SplineEphem cannot find body " + std::to_string(this->SPICE_ID) + " with respect to " + std::to_string(this->reference_body_SPICE_ID)
                    + " on epoch " + std::to_string(epoch / 86400.0) + ". Ephemeris window opens on MJD " + std::to_string(this->ephemeris_window_open / 86400.0) + " and closes on MJD "
                    + std::to_string(this->ephemeris_window_close / 86400.0)
                    + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            PositionArray[0] = gsl_spline_eval(this->Spline_x, epoch, this->SplineAccelerator);
            PositionArray[1] = gsl_spline_eval(this->Spline_y, epoch, this->SplineAccelerator);
            PositionArray[2] = gsl_spline_eval(this->Spline_z, epoch, this->SplineAccelerator);
        };

        inline void getVelocity(const double& epoch, double* VelocityArray)
        {
            if (epoch > this->ephemeris_window_close || epoch < this->ephemeris_window_open)
            {
                throw std::runtime_error("SplineEphem cannot find body " + std::to_string(this->SPICE_ID) + " with respect to " + std::to_string(this->reference_body_SPICE_ID)
                    + " on epoch " + std::to_string(epoch / 86400.0) + ". Ephemeris window opens on MJD " + std::to_string(this->ephemeris_window_open / 86400.0) + " and closes on MJD "
                    + std::to_string(this->ephemeris_window_close / 86400.0)
                    + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            VelocityArray[0] = gsl_spline_eval(this->Spline_xdot, epoch, this->SplineAccelerator);
            VelocityArray[1] = gsl_spline_eval(this->Spline_ydot, epoch, this->SplineAccelerator);
            VelocityArray[2] = gsl_spline_eval(this->Spline_zdot, epoch, this->SplineAccelerator);
        };

        void get6State(const double& epoch, double* StateArray);

        inline void getPositionDerivative(const double& epoch, double* PositionDerivativeArray)
        {
            if (epoch > this->ephemeris_window_close || epoch < this->ephemeris_window_open)
            {
                throw std::runtime_error("SplineEphem cannot find body " + std::to_string(this->SPICE_ID) + " with respect to " + std::to_string(this->reference_body_SPICE_ID)
                    + " on epoch " + std::to_string(epoch / 86400.0) + ". Ephemeris window opens on MJD " + std::to_string(this->ephemeris_window_open / 86400.0) + " and closes on MJD "
                    + std::to_string(this->ephemeris_window_close / 86400.0)
                    + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            PositionDerivativeArray[0] = gsl_spline_eval_deriv(this->Spline_x, epoch, this->SplineAccelerator);
            PositionDerivativeArray[1] = gsl_spline_eval_deriv(this->Spline_y, epoch, this->SplineAccelerator);
            PositionDerivativeArray[2] = gsl_spline_eval_deriv(this->Spline_z, epoch, this->SplineAccelerator);
        };

        inline void getVelocityDerivative(const double& epoch, double* VelocityDerivativeArray)
        {
            if (epoch > this->ephemeris_window_close || epoch < this->ephemeris_window_open)
            {
                throw std::runtime_error("SplineEphem cannot find body " + std::to_string(this->SPICE_ID) + " with respect to " + std::to_string(this->reference_body_SPICE_ID)
                    + " on epoch " + std::to_string(epoch / 86400.0) + ". Ephemeris window opens on MJD " + std::to_string(this->ephemeris_window_open / 86400.0) + " and closes on MJD "
                    + std::to_string(this->ephemeris_window_close / 86400.0)
                    + ". Place a breakpoint in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__));
            }

            VelocityDerivativeArray[0] = gsl_spline_eval_deriv(this->Spline_xdot, epoch, this->SplineAccelerator);
            VelocityDerivativeArray[1] = gsl_spline_eval_deriv(this->Spline_ydot, epoch, this->SplineAccelerator);
            VelocityDerivativeArray[2] = gsl_spline_eval_deriv(this->Spline_zdot, epoch, this->SplineAccelerator);
        };

        void get6StateDerivative(const double& epoch, double* StateDerivativeArray);
        void get6StateAndDerivative(const double& epoch, double* StateAndDerivativeArray);

        inline int getSPICE_ID() 
        {
            return this->SPICE_ID;
        };

        inline int getReferenceBody_SPICE_ID()
        {
            return this->reference_body_SPICE_ID;
        };

        inline double getEphemerisWindowOpen() const { return this->ephemeris_window_open; }
        inline double getEphemerisWindowClose() const { return this->ephemeris_window_close; }

    private:
        //private methods, if applicable

        void initialize(const int& SPICE_ID,
            const int& reference_body_SPICE_ID,
            const double& reference_body_mu,
            const size_t& number_of_steps_per_period = 100,
            const double& tLowerBound = 30000.0*86400.0,
            const double& tUpperBound = 100000.0*86400.0);

        void getCoverageWindow();

        void getPeriodFromSPICE();

        void free();

        //*********fields

        //setup information
        int SPICE_ID;
        int reference_body_SPICE_ID;
        double reference_body_mu;
        double ephemeris_window_open;//in MJD, defined as JD - 2400000.5
        double ephemeris_window_close;//in MJD, defined as JD - 2400000.5
        double ephemeris_window_width;//in days
        double Period;
        double time_step_width;
        size_t number_of_periods;
        size_t number_of_steps;
        size_t number_of_steps_per_period;
        
        //underlying data
        std::vector<double> t;
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;
        std::vector<double> xdot;
        std::vector<double> ydot;
        std::vector<double> zdot;

        //GSL spline structs
        gsl_interp_accel *SplineAccelerator;
        gsl_spline *Spline_x;
        gsl_spline *Spline_y;
        gsl_spline *Spline_z;
        gsl_spline *Spline_xdot;
        gsl_spline *Spline_ydot;
        gsl_spline *Spline_zdot;
    };//end class body
}//end namespace SplineEphem

#endif //SPLINEPHEM_BODY
#endif //#SPLINE_EPHEM