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

//Kepler propagator in the time domain

#include "KeplerPropagatorTimeDomain.h"

#include <exception>

namespace EMTG
{
    namespace Astrodynamics
    {

        KeplerPropagatorTimeDomain::KeplerPropagatorTimeDomain(const size_t& numStates) :
                                                               KeplerPropagator(numStates)
        {}

        //methods
        void KeplerPropagatorTimeDomain::propagate(const doubleType& PropagationTime, const bool& needSTM)
        {
			// In order to scale this propagator to canonical units, 
			// set LU/TU to their values from myUniverse (e.g. LU = this->myUniverse->LU) and set mu to 1.0
		    double LU = 1.0;
			double TU = 1.0;
			//double mu = 1.0;
			//double LU = this->myUniverse->LU;
			//double TU = this->myUniverse->TU;

            //Step 0: declare necessary variables
            int PropagationDirection = PropagationTime >= 0.0 ? 1 : -1;
			double sqmu = sqrt(this->mu);
            const double alphatol = 1.0e-12;
            const double Xtol = 1.0e-12;
            const int maximum_order = 10;
            const int maximum_iterations_per_order = 10;
            doubleType r, sigma, U0 = 0.0, U1 = 0.0, U2 = 0.0, U3 = 0.0, X = 1.0e+100, X_new = 0.0, dX = 1.0e+100;
            doubleType sqalpha = 0.0, sqmalpha = 0.0, sqmalphaX = 0.0;
            doubleType state0[6], state[6];
            double fstate0[6], fstate[6];
            math::Matrix<doubleType>& StateLeft = *this->StateLeftPointer;
            math::Matrix<doubleType>& StateRight = *this->StateRightPointer;
            math::Matrix<double>& STM = *this->STMpointer;
            math::Matrix<double>& dStatedIndependentVariable = *this->dStatedIndependentVariablePointer;

            //Step 1: compute r0 and v0 in LUs and TUs, and propTime in TUs
            for (int stateIndex = 0; stateIndex < 3; ++stateIndex)
            {
                state0[stateIndex] = StateLeft(stateIndex) / LU;
                fstate0[stateIndex] = state0[stateIndex] _GETVALUE;
            }
            for (int stateIndex = 3; stateIndex < 6; ++stateIndex)
            {
                state0[stateIndex] = StateLeft(stateIndex) / LU * TU;
                fstate0[stateIndex] = state0[stateIndex] _GETVALUE;
            }



            doubleType r0 = sqrt(state0[0] * state0[0] + state0[1] * state0[1] + state0[2] * state0[2]);
            doubleType v0 = sqrt(state0[3] * state0[3] + state0[4] * state0[4] + state0[5] * state0[5]);
            doubleType propTime = PropagationTime / TU;

            //Step 2: compute alpha, which determines whether we are on an ellipse, parabola, or hyperbola
            //and also sigma0 which we need in the universal Kepler iteration
            doubleType alpha = (2.0 / r0 - v0*v0 / this->mu);
            doubleType sigma0 = (state0[0] * state0[3] + state0[1] * state0[4] + state0[2] * state0[5]) / sqmu;
            if (alpha > alphatol) //ellipse
                sqalpha = sqrt(alpha);
            else if (alpha < -alphatol)
                sqmalpha = sqrt(-alpha);

            //Step 3: compute an initial guess
            if (alpha > alphatol)
                X_new = alpha * sqmu * propTime;
            else
                X_new = 0.1 * sqmu * propTime / r0;
            r = r0;

            //Step 4: Perform the Laguerre-Conway-Der iteration
            int N = 2; //current Laguerre iteration order
            int iteration_this_N = 0; //number of iterations for the current order
            int total_iterations = 0;
            while (fabs(X - X_new) > Xtol && N < maximum_order)
            {
                try
                {
                    //Step 4.1: increment the iteration count
                    //if we have maxed out the number of iterations for this order, increase the order
                    ++iteration_this_N;
                    ++total_iterations;
                    if (iteration_this_N >= maximum_iterations_per_order)
                    {
                        ++N; //increase the order
                        iteration_this_N = 0; //reset the number of iterations
                    }

                    //Step 4.2: update X
                    X = X_new;

                    //Step 4.3: compute U0, U1, U2, U3 for the candidate X
                    if (alpha > alphatol) //ellipse
                    {
                        //compute U0, U1, U2, U3 via Stumpff functions
                        doubleType y = alpha * X * X;
                        doubleType C = (1.0 - cos(sqrt(y))) / y;
                        doubleType S = (sqrt(y) - sin(sqrt(y))) / sqrt(y*y*y);
                        U1 = X * (1.0 - y * S);
                        U2 = X * X * C;
                        U3 = X * X * X * S;
                        U0 = 1.0 - alpha * U2;
                    }
                    else if (alpha < -alphatol) //hyperbola
                    {

                        sqmalphaX = sqmalpha * X;
                        if (sqmalphaX > 30.0 || sqmalphaX < -30.0)
                        {
                            throw std::runtime_error("Kepler solver failed, went too far out on a hyperbola");
                        }
                        U0 = cosh(sqmalphaX);
                        U1 = sinh(sqmalphaX) / sqmalpha;
                        U2 = (1.0 - U0) / alpha;
                        U3 = 1.0 / alpha * (X - U1);
                    }
                    else //parabola
                    {
                        U0 = 1.0;
                        U1 = X;
                        U2 = U1*X / 2.0;
                        U3 = U2*X / 3.0;
                    }



                    //Step 4.4: compute r and sigma
                    r = r0 * U0 + sigma0 * U1 + U2;
                    sigma = sigma0 * U0 + (1.0 - alpha * r0) * U1;

                    //Step 4.5 compute F, dF, and ddF
                    doubleType FX = r0 * U1 + sigma0 * U2 + U3 - sqmu * propTime;
                    doubleType dFX = r;
                    doubleType ddFX = sigma;

                    //Step 4.6 Laguerre-Conway or Newton iteration depending on the situation
                    int sgn = dFX >= 0 ? 1 : -1;
                    doubleType denom = fabs((N - 1)*(N - 1)*dFX*dFX - N * (N - 1) * FX * ddFX);
                    if (denom > 0.0) //Laguerre-Conway iteration
                    {
                        dX = N*FX / (dFX + sgn * sqrt(denom));
                    }
                    else //Newton iteration for very sensitive cases, i.e. when on a hyperbolic/parabolic asymptote
                        dX = FX / dFX;

                    X_new -= dX;
                } //end try-catch block on Kepler iteration
                catch (std::runtime_error error)
                {
                    if (N < maximum_order) //increase the order if we can
                    {
                        ++N;
                        iteration_this_N = 0;
                    }
                    else //otherwise break out of the while loop and return the partially converged solution
                        break;
                }
            }

            //Step 5: find F, G, Ft, Gt
            doubleType aF = 1.0 - U2 / r0;
            doubleType aG = (r0 * U1 + sigma0 * U2) / sqmu;
            doubleType aFt = -sqmu / (r0 * r) * U1;
            doubleType aGt = 1.0 - U2 / r;
            F = aF _GETVALUE;
            G = aG _GETVALUE;
            Ft = aFt _GETVALUE;
            Gt = aGt _GETVALUE;

            //Step 6: compute the final state as functions of F, G, Ft, Gt
            state[0] = (aF*state0[0] + aG*state0[3]);
            state[1] = (aF*state0[1] + aG*state0[4]);
            state[2] = (aF*state0[2] + aG*state0[5]);
            state[3] = (aFt*state0[0] + aGt*state0[3]);
            state[4] = (aFt*state0[1] + aGt*state0[4]);
            state[5] = (aFt*state0[2] + aGt*state0[5]);

            for (int stateIndex = 0; stateIndex < 3; ++stateIndex)
            {
                StateRight(stateIndex) = state[stateIndex] * LU;
                fstate[stateIndex] = state[stateIndex] _GETVALUE;
            }
            for (int stateIndex = 3; stateIndex < 6; ++stateIndex)
            {
                StateRight(stateIndex) = state[stateIndex] * LU / TU;
                fstate[stateIndex] = state[stateIndex] _GETVALUE;
            }

            //Step 7: Compute the state transition matrix if requested
            if (needSTM)
            {
                //STM calculation code
                //Step 7.1 compute the universal functions Ui and their derivatives
                //from the recursion relation 4.76
                double fr = r _GETVALUE;
                double fr0 = r0 _GETVALUE;
                double fv0 = v0 _GETVALUE;

                double a;
                if (fabs(alpha) < alphatol)
                    a = 1.0e+30;
                else
                    a = 1.0 / alpha _GETVALUE;

                double fX = X _GETVALUE;
                double U4 = a * (fX*fX / 2.0 - U2 _GETVALUE);
                double U5 = a * (fX*fX*fX / 6.0 - U3 _GETVALUE);

                //derivatives
                double dXdt = sqmu / fr;
                double U0dot = -alpha _GETVALUE * U1 _GETVALUE * dXdt;
                double U1dot = U0 _GETVALUE * dXdt;
                double U2dot = U1 _GETVALUE * dXdt;

                //Step 5.2 compute C, which along with F, G, Ft, and Gt can be used to compute the STM
                //Battin equation 9.74
                double C = 1.0 / sqmu * (3 * U5 - fX*U4) - propTime _GETVALUE * U2 _GETVALUE;

                //Step 5.3 compute R, V, R~ and V~, Battin equations 9.84 - 9.87
                double x0 = fstate0[0];
                double y0 = fstate0[1];
                double z0 = fstate0[2];
                double xdot0 = fstate0[3];
                double ydot0 = fstate0[4];
                double zdot0 = fstate0[5];
                double x = fstate[0];
                double y = fstate[1];
                double z = fstate[2];
                double xdot = fstate[3];
                double ydot = fstate[4];
                double zdot = fstate[5];
                double fr_2 = fr*fr;
                double fr0_2 = fr0*fr0;
                double fr_3 = fr_2 * fr;
                double fr0_3 = fr0_2 * fr0;
                double mu2 = this->mu * this->mu;

                //R~
                STM(0, 0) = fr / this->mu*(xdot - xdot0)*(xdot - xdot0) + (fr0*(1.0 - F)*(x * x0) + C*(x0 * xdot)) / (fr0_3)+F;
                STM(0, 1) = fr / this->mu*(xdot - xdot0)*(ydot - ydot0) + (fr0*(1.0 - F)*(x * y0) + C*(y0 * xdot)) / (fr0_3);
                STM(0, 2) = fr / this->mu*(xdot - xdot0)*(zdot - zdot0) + (fr0*(1.0 - F)*(x * z0) + C*(z0 * xdot)) / (fr0_3);
                STM(1, 0) = fr / this->mu*(ydot - ydot0)*(xdot - xdot0) + (fr0*(1.0 - F)*(y * x0) + C*(x0 * ydot)) / (fr0_3);
                STM(1, 1) = fr / this->mu*(ydot - ydot0)*(ydot - ydot0) + (fr0*(1.0 - F)*(y * y0) + C*(y0 * ydot)) / (fr0_3)+F;
                STM(1, 2) = fr / this->mu*(ydot - ydot0)*(zdot - zdot0) + (fr0*(1.0 - F)*(y * z0) + C*(z0 * ydot)) / (fr0_3);
                STM(2, 0) = fr / this->mu*(zdot - zdot0)*(xdot - xdot0) + (fr0*(1.0 - F)*(z * x0) + C*(x0 * zdot)) / (fr0_3);
                STM(2, 1) = fr / this->mu*(zdot - zdot0)*(ydot - ydot0) + (fr0*(1.0 - F)*(z * y0) + C*(y0 * zdot)) / (fr0_3);
                STM(2, 2) = fr / this->mu*(zdot - zdot0)*(zdot - zdot0) + (fr0*(1.0 - F)*(z * z0) + C*(z0 * zdot)) / (fr0_3)+F;

                //R
                STM(0, 3) = fr0 / this->mu*(1.0 - F)*(xdot0 * (x - x0) - x0 * (xdot - xdot0)) + C / this->mu*(xdot * xdot0) + G;
                STM(0, 4) = fr0 / this->mu*(1.0 - F)*(ydot0 * (x - x0) - y0 * (xdot - xdot0)) + C / this->mu*(xdot * ydot0);
                STM(0, 5) = fr0 / this->mu*(1.0 - F)*(zdot0 * (x - x0) - z0 * (xdot - xdot0)) + C / this->mu*(xdot * zdot0);
                STM(1, 3) = fr0 / this->mu*(1.0 - F)*(xdot0 * (y - y0) - x0 * (ydot - ydot0)) + C / this->mu*(ydot * xdot0);
                STM(1, 4) = fr0 / this->mu*(1.0 - F)*(ydot0 * (y - y0) - y0 * (ydot - ydot0)) + C / this->mu*(ydot * ydot0) + G;
                STM(1, 5) = fr0 / this->mu*(1.0 - F)*(zdot0 * (y - y0) - z0 * (ydot - ydot0)) + C / this->mu*(ydot * zdot0);
                STM(2, 3) = fr0 / this->mu*(1.0 - F)*(xdot0 * (z - z0) - x0 * (zdot - zdot0)) + C / this->mu*(zdot * xdot0);
                STM(2, 4) = fr0 / this->mu*(1.0 - F)*(ydot0 * (z - z0) - y0 * (zdot - zdot0)) + C / this->mu*(zdot * ydot0);
                STM(2, 5) = fr0 / this->mu*(1.0 - F)*(zdot0 * (z - z0) - z0 * (zdot - zdot0)) + C / this->mu*(zdot * zdot0) + G;

                //V~
                STM(3, 0) = (-C*mu2 * x*x0 + Ft*fr*fr0_3 * (this->mu*fr_2 - this->mu*x*x + fr*(xdot0 - xdot)*(y*(xdot*y - ydot*x) + z*(xdot*z - zdot*x))) + this->mu*fr_3 * fr0*x0*(xdot0 - xdot) + this->mu*fr*fr0_3 * x*(xdot0 - xdot)) / (this->mu*fr_3 * fr0_3);
                STM(3, 1) = (-C*mu2 * x*y0 + Ft*fr*fr0_3 * (-this->mu*x*y + fr*(ydot0 - ydot)*(y*(xdot*y - ydot*x) + z*(xdot*z - zdot*x))) + this->mu*fr_3 * fr0*y0*(xdot0 - xdot) + this->mu*fr*fr0_3 * x*(ydot0 - ydot)) / (this->mu*fr_3 * fr0_3);
                STM(3, 2) = (-C*mu2 * x*z0 + Ft*fr*fr0_3 * (-this->mu*x*z + fr*(zdot0 - zdot)*(y*(xdot*y - ydot*x) + z*(xdot*z - zdot*x))) + this->mu*fr_3 * fr0*z0*(xdot0 - xdot) + this->mu*fr*fr0_3 * x*(zdot0 - zdot)) / (this->mu*fr_3 * fr0_3);
                STM(4, 0) = (-C*mu2 * x0*y - Ft*fr*fr0_3 * (this->mu*x*y + fr*(xdot0 - xdot)*(x*(xdot*y - ydot*x) - z*(ydot*z - zdot*y))) + this->mu*fr_3 * fr0*x0*(ydot0 - ydot) + this->mu*fr*fr0_3 * y*(xdot0 - xdot)) / (this->mu*fr_3 * fr0_3);
                STM(4, 1) = (-C*mu2 * y*y0 - Ft*fr*fr0_3 * (-this->mu*fr_2 + this->mu*y*y + fr*(ydot0 - ydot)*(x*(xdot*y - ydot*x) - z*(ydot*z - zdot*y))) + this->mu*fr_3 * fr0*y0*(ydot0 - ydot) + this->mu*fr*fr0_3 * y*(ydot0 - ydot)) / (this->mu*fr_3 * fr0_3);
                STM(4, 2) = (-C*mu2 * y*z0 - Ft*fr*fr0_3 * (this->mu*y*z + fr*(zdot0 - zdot)*(x*(xdot*y - ydot*x) - z*(ydot*z - zdot*y))) + this->mu*fr_3 * fr0*z0*(ydot0 - ydot) + this->mu*fr*fr0_3 * y*(zdot0 - zdot)) / (this->mu*fr_3 * fr0_3);
                STM(5, 0) = (-C*mu2 * x0*z - Ft*fr*fr0_3 * (this->mu*x*z + fr*(xdot0 - xdot)*(x*(xdot*z - zdot*x) + y*(ydot*z - zdot*y))) + this->mu*fr_3 * fr0*x0*(zdot0 - zdot) + this->mu*fr*fr0_3 * z*(xdot0 - xdot)) / (this->mu*fr_3 * fr0_3);
                STM(5, 1) = (-C*mu2 * y0*z - Ft*fr*fr0_3 * (this->mu*y*z + fr*(ydot0 - ydot)*(x*(xdot*z - zdot*x) + y*(ydot*z - zdot*y))) + this->mu*fr_3 * fr0*y0*(zdot0 - zdot) + this->mu*fr*fr0_3 * z*(ydot0 - ydot)) / (this->mu*fr_3 * fr0_3);
                STM(5, 2) = (-C*mu2 * z*z0 - Ft*fr*fr0_3 * (-this->mu*fr_2 + this->mu*z*z + fr*(zdot0 - zdot)*(x*(xdot*z - zdot*x) + y*(ydot*z - zdot*y))) + this->mu*fr_3 * fr0*z0*(zdot0 - zdot) + this->mu*fr*fr0_3 * z*(zdot0 - zdot)) / (this->mu*fr_3 * fr0_3);

                //V
                STM(3, 3) = fr0 / this->mu*(xdot - xdot0)*(xdot - xdot0) + (fr0*(1.0 - F)*(x * x0) - C*(x * xdot0)) / (fr_3)+Gt;
                STM(3, 4) = fr0 / this->mu*(xdot - xdot0)*(ydot - ydot0) + (fr0*(1.0 - F)*(x * y0) - C*(x * ydot0)) / (fr_3);
                STM(3, 5) = fr0 / this->mu*(xdot - xdot0)*(zdot - zdot0) + (fr0*(1.0 - F)*(x * z0) - C*(x * zdot0)) / (fr_3);
                STM(4, 3) = fr0 / this->mu*(ydot - ydot0)*(xdot - xdot0) + (fr0*(1.0 - F)*(y * x0) - C*(y * xdot0)) / (fr_3);
                STM(4, 4) = fr0 / this->mu*(ydot - ydot0)*(ydot - ydot0) + (fr0*(1.0 - F)*(y * y0) - C*(y * ydot0)) / (fr_3)+Gt;
                STM(4, 5) = fr0 / this->mu*(ydot - ydot0)*(zdot - zdot0) + (fr0*(1.0 - F)*(y * z0) - C*(y * zdot0)) / (fr_3);
                STM(5, 3) = fr0 / this->mu*(zdot - zdot0)*(xdot - xdot0) + (fr0*(1.0 - F)*(z * x0) - C*(z * xdot0)) / (fr_3);
                STM(5, 4) = fr0 / this->mu*(zdot - zdot0)*(ydot - ydot0) + (fr0*(1.0 - F)*(z * y0) - C*(z * ydot0)) / (fr_3);
                STM(5, 5) = fr0 / this->mu*(zdot - zdot0)*(zdot - zdot0) + (fr0*(1.0 - F)*(z * z0) - C*(z * zdot0)) / (fr_3)+Gt;

                //scale the stm->operator()
                //lower left
                for (int i = 3; i < 6; ++i)
                    for (int j = 0; j < 3; ++j)
                        STM(i, j) /= TU;
                //upper right
                for (int i = 0; i < 3; ++i)
                    for (int j = 3; j < 6; ++j)
                        STM(i, j) *= TU;

                //Step 5.4 compute Ftt and Gtt
                double rdot = fr0 * U0dot + sigma0 _GETVALUE * U1dot + U2dot;
                double r2 = fr*fr;
                Ftt = -sqmu / fr0 * (U1dot / fr - U1 _GETVALUE * rdot / r2);

                Gtt = -(U2dot / fr - U2 _GETVALUE * rdot / r2);

                //Step 5.4 scale F, G, Ft, Gt, Ftt, Gtt
                //F converts km to km so there is no scalig
                //G converts km/s to km so it must be multiplied by TU
                G *= TU;
                //Ft converts km to km/s so it must be divided by TU
                Ft /= TU;
                //Gt converts km/s to km/s so there is no scaling
                //Ftt converts km to km/s^2 so it must be divided by TU twice
                Ftt /= (TU * TU);
                //Gtt converts km/s to km/s^2 so it must be divided by TU
                Gtt /= TU;

                for (size_t stateIndex = 0; stateIndex < 3; ++stateIndex)
                {
                    dStatedIndependentVariable(stateIndex) = PropagationDirection * *this->dPropagationTime_dIndependentVariable * (Ft * StateLeft(stateIndex) + Gt * StateLeft(stateIndex + 3)) _GETVALUE;
                    dStatedIndependentVariable(stateIndex + 3) = PropagationDirection * *this->dPropagationTime_dIndependentVariable * (Ftt * StateLeft(stateIndex) + Gtt * StateLeft(stateIndex + 3)) _GETVALUE;
                }
            }//end derivatives
        }//end propagate()
    }//end namespace Astrodynamics
}//end namespace EMTG 