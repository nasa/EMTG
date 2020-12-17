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

// aerodynamic drag term

#include "SpacecraftAccelerationModel.h"
#include "SpacecraftAccelerationModelTerm.h"
#include "doubleType.h"
#include "AerodynamicDragTerm.h"


namespace EMTG
{
	namespace Astrodynamics
	{
		// constructors
		AerodynamicDragTerm::AerodynamicDragTerm(SpacecraftAccelerationModel * acceleration_model_in) :
            SpacecraftAccelerationModelTerm::SpacecraftAccelerationModelTerm(acceleration_model_in)
		{
			this->sc_area = this->acceleration_model->my_journey_options->spacecraft_drag_area / (1000.0 * 1000.0); // convert from input (m^2) to internal units (km^2)
			this->Cd = this->acceleration_model->my_journey_options->coefficient_of_drag; // [0, inf)

			// resize class variables
			this->omega.resize(3, 1, 0.0);
			this->omega_BCI.resize(3, 1, 0.0);
			this->vRelative.resize(3, 1, 0.0);
			this->rBCF.resize(3, 1, 0.);
			this->dLLAdBCF.resize(3, 3, 0.);
			this->lla_output.resize(3, 1, 0.0);
			this->vAtmosphere.resize(3, 1, 0.0);

			// resize class variables
			this->DhDrBCF.resize(1, 3, 0.0); // derivative of altitude w.r.t. true-of-date body-fixed
			this->DrBCFDr.resize(3, 3, 0.0); // derivative of position of true-of-date body-fixed position w.r.t. ICRF position; same as transformation matrix from ICRF to true-of-date, body-fixed frame
			this->DhDr.resize(1, 3, 0.0); // derivative of altitude w.r.t. ICRF r
			this->Ddensity_kmDr.resize(1, 3, 0.0); // derivative of density w.r.t. ICRF r
			this->Ddensity_kmDt.resize(1, 1, 0.0); // derivative of density w.r.t. ICRF t
			this->DaccDstate.resize(3, 8, 0.0); // derivative of acceleration vector w.r.t. [r, v, m, t]
			this->DmInvDstate.resize(1, 8, 0.0); // derivative of 1/mass w.r.t. [r, v, m, t]
			this->Ddensity_kmDstate.resize(1, 8, 0.0); // derivative of density w.r.t. [r, v, m, t]
			this->DvRelativeMagDstate.resize(1, 8, 0.0); // derivative of norm(vRelative) w.r.t. [r, v, m, t]
			this->DvRelativeDstate.resize(3, 8, 0.0); // derivative of vRelative (ICRF) w.r.t. [r, v, m, t]
			this->DrDstate.resize(3, 8, 0.0); // derivative of r (ICRF) w.r.t. [r, v, m, t]
			this->DvDstate.resize(3, 8, 0.0); // derivative of v (ICRF) w.r.t. [r, v, m, t]
			this->DcrossTermDx.resize(3, 8, 0.0); // derivative of omega \cross r w.r.t. [r, v, m, t]
			this->DrBCFDt.resize(3, 1, 0.0); // derivative of r TOD BCF w.r.t t
			this->vRelativeDouble.resize(3, 1, 0.0); // double version of vRelative
			this->DaccDpropVars.resize(3, 2, 0.0); // derivative of acceleration w.r.t. prop vars
			this->DomegaDstate.resize(3, 8, 0.0); // derivative of omega w.r.t. [r, v, m, t]
			this->DrBCFDt_doubleType.resize(3, 1, 0.0); // derivative of r TOD BCF w.r.t t
			this->dICRFdBCI.resize(3, 3, 0.0); // rotation matrix  from BCI to ICRF
			this->referenceVector.resize(3, 1, 0.0); // don't care
			this->dreferenceVector_dt.resize(3, 1, 0.0); // don't care
			this->DomegaDt.resize(3, 1, 0.0); // derivative of rotated vector w.r.t. time
			this->dRstate_dreferenceVector.resize(3, 3, 0.0); // don't care
			this->DrBCFDr_doubleType.resize(3, 3, 0.0); // derivative of r in the BCF TOD frame w.r.t. r ICRF as a doubleType
			this->temp_matrix_out.resize(3, 1, 0.0);
		}

		AerodynamicDragTerm::~AerodynamicDragTerm() {}

		// methods
		void AerodynamicDragTerm::computeAccelerationTerm()
		{
			// zero-out class variables
			this->omega.assign_zeros();
			this->omega_BCI.assign_zeros();
			this->vRelative.assign_zeros();
			this->rBCF.assign_zeros();
			this->dLLAdBCF.assign_zeros();
			this->lla_output.assign_zeros();
			this->vAtmosphere.assign_zeros();

			// get velocity of atmosphere here using omega \times r
			// need to get omega
			// this comes from the universe class and is about the central body's pole
			std::vector<doubleType> central_body_reference_angles;
			central_body_reference_angles.resize(6);

			// get the reference angles in rad and rad/s
			doubleType alphadot = this->acceleration_model->my_universe->LocalFrame.getAlphadot();
			alphadot = alphadot / 36525. / 86400.; // convert from rad/century to rad/sec
			doubleType deltadot = this->acceleration_model->my_universe->LocalFrame.getDeltadot();
			deltadot = deltadot / 36525. / 86400.; // convert from rad/century to rad/sec
			doubleType Wdot = this->acceleration_model->my_universe->LocalFrame.getWdot();
			Wdot = Wdot / 86400.; // convert from rad/day to rad/sec
			doubleType delta = this->acceleration_model->my_universe->LocalFrame.getDelta(this->acceleration_model->current_epoch); // expects MJD seconds

			// construct angular velocity vector of true of date BCF w.r.t. ICRF, expressed in TOD BCI frame
			this->omega_BCI(0) = deltadot;
			this->omega_BCI(1) = -sin(delta - math::PIover2) * alphadot;
			this->omega_BCI(2) = Wdot + cos(delta - math::PIover2) * alphadot;

			// convert to icrf. we are starting with omega in BCI true of date
			this->acceleration_model->my_universe->LocalFrame.rotate_frame_to_frame(EMTG::ReferenceFrame::TrueOfDate_BCI,
				this->omega_BCI,
				EMTG::ReferenceFrame::ICRF,
				this->omega,
				this->acceleration_model->current_epoch);

			// cross product
			this->vAtmosphere = omega.cross(this->acceleration_model->r_cb2sc); // this is vAtmosphere in icrf

			// relative velocity is inertial velocity of spacecraft - velocity of atmosphere
			for (size_t k = 0; k < 3; ++k)
			{
				this->vRelative(k) = this->acceleration_model->v_cb2sc(k) - this->vAtmosphere(k); // icrf
			}
			this->vRelativeMag = vRelative.norm();

			// get the atmospheric density. first, need the body-fixed position vector 
			this->acceleration_model->my_universe->LocalFrame.rotate_frame_to_frame(EMTG::ReferenceFrame::ICRF,
				this->acceleration_model->r_cb2sc,
				EMTG::ReferenceFrame::TrueOfDate_BCF,
				this->rBCF,
				this->acceleration_model->current_epoch);

			// now, get altitude
			doubleType Re = this->acceleration_model->my_universe->central_body_radius;
			doubleType f = this->acceleration_model->my_universe->central_body.flattening_coefficient;
			BCF2LLA_oblate(this->rBCF, Re, f, this->lla_output, true, dLLAdBCF);
			doubleType altitude = this->lla_output(2);

			// now, call the current atmosphere
			this->acceleration_model->my_universe->TheAtmosphere->computeDensity(altitude, true);
			this->density = this->acceleration_model->my_universe->TheAtmosphere->getDensity(); // kg / m^3
			doubleType density_km = this->density * 1.0e9; // kg / km^3

			// compute the drag force and acceleration
			doubleType coef = -0.5 * (density_km) * this->Cd * this->sc_area * (this->vRelativeMag); // km/s
			doubleType dragForce[3];

			for (size_t k = 0; k < 3; ++k)
			{
				dragForce[k] = coef * (this->vRelative(k));
				this->term_acceleration(k) = (dragForce[k]) / this->acceleration_model->spacecraft_mass;
				this->acceleration_model->acceleration(k) += this->term_acceleration(k);
			}
		}// end computeAccelerationTerm()

		void AerodynamicDragTerm::computeAccelerationTerm(const bool & generate_derivatives)
		{
			// zero-out class variables
			this->DhDrBCF.assign_zeros(); // derivative of altitude w.r.t. true-of-date body-fixed
			this->DrBCFDr.assign_zeros(); // derivative of position of true-of-date body-fixed position w.r.t. ICRF position; same as transformation matrix from ICRF to true-of-date, body-fixed frame
			this->DhDr.assign_zeros(); // derivative of altitude w.r.t. ICRF r
			this->Ddensity_kmDr.assign_zeros(); // derivative of density w.r.t. ICRF r
			this->Ddensity_kmDt.assign_zeros(); // derivative of density w.r.t. ICRF t
			this->DaccDstate.assign_zeros(); // derivative of acceleration vector w.r.t. [r, v, m, t]
			this->DmInvDstate.assign_zeros(); // derivative of 1/mass w.r.t. [r, v, m, t]
			this->Ddensity_kmDstate.assign_zeros(); // derivative of density w.r.t. [r, v, m, t]
			this->DvRelativeMagDstate.assign_zeros(); // derivative of norm(vRelative) w.r.t. [r, v, m, t]
			this->DvRelativeDstate.assign_zeros(); // derivative of vRelative (ICRF) w.r.t. [r, v, m, t]
			this->DrDstate.assign_zeros(); // derivative of r (ICRF) w.r.t. [r, v, m, t]
			this->DvDstate.assign_zeros(); // derivative of v (ICRF) w.r.t. [r, v, m, t]
			this->DcrossTermDx.assign_zeros(); // derivative of omega \cross r w.r.t. [r, v, m, t]
			this->DrBCFDt.assign_zeros(); // derivative of r TOD BCF w.r.t t
			this->vRelativeDouble.assign_zeros(); // double version of vRelative
			this->DaccDpropVars.assign_zeros(); // derivative of acceleration w.r.t. prop vars
			this->DomegaDstate.assign_zeros(); // derivative of omega w.r.t. [r, v, m, t]
			this->DrBCFDt_doubleType.assign_zeros(); // derivative of r TOD BCF w.r.t t
			this->dICRFdBCI.assign_zeros(); // rotation matrix  from BCI to ICRF
			this->referenceVector.assign_zeros(); // don't care
			this->dreferenceVector_dt.assign_zeros(); // don't care
			this->DomegaDt.assign_zeros(); // derivative of omega in ICRF w.r.t. t
			this->dRstate_dreferenceVector.assign_zeros(); // don't care
			this->DrBCFDr_doubleType.assign_zeros(); // derivative of r in the BCF TOD frame w.r.t. r ICRF as a doubleType
			this->temp_matrix_out.assign_zeros();
			
			// preliminaries
			double spacecraft_mass = this->acceleration_model->spacecraft_mass _GETVALUE; // extract mass
			double mInv = 1. / spacecraft_mass;
			this->computeAccelerationTerm(); // call the version of the routine that gets the acceleration itself
			double density_km = (this->density * 1.0e9) _GETVALUE; // kg / km^3

			// extract derivative of density w.r.t. altitude
			this->DdensityDh = this->acceleration_model->my_universe->TheAtmosphere->getDdensityDh(); // d [density] / d [h] in [kg / m^3] / [km]
			double Ddensity_kmDh = (this->DdensityDh * 1.0e9) _GETVALUE; // [kg / km^3] / [km]

			// dr / dr is identity
			this->DrDstate(0, 0) = 1.;
			this->DrDstate(1, 1) = 1.;
			this->DrDstate(2, 2) = 1.;

			// dv / dv is identity
			this->DvDstate(0, 3) = 1.;
			this->DvDstate(1, 4) = 1.;
			this->DvDstate(2, 5) = 1.;

			// d omega / d state only has nonzero components for time
			double alphadot = this->acceleration_model->my_universe->LocalFrame.getAlphadot() _GETVALUE;
			alphadot = alphadot / 36525. / 86400.; // convert from rad/century to rad/sec
			double deltadot = this->acceleration_model->my_universe->LocalFrame.getDeltadot() _GETVALUE;
			deltadot = deltadot / 36525. / 86400.; // convert from rad/century to rad/sec
			double delta = this->acceleration_model->my_universe->LocalFrame.getDelta(this->acceleration_model->current_epoch) _GETVALUE; // expects seconds MJD
			this->DomegaDstate(1, 7) = -alphadot * deltadot * cos(delta - math::PIover2); // BCI true of date
			this->DomegaDstate(2, 7) = -alphadot * deltadot * sin(delta - math::PIover2); // BCI true of date

			// convert from BCI to ICRF, including time derivative
			this->acceleration_model->my_universe->LocalFrame.construct_rotation_matrices(this->acceleration_model->current_epoch, true); // initialize the rotation matrices
			
			this->acceleration_model->my_universe->LocalFrame.rotate_frame_to_frame(
				EMTG::ReferenceFrame::TrueOfDate_BCI,
				this->omega_BCI, // vector in BCI
				this->DomegaDstate.getcolumn(7), // d/dt of vector in BCI
				this->referenceVector,
				this->dreferenceVector_dt,
				EMTG::ReferenceFrame::ICRF,
				this->temp_matrix_out, // vector in ICRF (already contained in this->omega)
				this->dICRFdBCI, // R of BCI->ICRF
				this->DomegaDt, // d/dt of vector in ICRF
				this->dRstate_dreferenceVector,
				this->acceleration_model->current_epoch, // expects seconds MJD
				true);
			for (size_t k = 0; k < 3; ++k)
			{
				this->DomegaDstate(k, 7) = this->DomegaDt(k, 0); // DomegaDstate now contains ICRF time derivative
			}

			// derivative of d [omega X r] / d state
			for (size_t k = 0; k < 8; ++k)
			{
				this->DcrossTermDx(0, k) = (this->DomegaDstate(1, k) * this->acceleration_model->r_cb2sc(2) +
					this->omega(1) * this->DrDstate(2, k) -
					this->DomegaDstate(2, k) * this->acceleration_model->r_cb2sc(1) -
					this->omega(2) * this->DrDstate(1, k)) _GETVALUE;
				this->DcrossTermDx(1, k) = (this->DomegaDstate(2, k) * this->acceleration_model->r_cb2sc(0) +
					this->omega(2) * this->DrDstate(0, k) -
					this->DomegaDstate(0, k) * this->acceleration_model->r_cb2sc(2) -
					this->omega(0) * this->DrDstate(2, k)) _GETVALUE;
				this->DcrossTermDx(2, k) = (this->DomegaDstate(0, k) * this->acceleration_model->r_cb2sc(1) +
					this->omega(0) * this->DrDstate(1, k) -
					this->DomegaDstate(1, k) * this->acceleration_model->r_cb2sc(0) -
					this->omega(1) * this->DrDstate(0, k)) _GETVALUE;
			}

			// d [vRel] / d state
			this->DvRelativeDstate = this->DvDstate - this->DcrossTermDx; // assumes zero wind velocity

			// d [norm(vRel)] / d state
			for (size_t k = 0; k < 3; ++k)
			{
				this->vRelativeDouble(k) = this->vRelative(k) _GETVALUE;
			}
			this->DvRelativeMagDstate = this->vRelativeDouble.transpose() * this->DvRelativeDstate * (1. / this->vRelativeMag _GETVALUE);

			// DhDrBCF comes from derivatives of LLA state w.r.t. BCF true of date state
			for (size_t k = 0; k < 3; ++k)
			{
				this->DhDrBCF(k) = this->dLLAdBCF(2, k) _GETVALUE;
			}

			// d [rBCFTOD] / d [rICRF] is the rotation matrix
			this->acceleration_model->my_universe->LocalFrame.rotate_frame_to_frame(
				EMTG::ReferenceFrame::ICRF,
				this->acceleration_model->r_cb2sc, // vector in ICRF
				this->referenceVector, // d/dt of vector in ICRF is 0 (position state is not explicit function of time)
				this->referenceVector,
				this->dreferenceVector_dt,
				EMTG::ReferenceFrame::TrueOfDate_BCF,
				this->temp_matrix_out, // vector in BCF
				this->DrBCFDr_doubleType, // R of ICRF->BCF
				this->DrBCFDt_doubleType, // d/dt of vector in BCF
				this->dRstate_dreferenceVector,
				this->acceleration_model->current_epoch, // expects seconds MJD
				true);

			for (size_t k = 0; k < 3; ++k)
			{
				this->DrBCFDt(k, 0) = this->DrBCFDt_doubleType(k, 0) _GETVALUE;
				for (size_t j = 0; j < 3; ++j)
				{
					this->DrBCFDr(k, j) = this->DrBCFDr_doubleType(k, j) _GETVALUE;
				}
			}
			this->DhDr = this->DhDrBCF * this->DrBCFDr; // d [h] / d [rICRF]
			this->Ddensity_kmDr = this->DhDr * Ddensity_kmDh; // d [density] / d [rICRF]

			this->Ddensity_kmDt = this->DhDrBCF * this->DrBCFDt * Ddensity_kmDh; // d [density] / d [t]

			// d [density] / d state and d h / d state
			for (size_t k = 0; k < 3; ++k)
			{
				this->Ddensity_kmDstate(k) = this->Ddensity_kmDr(k); // the only state element that density depends on is rICRF
			}
			this->Ddensity_kmDstate(7) = this->Ddensity_kmDt(0);

			double constCoeff = (-0.5 * this->Cd * this->sc_area) _GETVALUE; // coefficient in km^2

			this->DmInvDstate(6) = -pow(spacecraft_mass, -2.); // d [m^-1] / d [m]
			// acc in km/s^2; state in km, km/s, kg
			this->DaccDstate = this->vRelativeDouble * this->DmInvDstate * density_km * this->vRelativeMag _GETVALUE +
				this->vRelativeDouble * this->Ddensity_kmDstate * mInv * this->vRelativeMag _GETVALUE +
				this->vRelativeDouble * this->DvRelativeMagDstate * mInv * density_km +
				this->DvRelativeDstate * this->vRelativeMag _GETVALUE * mInv * density_km;
			this->DaccDstate *= constCoeff; // multiply by the coefficient

			// populate the state propagation matrix:

			///////////////////////////
			// d a / d r: fx(3:5, 0:2)
			//////////////////////////
			for (size_t k = 0; k < 3; ++k)
			{
				for (size_t j = 0; j < 3; ++j)
				{
					this->acceleration_model->fx(k + 3, j) += this->DaccDstate(k, j);
				}
			}

			///////////////////////////
			// d a / d v: fx(3:5, 3:5)
			///////////////////////////
			for (size_t k = 0; k < 3; ++k)
			{
				for (size_t j = 0; j < 3; ++j)
				{
					this->acceleration_model->fx(k + 3, j + 3) += this->DaccDstate(k, j + 3);
				}
			}

			/////////////////////////
			// d a / d m: fx(3:5, 6)
			/////////////////////////
			for (size_t k = 0; k < 3; ++k)
			{
				this->acceleration_model->fx(k + 3, 6) += this->DaccDstate(k, 6);
			}

			/////////////////////////
			// d a / d t: fx(3:5, 7)
			/////////////////////////
			for (size_t k = 0; k < 3; ++k)
			{
				this->acceleration_model->fx(k + 3, 7) += this->DaccDstate(k, 7);
			}

			///////////////////////
			// d a / d [prop vars]
			///////////////////////
			
            /*
			for (size_t i = 0; i < 3; ++i)
			{
				this->DaccDpropVars(i, 0) = this->DaccDstate(i, 7) * this->acceleration_model->dcurrent_epochdProp_var_previous;
				this->DaccDpropVars(i, 1) = this->DaccDstate(i, 7) * this->acceleration_model->dcurrent_epochdProp_var;
			}

			for (size_t propVar = 0; propVar < 2; ++propVar)
			{
				this->acceleration_model->da_cb2scdPropVars(0, propVar) += this->DaccDpropVars(0, propVar);

				this->acceleration_model->da_cb2scdPropVars(1, propVar) += this->DaccDpropVars(1, propVar);

				this->acceleration_model->da_cb2scdPropVars(2, propVar) += this->DaccDpropVars(2, propVar);
			}
            */
		} // end computeAccelerationTerm(bool)

		void AerodynamicDragTerm::populateInstrumentationFile(std::ofstream & acceleration_model_file)
		{
			acceleration_model_file << "," << sqrt(this->term_acceleration(0)*this->term_acceleration(0)
				+ this->term_acceleration(1)*this->term_acceleration(1)
				+ this->term_acceleration(2)*this->term_acceleration(2));
			for (size_t k = 0; k < this->term_acceleration.get_n(); ++k)
			{
				acceleration_model_file << "," << this->term_acceleration(k);
			}
		}
	} // close namespace Astrodynamics
} // close namespace EMTG