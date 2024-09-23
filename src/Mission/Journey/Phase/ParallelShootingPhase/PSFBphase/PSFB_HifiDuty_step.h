//parallel shooting finite burn (PSFB) step for EMTGv9
//High-fidelity duty cycle
//Jacob Englander 3-19-2018

#pragma once

#include "doubleType.h"

#include "PSFBstep.h"

#include "IntegrationScheme.h"
#include "SpacecraftAccelerationModel.h"
#include "TimeDomainSpacecraftEOM.h"

namespace EMTG
{
    namespace Phases
    {
        //forward declare PSFBphase
        class PSFBphase;

        class PSFB_HifiDuty_step : virtual public PSFBstep
        {
        public:
            //constructor
            PSFB_HifiDuty_step();
            PSFB_HifiDuty_step(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stepIndex,
                const size_t& stageIndex,
                ParallelShootingStep* previousStep,
                ParallelShootingPhase* myPhase,
                Astrodynamics::universe* myUniverse,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            virtual void initialize(const std::string& name,
                const size_t& journeyIndex,
                const size_t& phaseIndex,
                const size_t& stepIndex,
                const size_t& stageIndex,
                ParallelShootingStep* previousStep,
                ParallelShootingPhase* myPhase,
                Astrodynamics::universe* myUniverse,
                HardwareModels::Spacecraft* mySpacecraft,
                missionoptions* myOptions);

            //clone
            virtual PSFB_HifiDuty_step* clone() const { return new PSFB_HifiDuty_step(*this); }

            //destructor
            virtual ~PSFB_HifiDuty_step();

            //output
            virtual void output(std::ofstream& outputfile,
                size_t& eventcount);

            virtual void output_ephemeris(std::ofstream& outputfile, std::ofstream& acceleration_model_file);
            
            virtual void output_STM(size_t& STMindex);

            virtual void output_maneuver_and_target_spec(std::ofstream& maneuver_spec_file, std::ofstream& target_spec_file, bool& haveManeuverNeedTarget);

            virtual void calcbounds_step();

            virtual void process_step(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

        protected:
            //configure propagator
            virtual void configure_propagator();

            //calcbounds
            virtual void calcbounds_deltav_contribution() {}; //empty for now

            //process

            virtual void process_step_main(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            virtual void process_deltav_contribution(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG) {};//empty for now
            
            void process_substep_times(const std::vector<doubleType>& X,
                size_t& Xindex,
                std::vector<doubleType>& F,
                size_t& Findex,
                std::vector<double>& G,
                const bool& needG);

            //fields
            math::Matrix<doubleType> StateAfterThrustInertial;
            Astrodynamics::PropagatorBase* myCoastPropagator;
            Integration::IntegrationScheme* myCoastIntegrationScheme;
            Astrodynamics::SpacecraftAccelerationModel* myCoastSpacecraftAccelerationModel;
            Astrodynamics::TimeDomainSpacecraftEOM myCoastEOM;
            std::vector< math::Matrix<double> > ThrustSTM;
            std::vector< math::Matrix<double> > ThrustAugmentedSTM;
            std::vector< math::Matrix<double> > ThrustCumulativeAugmentedSTM;
            std::vector< math::Matrix<double> > ThrustdPropagatedStatedIndependentVariable;
            math::Matrix<double> CoastSTM;
            math::Matrix<double> CoastAugmentedSTM;
            math::Matrix<double> CoastdPropagatedStatedIndependentVariable;
            doubleType StepThrustTime;
            doubleType subStepThrustTime;
            doubleType StepCoastTime;
            double dStepTime_dPhaseFlightTime_Coast;
            double dSubStepTime_dPhaseFlightTime;
        };

        inline PSFB_HifiDuty_step * new_clone(PSFB_HifiDuty_step const & other)
        {
            return other.clone();
        }
    }//close namespace Phases
}//close namespace EMTG