//PropulatorJourney

#include "missionoptions.h"
#include "universe.h"
#include "Spacecraft.h"

#include "IntegrationSchemeFactory.h"
#include "PropagatorFactory.h"

#include "PropulatorJourney.h"

namespace EMTG
{
    namespace Propulator
    {
        PropulatorJourney::PropulatorJourney() :
            StateBeforePropagation(math::Matrix<double>(10, 1, 0.0)),
            StateAfterPropagation(math::Matrix<double>(10, 1, 0.0)),
            STM(math::Matrix<double>(14, 14, math::identity)),
            dPropagatedStatedIndependentVariable(math::Matrix<double>(10, 2, 0.0)),
            total_number_of_states_to_integrate(10),
            dPropagationTime_dPropagationTime(1.0)
        {}

        PropulatorJourney::PropulatorJourney(const size_t& journeyIndex,
            size_t& stageIndex,
            missionoptions* options,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* spacecraft) :
            PropulatorJourney()
        {
            this->initialize(journeyIndex, stageIndex, options, Universe, spacecraft);
        }

        void PropulatorJourney::initialize(const size_t& journeyIndex,
            size_t& stageIndex,
            missionoptions* options,
            Astrodynamics::universe* Universe,
            HardwareModels::Spacecraft* spacecraft)
        {
            this->journeyIndex = journeyIndex;

            //options, universe, and spacecraft
            this->myOptions = options;
            this->myJourneyOptions = &this->myOptions->Journeys[this->journeyIndex];
            this->myUniverse = Universe;
            this->mySpacecraft = spacecraft;

            //staging
            if (this->myJourneyOptions->stage_after_departure)
                ++stageIndex;

            this->stageIndex = stageIndex;

            if (this->myJourneyOptions->stage_before_arrival)
                ++stageIndex;

            //dynamics and propagator
            //acceleration model object
            std::vector<std::string> dummy_Xdescriptions;
            this->mySpacecraftAccelerationModel = new Astrodynamics::SpacecraftAccelerationModel(this->myOptions,
                this->myJourneyOptions,
                this->myUniverse,
                &dummy_Xdescriptions,
                this->mySpacecraft,
                14); // STM dimension
                
            this->mySpacecraftAccelerationModel->setDutyCycle(1.0);

            //EOM
            this->myEOM.setSpacecraftAccelerationModel(this->mySpacecraftAccelerationModel);

            //integration scheme
            this->myIntegrationScheme = CreateIntegrationScheme(&this->myEOM, 10, 14); // EOM, state dimension, STM dimension

            //propagator
            this->myPropagator = CreatePropagator(this->myOptions,
                this->myUniverse,
                10,
                14,
                this->StateBeforePropagation,
                this->StateAfterPropagation,
                this->STM,
                this->dPropagatedStatedIndependentVariable,
                (Integration::Integrand*) &this->myEOM,
                this->myIntegrationScheme,
                &this->dPropagationTime_dPropagationTime,
                this->myJourneyOptions->override_integration_step_size
                ? this->myJourneyOptions->integration_step_size
                : this->myOptions->integration_time_step_size);
        }//end initialize()

        PropulatorJourney::~PropulatorJourney()
        {
            delete this->myPropagator;
            delete this->mySpacecraftAccelerationModel;
            delete this->myIntegrationScheme;

        }//end destructor

        void PropulatorJourney::propulate(const math::Matrix<double>& InputState,
            const math::Matrix<double>& ControlVector,
            const double& DutyCycle,
            const double& propagationTime,
            math::Matrix<double>& OutputState,
            const bool& needSTM)
        {
            this->StateBeforePropagation.shallow_copy(InputState);

            //rotate to ICRF

            //propagate
            if (this->myOptions->duty_cycle_type == DutyCycleType::Averaged || DutyCycle > 0.999)
            {
                this->mySpacecraftAccelerationModel->setDutyCycle(DutyCycle);

                this->myPropagator->setCurrentEpoch(InputState(7));
                this->myPropagator->setIndexOfEpochInStateVec(7);
                this->myPropagator->setCurrentIndependentVariable(InputState(7));

                this->myPropagator->propagate(propagationTime, ControlVector, needSTM);
            }
            else
            {
                //propagate with thrust
                this->mySpacecraftAccelerationModel->setDutyCycle(1.0);

                this->myPropagator->setCurrentEpoch(InputState(7));
                this->myPropagator->setIndexOfEpochInStateVec(7);
                this->myPropagator->setCurrentIndependentVariable(InputState(7));

                this->myPropagator->propagate(propagationTime * DutyCycle, ControlVector, needSTM);

                math::Matrix<double> ThrustSTM = this->STM;

                //propagate without thrust
                this->StateBeforePropagation.shallow_copy(this->StateAfterPropagation);

                this->myPropagator->setCurrentEpoch(StateBeforePropagation(7));
                this->myPropagator->setIndexOfEpochInStateVec(7);
                this->myPropagator->setCurrentIndependentVariable(StateBeforePropagation(7));

                this->myPropagator->propagate(propagationTime * (1.0 - DutyCycle), needSTM);

                this->STM = ThrustSTM * this->STM;
            }

            //rotate from ICRF

            OutputState.shallow_copy(this->StateAfterPropagation);
        }

        void PropulatorJourney::propulate(const math::Matrix<double>& InputState,
            const double& propagationTime,
            math::Matrix<double>& OutputState,
            const bool& needSTM)
        {
            this->StateBeforePropagation = InputState;

            //propagate

            this->myPropagator->setCurrentEpoch(InputState(7));
            this->myPropagator->setIndexOfEpochInStateVec(7);
            this->myPropagator->setCurrentIndependentVariable(InputState(7));

            this->myPropagator->propagate(propagationTime, needSTM);

            OutputState = this->StateAfterPropagation;
        }
    }//close namespace Propulator
}//close namespace EMTG