//PropulatorJourney

#include "missionoptions.h"
#include "journeyoptions.h"
#include "universe.h"
#include "Spacecraft.h"

#include "IntegrationScheme.h"
#include "SpacecraftAccelerationModel.h"
#include "TimeDomainSpacecraftEOM.h"
#include "PropagatorBase.h"

namespace EMTG
{
    namespace Propulator
    {
        class PropulatorJourney
        {
        public:
            PropulatorJourney();
            PropulatorJourney(const size_t& journeyIndex,
                size_t& stageIndex,
                missionoptions* options,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* spacecraft);

            ~PropulatorJourney();

            void initialize(const size_t& journeyIndex,
                size_t& stageIndex,
                missionoptions* options,
                Astrodynamics::universe* Universe,
                HardwareModels::Spacecraft* spacecraft);

            void setStepSize(const double& stepSize) { this->myPropagator->setPropagationStepSize(stepSize); }

            math::Matrix<double> getSTM() { return this->STM; }

            void propulate(const math::Matrix<double>& InputState,
                const math::Matrix<double>& ControlVector,
                const double& DutyCycle,
                const double& propagationTime,
                math::Matrix<double>& OutputState,
                const bool& needSTM);

            void propulate(const math::Matrix<double>& InputState,
                const double& propagationTime,
                math::Matrix<double>& OutputState,
                const bool& needSTM);

        protected:
            //structure
            size_t journeyIndex;
            size_t stageIndex;
            JourneyOptions* myJourneyOptions;
            Astrodynamics::universe* myUniverse;
            missionoptions* myOptions;
            HardwareModels::Spacecraft* mySpacecraft;

            //states
            math::Matrix<double> StateBeforePropagation;
            math::Matrix<double> StateAfterPropagation;
            math::Matrix<double> STM;
            math::Matrix<double> dPropagatedStatedIndependentVariable;

            //propagator and dynamics
            Astrodynamics::PropagatorBase* myPropagator;
            Astrodynamics::SpacecraftAccelerationModel* mySpacecraftAccelerationModel;
            Astrodynamics::TimeDomainSpacecraftEOM myEOM;
            Integration::IntegrationScheme* myIntegrationScheme;
            size_t total_number_of_states_to_integrate;
            double dPropagationTime_dPropagationTime;
        };
    }//close namespace Propulator
}//close namespace EMTG