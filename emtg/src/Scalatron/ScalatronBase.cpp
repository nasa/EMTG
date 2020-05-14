//base Scalatron class
//for use with EMTG and CSALT
//Jacob Englander 8/15/2018

#include <cmath>
#include "ScalatronBase.h"

namespace Scalatron
{
    ScalatronBase::ScalatronBase(const std::vector<double>& X0,
        const std::vector<double>& Xlowerbounds,
        const std::vector<double>& Xupperbounds,
        const std::vector<double>& F0,
        const std::vector<double>& Flowerbounds,
        const std::vector<double>& Fupperbounds,
        const std::vector<size_t>& iGfun,
        const std::vector<size_t>& jGvar,
        const std::vector<double>& G0,
        size_t ObjectiveIndex)
        : ScalatronBase()
    {
        this->initialize(X0,
            Xlowerbounds,
            Xupperbounds,
            F0,
            Flowerbounds,
            Fupperbounds,
            iGfun,
            jGvar,
            G0, 
            ObjectiveIndex);
    }//end constructor

    void ScalatronBase::initialize(const std::vector<double>& X0,
        const std::vector<double>& Xlowerbounds,
        const std::vector<double>& Xupperbounds,
        const std::vector<double>& F0,
        const std::vector<double>& Flowerbounds,
        const std::vector<double>& Fupperbounds,
        const std::vector<size_t>& iGfun,
        const std::vector<size_t>& jGvar,
        const std::vector<double>& G0,
        size_t ObjectiveIndex)
    {
        this->X0 = X0;
        this->Xlowerbounds = Xlowerbounds;
        this->Xupperbounds = Xupperbounds;
        this->F0 = F0;
        this->Flowerbounds = Flowerbounds;
        this->Fupperbounds = Fupperbounds;
        this->iGfun = iGfun;
        this->jGvar = jGvar;
        this->G0 = G0;

        this->nX = this->X0.size();
        this->nF = this->F0.size();
        this->nG = this->G0.size();

        this->FTypes.resize(nF);
        this->Kx.resize(this->nX, 1.0);
        this->bx.resize(this->nX, 0.0);
        this->Kf.resize(this->nF, 1.0);
        
        this->classify_constraints();

        this->FTypes[ObjectiveIndex] = F_EntryType::ObjectiveFunction;
    }//end initialize()

    void ScalatronBase::classify_constraints()
    {
        for (size_t Findex = 0; Findex < this->nF; ++Findex)
        {
            if (std::fabs(this->Fupperbounds[Findex] - this->Flowerbounds[Findex]) < 1.0e-8)
            {
                this->FTypes[Findex] = F_EntryType::EqualityConstraint;
                this->EqualityConstraintIndices.push_back(Findex);
            }
            else
            {
                this->FTypes[Findex] = F_EntryType::InequalityConstraint;
                this->InequalityConstraintIndices.push_back(Findex);
            }
        }//end loop over constraints
    }//end classify_constraints()

    void ScalatronBase::defineDefectIndices(const std::vector<size_t>& DefectIndices)
    {
        for (size_t DefectIndex : DefectIndices)
        {
            this->FTypes[DefectIndex] = F_EntryType::DefectConstraint;
        }
    }//end defineDefectIndices()

    void ScalatronBase::compute_scale_factors()
    {
        //This method is abstract, so derived classes MUST define it. But we also define some default behavior here that derived classes will want to be able to access.
        
        this->compute_X_scaling();
    }//end compute_scale_factors

    void ScalatronBase::compute_X_scaling()
    {
        //from Sagliano, "Performance analysis of linear and nonlinear techniques for automatic scaling of discretized control problems"
        //Operations Research Letters 42 (2014) 213-216
        //https://core.ac.uk/download/pdf/31014768.pdf
        //Equation 4

        for (size_t Xindex = 0; Xindex < this->nX; ++Xindex)
        {
            this->Kx[Xindex] = 1.0 / (this->Xupperbounds[Xindex] - this->Xlowerbounds[Xindex]);
            
            this->bx[Xindex] = -this->Xlowerbounds[Xindex] * this->Kx[Xindex];
        }//end loop over states
    }//compute_X_scaling()

    void ScalatronBase::compute_equality_constraint_scaling()
    {
        //PJRN method
        //from Sagliano, "Performance analysis of linear and nonlinear techniques for automatic scaling of discretized control problems"
        //Operations Research Letters 42 (2014) 213-216
        //https://core.ac.uk/download/pdf/31014768.pdf
        //Equation 7

        for (size_t Findex : this->EqualityConstraintIndices)
        {
			double Kf_thisConstraint = 0.0;

			for (size_t Gindex = 0; Gindex < this->nG; ++Gindex)
			{
				if (this->iGfun[Gindex] == Findex)
				{
					size_t Xindex = this->jGvar[Gindex];

					double candidate_Kf = std::fabs(this->G0[Gindex] / this->Kx[Xindex]);

					Kf_thisConstraint = (candidate_Kf > Kf_thisConstraint) ? candidate_Kf : Kf_thisConstraint;
				}
			}//end loop over Jacobian entries

			this->Kf[Findex] = Kf_thisConstraint;
        }//end loop over constraints
    }//end compute_equality_constraint_scaling
}//end namespace Scalatron