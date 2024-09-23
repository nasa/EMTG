//base Scalatron class
//for use with EMTG and CSALT
//Jacob Englander 8/15/2018

#pragma once

#include <vector>
#include <cmath>

namespace Scalatron
{
    enum F_EntryType { ObjectiveFunction, InequalityConstraint, EqualityConstraint, DefectConstraint };

    class ScalatronBase
    {
    public:
        //default constructor
        ScalatronBase() {};

        //constructor with arguments
        ScalatronBase(const std::vector<double>& X0,
                      const std::vector<double>& Xlowerbounds,
                      const std::vector<double>& Xupperbounds,
                      const std::vector<double>& F0,
                      const std::vector<double>& Flowerbounds,
                      const std::vector<double>& Fupperbounds,
                      const std::vector<size_t>& iGfun,
                      const std::vector<size_t>& jGvar,
                      const std::vector<double>& G0,
                      size_t ObjectiveIndex);

        //destructor
        virtual ~ScalatronBase() {}; //doesn't need to do anything right now, all members clean themselves up

        //initialize
        void initialize(const std::vector<double>& X0,
                        const std::vector<double>& Xlowerbounds,
                        const std::vector<double>& Xupperbounds,
                        const std::vector<double>& F0,
                        const std::vector<double>& Flowerbounds,
                        const std::vector<double>& Fupperbounds,
                        const std::vector<size_t>& iGfun,
                        const std::vector<size_t>& jGvar,
                        const std::vector<double>& G0,
                        size_t ObjectiveIndex);

        //define defect indices
        void defineDefectIndices(const std::vector<size_t>& DefectIndices);

        //compute scale factors
        virtual void compute_scale_factors() = 0;        

        //set/get
        std::vector<double> get_Kx() { return this->Kx; }
        std::vector<double> get_bx() { return this->bx; }
        std::vector<double> get_Kf() { return this->Kf; }
        std::vector<double> get_KF() { return this->Kf; }

    protected:
        //protected methods
        void classify_constraints();

        void compute_X_scaling();

		virtual void compute_objective_scaling() = 0;

        virtual void compute_equality_constraint_scaling();

		virtual void compute_inequality_constraint_scaling() = 0;

        //protected fields

        //copies of the original input data
        std::vector<double> X0;
        std::vector<double> Xlowerbounds, Xupperbounds, Flowerbounds, Fupperbounds;
        std::vector<double> F0;
        std::vector<size_t> iGfun;
        std::vector<size_t> jGvar;
        std::vector<double> G0;

        //book-keeping
        size_t nX;
        size_t nF;
        size_t nG;
        std::vector<F_EntryType> FTypes;
        size_t ObjectiveIndex;
        std::vector<size_t> EqualityConstraintIndices;
        std::vector<size_t> InequalityConstraintIndices;

        //scale factors to be written out
        std::vector<double> Kx;
        std::vector<double> bx;
        std::vector<double> Kf;
		std::vector<double> KG;
    };
}//end namespace Scalatron