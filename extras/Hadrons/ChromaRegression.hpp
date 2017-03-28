
#include <Grid/Grid.h>

#ifndef Hadrons_ChromaRegression_hpp_
#define Hadrons_ChromaRegression_hpp_
//#define CHROMA_REGRESSION
enum ChromaAction {
                 DWF,           // CPS style preconditioning
		 WilsonFermion, // Wilson
		 HwPartFracZolo, // KEK's approach
		 HwContFracZolo, // Edwards, Kennedy et al prefer this
		 HwPartFracTanh, // 
		 HwContFracTanh, // 
		 HwCayleyZolo, // Chiu Optimal
		 HtCayleyZolo, // 
		 HmCayleyZolo, // Scaled shamir 13
		 HwCayleyTanh, // Scaled shamir
		 HtCayleyTanh, // Plain old DWF.
		 HmCayleyTanh, // Scaled shamir 13
		 HtContFracTanh,
		 HtContFracZolo
};

struct DWF_parms
{
    int    Ls;
    double M5;
    double mq;
    double zolo_lo;
    double zolo_hi;
    double mobius_scale;
};

void init_chroma(const std::vector<int> &dimensions, int *argc, char ***argv);
void calc_chroma(ChromaAction action,Grid::QCD::LatticeGaugeField & lat, Grid::QCD::LatticeFermion &src, Grid::QCD::LatticeFermion &res, Grid::QCD::LatticeFermion &chroma_res, int dag, DWF_parms &dwf_par);

#endif
