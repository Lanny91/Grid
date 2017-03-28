
#include "ChromaRegression.hpp"
#include <chroma.h>
#include <actions/ferm/invert/syssolver_linop_cg_array.h>
#include <actions/ferm/invert/syssolver_linop_aggregate.h>

namespace Chroma { 

class ChromaWrapper {
public:

  
  typedef multi1d<LatticeColorMatrix> U;
  typedef LatticeFermion T4;
  typedef multi1d<LatticeFermion> T5;
  
  static void ImportGauge(Grid::QCD::LatticeGaugeField & gr,
			  QDP::multi1d<QDP::LatticeColorMatrix> & ch) 
  {
    Grid::QCD::LorentzColourMatrix LCM;
    Grid::Complex cc;
    QDP::ColorMatrix cm;
    QDP::Complex c;

    std::vector<int> x(4);
    QDP::multi1d<int> cx(4);
    std::vector<int> gd= gr._grid->GlobalDimensions();

    for (x[0]=0;x[0]<gd[0];x[0]++){
    for (x[1]=0;x[1]<gd[1];x[1]++){
    for (x[2]=0;x[2]<gd[2];x[2]++){
    for (x[3]=0;x[3]<gd[3];x[3]++){
      cx[0] = x[0];
      cx[1] = x[1];
      cx[2] = x[2];
      cx[3] = x[3];
      Grid::peekSite(LCM,gr,x);

      for(int mu=0;mu<4;mu++){
	for(int i=0;i<3;i++){
	for(int j=0;j<3;j++){
	  cc = LCM(mu)()(i,j);
	  c = QDP::cmplx(QDP::Real(real(cc)),QDP::Real(imag(cc)));
	  QDP::pokeColor(cm,c,i,j);
	}}
	QDP::pokeSite(ch[mu],cm,cx);
	/*
	std::cout << "("<<x[0]<<",";
	std::cout << x[1]<<",";
	std::cout << x[2]<<",";
	std::cout << x[3]<<") "<< Grid::norm2(LCM(mu)) << " " <<QDP::norm2(cm)<<std::endl ;
	*/
      }

    }}}}
  }
  
  static void ImportFermion(Grid::QCD::LatticeFermion & gr,
			    QDP::multi1d<QDP::LatticeFermion> & ch  ) 
  {
    Grid::QCD::SpinColourVector F;
    Grid::Complex c;

    QDP::Fermion cF;
    QDP::SpinVector cS;
    QDP::Complex cc;

    std::vector<int> x(5);
    QDP::multi1d<int> cx(4);
    std::vector<int> gd= gr._grid->GlobalDimensions();

    for (x[0]=0;x[0]<gd[0];x[0]++){
    for (x[1]=0;x[1]<gd[1];x[1]++){
    for (x[2]=0;x[2]<gd[2];x[2]++){
    for (x[3]=0;x[3]<gd[3];x[3]++){
    for (x[4]=0;x[4]<gd[4];x[4]++){
      int s = x[0];
      cx[0] = x[1];
      cx[1] = x[2];
      cx[2] = x[3];
      cx[3] = x[4];

      Grid::peekSite(F,gr,x);

      for(int j=0;j<3;j++){
	for(int sp=0;sp<4;sp++){

	  c= F()(sp)(j) ;

	  cc = QDP::cmplx(QDP::Real(real(c)),QDP::Real(imag(c)));

	  QDP::pokeSpin(cS,cc,sp);

	}
	QDP::pokeColor(cF,cS,j);
      }
      QDP::pokeSite(ch[s],cF,cx);
    }}}}}
  }
  static void ExportFermion(Grid::QCD::LatticeFermion & gr,
			    QDP::multi1d<QDP::LatticeFermion> & ch  ) 
  {
    Grid::QCD::SpinColourVector F;
    Grid::Complex c;

    QDP::Fermion cF;
    QDP::SpinVector cS;
    QDP::Complex cc;

    std::vector<int> x(5);
    QDP::multi1d<int> cx(4);
    std::vector<int> gd= gr._grid->GlobalDimensions();

    for (x[0]=0;x[0]<gd[0];x[0]++){
    for (x[1]=0;x[1]<gd[1];x[1]++){
    for (x[2]=0;x[2]<gd[2];x[2]++){
    for (x[3]=0;x[3]<gd[3];x[3]++){
    for (x[4]=0;x[4]<gd[4];x[4]++){
      int s = x[0];
      cx[0] = x[1];
      cx[1] = x[2];
      cx[2] = x[3];
      cx[3] = x[4];

      cF = QDP::peekSite(ch[s],cx);
      for(int sp=0;sp<4;sp++){
	for(int j=0;j<3;j++){
	  cS =QDP::peekColor(cF,j);
	  cc =QDP::peekSpin(cS,sp);
	  c = Grid::Complex(QDP::toDouble(QDP::real(cc)), 
			    QDP::toDouble(QDP::imag(cc)));
	  F()(sp)(j) = c;
	}
      }
      Grid::pokeSite(F,gr,x);
    }}}}}
  }

  static Handle< LinearOperatorArray<T4> >  GetLinOp (U &u, ChromaAction parms, DWF_parms & dwf_par)
  {
    QDP::Real _mq(dwf_par.mq);
    QDP::Real eps_lo(dwf_par.zolo_lo);
    QDP::Real eps_hi(dwf_par.zolo_hi);
    QDP::Real scale(dwf_par.mobius_scale);
    QDP::Real M5(dwf_par.M5);
    int Ls = dwf_par.Ls;

    QDP::multi1d<int> bcs(QDP::Nd);

    bcs[0] = bcs[1] = bcs[2] = bcs[3] = 1;

    Chroma::Handle<Chroma::FermBC<T4,U,U> > fbc(new Chroma::SimpleFermBC< T4, U, U >(bcs));
    Chroma::Handle<Chroma::CreateFermState<T4,U,U> > cfs( new Chroma::CreateSimpleFermState<T4,U,U>(fbc));


    Chroma::GroupXML_t invparm;
    invparm.xml=std::string(
"   <InvertParam>\n"
"   <invType>CG_INVERTER</invType>\n"
"   <RsdCG>1.0e-9</RsdCG>\n"
"   <MaxCG>3000</MaxCG>\n"
"   </InvertParam>"
);

    invparm.id=std::string("CG_INVERTER");
    invparm.path=std::string("/InvertParam");

    if ( (parms == HtCayleyTanh)|| (parms==DWF) ) {
      Chroma::UnprecDWFermActArray  S_f(cfs, M5, _mq, Ls);
      Chroma::Handle< Chroma::FermState<T4,U,U> > fs( S_f.createState(u) );
      Chroma::Handle< Chroma::LinearOperatorArray<T4> > M(S_f.unprecLinOp(fs,_mq));
     return  M;
    }
    if ( parms == HwCayleyTanh ) {
      QDP::Real b5 = 0.5;
      QDP::Real c5 = 0.5;
      Chroma::UnprecNEFFermActArray  S_f(cfs, M5,b5,c5, _mq, Ls);
      Chroma::Handle< Chroma::FermState<T4,U,U> > fs( S_f.createState(u) );
      Chroma::Handle< Chroma::LinearOperatorArray<T4> > M(S_f.unprecLinOp(fs,_mq));
      return  M;
    }
    if ( parms == HmCayleyTanh ) {
      Real b5 = 0.5*(scale +1.0);
      Real c5 = 0.5*(scale -1.0);
      UnprecNEFFermActArray  S_f(cfs, M5,b5,c5, _mq, Ls);
      Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
      Handle< LinearOperatorArray<T4> > M(S_f.unprecLinOp(fs,_mq));
      return  M;
    }
   if ( parms == HwCayleyZolo ) {
     UnprecZoloNEFFermActArrayParams params;
     params.OverMass=M5;
     params.Mass=_mq;
     params.b5=0.5;
     params.c5=0.5;
     params.N5=Ls;
     params.approximation_type = COEFF_TYPE_ZOLOTAREV;
     params.ApproxMin=eps_lo;
     params.ApproxMax=eps_hi;
     UnprecZoloNEFFermActArray  S_f(cfs, params);
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperatorArray<T4> > M(S_f.unprecLinOp(fs,_mq));
     return  M;
   }
   if ( parms == HtCayleyZolo ) {
     UnprecZoloNEFFermActArrayParams params;
     params.OverMass=M5;
     params.Mass=_mq;
     params.b5=1.0;
     params.c5=0.0;
     params.N5=Ls;
     params.approximation_type = COEFF_TYPE_ZOLOTAREV;
     params.ApproxMin=eps_lo;
     params.ApproxMax=eps_hi;
     UnprecZoloNEFFermActArray  S_f(cfs, params);
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperatorArray<T4> > M(S_f.unprecLinOp(fs,_mq));
     return M;
   }
   if ( parms == HmCayleyZolo ) {
     UnprecZoloNEFFermActArrayParams params;
     params.OverMass=M5;
     params.Mass=_mq;
     params.b5= 0.5*(scale +1.0);
     params.c5= 0.5*(scale -1.0);
     params.N5=Ls;
     params.approximation_type = COEFF_TYPE_ZOLOTAREV;
     params.ApproxMin=eps_lo;
     params.ApproxMax=eps_hi;
     UnprecZoloNEFFermActArray  S_f(cfs, params);
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperatorArray<T4> > M(S_f.unprecLinOp(fs,_mq));
     return M;
   }
   if ( parms == HwPartFracZolo ) {
     if ( Ls%2 == 0 ) { 
       printf("Ls is not odd\n");
       exit(-1);
     }
     UnprecOvExtFermActArrayParams param;
     param.OverMass=M5; 
     param.Mass=_mq;
     param.RatPolyDeg = Ls;
     param.ApproxMin =eps_lo;
     param.ApproxMax =eps_hi;
     param.b5 =1.0;
     param.c5 =1.0;
     param.approximation_type=COEFF_TYPE_ZOLOTAREV;
     //     param.approximation_type=COEFF_TYPE_TANH_UNSCALED;
     //     param.approximation_type=COEFF_TYPE_TANH;
     param.tuning_strategy_xml=
"<TuningStrategy><Name>OVEXT_CONSTANT_STRATEGY</Name></TuningStrategy>\n";
     UnprecOvExtFermActArray S_f(cfs,param);
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperatorArray<T4> > M(S_f.linOp(fs));
     return M;
   }
   if ( parms == HwContFracZolo ) {
     UnprecOvlapContFrac5DFermActParams param;
     param.Mass=_mq; // How is M5 set? Wilson mass In AuxFermAct
     param.ApproxMin=eps_lo;
     param.ApproxMax=eps_hi;
     param.approximation_type=COEFF_TYPE_ZOLOTAREV;
     param.RatPolyDeg=Ls;
     // The following is why I think Chroma made some directional errors:
     param.AuxFermAct= std::string(
"<AuxFermAct>\n"
"  <FermAct>UNPRECONDITIONED_WILSON</FermAct>\n"
"  <Mass>-1.8</Mass>\n"
"  <b5>1</b5>\n"
"  <c5>0</c5>\n"
"  <MaxCG>1000</MaxCG>\n"
"  <RsdCG>1.0e-9</RsdCG>\n"
"  <FermionBC>\n"
"      <FermBC>SIMPLE_FERMBC</FermBC>\n"
"      <boundary>1 1 1 1</boundary>\n"
"   </FermionBC> \n"
"</AuxFermAct>"
);
     param.AuxFermActGrp= std::string("");
     UnprecOvlapContFrac5DFermActArray S_f(fbc,param);
     Handle< FermState<T4,U,U> > fs( S_f.createState(u) );
     Handle< LinearOperatorArray<T4> > M(S_f.linOp(fs));
     return  M;
   }
   assert(0);
  }

  static Chroma::Handle< Chroma::SystemSolver<QDP::LatticeFermion> > GetSolver(QDP::multi1d<QDP::LatticeColorMatrix> &u, ChromaAction parms, DWF_parms &dwf_par)
  {
    QDP::multi1d<int> bcs(Nd);
    bcs[0] = bcs[1] = bcs[2] = bcs[3] = 1;
    Chroma::Handle<Chroma::FermBC<T4,U,U> > fbc(new Chroma::SimpleFermBC< T4, U, U >(bcs));
    Chroma::Handle<Chroma::CreateFermState<T4,U,U> > cfs( new Chroma::CreateSimpleFermState<T4,U,U>(fbc));
    int Ls = dwf_par.Ls;
    QDP::Real M5(dwf_par.M5);
    QDP::Real mq(dwf_par.mq);
    
    Chroma::GroupXML_t invparm;
    invparm.xml=std::string(
			    "   <InvertParam>\n"
			    "   <invType>CG_INVERTER</invType>\n"
			    "   <RsdCG>1.0e-10</RsdCG>\n"
			    "   <MaxCG>3000</MaxCG>\n"
			    "   </InvertParam>"
			    );
    invparm.id=std::string("CG_INVERTER");
    invparm.path=std::string("/InvertParam");
    
    Chroma::UnprecDWFermActArray  S_f(cfs, M5, mq, Ls);
    std::cout << "GetSolver: DWF 4d prec "<<std::endl;
    std::cout << "GetSolver: M5 "<<M5<<std::endl;
    std::cout << "GetSolver: mq "<<mq<<std::endl;
    std::cout << "GetSolver: Ls "<<Ls<<std::endl;
    Chroma::Handle< Chroma::FermState<T4,U,U> > fs( S_f.createState(u) );
    Chroma::Handle< LinearOperatorArray<T4> > M(S_f.unprecLinOp(fs,mq));
    return  S_f.qprop(fs,invparm);
  }
};
}

void init_chroma(const std::vector<int> &dimensions, int *argc, char ***argv)
{
  Chroma::initialize(argc, argv);
  Chroma::WilsonTypeFermActs4DEnv::registerAll();

  QDP::multi1d<int> nrow(QDP::Nd);
  for(int mu=0;mu<QDP::Nd;mu++) nrow[mu] = dimensions[mu];

  QDP::Layout::setLattSize(nrow);
  QDP::Layout::create();
}

void calc_chroma(ChromaAction action,Grid::QCD::LatticeGaugeField & lat, Grid::QCD::LatticeFermion &src, Grid::QCD::LatticeFermion &res, Grid::QCD::LatticeFermion &chroma_res, int dag, DWF_parms &dwf_par)
{
  QDP::multi1d<QDP::LatticeColorMatrix> u(4);
  //  Chroma::HotSt(u);
  Chroma::ChromaWrapper::ImportGauge(lat,u) ;
  QDP::multi1d<QDP::LatticeFermion>  check(dwf_par.Ls);
  QDP::multi1d<QDP::LatticeFermion> result(dwf_par.Ls);
  QDP::multi1d<QDP::LatticeFermion>  psi(dwf_par.Ls);
  Chroma::ChromaWrapper::ImportFermion(res,psi);
  for(int mu=0;mu<4;mu++){
    std::cout <<"Imported Gauge norm ["<<mu<<"] "<< QDP::norm2(u[mu])<<std::endl;
  }
  std::cout <<"Imported Fermion norm "<< QDP::norm2(psi)<<std::endl;

  typedef QDP::LatticeFermion T;
  typedef QDP::multi1d<QDP::LatticeColorMatrix> U;
  
  auto linop =Chroma::ChromaWrapper::GetLinOp(u, action, dwf_par);

  printf("Calling Chroma Linop\n"); fflush(stdout);

  if ( dag == 1 ) 
    (*linop)(check,psi,Chroma::MINUS);
  else
    (*linop)(check,psi,Chroma::PLUS);

  std::cout << "Called Chroma Linop" << std::endl;

  Chroma::ChromaWrapper::ExportFermion(chroma_res, check);
  
  chroma_res -= src;
  std::cout << "Norm2 of difference = " << (Grid::norm2(chroma_res) / Grid::norm2(src)) << std::endl;
}
