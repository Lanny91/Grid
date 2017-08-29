/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/fermion/MobiusFermion.cc

Copyright (C) 2017

Author: Andrew Lawson <andrew.lawson1991@gmail.com>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

#include <Grid/qcd/action/fermion/MobiusFermion.h>

namespace Grid {
namespace QCD {

// Constructor
template <class Impl>
MobiusFermion<Impl>::MobiusFermion(GaugeField &_Umu,
        GridCartesian         &FiveDimGrid,
        GridRedBlackCartesian &FiveDimRedBlackGrid,
        GridCartesian         &FourDimGrid,
        GridRedBlackCartesian &FourDimRedBlackGrid,
        RealD _mass,RealD _M5,
        RealD b, RealD c,const ImplParams &p) : 
    CayleyFermion5D<Impl>(_Umu,
                          FiveDimGrid,
                          FiveDimRedBlackGrid,
                          FourDimGrid,
                          FourDimRedBlackGrid,_mass,_M5,p)
{
    RealD eps = 1.0;

    std::cout<<GridLogMessage << "MobiusFermion (b="<<b<<",c="<<c<<") with Ls= "<<this->Ls<<" Tanh approx"<<std::endl;
    Approx::zolotarev_data *zdata = Approx::higham(eps,this->Ls);// eps is ignored for higham
    assert(zdata->n==this->Ls);

    // Call base setter
    this->SetCoefficientsTanh(zdata,b,c);

    Approx::zolotarev_free(zdata);
}

/*******************************************************************************
 * Conserved current implementation for Mobius domain wall fermions, for
 * contracting propagators to make a conserved current sink or inserting the
 * conserved current sequentially.
 ******************************************************************************/
// Must pass in 4D source (src) used to generate q_in_1.
template <class Impl>
void MobiusFermion<Impl>::ContractConservedCurrent(PropagatorField &q_in_1,
                                                   PropagatorField &q_in_2,
                                                   PropagatorField &q_out,
                                                   Current curr_type,
                                                   unsigned int mu,
                                                   PropagatorField *src)
{
    conformable(q_in_1._grid, this->FermionGrid());
    conformable(q_in_1._grid, q_in_2._grid);
    conformable(this->_FourDimGrid, q_out._grid);
    conformable(this->_FourDimGrid, src->_grid);
    PropagatorField qMob1(this->FermionGrid()), qMob2(this->FermionGrid());

    ConservedCurrentSetup(q_in_1, qMob1, *src);
    ConservedCurrentSetup(q_in_2, qMob2, *src);

    this->ContractConservedCurrentHt(qMob1, qMob2, q_out, curr_type, mu);
}


// Must pass in 4D source (src) used to generate q_in.
template <class Impl>
void MobiusFermion<Impl>::SeqConservedCurrent(PropagatorField &q_in,
                                              PropagatorField &q_out,
                                              Current curr_type,
                                              unsigned int mu,
                                              std::vector<Real> mom,
                                              unsigned int tmin,
                                              unsigned int tmax,
                                              PropagatorField *src)
{
    conformable(q_in._grid, this->FermionGrid());
    conformable(q_in._grid, q_out._grid);
    PropagatorField qMob(this->FermionGrid());

    ConservedCurrentSetup(q_in, qMob, *src);
    this->SeqConservedCurrentHt(qMob, q_out, curr_type, mu, mom, tmin, tmax);
}

template <class Impl>
void MobiusFermion<Impl>::ConservedCurrentSetup(PropagatorField &q_in,
                                                PropagatorField &q_out,
                                                PropagatorField &src)
{
    PropagatorField tmpBwd(this->FermionGrid()), tmpFwd(this->FermionGrid());
    PropagatorField tmps1(this->_FourDimGrid), tmps2(this->_FourDimGrid);
    unsigned int LLs = q_in._grid->_rdimensions[0];
    Gamma g5(Gamma::Algebra::Gamma5), Id(Gamma::Algebra::Identity);
    Coeff_t bdiv = this->bs[0] / (this->bs[0] + this->cs[0]);
    Coeff_t cdiv = this->cs[0] / (2.*(this->bs[0] + this->cs[0])); // Factor of 1/2 from +/- projection
    unsigned int Ls   = this->Ls;
    RealD mass        = this->mass;

    // Construct object mixing adjacent sites along 5th dimension.
    //              (b * Pm(q(s)) + c*Pp(q(s-1)) + 
    //               b * Pp(q(s)) + c*Pm(q(s+1))) / (b + c)
    tmpBwd = Cshift(q_in, 0, -1);
    tmpFwd = Cshift(q_in, 0,  1);
    tmpBwd = (tmpBwd + g5*tmpBwd);
    tmpFwd = (tmpFwd - g5*tmpFwd);
    axpby(q_out, bdiv, cdiv, q_in, tmpBwd);
    axpy(q_out, cdiv, tmpFwd, q_out);

    // Deal with edge cases at s = 0 and s = Ls - 1.
    // 
    // s = 0:       (b * Pm(q(0)) + c*Pp(src - mass*q(Ls - 1))
    //               b * Pp(q(0)) + c*Pm(q(1))) / (b + c)
    // Hence make correction for 2nd term.
    ExtractSlice(tmps1, q_in,  Ls - 1, 0);
    ExtractSlice(tmps2, q_out, 0, 0);
    axpby(tmps1, cdiv, -cdiv*(1. + mass), src, tmps1);
    tmps2 += tmps1 + g5*tmps1;
    InsertSlice(tmps2, q_out, 0, 0);

    // s = Ls - 1:  (b * Pm(q(Ls - 1)) + c*Pp(q(Ls - 2)) + 
    //               b * Pp(q(Ls - 1)) + c*Pm(src - mass*q(0))) / (b + c)
    // Hence make correction for 4th term.
    ExtractSlice(tmps1, q_in,  0, 0);
    ExtractSlice(tmps2, q_out, Ls - 1, 0);
    axpby(tmps1, cdiv, -cdiv*(1. + mass), src, tmps1);
    tmps2 += tmps1 - g5*tmps1;
    InsertSlice(tmps2, q_out, Ls - 1, 0);
}

FermOpTemplateInstantiate(MobiusFermion);
GparityFermOpTemplateInstantiate(MobiusFermion);

}}