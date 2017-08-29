    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/fermion/MobiusFermion.h

    Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
#ifndef  GRID_QCD_MOBIUS_FERMION_H
#define  GRID_QCD_MOBIUS_FERMION_H

#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/CayleyFermion5D.h>

namespace Grid {

  namespace QCD {

    template<class Impl>
    class MobiusFermion : public CayleyFermion5D<Impl>
    {
    public:
     INHERIT_IMPL_TYPES(Impl);
    public:
      virtual void   Instantiatable(void) {};
      // Constructors
      MobiusFermion(GaugeField &_Umu,
		    GridCartesian         &FiveDimGrid,
		    GridRedBlackCartesian &FiveDimRedBlackGrid,
		    GridCartesian         &FourDimGrid,
		    GridRedBlackCartesian &FourDimRedBlackGrid,
		    RealD _mass,RealD _M5,
		    RealD b, RealD c,const ImplParams &p= ImplParams());
      
      // Conserved Current utilities
      void ContractConservedCurrent(PropagatorField &q_in_1,
                                    PropagatorField &q_in_2,
                                    PropagatorField &q_out,
                                    Current curr_type,
                                    unsigned int mu,
                                    PropagatorField *src);
      void SeqConservedCurrent(PropagatorField &q_in,
                               PropagatorField &q_out,
                               Current curr_type,
                               unsigned int mu, 
                               std::vector<Real> mom,
                               unsigned int tmin,
                               unsigned int tmax,
                               PropagatorField *src);
    private:
      void ConservedCurrentSetup(PropagatorField &q_in,
                                 PropagatorField &q_out,
                                 PropagatorField &src);
    };

  }
}

#endif
