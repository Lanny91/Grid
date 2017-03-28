/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/FreeQuark.hpp

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

#ifndef Hadrons_FreeQuark_hpp_
#define Hadrons_FreeQuark_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                               TFreeQuark                                   *
 ******************************************************************************/
class FreeQuarkPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FreeQuarkPar,
                                    std::string, source,
                                    std::string, action,
                                    double     , mass);
};

template <typename FImpl>
class TFreeQuark: public Module<FreeQuarkPar>
{
public:
    TYPE_ALIASES(FImpl,);
public:
    // constructor
    TFreeQuark(const std::string name);
    // destructor
    virtual ~TFreeQuark(void) = default;
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Ls_;
};

MODULE_REGISTER(FreeQuark, TFreeQuark<FIMPL>);

/******************************************************************************
 *                          TFreeQuark implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TFreeQuark<FImpl>::TFreeQuark(const std::string name)
: Module(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TFreeQuark<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().source, par().action};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TFreeQuark<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TFreeQuark<FImpl>::setup(void)
{
    env().template registerLattice<PropagatorField>(getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TFreeQuark<FImpl>::execute(void)
{
    LOG(Message) << "Computing free quark propagator '" << getName() << "'"
                 << std::endl;

    PropagatorField &prop   = *env().template createLattice<PropagatorField>(getName());
    PropagatorField &source = *env().template getObject<PropagatorField>(par().source);
    FermionField ferm(env().getGrid()), src(env().getGrid());
    auto &mat = *(env().template getObject<FMat>(par().action));

    for (unsigned int s = 0; s < Ns; ++s)
    for (unsigned int c = 0; c < Nc; ++c)
    {
        PropToFerm(src, source, s, c);
        mat.FreePropagator(src, ferm, par().mass);
        FermToProp(prop, ferm, s, c);
    }

    std::cout << "Norm of propagator  '" << getName() << "' is " 
              << sqrt(norm2(p4)) << std::endl;
}

END_HADRONS_NAMESPACE

#endif // Hadrons_FreeQuark_hpp_
