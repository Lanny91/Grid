/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MUtilities/OutputProp.hpp

Copyright (C) 2017

Author: Andrew Lawson    <andrew.lawson1991@gmail.com>

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

#ifndef Hadrons_OutputProp_hpp_
#define Hadrons_OutputProp_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         OutputProp                                         *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

GRID_SERIALIZABLE_ENUM(OutputOpt, undef,
                       coordinate, 1,
                       timeslice, 2,
                       full, 3);

class OutputPropPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(OutputPropPar,
                                    std::string, q,
                                    std::string, coords,
                                    OutputOpt,   opt,
                                    std::string, output);
};

template <typename FImpl>
class TOutputProp: public Module<OutputPropPar>
{
public:
    TYPE_ALIASES(FImpl,);
    typedef typename SitePropagator::scalar_object SiteProp;
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::string, name,
                                        std::vector<SiteProp>, out);
    };
public:
    // constructor
    TOutputProp(const std::string name);
    // destructor
    virtual ~TOutputProp(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(OutputProp, TOutputProp<FIMPL>, MUtilities);

/******************************************************************************
 *                 TOutputProp implementation                                 *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TOutputProp<FImpl>::TOutputProp(const std::string name)
: Module<OutputPropPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TOutputProp<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TOutputProp<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TOutputProp<FImpl>::setup(void)
{
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TOutputProp<FImpl>::execute(void)
{
    CorrWriter      writer(par().output);
    PropagatorField &q = *env().template getObject<PropagatorField>(par().q);
    Result          result;

    switch(par().opt)
    {
        case OutputOpt::coordinate:
        {
            SiteProp s;
            std::vector<int> siteCoord;
            siteCoord = strToVec<int>(par().coords);
            peekSite(s, q, siteCoord);
            result.name = par().coords;
            result.out.push_back(s);
            break;
        }
        case OutputOpt::timeslice:
            HADRON_ERROR("Timeslice output not yet implemented");
        case OutputOpt::full:
            HADRON_ERROR("Full propagator output not yet implemented");
        default:
            HADRON_ERROR("No output option specified.");
    }

    write(writer, "propagator", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_OutputProp_hpp_
