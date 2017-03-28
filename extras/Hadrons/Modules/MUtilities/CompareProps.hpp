#ifndef Hadrons_CompareProps_hpp_
#define Hadrons_CompareProps_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                             CompareProps                                   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class ComparePropsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ComparePropsPar,
                                    std::string, q1,
                                    std::string, q2,
                                    double, precision);
};

template <typename FImpl>
class TCompareProps: public Module<ComparePropsPar>
{
public:
    TYPE_ALIASES(FImpl,);
public:
    // constructor
    TCompareProps(const std::string name);
    // destructor
    virtual ~TCompareProps(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(CompareProps, TCompareProps<FIMPL>, MUtilities);

/******************************************************************************
 *                   TCompareProps implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TCompareProps<FImpl>::TCompareProps(const std::string name)
: Module<ComparePropsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TCompareProps<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q1, par().q2};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TCompareProps<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TCompareProps<FImpl>::setup(void)
{
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TCompareProps<FImpl>::execute(void)
{
    PropagatorField &q1  = *env().template getObject<PropagatorField>(par().q1);
    PropagatorField &q2  = *env().template getObject<PropagatorField>(par().q2);
    PropagatorField diff(env().getGrid());
    diff = q1 - q2;
    if (norm2(diff) < par().precision)
    {
        LOG(Message) << "Propagators '" << par().q1 << "' and '" << par().q2 
                     << "' are equivalent." << std::endl;
    }
    else
    {
        LOG(Message) << "Propagators '" << par().q1 << "' and '" << par().q2 
                     << "' are not equivalent." << std::endl;
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_CompareProps_hpp_
