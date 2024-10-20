#pragma once

#include "GenericMoltresMaterial.h"
#include "ScalarTransportBase.h"

class PoisonAbsorptionXSUpdate : public GenericMoltresMaterial, public ScalarTransportBase
{
public:
    static InputParameters validParams();
    PoisonAbsorptionXSUpdate(const InputParameters & parameters);

protected:
    virtual void computeQpProperties() override;
    virtual void computeSplineAbsorbingQpProperties();
    
    const VariableValue & _poison;
    const std::string _nuclide;
    const MaterialProperty<std::vector<Real>> & _abs_xs;
};