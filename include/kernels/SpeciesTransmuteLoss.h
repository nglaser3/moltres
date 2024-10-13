#pragma once

#include "Kernel.h"
#include "ScalarTransportBase.h"

class SpeciesTransmuteLoss : public Kernel, public ScalarTransportBase
{
public:
    SpeciesTransmuteLoss(const InputParameters & parameters);

    static InputParameters validParams();

protected:
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
    
    unsigned int _num_groups;
    const MaterialProperty<std::vector<Real>> & _abs_xs;

    std::vector<const VariableValue*> _group_fluxes;
    std::vector<unsigned int> _flux_ids;
};