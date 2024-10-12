#pragma once

#include "Kernel.h"
#include "ScalarTransportBase.h"

class SpeciesFissionSource : public Kernel, public ScalarTransportBase
{
public:
    SpeciesFissionSource(const InputParameters & parameters);
    static InputParameters validParams();

protected: 
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

    unsigned int _num_groups;

    const VariableValue & _temp;
    unsigned int _temp_id;
    const MaterialProperty<Real> & _fis_yield;
    const MaterialProperty<std::vector<Real>> & _nsf;
    const MaterialProperty<std::vector<Real>> & _d_nsf_d_temp;
    
    std::vector<const VariableValue*> _group_fluxs;
    std::vector<unsigned int> _flux_ids;
};

