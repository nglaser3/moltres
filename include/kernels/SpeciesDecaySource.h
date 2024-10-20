#pragma once

#include "ScalarTransportBase.h"
#include "Kernel.h"

class SpeciesDecaySource : public Kernel, public ScalarTransportBase
{
public:
    SpeciesDecaySource(const InputParameters & parameters);
    static InputParameters validParams();

protected:
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

    std::vector<std::string> _parents;
    unsigned int _num_parents;

    const MaterialProperty<std::unordered_map<std::string,Real>> & _dp_map;

    std::vector<const VariableValue*> _parent_concs;
    std::vector<unsigned int> _parent_ids;
};