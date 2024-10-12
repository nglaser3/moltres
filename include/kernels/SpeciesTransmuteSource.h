#pragma once

#include "ScalarTransportBase.h"
#include "Kernel.h"

class SpeciesTransmuteSource : public Kernel, public ScalarTransportBase
{
public:
    SpeciesTransmuteSource(const InputParameters & parameters);
    static InputParameters validParams();

protected:
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

    unsigned int _num_groups;
    std::vector<std::string> _parents;
    unsigned int _num_parents;

    const MaterialProperty<std::unordered_map<std::string,std::vector<Real>>> & _tp_map;

    std::vector<const VariableValue*> _parent_concs;
    std::vector<unsigned int> _parent_ids;


    std::vector<const VariableValue*> _group_fluxes;
    std::vector<unsigned int> _flux_ids;

};

