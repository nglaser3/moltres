#pragma once

#include "ScalarTransportBase.h"

class FPSpeciesTransmutationSource : public Kernel, public ScalarTransportBase
{

public:
    static InputParameters validParams();
    FPSpeciesTransmutationSource(const InputParameters & parameters);

protected:
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

    std::vector<std::string> _parent_zaids;
    int _num_parents;
    int _num_groups;

    bool _transmute;
    bool _decay;

    std::vector<const VariableValue*> _parent_concs;
    std::vector<unsigned int> _parent_ids;
    std::vector<const VariableValue*> _group_fluxs;
    std::vector<unsigned int> _flux_ids;
};
