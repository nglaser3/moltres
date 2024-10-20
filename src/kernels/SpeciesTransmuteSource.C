#pragma once

#include "SpeciesTransmuteSource.h"

registerMooseObject("MoltresApp",SpeciesTransmuteSource);

InputParameters
SpeciesTransmuteSource::validParams()
{
    InputParameters params = Kernel::validParams();
    params += ScalarTransportBase::validParams();
    params.addRequiredParam<unsigned>("num_groups",
    "Number of group fluxes in simulation");
    params.addRequiredCoupledVar("group_fluxes", 
    "All the variables that hold the group fluxes."
    "These MUST be listed by decreasing" 
    "energy/increasing group number.");
    params.addRequiredCoupledVar("parent_concs",
    "All the variables that hold the flux-induced" 
    "transmutation parent concentrations. These MUST"
    "be listed in the same order as in parent_names.");
    params.addRequiredParam<std::vector<std::string>>("parent_names",
    "All flux-induced transmutation parent names, MUST" 
    "be listed in the same order as in parent_concs.");
    params.addRequiredParam<MaterialPropertyName>("transmute_material",
    "Material holding all parent transmutation properties.");
    return params;
}

SpeciesTransmuteSource::SpeciesTransmuteSource(const InputParameters & parameters)
    : Kernel(parameters),
      ScalarTransportBase(parameters),
      _num_groups(getParam<int>("num_groups")),
      _parents(getParam<std::vector<std::string>>("parent_names")),
      _num_parents(_parents.size()),
      _tp_map(getMaterialProperty<std::unordered_map<std::string,std::vector<Real>>>("transmute_material"))
{

    _parent_concs.resize(_num_parents);
    _parent_ids.resize(_num_parents);
    for (unsigned i = 0; i < _num_parents; ++i)
    {
        _parent_concs[i] = &coupledValue("parent_concs",i);
        _parent_ids[i] = coupled("parent_concs",i);
    }
    

    _group_fluxes.resize(_num_groups);
    _flux_ids.resize(_num_groups);
    for (unsigned int i = 0; i < _num_groups; ++i)
    {
        _group_fluxes[i] = &coupledValue("group_fluxes", i);
        _flux_ids[i] = coupled("group_fluxes", i);
    }

}

Real
SpeciesTransmuteSource::computeQpResidual()
{
    Real res = 0.0;
    for (int i = 0; i < _num_groups; i++)
    {
        for (int j = 0; j < _num_parents; j++)
        {
            res += _test[_i][_qp] * _tp_map[_qp].at(_parents[j])[i]
            * computeConcentration((*_parent_concs[j]),_qp)
            * computeConcentration((*_group_fluxes[i]),_qp);
        }
    }
    return res;
}

Real 
SpeciesTransmuteSource::computeQpJacobian()
{
    return 0.0;
}

Real
SpeciesTransmuteSource::computeQpOffDiagJacobian(unsigned int jvar)
{
    Real jac = 0.0;
    for (int i = 0; i < _num_groups; ++i)
    {
        if (jvar == _flux_ids[i])
        {
            for (int j = 0; j < _num_parents; ++j)
            {
                jac -= _test[_i][_qp] * _tp_map[_qp].at(_parents[j])[i]
                * computeConcentrationDerivative((*_group_fluxes[i]),_phi, _j, _qp) 
                * computeConcentration((*_parent_concs[j]),_qp);
            }
            break;
        }
        
    }
    for (int i = 0; i < _num_parents; ++i)
    {
        if (jvar == _parent_ids[i])
        {
            for (int j = 0; j < _num_groups; ++j)
            {
                jac -= _test[_i][_qp] * _tp_map[_qp].at(_parents[i])[j]
                * computeConcentration((*_group_fluxes[j]),_qp)
                * computeConcentrationDerivative((*_parent_concs[i]),_phi, _j, _qp);
            }
            break;
        }
    }
    
    return jac;
}
