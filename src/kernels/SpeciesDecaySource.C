#pragma once

#include "SpeciesDecaySource.h"

registerMooseObject("MoltresApp",SpeciesDecaySource);

InputParameters
SpeciesDecaySource::validParams()
{
    InputParameters params = Kernel::validParams();
    params += ScalarTransportBase::validParams();
    params.addRequiredCoupledVar("parent_concs",
    "All the variables that hold the" 
    "decay parent concentrations. These MUST"
    "be listed in the same order as in parent_names.");
    params.addRequiredParam<std::vector<std::string>>("parent_names",
    "All decay parent names, MUST be listed in"
    "the same order as parent_concs");
    params.addRequiredParam<MaterialPropertyName>("decay_parent_material",
    "Material holding all parent decay properties");
    return params;
}

SpeciesDecaySource::SpeciesDecaySource(const InputParameters & parameters)
    : Kernel(parameters),
      ScalarTransportBase(parameters),
      _parents(getParam<std::vector<std::string>>("parent_names")),
      _num_parents(_parents.size()),
      _dp_map(getMaterialProperty<std::unordered_map<std::string,Real>>("decay_parent_material"))
{
    _parent_concs.resize(_num_parents);
    _parent_ids.resize(_num_parents);
    for (unsigned i = 0; i < _num_parents; ++i)
    {
        _parent_concs[i] = &coupledValue("parent_concs",i);
        _parent_ids[i] = coupled("parent_concs",i);
    }

}

Real
SpeciesDecaySource::computeQpResidual()
{
    Real res = 0.0;

    for (int i = 0; i < _num_parents; i++)
    {
        res += _test[_i][_qp] * _dp_map[_qp].at(_parents[i])
        * computeConcentration((*_parent_concs[i]),_qp);
    }
    return res;
}

Real 
SpeciesDecaySource::computeQpJacobian()
{
    return 0.0;
}

Real 
SpeciesDecaySource::computeQpOffDiagJacobian(unsigned int jvar)
{
    Real jac = 0.0;

    for (int i = 0; i < _num_parents; i++)
    {
        if (jvar == _parent_ids[i])
        {
            jac -= _test[_i][_qp] * _dp_map[_qp].at(_parents[i])
            * computeConcentrationDerivative((*_parent_concs[i]),_phi,_j,_qp);
            break;
        }
    }
    return jac;
}