#pragma once

#include "SpeciesTransmuteLoss.h"

registerMooseObject("MoltresApp",SpeciesTransmuteLoss);

InputParameters
SpeciesTransmuteLoss::validParams()
{
    InputParameters params = Kernel::validParams();
    params += ScalarTransportBase::validParams();
    params.addRequiredParam<unsigned>("num_groups",
    "Number of group fluxes in simulation");
    params.addRequiredCoupledVar("group_fluxes", 
    "All the variables that hold the group fluxes."
    "These MUST be listed by decreasing" 
    "energy/increasing group number.");
    params.addRequiredParam<MaterialPropertyName>("micro_abs_xs_material",
    "Name of material property holding microscopic" 
    "absorbtion cross section of species.");
    return params;
}

SpeciesTransmuteLoss::SpeciesTransmuteLoss(const InputParameters & parameters)
    : Kernel(parameters),
      ScalarTransportBase(parameters),
      _num_groups(getParam<unsigned>("num_groups")),
      _abs_xs(getMaterialProperty<std::vector<Real>>("micro_abs_xs_material"))
{
    _group_fluxes.resize(_num_groups);
    _flux_ids.resize(_num_groups);
    for (unsigned int i = 0; i < _num_groups; ++i)
    {
        _group_fluxes[i] = &coupledValue("group_fluxes", i);
        _flux_ids[i] = coupled("group_fluxes", i);
    }
}

Real 
SpeciesTransmuteLoss::computeQpResidual()
{
    Real res = 0.0;
    for (int i = 0; i < _num_groups; i++)
    {
        res += _test[_i][_qp] * _abs_xs[_qp][i]
        * computeConcentration(_u,_qp)
        * computeConcentration((*_group_fluxes[i]),_qp);
    }
    return res;
}

Real
SpeciesTransmuteLoss::computeQpJacobian()
{
    Real jac = 0.0;
    for (int i = 0; i < _num_groups; i++)
    {
        jac -=_test[_i][_qp] * _abs_xs[_qp][i] 
        * computeConcentration((*_group_fluxes[i]),_qp)
        * computeConcentrationDerivative(_u, _phi, _j, _qp);
    }
    return jac;
}

Real 
SpeciesTransmuteLoss::computeQpOffDiagJacobian(unsigned int jvar)
{
    Real jac = 0.0;
    for (int i = 0; i < _num_groups; i++)
    {
        if (jvar == _flux_ids[i])
        {
            jac -= _test[_i][_qp] * _abs_xs[_qp][i]
            * computeConcentrationDerivative((*_group_fluxes[i]),_phi, _j, _qp)
            * computeConcentration(_u,_qp);
            break;
        }
        
    }
    return jac;
}