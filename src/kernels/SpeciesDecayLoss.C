#pragma once

#include "SpeciesDecayLoss.h"

registerMooseObject("MoltresApp",SpeciesDecayLoss);

InputParameters
SpeciesDecayLoss::validParams()
{
    InputParameters params = Kernel::validParams();
    params += ScalarTransportBase::validParams();
    params.addRequiredParam<MaterialPropertyName>("decay_constant",
    "Name of material property holding decay constant on species.");
    return params;
}

SpeciesDecayLoss::SpeciesDecayLoss(const InputParameters & parameters)
    : Kernel(parameters),
      ScalarTransportBase(parameters),
      _lambda(getMaterialProperty<Real>("decay_constant"))
{
}

Real 
SpeciesDecayLoss::computeQpResidual()
{
    return _test[_i][_qp] * _lambda[_qp] * computeConcentration(_u,_qp);
}

Real
SpeciesDecayLoss::computeQpJacobian()
{
    return _test[_i][_qp] * _lambda[_qp] * computeConcentrationDerivative(_u, _phi, _j, _qp);
}