#pragma once

#include "SpeciesTransmuteLoss.h"

registerMooseObject("MoltresApp",SpeciesTransmuteLoss);

InputParameters
SpeciesTransmuteLoss::validParams()
{
    InputParameters params = Kernel::validParams();
    params += ScalarTransportBase::validParams();
    params.addRequiredParam<MaterialPropertyName>("micro_abs_xs_material",
    "Name of material property holding microscopic" 
    "absorbtion cross section of species.");
    return params;
}

SpeciesTransmuteLoss::SpeciesTransmuteLoss(const InputParameters & parameters)
    : Kernel(parameters),
      ScalarTransportBase(parameters),
      _abs_xs(getMaterialProperty<Real>("micro_abs_xs_material"))
{
}

Real 
SpeciesTransmuteLoss::computeQpResidual()
{
    return _test[_i][_qp] * _abs_xs[_qp] * computeConcentration(_u,_qp);
}

Real
SpeciesTransmuteLoss::computeQpJacobian()
{
    return _test[_i][_qp] * _abs_xs[_qp] * computeConcentrationDerivative(_u, _phi, _j, _qp);
}