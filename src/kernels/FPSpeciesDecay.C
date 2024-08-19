#include "FPSpeciesDecay.h"

registerMooseObject("MoltresApp", FPSpeciesDecay);

InputParameters 
FPSpeciesDecay::validParams()
{
    InputParameters params = Kernel::validParams();
    params += ScalarTransportBase::validParams();
    params.addRequiredParam<Real>("decay_constant",
    "Decay constant of fission product species");
    return params;
}

FPSpeciesDecay::FPSpeciesDecay(const InputParameters & parameters)
  : Kernel(parameters),
    ScalarTransportBase(parameters),
    _lambda(getParam<Real>("decay_constant"))
{ 
}

Real
FPSpeciesDecay::computeQpResidual()
{
  return _test[_i][_qp] * _lambda * computeConcentration(_u,_qp);
}

Real
FPSpeciesDecay::computeQpJacobian()
{
  return _test[_i][_qp] * _lambda * computeConcentrationDerivative(_u, _phi, _j, _qp);
}

