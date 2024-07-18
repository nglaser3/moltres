#include "ConvectiveHTBC.h"

registerMooseObject("MoltresApp",ConvectiveHTBC);

InputParameters
ConvectiveHTBC::validParams()
{
    InputParameters params = IntegratedBC::validParams();
    params.addRequiredParam<std::string>("h_coeff_material",
    "Name of material holding the heat transfer coefficient");
    params.addRequiredParam<Real>("bulk_temperature",
    "Constant temperature of the bulk fluid");
    return params;
}

ConvectiveHTBC::ConvectiveHTBC(const InputParameters & parameters)
    : IntegratedBC(parameters),
      _h(getMaterialProperty<Real>("h_coeff_material")),
      _bulk_temp(getParam<Real>("bulk_temperature"))
{
}

Real 
ConvectiveHTBC::computeQpResidual()
{
    return (_test[_i][_qp] * (_h[_qp]) * (_u[_qp] - _bulk_temp));
}

Real 
ConvectiveHTBC::computeQpJacobian()
{
    return (_test[_i][_qp] * (_h[_qp])*_phi[_j][_qp]);
}
