#include "ConvectiveHTBC.h"

registerMooseObject("MoltresApp",ConvectiveHTBC);

InputParameters
ConvectiveHTBC::validParams()
{
    InputParameters params = IntegratedBC::validParams();
    params.addRequiredParam<MaterialPropertyName>("h_coeff_material",
    "Name of material holding the heat transfer coefficient");
    params.addRequiredCoupledVar("bulk_temperature",
    "Variable holding the bulk temperature of the fluid");
    return params;
}

ConvectiveHTBC::ConvectiveHTBC(const InputParameters & parameters)
    : IntegratedBC(parameters),
      _h(getMaterialProperty<Real>("h_coeff_material")),
      _bulk_temp(coupledValue("bulk_temperature"))
{
}

Real 
ConvectiveHTBC::computeQpResidual()
{
    return (_test[_i][_qp] * (_h[_qp]) * (_u[_qp] - _bulk_temp[_qp]));
}

Real 
ConvectiveHTBC::computeQpJacobian()
{
    return (_test[_i][_qp] * (_h[_qp])*_phi[_j][_qp]);
}
