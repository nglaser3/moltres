#include "Nusselt.h"

registerMooseObject("MoltresApp", Nusselt);

InputParameters
Nusselt::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<std::string>("thermal_conductivity", 
                                      "Thermal conductivity of fluid");
  params.addRequiredParam<Real>("L", "The channel radius");
  params.addParam<std::string>("prandtl_number","Prandtl number of fluid");
  params.addParam<std::string>("kinematic_viscosity",
  "The kinematic viscosity of the fluid, dynamic viscosity divided by the density");
  MooseEnum nusselt_correlation("Krepel Dittus-Boelter", "Krepel");
  params.addParam<MooseEnum>("nusselt_correlation",nusselt_correlation,
          "Which correlation to use to approximate the Nusselt number");
  params.addCoupledVar("u_velocity","scalar variable holding u component of velocity");
  params.addCoupledVar("v_velocity","scalar variable holding v component of velocity");
  params.addCoupledVar("w_velocity","scalar variable holding w component of velocity");
  return params;
}

Nusselt::Nusselt(const InputParameters & parameters)
  : Material(parameters),
    _k(getMaterialProperty<Real>(getParam<std::string>("thermal_conductivity"))),
    _l_value(getParam<Real>("L")),
    _pr(getMaterialProperty<Real>(getParam<std::string>("prandtl_number"))),
    _kin_visc(getMaterialProperty<Real>(getParam<std::string>("kinematic_viscosity"))),
    _correlation(getParam<MooseEnum>("nusselt_correlation")),
    _re(declareProperty<Real>("ReynoldsNumber")),
    _nu(declareProperty<Real>("NusseltNumber")),
    _h(declareProperty<Real>("HeatTransferCoeff"))
{
  if (_correlation != "Krepel")
  {
    _u_vel(&coupledValue("u_velocity"));
    _v_vel(&coupledValue("v_velocity"));
    _w_vel(&coupledValue("w_velocity"));
  }
}

/* assumes velocity is an elemental variable
*See moose/framework/<include or src>/auxkernels/VectorVariableComponentAux
*for finding components of elemantal vs nodal vector variables
*/
Real Nusselt::computeSpeed()
{
    return std::pow(std::pow((*_u_vel)[_qp],2) + std::pow((*_v_vel)[_qp],2) + std::pow((*_w_vel)[_qp],2),0.5);
}

void
Nusselt::computeQpProperties()
{
  if (_correlation=="Dittus-Boelter")
  {
    _re[_qp] = computeSpeed() * _l_value / _kin_visc[_qp];
    _nu[_qp] = 0.023 * std::pow(_re[_qp],0.8) * std::pow(_pr[_qp],0.4);
  };
  if (_correlation=="Krepel")
  {
    /* Forced laminar flow; J. Křepel et al. / Annals of Nuclear Energy 34 (2007) 449–462 */
    _re[_qp] = 0.; //place holder for no reynolds
    _nu[_qp] = 4.;
  };
  _h[_qp] = _nu[_qp] * _k[_qp] / _l_value;
}
