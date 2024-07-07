#include "Nusselt.h"

registerMooseObject("MoltresApp", Nusselt);

InputParameters
Nusselt::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<Real>("thermal_conductivity", 
                        "Thermal conductivity of fluid");
  params.addRequiredParam<Real>("L", "The channel radius");
  MooseEnum nusselt_correlation("Krepel Dittus-Boelter", "Krepel");
  params.addParam<MooseEnum>("nusselt_correlation",nusselt_correlation,
                             "Which correlation to use to approximate the Nusselt number");
  params.addParam<Real>("prandtl_number",10.0,"Prandtl number of fluid");
  params.addParam<Real>("kinematic_viscosity",1.0e-6,
                        "The kinematic viscosity of the fluid, dynamic viscosity divided by the density");
  params.addParam<bool>("vector_velocity",false,
                        "Whether the velocity variable is a vector or in components, true for vector");
  params.addCoupledVar("velocity","vector variable holding velocity");
  params.addCoupledVar("u_velocity","scalar variable holding u component of velocity");
  params.addCoupledVar("v_velocity","scalar variable holding v component of velocity");
  params.addCoupledVar("w_velocity","scalar variable holding w component of velocity");
  return params;
}

Nusselt::Nusselt(const InputParameters & parameters)
  : Material(parameters),
    _k(getParam<Real>("thermal_conductivity")),
    _l_value(getParam<Real>("L")),
    _correlation(getParam<MooseEnum>("nusselt_correlation")),
    _pr(getParam<Real>("prandtl_number")),
    _kin_visc(getParam<Real>("kinematic_viscosity")),
    _vec_velocity(getParam<bool>("vector_velocity")),
    _velocity(((getParam<MooseEnum>("nusselt_correlation") != "Krepel") && (getParam<bool>("vector_velocity"))) ? &coupledVectorValue("velocity") : nullptr),
    _u_vel(((getParam<MooseEnum>("nusselt_correlation") != "Krepel") && !(getParam<bool>("vector_velocity"))) ? &coupledValue("u_velocity") : nullptr),
    _v_vel(((getParam<MooseEnum>("nusselt_correlation") != "Krepel") && !(getParam<bool>("vector_velocity"))) ? &coupledValue("v_velocity") : nullptr),
    _w_vel(((getParam<MooseEnum>("nusselt_correlation") != "Krepel") && !(getParam<bool>("vector_velocity"))) ? &coupledValue("w_velocity") : nullptr),
    _re(declareProperty<Real>("ReynoldsNumber")),
    _nu(declareProperty<Real>("NusseltNumber")),
    _h(declareProperty<Real>("HeatTransferCoeff"))
{
}

/* assumes velocity is an elemental variable
*See moose/framework/<include or src>/auxkernels/VectorVariableComponentAux
*for finding components of elemantal vs nodal vector variables
*/
Real Nusselt::computeSpeed()
{
  if (_vec_velocity)
  {
    return std::pow(std::pow((*_velocity)[_qp](0),2)+std::pow((*_velocity)[_qp](1),2)+std::pow((*_velocity)[_qp](2),2),.5);
  }
  else
  {
    return std::pow(std::pow((*_u_vel)[_qp],2) + std::pow((*_v_vel)[_qp],2) + std::pow((*_w_vel)[_qp],2),0.5);
  }
}

void
Nusselt::computeQpProperties()
{
  if (_correlation=="Dittus-Boelter")
  {
    _re[_qp] = Nusselt::computeSpeed() * _l_value / _kin_visc;
    _nu[_qp] = 0.023 * std::pow(_re[_qp],0.8) * std::pow(_pr,0.4);
  };
  if (_correlation=="Krepel")
  {
    /* Forced laminar flow; J. Křepel et al. / Annals of Nuclear Energy 34 (2007) 449–462 */
    _re[_qp] = 0.; //place holder for no reynolds
    _nu[_qp] = 4.;
  };
  _h[_qp] = _nu[_qp] * _k / _l_value;
}
