#include "Nusselt.h"

registerMooseObject("MoltresApp", Nusselt);

InputParameters
Nusselt::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<MaterialPropertyName>("thermal_conductivity", 
        "Name of material property holding thermal conductivity of fluid");
  params.addRequiredParam<Real>("L", "The channel radius");
  params.addParam<MaterialPropertyName>("prandtl_number",
  "Name of material property holding prandtl number of fluid");
  params.addParam<MaterialPropertyName>("kinematic_viscosity",
  "Name of material property holding the kinematic viscosity of the fluid, dynamic viscosity divided by the density");
  MooseEnum nusselt_correlation("Krepel Dittus-Boelter", "Krepel");
  params.addParam<MooseEnum>("nusselt_correlation",nusselt_correlation,
          "Which correlation to use to approximate the Nusselt number");
  params.addCoupledVar("u_velocity","scalar variable holding u component of velocity");
  params.addCoupledVar("v_velocity","scalar variable holding v component of velocity");
  params.addCoupledVar("w_velocity","scalar variable holding w component of velocity");
  params.addParam<std::string>("re_material_name","ReynoldsNumber",
  "Name of generated reynolds number material property");
  params.addParam<std::string>("nu_material_name","NusseltNumber",
  "Name of generated nusselt number material property");
  params.addParam<std::string>("h_material_name","HeatTransferCoeff",
  "Name of generated heat transfer coefficient material property");
  return params;
}

Nusselt::Nusselt(const InputParameters & parameters)
  : Material(parameters),
    _k(getMaterialProperty<Real>("thermal_conductivity")),
    _l_value(getParam<Real>("L")),
    _pr(getMaterialProperty<Real>("prandtl_number")),
    _kin_visc(getMaterialProperty<Real>("kinematic_viscosity")),
    _correlation(getParam<MooseEnum>("nusselt_correlation")),
    _u_vel(_correlation!="Krepel" && isParamValid("u_velocity") ? 
          &coupledValue("u_velocity") : nullptr),
    _v_vel(_correlation!="Krepel" && isParamValid("v_velocity") ? 
          &coupledValue("v_velocity") : nullptr),
    _w_vel(_correlation!="Krepel" && isParamValid("w_velocity") ?  
          &coupledValue("w_velocity") : nullptr),
    _re(declareProperty<Real>(getParam<std::string>("re_material_name"))),
    _nu(declareProperty<Real>(getParam<std::string>("nu_material_name"))),
    _h(declareProperty<Real>(getParam<std::string>("h_material_name")))
{
    _velocity = {_u_vel, _v_vel, _w_vel};
    for (int i = 2; i >= 0; i--)
    {
      if(_velocity[i] == nullptr)
        {
          _velocity.erase(_velocity.begin() + i);
        }
    }
    
}
 
/* assumes velocity is an elemental variable
*See moose/framework/<include or src>/auxkernels/VectorVariableComponentAux
*for finding components of elemantal vs nodal vector variables
*/
Real Nusselt::computeVelocityMagnitude(unsigned int qp)
{
    Real speed = 0.0;
    for (int i = 0; i < _velocity.size(); i++)
    {
      speed += std::pow((*_velocity[i])[qp],2);
    }
    return std::pow(speed,0.5);
    //return std::pow(std::pow((*_u_vel)[qp],2) + std::pow((*_v_vel)[qp],2) + std::pow((*_w_vel)[qp],2),0.5);
}

void
Nusselt::computeQpProperties()
{
  if (_correlation=="Dittus-Boelter")
  {
    _re[_qp] = computeVelocityMagnitude(_qp) * _l_value / _kin_visc[_qp];
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
