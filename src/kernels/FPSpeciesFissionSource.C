#include "FPSpeciesFissionSource.h"

 registerMooseObject("MoltresApp",FPSpeciesFissionSource);

 InputParameters
 FPSpeciesFissionSource::validParams()
 {
    InputParameters params = Kernel::validParams();
    params += ScalarTransportBase::validParams();
    params.addRequiredParam<unsigned>("num_groups",
    "Number of group fluxes in simulation");
    params.addRequiredCoupledVar("group_fluxes", 
    "All the variables that hold the group fluxes."
    "These MUST be listed by decreasing" 
    "energy/increasing group number.");
    params.addRequiredParam<std::vector<Real>>("fission_yield",
    "Fission yield, gamma, of species for each energy group.");
    params.addCoupledVar("temperature", 800, 
    "The temperature used to interpolate material properties.");
    return params;
 }

 FPSpeciesFissionSource::FPSpeciesFissionSource(const InputParameters & parameters)
   : Kernel(parameters),
     ScalarTransportBase(parameters),
     _num_groups(getParam<unsigned>("num_groups")),
     _gamma(getParam<std::vector<Real>>("fission_yield")),
     _temp(coupledValue("temperature")),
     _temp_id(coupled("temperature")),
     _nsf(getMaterialProperty<std::vector<Real>>("nsf")),
     _d_nsf_d_temp(getMaterialProperty<std::vector<Real>>("d_nsf_d_temp"))
     

{
    
    _group_fluxs.resize(_num_groups);
    _flux_ids.resize(_num_groups);
    for (unsigned int i = 0; i < _num_groups; ++i)
    {
        _group_fluxs[i] = &coupledValue("group_fluxes", i);
        _flux_ids[i] = coupled("group_fluxes", i);
    }
}

Real
FPSpeciesFissionSource::computeQpResidual()
{
  Real r = 0;
  for (unsigned int i = 0; i < _num_groups; ++i)
  {
    r += -_test[_i][_qp] *_gamma[i] * _nsf[_qp][i] * computeConcentration((*_group_fluxs[i]),_qp);
  }
    return r;
}

Real
FPSpeciesFissionSource::computeQpJacobian()
{
  return 0.0; //decay has no temperature dependence, neither does fission yield
}


Real
FPSpeciesSource::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real jac = 0.0;
  for (unsigned int i = 0; i < _num_groups; ++i)
    if (jvar == _flux_ids[i])
    {
      jac = -_test[_i][_qp] * _gamma[i] * _nsf[_qp][i] *
            computeConcentrationDerivative((*_group_fluxs[i]), _phi, _j, _qp);
      break;
    }

  if (jvar == _temp_id)
  {
    for (unsigned int i = 0; i < _num_groups; ++i)
      jac += -_test[_i][_qp] *_gamma[i] * _d_nsf_d_temp[_qp][i] * 
      _phi[_j][_qp] * computeConcentration((*_group_fluxs[i]), _qp);
  }
  return jac;
}
