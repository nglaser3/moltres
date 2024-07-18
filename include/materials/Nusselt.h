#pragma once

#include "Material.h"
#include "SplineInterpolation.h"
#include "BicubicSplineInterpolation.h"
#include <cmath>
#include "INSADTauMaterial.h"

class Nusselt : public Material
{
public:
  Nusselt(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual void computeQpProperties() override;
  Real computeSpeed();

  MaterialProperty<Real> _k;
  MaterialProperty<Real> _pr;
  MaterialProperty<Real> _kin_visc;

  Real _l_value;
  MooseEnum _correlation; 

  const VariableValue * const _u_vel;
  const VariableValue * const _v_vel;
  const VariableValue * const _w_vel;
  
  MaterialProperty<Real> & _re;
  MaterialProperty<Real> & _nu;
  MaterialProperty<Real> & _h;

};
