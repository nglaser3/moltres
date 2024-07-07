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
  Real _k;
  Real _l_value;
  MooseEnum _correlation; 
  Real _pr;
  Real _kin_visc;
  bool _vec_velocity;
  const VectorVariableValue * const _velocity;
  const VariableValue * const _u_vel;
  const VariableValue * const _v_vel;
  const VariableValue * const _w_vel;
  MaterialProperty<Real> & _re;
  MaterialProperty<Real> & _nu;
  MaterialProperty<Real> & _h;
  //using T::_velocity;
};
