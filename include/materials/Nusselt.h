#pragma once

#include "Material.h"
#include "SplineInterpolation.h"
#include "BicubicSplineInterpolation.h"
#include <cmath>
#include <vector>
#include "INSADTauMaterial.h"

class Nusselt : public Material
{
public:
  Nusselt(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual void computeQpProperties() override;
  Real computeVelocityMagnitude(unsigned int qp);

  const MaterialProperty<Real> & _k;
  Real _l_value;

  const MaterialProperty<Real> & _pr;
  const MaterialProperty<Real> & _kin_visc;
  
  MooseEnum _correlation; 
  
  const VariableValue * const _u_vel;
  const VariableValue * const _v_vel;
  const VariableValue * const _w_vel;

  MaterialProperty<Real> & _re;
  MaterialProperty<Real> & _nu;
  MaterialProperty<Real> & _h;

  std::vector< const VariableValue *> _velocity;
};
