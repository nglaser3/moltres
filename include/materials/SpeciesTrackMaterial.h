#pragma once

#include "GenericConstantMaterial.h"
#include "nlohmann/json.h"
#include <fstream>
#include <vector>
#include <map>

class SpeciesTrackMaterial : public GenericConstantMaterial
{
public:
    SpeciesTrackMaterial(const InputParameters & parameters);
    static InputParameters validParams();

protected:
    virtual void ComputeQpProperties();

    std::string _filename;
    std::string _nuclide;
    MaterialProperty<Real> & _lambda;
    MaterialProperty<Real> & _fis_yield;
    MaterialProperty<std::vector<Real>> & _abs_xs;
    MaterialProperty<std::map<std::string,std::vector<Real>>> & _dp_props;
    MaterialProperty<std::map<std::string,std::vector<Real>>> & _tp_props;
};
