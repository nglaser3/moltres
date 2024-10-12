#pragma once

#include "SpeciesTrackMaterial.h"

registerMooseObject("MoltresApp", SpeciesTrackMaterial);

using json = nlohmann::json;
using namespace std;

InputParameters
SpeciesTrackMaterial::validParams()
{
    InputParameters params = GenericConstantMaterial::validParams();
    params.addRequiredParam<string>("filepath",
        "The file path of the json file" 
        "constructed with moltres_xs.speciestracking.");
    params.addRequiredParam<string>("Isotope",
        "The isotope this material will store properties for.");
    return params;
}

SpeciesTrackMaterial::SpeciesTrackMaterial(const InputParameters & parameters)
    : GenericConstantMaterial(parameters), 
    _filename(getParam<string>("filepath")),
    _nuclide(getParam<string>("isotope")),
    _lambda(declareProperty<Real>(_nuclide+"_lambda")),
    _fis_yield(declareProperty<Real>(_nuclide+"_fisyield")),
    _abs_xs(declareProperty<vector<Real>>(_nuclide+"_absxs")),
    _dp_props(declareProperty<map<string, vector<Real>>>(_nuclide+"_dp_data")),
    _tp_props(declareProperty<map<string,vector<Real>>>(_nuclide+"_tp_data"))
{       
}

void 
SpeciesTrackMaterial::ComputeQpProperties()
{
    std::ifstream file(_filename);
    json data = json::parse(file);

    json decay = data[_nuclide]["decay"];
    auto dpars = decay["parents"].template get<vector<string>>();
    auto dp_lam = decay["parent_lambdas"].template get<vector<Real>>();
    auto dp_bra = decay["parent_branching"].template get<vector<Real>>();
    for (int i =0; i < dpars.size(); i++)
    {
        vector _p_data(dp_lam[i],dp_bra[i]);
        _dp_props[_qp][dpars[i]] = _p_data;
    }
    _lambda[_qp] = decay["lambda"].template get<Real>();

    json transmute = data[_nuclide]["transmute"];
    auto tpars = transmute["parents"].template get<map<string,vector<Real>>>();
    for (auto it = tpars.begin(); it != tpars.end(); it++)
    {
        _tp_props[_qp][it->first] = it->second;
    }
    _abs_xs[_qp] = transmute["abs_xs"].template get<vector<Real>>();
    _fis_yield[_qp] = transmute["fission_yield"].template get<Real>();

}