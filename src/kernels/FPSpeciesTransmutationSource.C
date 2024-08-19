#pragma once

#include "FPSpeciesTransmutationSource.h"

registerMooseObject("MoltresApp", FPSpeciesTransmutationSource);

InputParameters 
FPSpeciesTransmutationSource::validParams()
{
    params = Kernel::validParams();
    params += ScalarTransportBase::validParams();
    params.addRequiredParam<std::vector<std::string>>(
    "parents_identifiers", "The parent identifier, in ZAID format."
    "Multiply the atomic number, Z, by 1000; and then add the atomic mass, A."
    "For example, Xe-135 : 54135; C-12 : 06012.");
    params.addParam<std::vector<Real>>("decay_constant", 
    "The decay constant of the parents that yields the species of interest."
    "This is the decay constant of the parents multiplied by the branching"
    "fraction and factor that yield the species of interest.")
    params.addCoupledVar("group_fluxes", 
    "All the variables that hold the group fluxes."
    "These MUST be listed by decreasing" 
    "energy/increasing group number.");
    params.addParam<unsigned>("num_groups",
    "Number of group fluxes in simulation");
    params.addParam<int>("num_groups",
    "The number of energy groups in the simulation")
    params.addCoupledVar("parent_concentrations",
    "The variable(s) holding the concentrations of the"
    "DIRECT parents of the species of interest.");
    params.addCoupledVar("temperature", 800, 
    "The temperature used to interpolate material properties.");
}

FPSpeciesTransmutationSource::FPSpeciesTransmutationSource(const InputParameters & parameters)
  : Kernel(parameters),
    ScalarTransportBase(parameters),
    _parent_zaids(getParam<std::vector<std::string>>("parent_identifiers")),
    _num_parents(_parent_zaids.size()),
    _num_groups(getParam<int>("num_groups"))
    _transmute(isValidParam("group_fluxes")),
    _decay(isValidParam("decay_constant"))
{
    _parent_concs.resize(_num_parents);
    _parent_ids.resize(_num_parents);
    for (unsigned int i = 0; i < _num_parents; ++i)
    {
        _parent_concs[i] = &coupledValue("parent_concentrations", i);
        _parent_ids[i] = coupled("parent_concentrations", i);
    }

    if (_transmute)
    {
        _group_fluxs.resize(_num_groups);
        _flux_ids.resize(_num_groups);
        for (unsigned int i = 0; i < _num_groups; ++i)
        {
            _group_fluxs[i] = &coupledValue("group_fluxes", i);
            _flux_ids[i] = coupled("group_fluxes", i);
        }
    }
    if (!_decay)
    {
        _lambda = getParam<std::vector<Real>>("decay_constant");
    }
    if (!_decay && !_transmute)
    {
        mooseError("Parent must either decay or transmute into the species of interest. 
        If decay is desired, please provide the decay_constant parameter. If transmutation
        is desired, please provide the group_fluxes and num_groups parameter, ASWELL as
        insure the material block 'FileSpeciesTransmute' is in the simulation!");
    }
}

Real
FPSpeciesTransmutationSource::computeQpResidual()
{
    Real res = 0.0;
    for (int i = 0; i < num_parents; i++)
    {
        if (_decay)
        {
            res += _test[_i][_qp]*_lambda * computeConcentration((*_parent_concs[i]),_qp);
        }
        if (_transmute)
        {
            for (int j = 0; j < num_groups; j++)
            {
                res +=  _test[_i][_qp] * _microxs[i] * computeConcentration((*parent_concs[i]),_qp) 
            * computeConcentration((*_group_fluxes[j]),_qp);
            }
        }
    }
}

Real 
FPSpeciesTransmutationSource::computeQpJacobian()
{
    return 0.0;
}


