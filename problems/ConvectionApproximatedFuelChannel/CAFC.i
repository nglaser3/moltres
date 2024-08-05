[Mesh]
    [cartesianmesh]
        type = CartesianMeshGenerator
        dim = 2
        dx = 1
        dy = 1
    []
[]

[Variables]
    [temp]
    []
[]

[Kernels]
    [conduction]
        type = HeatConduction
        variable = temp
        diffusion_coefficient = graphite_therm_cond_material
    []
[]

[AuxVariables]
    [v_velocity]
    []
    [bulk_temp]
        initial_condition = 800
    []
[]

[Functions]
    [velocity_func]
        type = ParsedFunction
        expression = '5 * y'
    []
[]

[AuxKernels]
    [./velocity_aux]
        type = FunctionAux
        variable = v_velocity
        function = velocity_func
    [../]
[]

[Materials]
    [thermal_conductivity_graphite]
        type = GenericConstantMaterial
        prop_names = 'graphite_therm_cond_material'
        prop_values = '10.0'
    []

    
    [thermal_conductivity_fuel]
        type = GenericConstantMaterial
        prop_names = 'fuel_therm_cond_material'
        prop_values = '10.0'
    []
    [prandtl]
        type = GenericConstantMaterial
        prop_names = 'pr_material'
        prop_values = '10.0'
    []
    [kinematic_viscosity]
        type = GenericConstantMaterial
        prop_names = 'kin_visc_material'
        prop_values = '10.0'
    []
    [nusselt]
        type = Nusselt 
        L = 1.0
        nusselt_correlation = 'Dittus-Boelter'
        thermal_conductivity = fuel_therm_cond_material
        prandtl_number = pr_material
        kinematic_viscosity = kin_visc_material
        v_velocity = v_velocity

    []

[]

[BCs]
    [right]
        type = ConvectiveHTBC
        variable = temp
        boundary = right
        bulk_temperature = bulk_temp
        h_coeff_material = HeatTransferCoeff
    []
    [left]
        type = NeumannBC
        variable = temp
        boundary = left
        value = 0
    []
[]

[Executioner]
    type = Steady
[]


[Outputs]
  exodus = true
[]
