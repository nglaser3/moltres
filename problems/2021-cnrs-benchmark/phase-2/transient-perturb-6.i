density = .002  # kg cm-3
cp = 3075       # J kg-1 K-1, 6.15 / 2.0e-3
k = .005        # W cm-1 K-1
gamma = 1       # W cm-3 K-1, Volumetric heat transfer coefficient
viscosity = .5  # dynamic viscosity
alpha = 1       # SUPG stabilization parameter
t_alpha = 2e-4  # K-1, Thermal expansion coefficient

[GlobalParams]
  num_groups = 6
  num_precursor_groups = 8
  use_exp_form = false
  group_fluxes = 'group1 group2 group3 group4 group5 group6'
  pre_concs = 'pre1 pre2 pre3 pre4 pre5 pre6 pre7 pre8'
  temperature = temp
  sss2_input = true
  account_delayed = true
  integrate_p_by_parts = true
[../]

[Mesh]
  [./file_mesh]
    type = FileMeshGenerator
    file = '../full-coupling/u5p10_out_ntsApp0_exodus.e'
    use_for_exodus_restart = true
  [../]
  [./corner_node]
    type = ExtraNodesetGenerator
    input = file_mesh
    new_boundary = 'pinned_node'
    coord = '200 200'
  [../]
[]

[Problem]
  type = FEProblem
[]

[Variables]
  [./temp]
    family = LAGRANGE
    order = FIRST
    scaling = 1e-2
  [../]
  [./vel]
    family = LAGRANGE_VEC
    order = FIRST
    scaling = 1e2
  [../]
  [./p]
    family = LAGRANGE
    order = FIRST
    scaling = 1e2
  [../]
[]

[Nt]
  var_name_base = group
  vacuum_boundaries = 'bottom left right top'
  create_temperature_var = false
  init_nts_from_file = true
  scaling = 1e-8
  eigenvalue_scaling = 0.9927821803
[]

[Precursors]
  [./pres]
    var_name_base = pre
    outlet_boundaries = ''
    constant_velocity_values = false
    uvel = vel_x
    vvel = vel_y
    nt_exp_form = false
    family = MONOMIAL
    order = CONSTANT
    transient = true
    scaling = 1e-8
  [../]
[]

[AuxVariables]
  [./vel_x]
    family = LAGRANGE
    order = FIRST
  [../]
  [./vel_y]
    family = LAGRANGE
    order = FIRST
  [../]
  [./heat]
    family = MONOMIAL
    order = FIRST
  [../]
[]

[Kernels]
  [./mass]
    type = INSADMass
    variable = p
  [../]
  [./mass_pspg]
    type = INSADMassPSPG
    variable = p
  [../]
  [./momentum_time]
    type = INSADMomentumTimeDerivative
    variable = vel
  [../]
  [./momentum_advection]
    type = INSADMomentumAdvection
    variable = vel
  [../]
  [./momentum_viscous]
    type = INSADMomentumViscous
    variable = vel
  [../]
  [./momentum_pressure]
    type = INSADMomentumPressure
    variable = vel
    p = p
  [../]
  [./momentum_supg]
    type = INSADMomentumSUPG
    variable = vel
    velocity = vel
  [../]
  [./buoyancy]
    type = INSADBoussinesqBodyForce
    variable = vel
    # gravity vector, cm s-2
    gravity = '0 -981 0'
    alpha_name = 't_alpha'
    ref_temp = 'temp_ref'
  [../]
  [./gravity]
    type = INSADGravityForce
    variable = vel
    # gravity vector, cm s-2
    gravity = '0 -981 0'
  [../]

  [./temp_time]
    type = INSADHeatConductionTimeDerivative
    variable = temp
  [../]
  [./temp_source]
    type = INSADEnergySource
    variable = temp
    source_variable = heat
  [../]
  [./temp_advection]
    type = INSADEnergyAdvection
    variable = temp
  [../]
  [./temp_conduction]
    type = ADHeatConduction
    variable = temp
    thermal_conductivity = 'k'
  [../]
  [./temp_supg]
    type = INSADEnergySUPG
    variable = temp
    velocity = vel
  [../]
  [./temp_sink]
    type = ManuHX
    variable = temp
    htc = ${gamma}
    tref = 900
  [../]
[]

[AuxKernels]
  [./vel_x]
    type = VectorVariableComponentAux
    variable = vel_x
    vector_variable = vel
    component = 'x'
  [../]
  [./vel_y]
    type = VectorVariableComponentAux
    variable = vel_y
    vector_variable = vel
    component = 'y'
  [../]
  [./heat_source]
    type = FissionHeatSourceTransientAux
    variable = heat
  [../]
[]

[UserObjects]
  [./initial_th]
    type = SolutionUserObject
    mesh = '../full-coupling/u5p10_exodus.e'
    system_variables = 'vel_x vel_y p temp'
    timestep = LATEST
    execute_on = INITIAL
  [../]
  [./initial_pre]
    type = SolutionUserObject
    mesh = '../full-coupling/u5p10_out_ntsApp0_exodus.e'
    system_variables = 'pre1 pre2 pre3 pre4 pre5 pre6 pre7 pre8 heat'
    timestep = LATEST
    execute_on = INITIAL
  [../]
[]

[Controls]
  [./func_control]
    type = RealFunctionControl
    parameter = 'Kernels/temp_sink/htc'
    function = func_alpha
    execute_on = 'initial timestep_begin'
  [../]
[]

[Functions]
  [./func_alpha]
    type = ParsedFunction
    value = '1 + 0.1 * sin(2*pi*t*.4)'  # Perturbation frequency = 0.4Hz
  [../]
  [./pre1f]
    type = SolutionFunction
    from_variable = pre1
    solution = initial_pre
    scale_factor = 7.66908e+17
  [../]
  [./pre2f]
    type = SolutionFunction
    from_variable = pre2
    solution = initial_pre
    scale_factor = 7.66908e+17
  [../]
  [./pre3f]
    type = SolutionFunction
    from_variable = pre3
    solution = initial_pre
    scale_factor = 7.66908e+17
  [../]
  [./pre4f]
    type = SolutionFunction
    from_variable = pre4
    solution = initial_pre
    scale_factor = 7.66908e+17
  [../]
  [./pre5f]
    type = SolutionFunction
    from_variable = pre5
    solution = initial_pre
    scale_factor = 7.66908e+17
  [../]
  [./pre6f]
    type = SolutionFunction
    from_variable = pre6
    solution = initial_pre
    scale_factor = 7.66908e+17
  [../]
  [./pre7f]
    type = SolutionFunction
    from_variable = pre7
    solution = initial_pre
    scale_factor = 7.66908e+17
  [../]
  [./pre8f]
    type = SolutionFunction
    from_variable = pre8
    solution = initial_pre
    scale_factor = 7.66908e+17
  [../]
  [./velxf]
    type = SolutionFunction
    from_variable = vel_x
    solution = initial_th
  [../]
  [./velyf]
    type = SolutionFunction
    from_variable = vel_y
    solution = initial_th
  [../]
  [./pf]
    type = SolutionFunction
    from_variable = p
    solution = initial_th
  [../]
  [./tempf]
    type = SolutionFunction
    from_variable = temp
    solution = initial_th
  [../]
  [./heatf]
    type = SolutionFunction
    from_variable = heat
    solution = initial_pre
  [../]
[]

[ICs]
  [./vel_ic]
    type = VectorFunctionIC
    variable = vel
    function_x = velxf
    function_y = velyf
  [../]
  [./pre1_ic]
    type = FunctionIC
    variable = pre1
    function = pre1f
  [../]
  [./pre2_ic]
    type = FunctionIC
    variable = pre2
    function = pre2f
  [../]
  [./pre3_ic]
    type = FunctionIC
    variable = pre3
    function = pre3f
  [../]
  [./pre4_ic]
    type = FunctionIC
    variable = pre4
    function = pre4f
  [../]
  [./pre5_ic]
    type = FunctionIC
    variable = pre5
    function = pre5f
  [../]
  [./pre6_ic]
    type = FunctionIC
    variable = pre6
    function = pre6f
  [../]
  [./pre7_ic]
    type = FunctionIC
    variable = pre7
    function = pre7f
  [../]
  [./pre8_ic]
    type = FunctionIC
    variable = pre8
    function = pre8f
  [../]
  [./p_ic]
    type = FunctionIC
    variable = p
    function = pf
  [../]
  [./temp_ic]
    type = FunctionIC
    variable = temp
    function = tempf
  [../]
  [./heat_ic]
    type = FunctionIC
    variable = heat
    function = heatf
  [../]
[]

[BCs]
  [./no_slip]
    type = VectorDirichletBC
    variable = vel
    boundary = 'bottom left right'
    values = '0 0 0'
  [../]
  [./lid]
    type = VectorDirichletBC
    variable = vel
    boundary = 'top'
    values = '50 0 0'
  [../]
  [./pressure_pin]
    type = DirichletBC
    variable = p
    boundary = 'pinned_node'
    value = 0
  [../]
[]

[Materials]
  [./fuel]
    type = GenericMoltresMaterial
    property_tables_root = '../neutron-data/benchmark_'
    interp_type = 'linear'
    prop_names = 'temp_ref'
    prop_values = '900'
  [../]
  [./ad_mat]
    type = ADGenericConstantMaterial
    prop_names = 'k rho cp mu t_alpha'
    prop_values = '${k} ${density} ${cp} ${viscosity} ${t_alpha}'
  [../]
  [./ins_temp]
    type = INSADStabilized3Eqn
    alpha = ${alpha}
    velocity = vel
    pressure = p
    temperature = temp
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  num_steps = 5000

  solve_type = 'NEWTON'
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_linesearch_monitor'
  petsc_options_iname = '-pc_type -sub_pc_type -ksp_gmres_restart -pc_gasm_overlap -sub_pc_factor_shift_type -pc_gasm_blocks -sub_pc_factor_mat_solver_type'
  petsc_options_value = 'gasm     lu           200                1                NONZERO                   64             superlu_dist'
  line_search = 'none'

  nl_abs_tol = 1e-6
  nl_rel_tol = 1e-6
  nl_max_its = 25
  l_max_its = 400
  l_tol = 1e-4

  dt = 0.0125
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./memory]
    type = MemoryUsage
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./tot_fissions]
    type = ElmIntegTotFissPostprocessor
    execute_on = linear
  [../]
  [./powernorm]
    type = ElmIntegTotFissHeatPostprocessor
    execute_on = linear
  [../]
  [./group1norm]
    type = ElementIntegralVariablePostprocessor
    variable = group1
    execute_on = linear
  [../]
  [./group1max]
    type = NodalMaxValue
    variable = group1
    execute_on = timestep_end
  [../]
  [./group2norm]
    type = ElementIntegralVariablePostprocessor
    variable = group2
    execute_on = linear
  [../]
  [./group2max]
    type = NodalMaxValue
    variable = group2
    execute_on = timestep_end
  [../]
[]

[VectorPostprocessors]
[]

[Outputs]
  perf_graph = true
  [./exodus]
    type = Exodus
    interval = 20
  [../]
  [./csv]
    type = CSV
    append_restart = true
  [../]
  [./checkpoint]
    type = Checkpoint
    num_files = 3
    interval = 20
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]
