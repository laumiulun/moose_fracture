[GlobalParams]
  order = SECOND
  family = LAGRANGE
[]

[Problem]
  coord_type = RZ
[]

[Mesh]
  file = recover_in.e
[]

[Variables]
  [./temp]
    initial_condition = 580.0
  [../]
[]

[AuxVariables]
  [./gap_cond]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./heat]
    type = HeatConduction
    variable = temp
  [../]
  [./heat_source]
    type = BodyForce
    variable = temp
    block = pellet_type_1
    value = 1e-2
    function = 't'
  [../]
[]

[ThermalContact]
  [./thermal_contact]
    type = GapHeatTransfer
    variable = temp
    master = 5
    slave = 10
    quadrature = true
  [../]
[]

[BCs]
  [./outside]
    type = DirichletBC
    value = 580
    boundary = '1 2 3'
    variable = temp
  [../]
  [./edge]
    type = PresetBC
    value = 700
    boundary = 10
    variable = temp
  [../]
[]

[Materials]
  [./thermal_3]
    type = HeatConductionMaterial
    block = 3
    thermal_conductivity = 5
    specific_heat = 12
  [../]
  [./thermal_1]
    type = HeatConductionMaterial
    block = 1
    thermal_conductivity = 16.0
    specific_heat = 330.0
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'

  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = ' lu       superlu_dist'

  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-11

  start_time = -200
  n_startup_steps = 1
  end_time = 1.02e5
  num_steps = 10

  dtmax = 2e6
  dtmin = 1

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 2.0e2
    optimal_iterations = 15
    iteration_window = 2
  [../]

  [./Quadrature]
    order = FIFTH
    side_order = SEVENTH
  [../]
[]

[Postprocessors]
  [./ave_temp_interior]
     type = SideAverageValue
     boundary = 9
     variable = temp
     execute_on = 'initial linear'
  [../]
  [./avg_clad_temp]
    type = SideAverageValue
    boundary = 7
    variable = temp
    execute_on = 'initial timestep_end'
  [../]
  [./flux_from_clad]
    type = SideFluxIntegral
    variable = temp
    boundary = 5
    diffusivity = thermal_conductivity
    execute_on = timestep_end
  [../]
  [./_dt]
    type = TimestepSize
    execute_on = timestep_end
  [../]
[]

[Outputs]
  exodus = true
  [./console]
    type = Console
    perf_log = true
    max_rows = 25
  [../]
[]
