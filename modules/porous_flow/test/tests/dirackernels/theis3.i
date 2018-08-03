# Two phase Theis problem: Flow from single source
# Constant rate injection 0.5 kg/s
# 1D cylindrical mesh

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 100
  xmax = 2000
  bias_x = 1.05
[]

[Problem]
  type = FEProblem
  coord_type = RZ
  rz_coord_axis = Y
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]

[Variables]
  [./ppwater]
    initial_condition = 20e6
  [../]
  [./sgas]
    initial_condition = 0
  [../]
[]

[AuxVariables]
  [./massfrac_ph0_sp0]
    initial_condition = 1
  [../]
  [./massfrac_ph1_sp0]
    initial_condition = 0
  [../]
[]

[Kernels]
  [./mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = ppwater
  [../]
  [./flux0]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = ppwater
  [../]
  [./mass1]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = sgas
  [../]
  [./flux1]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    variable = sgas
  [../]
[]

[UserObjects]
  [./dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'ppwater sgas'
    number_fluid_phases = 2
    number_fluid_components = 2
  [../]
  [./pc]
    type = PorousFlowCapillaryPressureConst
    pc = 1e5
  [../]
[]

[Modules]
  [./FluidProperties]
    [./simple_fluid0]
      type = SimpleFluidProperties
      bulk_modulus = 2e9
      density0 = 1000
      viscosity = 1e-3
      thermal_expansion = 0
    [../]
    [./simple_fluid1]
      type = SimpleFluidProperties
      bulk_modulus = 2e9
      density0 = 10
      viscosity = 1e-4
      thermal_expansion = 0
    [../]
  [../]
[]

[Materials]
  [./temperature]
    type = PorousFlowTemperature
    at_nodes = true
  [../]
  [./temperature_qp]
    type = PorousFlowTemperature
  [../]
  [./ppss_qp]
    type = PorousFlow2PhasePS
    phase0_porepressure = ppwater
    phase1_saturation = sgas
    capillary_pressure = pc
  [../]
  [./ppss]
    type = PorousFlow2PhasePS
    at_nodes = true
    phase0_porepressure = ppwater
    phase1_saturation = sgas
    capillary_pressure = pc
  [../]
  [./massfrac]
    type = PorousFlowMassFraction
    at_nodes = true
    mass_fraction_vars = 'massfrac_ph0_sp0 massfrac_ph1_sp0'
  [../]
  [./simple_fluid0]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid0
    phase = 0
    at_nodes = true
    compute_enthalpy = false
    compute_internal_energy = false
  [../]
  [./simple_fluid0_qp]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid0
    phase = 0
    compute_enthalpy = false
    compute_internal_energy = false
  [../]
  [./simple_fluid1]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid1
    phase = 1
    at_nodes = true
    compute_enthalpy = false
    compute_internal_energy = false
  [../]
  [./simple_fluid1_qp]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid1
    phase = 1
    compute_enthalpy = false
    compute_internal_energy = false
  [../]
  [./dens_all]
    type = PorousFlowJoiner
    include_old = true
    at_nodes = true
    material_property = PorousFlow_fluid_phase_density_nodal
  [../]
  [./dens_all_at_quadpoints]
    type = PorousFlowJoiner
    material_property = PorousFlow_fluid_phase_density_qp
    at_nodes = false
  [../]
  [./visc_all]
    type = PorousFlowJoiner
    at_nodes = true
    material_property = PorousFlow_viscosity_nodal
  [../]
  [./porosity]
    type = PorousFlowPorosityConst
    at_nodes = true
    porosity = 0.2
  [../]
  [./permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1e-12 0 0 0 1e-12 0 0 0 1e-12'
  [../]
  [./relperm_water]
    type = PorousFlowRelativePermeabilityCorey
    at_nodes = true
    n = 1
    phase = 0
  [../]
  [./relperm_gas]
    type = PorousFlowRelativePermeabilityCorey
    at_nodes = true
    n = 1
    phase = 1
  [../]
  [./relperm_all]
    type = PorousFlowJoiner
    at_nodes = true
    material_property = PorousFlow_relative_permeability_nodal
  [../]
[]

[BCs]
  [./rightwater]
    type = DirichletBC
    boundary = right
    value = 20e6
    variable = ppwater
  [../]
[]

[DiracKernels]
  [./source]
    type = PorousFlowSquarePulsePointSource
    point = '0 0 0'
    mass_flux = 0.5
    variable = sgas
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol -snes_max_it'
    petsc_options_value = 'gmres      asm      lu           NONZERO                   2               1E-8       1E-10 20'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = Newton
  end_time = 2e4
  dtmax = 1e4
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    growth_factor = 2
  [../]
[]

[VectorPostprocessors]
  [./line]
    type = NodalValueSampler
    sort_by = x
    variable = 'ppwater sgas'
    execute_on = 'timestep_end'
  [../]
[]

[Postprocessors]
  [./ppwater]
    type = PointValue
    point =  '4 0 0'
    variable = ppwater
  [../]
  [./sgas]
    type = PointValue
    point = '4 0 0'
    variable = sgas
  [../]
  [./massgas]
    type = PorousFlowFluidMass
    fluid_component = 1
  [../]
[]

[Outputs]
  file_base = theis3
  print_linear_residuals = false
  print_perf_log = true
  exodus = true
  [./csv]
    type = CSV
    execute_on = timestep_end
    execute_vector_postprocessors_on = final
  [../]
[]
