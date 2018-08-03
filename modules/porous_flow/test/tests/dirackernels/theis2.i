# Theis problem: Flow to single sink
# Constant rate injection between 700 and 1700 s.
# Final state matches theis1.i final state.
# Cartesian mesh with logarithmic distribution in x and y.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 40
  ny = 40
  bias_x = 1.1
  bias_y = 1.1
  ymax = 100
  xmax = 100
[]

[GlobalParams]
  PorousFlowDictator = dictator
  compute_enthalpy = false
  compute_internal_energy = false
[]

[Variables]
  [./pp]
  [../]
[]

[Kernels]
  [./mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pp
  [../]
  [./flux]
    type = PorousFlowAdvectiveFlux
    variable = pp
    gravity = '0 0 0'
    fluid_component = 0
  [../]
[]

[UserObjects]
  [./dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pp'
    number_fluid_phases = 1
    number_fluid_components = 1
  [../]
  [./pc]
    type = PorousFlowCapillaryPressureVG
    m = 0.5
    alpha = 1e-7
  [../]
[]

[Modules]
  [./FluidProperties]
    [./simple_fluid]
      type = SimpleFluidProperties
      bulk_modulus = 2e9
      density0 = 1000
      viscosity = 0.001
      thermal_expansion = 0
    [../]
  [../]
[]

[Materials]
  [./temperature_nodal]
    type = PorousFlowTemperature
    at_nodes = true
  [../]
  [./temperature]
    type = PorousFlowTemperature
  [../]
  [./ppss]
    type = PorousFlow1PhaseP
    porepressure = pp
    capillary_pressure = pc
  [../]
  [./ppss_nodal]
    type = PorousFlow1PhaseP
    at_nodes = true
    porepressure = pp
    capillary_pressure = pc
  [../]
  [./massfrac]
    type = PorousFlowMassFraction
    at_nodes = true
  [../]
  [./simple_fluid]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid
    phase = 0
    at_nodes = true
  [../]
  [./simple_fluid_qp]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid
    phase = 0
  [../]
  [./dens_all]
    type = PorousFlowJoiner
    at_nodes = true
    material_property = PorousFlow_fluid_phase_density_nodal
    include_old = true
  [../]
  [./dens_all_at_quadpoints]
    type = PorousFlowJoiner
    material_property = PorousFlow_fluid_phase_density_qp
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
    permeability = '1E-14 0 0 0 1E-14 0 0 0 1E-14'
  [../]
  [./relperm]
    type = PorousFlowRelativePermeabilityCorey
    at_nodes = true
    n = 0
    phase = 0
  [../]
  [./relperm_all]
    type = PorousFlowJoiner
    at_nodes = true
    material_property = PorousFlow_relative_permeability_nodal
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = Newton
  dt = 100
  end_time = 1700
  nl_abs_tol = 1e-10
[]

[Outputs]
  print_perf_log = true
  file_base = theis2
  csv = true
  execute_on = 'final'
  [./con]
    output_linear = true
    type = Console
  [../]
[]

[ICs]
  [./PressureIC]
    variable = pp
    type = ConstantIC
    value = 20e6
  [../]
[]

[DiracKernels]
  [./sink]
    type = PorousFlowSquarePulsePointSource
    start_time = 700
    end_time = 1700
    point = '0 0 0'
    mass_flux = -0.04
    variable = pp
  [../]
[]

[BCs]
  [./right]
    type = DirichletBC
    variable = pp
    value = 20e6
    boundary = right
  [../]
  [./top]
    type = DirichletBC
    variable = pp
    value = 20e6
    boundary = top
  [../]
[]

[VectorPostprocessors]
  [./pressure]
    type = SideValueSampler
    variable = pp
    sort_by = x
    execute_on = timestep_end
    boundary = bottom
  [../]
[]
