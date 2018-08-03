# checking that the mass postprocessor correctly calculates the mass
# 1phase, 1component, constant porosity
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 3
  xmin = -1
  xmax = 1
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [./pp]
  [../]
[]

[ICs]
  [./pinit]
    type = FunctionIC
    function = x
    variable = pp
  [../]
[]

[Kernels]
  [./mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pp
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
    alpha = 1
  [../]
[]

[Modules]
  [./FluidProperties]
    [./simple_fluid]
      type = SimpleFluidProperties
      bulk_modulus = 1
      density0 = 1
      thermal_expansion = 0
    [../]
  [../]
[]

[Materials]
  [./temperature]
    type = PorousFlowTemperature
    at_nodes = true
  [../]
  [./ppss]
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
    at_nodes = true
    phase = 0
  [../]
  [./dens_all]
    type = PorousFlowJoiner
    include_old = true
    at_nodes = true
    material_property = PorousFlow_fluid_phase_density_nodal
  [../]
  [./porosity]
    type = PorousFlowPorosityConst
    at_nodes = true
    porosity = 0.1
  [../]
[]

[Postprocessors]
  [./total_mass]
    type = PorousFlowFluidMass
  [../]
[]

[Preconditioning]
  [./andy]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -pc_type -snes_atol -snes_rtol -snes_max_it'
    petsc_options_value = 'bcgs bjacobi 1 1 10000'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = Newton
  dt = 1
  end_time = 1
[]

[Outputs]
  execute_on = 'timestep_end'
  file_base = mass01
  csv = true
[]
