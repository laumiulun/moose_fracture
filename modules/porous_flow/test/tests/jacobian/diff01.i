# Test the Jacobian of the diffusive component of the PorousFlowDisperiveFlux kernel.
# By setting disp_long and disp_trans to zero, the purely diffusive component of the flux
# can be isolated.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 3
  xmin = 0
  xmax = 1
  ny = 1
  ymin = 0
  ymax = 1
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [./pp]
  [../]
  [./massfrac0]
  [../]
[]

[ICs]
  [./pp]
    type = RandomIC
    variable = pp
    max = 2e1
    min = 1e1
  [../]
  [./massfrac0]
    type = RandomIC
    variable = massfrac0
    min = 0
    max = 1
  [../]
[]

[Kernels]
  [./diff0]
    type = PorousFlowDispersiveFlux
    fluid_component = 0
    variable = pp
    gravity = '1 0 0'
    disp_long = 0
    disp_trans = 0
  [../]
  [./diff1]
    type = PorousFlowDispersiveFlux
    fluid_component = 1
    variable = massfrac0
    gravity = '1 0 0'
    disp_long = 0
    disp_trans = 0
  [../]
[]

[UserObjects]
  [./dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pp massfrac0'
    number_fluid_phases = 1
    number_fluid_components = 2
  [../]
[]

[Modules]
  [./FluidProperties]
    [./simple_fluid]
      type = SimpleFluidProperties
      bulk_modulus = 1e7
      density0 = 10
      thermal_expansion = 0
      viscosity = 1
    [../]
  [../]
[]

[Materials]
  [./temp]
    type = PorousFlowTemperature
    at_nodes = false
  [../]
  [./ppss]
    type = PorousFlow1PhaseFullySaturated
    porepressure = pp
    at_nodes = false
  [../]
  [./massfrac]
    type = PorousFlowMassFraction
    mass_fraction_vars = 'massfrac0'
    at_nodes = false
  [../]
  [./simple_fluid]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid
    at_nodes = false
    phase = 0
  [../]
  [./dens_qp_all]
    type = PorousFlowJoiner
    at_nodes = false
    material_property = PorousFlow_fluid_phase_density_qp
  [../]
  [./visc_all]
    type = PorousFlowJoiner
    at_nodes = false
    material_property = PorousFlow_viscosity_qp
  [../]
  [./poro]
    type = PorousFlowPorosityConst
    porosity = 0.1
    at_nodes = false
  [../]
  [./diff]
    type = PorousFlowDiffusivityConst
    diffusion_coeff = '1e-2 1e-1'
    tortuosity = '0.1'
    at_nodes = false
  [../]
  [./permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1 0 0 0 2 0 0 0 3'
    at_nodes = false
  [../]
  [./relperm]
    type = PorousFlowRelativePermeabilityConst
    at_nodes = false
    phase = 0
  [../]
  [./relperm_all]
    type = PorousFlowJoiner
    at_nodes = false
    material_property = PorousFlow_relative_permeability_qp
  [../]
[]

[Preconditioning]
  active = smp
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = Newton
  dt = 1
  end_time = 1
[]

[Outputs]
  exodus = false
[]
