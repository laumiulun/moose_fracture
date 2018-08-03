# This test verifies that if an Initial Condition depends on a UO and is
# set to "initial" that it will be executed _BEFORE_ the initial condition

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 8
  ny = 8

  # We are testing geometric ghosted functors
  # so we have to use distributed mesh
  parallel_type = distributed

  # We are testing these two parameters
  num_ghosted_layers = 1
  ghost_point_neighbors = false
[]

[Variables]
  [./ghost]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[ICs]
  [./ghost_ic]
    type = GhostUserObjectIC
    variable = ghost
    ghost_uo = ghost_uo
  [../]
[]

[UserObjects]
  [./ghost_uo]
    type = GhostUserObject
    execute_on = initial
  [../]
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
[]

[Problem]
  solve = false
  kernel_coverage_check = false
[]
