[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 2
  ny = 2
[]

[Debug]
  show_material_props = true
[]

[MeshModifiers]
  [./subdomain1]
    type = SubdomainBoundingBox
    bottom_left = '0 0 0'
    top_right = '0.5 0.5 0'
    block_id = 1
  [../]
  [./subdomain2]
    type = SubdomainBoundingBox
    bottom_left = '0.5 0 0'
    top_right = '1 0.5 0'
    block_id = 2
  [../]
  [./subdomain3]
    type = SubdomainBoundingBox
    bottom_left = '0 0.5 0'
    top_right = '0.5 1 0'
    block_id = 3
  [../]
  [./subdomain4]
    type = SubdomainBoundingBox
    bottom_left = '0.5 0.5 0'
    top_right = '1 1 0'
    block_id = 4
  [../]
[]

[Variables]
  [./dummy]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./dummy]
    type = Diffusion
    variable = dummy
  [../]
[]

[BCs]
  [./dummy_left]
    type = DirichletBC
    variable = dummy
    boundary = left
    value = 0
  [../]

  [./dummy_right]
    type = DirichletBC
    variable = dummy
    boundary = right
    value = 1
  [../]
[]

[AuxVariables]
  [./var1]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./var2]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./var3]
    family = MONOMIAL
    order = CONSTANT
  [../]
[../]

[AuxKernels]
  [./var1]
    variable = var1
    type = MaterialPropertyBlockAux
    mat_prop_name = prop1
  [../]
  [./var2]
    variable = var2
    type = MaterialPropertyBlockAux
    mat_prop_name = prop2
  [../]
  [./var3]
    variable = var3
    type = MaterialRealAux
    property = prop3
    block = '2 3 4'
  [../]
[]

[Materials]
  [./mat1]
    type = GenericConstantMaterial
    block = '1 2 4'
    prop_names = 'prop1'
    prop_values = '0'
  [../]
  [./mat2]
    type = GenericConstantMaterial
    block = '2 3 4'
    prop_names = 'prop2'
    prop_values = '0'
  [../]
  [./mat3]
    type = SubdomainConstantMaterial
    block = '2 3 4'
    mat_prop_name = 'prop3'
    values = '4 2 1'
  [../]
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
[]
