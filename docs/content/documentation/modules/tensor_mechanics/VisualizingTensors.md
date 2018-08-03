#Visualization of Tensor Components

To visualize stresses, strains, and elasticity tensor components, the material objects must be outputted to auxiliary variables using auxiliary kernels available. RankTwoAux is used to output a RankTwoTensor component, and RankFourAux is used to output a RankFourTensor component.  For example, $\sigma_{11}, \epsilon_{22}$, and $C_{1122}$ can be visualized by first declaring auxiliary variables for them in the input file.

Stress $\sigma_{xx}$ component:
!listing modules/combined/test/tests/eigenstrain/inclusion.i block=AuxVariables/s11_aux

Strain $\epsilon_{yy}$ component:
!listing modules/combined/test/tests/eigenstrain/variable.i block=AuxVariables/e22_aux

Elasticity Tensor $C_{1133}$ component:
!listing modules/combined/test/tests/elasticitytensor/composite.i block=AuxVariables/C1133_aux

Next, the appropriate auxiliary kernels are used.  These elements require the name of the material property of which you wish to see the field value and the indices of the tensor value (either 0, 1, or 2) in addition to the name of the output AuxVariable.  The corresponding AuxKernels blocks for each of the AuxVariables are given below.

Stress $\sigma_{xx}$ component:
!listing modules/combined/test/tests/eigenstrain/inclusion.i block=AuxKernels/matl_s11

Strain $\epsilon_{yy}$ component:
!listing modules/combined/test/tests/eigenstrain/variable.i block=AuxKernels/matl_e22

Elasticity Tensor $C_{1133}$ component:
!listing modules/combined/test/tests/elasticitytensor/composite.i block=AuxKernels/matl_C1133