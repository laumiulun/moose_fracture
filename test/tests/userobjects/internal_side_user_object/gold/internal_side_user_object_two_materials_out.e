CDF      
      
len_string     !   len_line   Q   four      	time_step          len_name   !   num_dim       	num_nodes      	   num_elem      
num_el_blk        num_node_sets         num_side_sets         num_el_in_blk1        num_nod_per_el1       num_el_in_blk2        num_nod_per_el2       num_side_ss1      num_side_ss2      num_side_ss3      num_side_ss4      num_nod_ns1       num_nod_ns2       num_nod_ns3       num_nod_ns4       num_nod_var       num_glo_var       num_info  �         api_version       @�
=   version       @�
=   floating_point_word_size            	file_size               int64_status             title         .internal_side_user_object_two_materials_out.e      maximum_name_length                 "   
time_whole                            ��   	eb_status                             	�   eb_prop1               name      ID              	�   	ns_status         	                    	�   ns_prop1      	         name      ID              	�   	ss_status         
                    
   ss_prop1      
         name      ID              
   coordx                      H      
$   coordy                      H      
l   eb_names                       D      
�   ns_names      	                 �      
�   ss_names      
                 �      |   
coor_names                         D          node_num_map                    $      D   connect1                  	elem_type         QUAD4               h   connect2                  	elem_type         QUAD4         0      x   elem_num_map                          �   elem_ss1                          �   side_ss1                          �   elem_ss2                          �   side_ss2                          �   elem_ss3                          �   side_ss3                          �   elem_ss4                          �   side_ss4                          �   node_ns1                          �   node_ns2                             node_ns3                             node_ns4                             vals_nod_var1                          H      ��   name_nod_var                       $      (   name_glo_var                       $      L   vals_glo_var                             �   info_records                      |\      p                                                               ��                      ��      ?�      ?�              ��      ?�      ��      ��                      ��              ?�      ?�      ?�                                                                          left                             bottom                           right                            top                              bottom                           left                             right                            top                                                                                                                             	                                             	                                                                                          	         	u                                   value                               ####################                                                             # Created by MOOSE #                                                             ####################                                                             ### Command Line Arguments ###                                                   -i                                                                               internal_side_user_object_two_materials.i                                                                                                                         ### Version Info ###                                                             Framework Information:                                                           MOOSE version:           git commit 6b35ab8 on 2015-09-08                        PETSc Version:           3.6.0                                                   Current Time:            Wed Oct  7 15:47:05 2015                                Executable Timestamp:    Wed Oct  7 15:47:00 2015                                                                                                                                                                                                  ### Input File ###                                                                                                                                                []                                                                                 name                           =                                                 initial_from_file_timestep     = 2                                               initial_from_file_var          = INVALID                                         block                          = INVALID                                         coord_type                     = XYZ                                             fe_cache                       = 0                                               kernel_coverage_check          = 1                                               material_coverage_check        = 1                                               rz_coord_axis                  = Y                                               type                           = FEProblem                                       use_legacy_uo_aux_computation  = INVALID                                         use_legacy_uo_initialization   = INVALID                                         element_order                  = AUTO                                            order                          = AUTO                                            side_order                     = AUTO                                            active_bcs                     = INVALID                                         active_kernels                 = INVALID                                         inactive_bcs                   = INVALID                                         inactive_kernels               = INVALID                                         start                          = 0                                               dimNearNullSpace               = 0                                               dimNullSpace                   = 0                                               error_on_jacobian_nonzero_reallocation = 0                                       petsc_inames                   =                                                 petsc_options                  = INVALID                                         petsc_values                   =                                                 solve                          = 1                                               use_nonlinear                  = 1                                             []                                                                                                                                                                [BCs]                                                                                                                                                               [./all]                                                                            boundary                     = '0 1 2 3'                                         implicit                     = 1                                                 name                         = BCs/all                                           type                         = FunctionDirichletBC                               use_displaced_mesh           = 0                                                 variable                     = u                                                 diag_save_in                 = INVALID                                           function                     = fn_exact                                          save_in                      = INVALID                                           seed                         = 0                                               [../]                                                                          []                                                                                                                                                                [Executioner]                                                                      name                           = Executioner                                     type                           = Steady                                          compute_initial_residual_before_preset_bcs = 0                                   l_abs_step_tol                 = -1                                              l_max_its                      = 10000                                           l_tol                          = 1e-05                                           line_search                    = default                                         nl_abs_step_tol                = 1e-50                                           nl_abs_tol                     = 1e-50                                           nl_max_funcs                   = 10000                                           nl_max_its                     = 50                                              nl_rel_step_tol                = 1e-50                                           nl_rel_tol                     = 1e-08                                           no_fe_reinit                   = 0                                               petsc_options                  = INVALID                                         petsc_options_iname            = INVALID                                         petsc_options_value            = INVALID                                         solve_type                     = INVALID                                         restart_file_base              =                                                 splitting                      = INVALID                                       []                                                                                                                                                                [Executioner]                                                                      _fe_problem                    = 0x7fc0ca813a00                                []                                                                                                                                                                [Functions]                                                                                                                                                         [./ffn]                                                                            name                         = Functions/ffn                                     type                         = ParsedFunction                                    vals                         = INVALID                                           value                        = -4                                                vars                         = INVALID                                         [../]                                                                                                                                                             [./fn_exact]                                                                       name                         = Functions/fn_exact                                type                         = ParsedFunction                                    vals                         = INVALID                                           value                        = x*x+y*y                                           vars                         = INVALID                                         [../]                                                                          []                                                                                                                                                                [Kernels]                                                                                                                                                           [./diff]                                                                           name                         = Kernels/diff                                      type                         = Diffusion                                         block                        = INVALID                                           diag_save_in                 = INVALID                                           implicit                     = 1                                                 save_in                      = INVALID                                           seed                         = 0                                                 use_displaced_mesh           = 0                                                 variable                     = u                                               [../]                                                                                                                                                             [./ffn]                                                                            name                         = Kernels/ffn                                       type                         = UserForcingFunction                               block                        = INVALID                                           diag_save_in                 = INVALID                                           function                     = ffn                                               implicit                     = 1                                                 save_in                      = INVALID                                           seed                         = 0                                                 use_displaced_mesh           = 0                                                 variable                     = u                                               [../]                                                                          []                                                                                                                                                                [Materials]                                                                                                                                                         [./stateful1]                                                                      name                         = Materials/stateful1                               type                         = StatefulMaterial                                  block                        = 0                                                 boundary                     = INVALID                                           initial_diffusivity          = 1                                                 output_properties            = INVALID                                           outputs                      = none                                              seed                         = 0                                                 use_displaced_mesh           = 0                                               [../]                                                                                                                                                             [./stateful2]                                                                      name                         = Materials/stateful2                               type                         = StatefulMaterial                                  block                        = 1                                                 boundary                     = INVALID                                           initial_diffusivity          = 2                                                 output_properties            = INVALID                                           outputs                      = none                                              seed                         = 0                                                 use_displaced_mesh           = 0                                               [../]                                                                          []                                                                                                                                                                [Mesh]                                                                             name                           = Mesh                                            displacements                  = INVALID                                         block_id                       = INVALID                                         block_name                     = INVALID                                         boundary_id                    = INVALID                                         boundary_name                  = INVALID                                         construct_side_list_from_node_list = 0                                           ghosted_boundaries             = INVALID                                         ghosted_boundaries_inflation   = INVALID                                         patch_size                     = 40                                              second_order                   = 0                                               skip_partitioning              = 0                                               type                           = FileMesh                                        uniform_refine                 = 0                                               centroid_partitioner_direction = INVALID                                         dim                            = 3                                               distribution                   = DEFAULT                                         file                           = TwoBlockMesh.e                                  nemesis                        = 0                                               partitioner                    = default                                         patch_update_strategy          = never                                         []                                                                                                                                                                [Outputs]                                                                          additional_output_on           = INVALID                                         checkpoint                     = 0                                               color                          = 1                                               console                        = 1                                               csv                            = 0                                               dofmap                         = 0                                               execute_on                     = 'INITIAL TIMESTEP_END'                          exodus                         = 1                                               file_base                      = INVALID                                         gmv                            = 0                                               gnuplot                        = 0                                               hide                           = INVALID                                         interval                       = 1                                               name                           = Outputs                                         nemesis                        = 0                                               output_failed                  = 0                                               output_final                   = 0                                               output_if_base_contains        = INVALID                                         output_initial                 = 0                                               output_intermediate            = 1                                               output_on                      = TIMESTEP_END                                    output_timestep_end            = 1                                               print_linear_residuals         = 1                                               print_mesh_changed_info        = 0                                               print_perf_log                 = 0                                               show                           = INVALID                                         solution_history               = 0                                               sync_times                     =                                                 tecplot                        = 0                                               vtk                            = 0                                               xda                            = 0                                               xdr                            = 0                                                                                                                                [./console]                                                                        name                         = Outputs/console                                   type                         = Console                                           additional_execute_on        = INVALID                                           additional_output_on         = INVALID                                           all_variable_norms           = 0                                                 append_displaced             = 0                                                 append_restart               = 0                                                 end_time                     = INVALID                                           execute_elemental_variables  = 1                                                 execute_input                = 1                                                 execute_input_on             = INVALID                                           execute_nodal_variables      = 1                                                 execute_on                   = 'FAILED INITIAL LINEAR NONLINEAR TIMESTEP_... BEGIN TIMESTEP_END'                                                                  execute_postprocessors_on    = TIMESTEP_END                                      execute_scalar_variables     = 1                                                 execute_scalars_on           = TIMESTEP_END                                      execute_system_information   = 1                                                 execute_system_information_on = INITIAL                                          execute_vector_postprocessors = 1                                                execute_vector_postprocessors_on = TIMESTEP_END                                  file_base                    = INVALID                                           fit_mode                     = ENVIRONMENT                                       hide                         = INVALID                                           interval                     = 1                                                 linear_residual_dt_divisor   = 1000                                              linear_residual_end_time     = INVALID                                           linear_residual_start_time   = INVALID                                           linear_residuals             = 0                                                 max_rows                     = 15                                                nonlinear_residual_dt_divisor = 1000                                             nonlinear_residual_end_time  = INVALID                                           nonlinear_residual_start_time = INVALID                                          nonlinear_residuals          = 0                                                 outlier_multiplier           = '0.8 2'                                           outlier_variable_norms       = 1                                                 output_failed                = 0                                                 output_file                  = 0                                                 output_final                 = 0                                                 output_if_base_contains      =                                                   output_initial               = 0                                                 output_intermediate          = 1                                                 output_linear                = 0                                                 output_nonlinear             = 0                                                 output_on                    = TIMESTEP_END                                      output_postprocessors        = 1                                                 output_screen                = 1                                                 output_timestep_end          = 1                                                 padding                      = 4                                                 perf_header                  = INVALID                                           perf_log                     = 0                                                 print_mesh_changed_info      = 0                                                 scientific_time              = 0                                                 setup_log                    = INVALID                                           setup_log_early              = 0                                                 show                         = INVALID                                           show_multiapp_name           = 0                                                 solve_log                    = INVALID                                           start_time                   = INVALID                                           sync_only                    = 0                                                 sync_times                   =                                                   system_info                  = 'AUX EXECUTION FRAMEWORK MESH NONLINEAR'          time_precision               = INVALID                                           time_tolerance               = 1e-14                                             use_displaced                = 0                                                 verbose                      = 0                                               [../]                                                                                                                                                             [./exodus]                                                                         name                         = Outputs/exodus                                    type                         = Exodus                                            additional_execute_on        = INVALID                                           additional_output_on         = INVALID                                           append_displaced             = 0                                                 append_oversample            = 0                                                 elemental_as_nodal           = 0                                                 end_time                     = INVALID                                           execute_elemental_on         = INVALID                                           execute_elemental_variables  = 1                                                 execute_input                = 1                                                 execute_input_on             = INITIAL                                           execute_nodal_on             = INVALID                                           execute_nodal_variables      = 1                                                 execute_on                   = 'INITIAL TIMESTEP_END'                            execute_postprocessors_on    = INVALID                                           execute_scalar_variables     = 1                                                 execute_scalars_on           = INVALID                                           execute_system_information   = 1                                                 execute_vector_postprocessors = 1                                                file                         = INVALID                                           file_base                    = INVALID                                           hide                         = INVALID                                           interval                     = 1                                                 linear_residual_dt_divisor   = 1000                                              linear_residual_end_time     = INVALID                                           linear_residual_start_time   = INVALID                                           linear_residuals             = 0                                                 nonlinear_residual_dt_divisor = 1000                                             nonlinear_residual_end_time  = INVALID                                           nonlinear_residual_start_time = INVALID                                          nonlinear_residuals          = 0                                                 output_failed                = 0                                                 output_final                 = 0                                                 output_if_base_contains      =                                                   output_initial               = 0                                                 output_intermediate          = 1                                                 output_linear                = 0                                                 output_material_properties   = 0                                                 output_nonlinear             = 0                                                 output_on                    = TIMESTEP_END                                      output_postprocessors        = 1                                                 output_timestep_end          = 1                                                 oversample                   = 0                                                 padding                      = 3                                                 position                     = INVALID                                           refinements                  = 0                                                 scalar_as_nodal              = 0                                                 sequence                     = INVALID                                           show                         = INVALID                                           show_material_properties     = INVALID                                           start_time                   = INVALID                                           sync_only                    = 0                                                 sync_times                   =                                                   time_tolerance               = 1e-14                                             use_displaced                = 0                                               [../]                                                                          []                                                                                                                                                                [Postprocessors]                                                                                                                                                    [./value]                                                                          name                         = Postprocessors/value                              type                         = InsideValuePPS                                    execute_on                   = TIMESTEP_END                                      outputs                      = INVALID                                           use_displaced_mesh           = 0                                                 user_object                  = isuo                                            [../]                                                                          []                                                                                                                                                                [UserObjects]                                                                                                                                                       [./isuo]                                                                           name                         = UserObjects/isuo                                  type                         = InsideUserObject                                  block                        = INVALID                                           diffusivity                  = diffusivity                                       execute_on                   = 'INITIAL TIMESTEP_END'                            use_displaced_mesh           = 0                                                 use_old_prop                 = 0                                                 variable                     = u                                               [../]                                                                          []                                                                                                                                                                [Variables]                                                                                                                                                         [./u]                                                                              block                        = INVALID                                           eigen                        = 0                                                 family                       = LAGRANGE                                          initial_condition            = INVALID                                           name                         = Variables/u                                       order                        = FIRST                                             outputs                      = INVALID                                           scaling                      = 1                                                 initial_from_file_timestep   = 2                                                 initial_from_file_var        = INVALID                                         [../]                                                                          []                                                                                                                                                                          ?�      @       ?�      �      ?�      @       ?�      ?�      @       @       @*����