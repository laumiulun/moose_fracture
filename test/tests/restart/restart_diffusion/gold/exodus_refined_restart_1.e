CDF      
      
len_string     !   len_line   Q   four      	time_step          len_name   !   num_dim       	num_nodes      Q   num_elem   @   
num_el_blk        num_node_sets         num_side_sets         num_el_in_blk1     @   num_nod_per_el1       num_side_ss1      num_side_ss2      num_nod_ns1    	   num_nod_ns2    	   num_nod_var       num_info   �         api_version       @��H   version       @��H   floating_point_word_size            	file_size               title         exodus_restart_1.e     maximum_name_length                    
time_whole                            X�   	eb_status                                eb_prop1               name      ID                 	ns_status         	                        ns_prop1      	         name      ID              (   	ss_status         
                    0   ss_prop1      
         name      ID              8   coordx                     �      @   coordy                     �      	�   eb_names                       $      P   ns_names      	                 D      t   ss_names      
                 D      �   
coor_names                         D      �   connect1                  	elem_type         QUAD4               @   elem_num_map                          @   elem_ss1                           @   side_ss1                           `   elem_ss2                           �   side_ss2                           �   node_ns1                    $      �   node_ns2                    $      �   vals_nod_var1                         �      X�   name_nod_var                       $         info_records                      E�      ,                                      ?�      ?�              ?�      ?�      ?�              ?�      ?�      ?�              ?�      ?�      ?�      ?�      ?�      ?�      ?�              ?�      ?�      ?�      ?�      ?�      ?�      ?�              ?�      ?�      ?�      ?�      ?�      ?�      ?�              ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�              ?�      ?�      ?�      ?�      ?�      ?�              ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�                      ?�      ?�              ?�      ?�      ?�      ?�              ?�      ?�      ?�      ?�              ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�              ?�      ?�      ?�      ?�              ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�              ?�      ?�      ?�              ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�      ?�                                                                                                                                                                                                                                                                        
                                 
                   !              "          !      "         #   $         %   #   $   #   &      #   %      &      "   '   %   "      (   '   %   '   )      '   (      )      *   +   !   *      ,   +   !   +   -      +   ,      -      .   /   ,   .      0   /   ,   /   1      /   0      1      -   2   (   -      3   2   (   2   4      2   3      4      1   5   3   1      6   5   3   5   7      5   6      7      &   8   9   &      :   8   9   8   ;      8   :      ;      )   <   :   )      =   <   :   <   >      <   =      >      ;   ?   @   ;      A   ?   @   ?   B      ?   A      B      >   C   A   >      D   C   A   C   E      C   D      E      4   F   =   4      G   F   =   F   H      F   G      H      7   I   G   7      J   I   G   I   K      I   J      K      H   L   D   H      M   L   D   L   N      L   M      N      K   O   M   K      P   O   M   O   Q      O   P   	   Q                                        !   "   #   $   %   &   '   (   )   *   +   ,   -   .   /   0   1   2   3   4   5   6   7   8   9   :   ;   <   =   >   ?   @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R   S   T         	      !   #   )   +                                        6   8   >   @                                    $            9   @         	   0   J   P         6u                                   ####################                                                             # Created by MOOSE #                                                             ####################                                                             ### Command Line Arguments ###                                                   -i                                                                               exodus_refined_restart_1_test.i                                                                                                                                   ### Input File ###                                                                                                                                                                                                                                 [Mesh]                                                                             action                         = setup_mesh                                      construct_side_list_from_node_list = 0                                           name                           = Mesh                                            parser_handle                  = 0x7fff60d330c8                                  second_order                   = 0                                               unique_id                      = 1                                               file                           = square.e                                        isObjectAction                 = 1                                               nemesis                        = 0                                               patch_size                     = 40                                              skip_partitioning              = 0                                               type                           = MooseMesh                                       _dimension                     = 1                                               uniform_refine                 = 2                                             []                                                                                                                                                                [Variables]                                                                        [./u]                                                                              action                       = add_variable                                      family                       = LAGRANGE                                          initial_condition            = 0                                                 name                         = Variables/u                                       order                        = FIRST                                             parser_handle                = 0x7fff60d330c8                                    scaling                      = 1                                                 unique_id                    = 10                                                initial_from_file_timestep   = 2                                               [../]                                                                                                                                                           [Variables]                                                                        action                         = no_action                                       active                         = u                                               name                           = Variables                                       parser_handle                  = 0x7fff60d330c8                                  unique_id                      = 21                                            []                                                                                                                                                                [Kernels]                                                                          [./diff]                                                                           action                       = add_kernel                                        isObjectAction               = 1                                                 name                         = Kernels/diff                                      parser_handle                = 0x7fff60d330c8                                    type                         = Diffusion                                         unique_id                    = 15                                                built_by_action              = add_kernel                                        execute_on                   = residual                                          start_time                   = -1.79769e+308                                     stop_time                    = 1.79769e+308                                      use_displaced_mesh           = 0                                                 variable                     = u                                               [../]                                                                                                                                                           [Kernels]                                                                          action                         = no_action                                       active                         = diff                                            name                           = Kernels                                         parser_handle                  = 0x7fff60d330c8                                  unique_id                      = 21                                            []                                                                                                                                                                [BCs]                                                                              action                         = no_action                                       active                         = 'left right'                                    name                           = BCs                                             parser_handle                  = 0x7fff60d330c8                                  unique_id                      = 21                                                                                                                               [./right]                                                                          action                       = add_bc                                            isObjectAction               = 1                                                 name                         = BCs/right                                         parser_handle                = 0x7fff60d330c8                                    type                         = DirichletBC                                       unique_id                    = 20                                                boundary                     = 2                                                 built_by_action              = add_bc                                            execute_on                   = residual                                          use_displaced_mesh           = 0                                                 value                        = 1                                                 variable                     = u                                               [../]                                                                                                                                                             [./left]                                                                           action                       = add_bc                                            isObjectAction               = 1                                                 name                         = BCs/left                                          parser_handle                = 0x7fff60d330c8                                    type                         = DirichletBC                                       unique_id                    = 20                                                boundary                     = 1                                                 built_by_action              = add_bc                                            execute_on                   = residual                                          use_displaced_mesh           = 0                                                 value                        = 0                                                 variable                     = u                                               [../]                                                                          []                                                                                                                                                                [Executioner]                                                                      action                         = setup_executioner                               isObjectAction                 = 1                                               l_abs_step_tol                 = -1                                              l_max_its                      = 10000                                           l_tol                          = 1e-05                                           name                           = Executioner                                     nl_abs_step_tol                = 1e-50                                           nl_abs_tol                     = 1e-50                                           nl_max_funcs                   = 10000                                           nl_max_its                     = 50                                              nl_rel_step_tol                = 1e-50                                           nl_rel_tol                     = 1e-08                                           no_fe_reinit                   = 0                                               parser_handle                  = 0x7fff60d330c8                                  petsc_options                  = -snes_mf_operator                               scheme                         = backward-euler                                  type                           = Steady                                          unique_id                      = 4                                               _mesh                          = 0x7f95284254e0                                  built_by_action                = setup_executioner                             []                                                                                                                                                                [Output]                                                                           action                         = setup_output                                    elemental_as_nodal             = 0                                               exodus                         = 1                                               exodus_inputfile_output        = 1                                               file_base                      = exodus_restart_1                                gmv                            = 0                                               gnuplot_format                 = ps                                              interval                       = 1                                               iteration_plot_start_time      = 1.79769e+308                                    max_pps_rows_screen            = 0                                               name                           = Output                                          nemesis                        = 0                                               output_displaced               = 0                                               output_initial                 = 1                                               output_solution_history        = 0                                               parser_handle                  = 0x7fff60d330c8                                  perf_log                       = 1                                               postprocessor_csv              = 0                                               postprocessor_gnuplot          = 0                                               postprocessor_screen           = 1                                               print_linear_residuals         = 0                                               screen_interval                = 1                                               show_setup_log_early           = 0                                               tecplot                        = 0                                               tecplot_binary                 = 0                                               unique_id                      = 7                                               xda                            = 0                                             []                                                                                                                                                                [init_problem]                                                                     action                         = init_problem                                    name                           = init_problem                                    parser_handle                  = 0x7fff60d330c8                                  unique_id                      = 43                                            []                                                                                                                                                                [copy_nodal_vars]                                                                  action                         = copy_nodal_vars                                 initial_from_file_timestep     = 2                                               name                           = copy_nodal_vars                                 parser_handle                  = 0x7fff60d330c8                                  unique_id                      = 11                                            []                                                                                                                                                                [check_integrity]                                                                  action                         = check_integrity                                 name                           = check_integrity                                 parser_handle                  = 0x7fff60d330c8                                  unique_id                      = 44                                            []                                                                                                                                                                [no_action]                                                                        action                         = no_action                                       name                           = no_action                                       parser_handle                  = 0x7fff60d330c8                                  unique_id                      = 21                                            []                                                                                                                                                                [setup_dampers]                                                                    action                         = setup_dampers                                   name                           = setup_dampers                                   parser_handle                  = 0x7fff60d330c8                                  unique_id                      = 41                                            []                                                                                                                                                                [setup_quadrature]                                                                 action                         = setup_quadrature                                name                           = setup_quadrature                                order                          = AUTO                                            parser_handle                  = 0x7fff60d330c8                                  type                           = GAUSS                                           unique_id                      = 31                                            []                                                                                                                                                                [setup_subproblem]                                                                 action                         = setup_subproblem                                coord_type                     = XYZ                                             name                           = setup_subproblem                                parser_handle                  = 0x7fff60d330c8                                  unique_id                      = 6                                             []                                                                                                                                                                [no_action]                                                                        action                         = no_action                                       name                           = no_action                                       unique_id                      = 21                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ?�              ?�    -2?�    M�        ?�     U?�     U?�            ?�     U?�   ?�   B        ?�    Q�?�   �'?������?������?�     U?������?�   &~        ?�������?�    ��?�����{`?�     U?�����=?�   �?�   ,        ?�   ��?�   f^?�   �??�   ��?�    A�?�   ��?�   v�        ?�   �)?�   �6?�   �f?�    4I?�   ۰?�������?�������?�����Si?�����
J?������?�����+-?�     U?�����&>?�����t�?�����
�?������H?������ ?�     U?�    �?�   >�        ?�   �?�   �@?�   f�?�    >8?�   (T?�   6        ?�   =?�    ݾ?�    ��?������^?�    ?������b?�����}�?�����_�?������1?�     U?�����ԇ?�����78?�������?������r?������?�     U?�������