### By Upadesh Subedi 19th July 2024  
## I150 v5


[GlobalParams]
    gravity = '0 -9.8e-5 0'     # gravity is scaled to length_scale = 1e9, time_scale = 1e7
    integrate_p_by_parts = true
[]

Irradiance = 150.0  # kW/cm^2
scan_speed = 5.0    # nm/ms


[Mesh] 
  type = GeneratedMesh  # MOOSE Mesh generator
  dim = 2               # 2 Dimensional Simulation
  nx = 185              # Mesh Elements in x-direction
  ny = 80               # Mesh Elements in y-direction
  xmin = 0              # Starting of Domain in x-direction  (Actual Scale will be defined in MATERIALS Block)
  xmax = 1300           # End of Domain in x-driection       (Actual Scale will be defined in MATERIALS Block)
  ymin = 0              # Starting of Domain in y-direction  (Actual Scale will be defined in MATERIALS Block)
  ymax = 500            # End of Domain in y-driection       (Actual Scale will be defined in MATERIALS Block)
  elem_type = QUAD9     # FEM element type is Quadrilateral with 9 nodal points
[]

################ VARIABLE BLOCK ################
#########################################################

[Variables]

    [vel_x] # x-component of meltpool velocity variable in N-S equation
        order = SECOND
        family = LAGRANGE
    []

    [vel_y] # y-component of meltpool velocity variable in N-S equation
        order = SECOND
        family = LAGRANGE
    []

    [p] # pressure velocity variable in N-S equation
        order = FIRST
        family = LAGRANGE
    []

    [temp] # temperature variable
        initial_condition = 300
        scaling = 1.0e-10
    []

    [w] # potential variable used in SplitCHCRes and kkssplitchcres (global)
        order = FIRST
        family = LAGRANGE
        scaling = 1e-10
    []

    [c] # concentration (global)
        order = FIRST
        family = LAGRANGE
        scaling = 1e-7
    []

    [c1] # local concentration for Phase 1 (Ti-rich HCP phase)
        order = FIRST
        family = LAGRANGE
        scaling =1e-8
    []

    [c2] # local concentration for Phase 2 (Ti3Au IMC Phase)
        order = FIRST
        family = LAGRANGE
        scaling =1e-10
    []

    [c3] # local concentration for Phase 3 (Au-rich Meltpool/LIQUID phase)
        order = FIRST
        family = LAGRANGE
        scaling = 1.0e-10
    []

    [c4] # local concentration for Phase 4 (Au-rich FCC phase)
        order = FIRST
        family = LAGRANGE
        scaling = 1.0e-7
    []

    [eta1] # Non-conserved Order parameter for Ti-rich HCP phase
        order = FIRST
        family = LAGRANGE
        scaling = 1e-2
    []

    [eta2] # Non-conserved Order parameter for Ti3Au IMC phase Grain 1
        order = FIRST
        family = LAGRANGE
        scaling = 1e-2
    []

    [eta3] # Non-conserved Order parameter for Ti3Au IMC phase Grain 2
        order = FIRST
        family = LAGRANGE
        scaling = 1e-2
    []

    [eta4] # Non-conserved Order parameter for Ti3Au IMC phase Grain 3
        order = FIRST
        family = LAGRANGE
        scaling = 1e-2
    []

    [eta5] # Non-conserved Order parameter for Ti3Au IMC phase Grain 4
        order = FIRST
        family = LAGRANGE
        scaling = 1e-2
    []

    [eta6] # Non-conserved Order parameter for Au-rich LIQUID phase
        order = FIRST
        family = LAGRANGE
        scaling = 1e-4
    []

    [eta7] # Non-conserved Order parameter for Au-rich FCC phase
        order = FIRST
        family = LAGRANGE
        scaling = 1e-4
    []

[]


################ INITIAL CONDITION BLOCK ################
#########################################################
[ICs]

    [eta1] # Initial Geometry condition setup for Order Parameter 1 (HCP Phase)
        variable = eta1
        type = FunctionIC
        function = 'if(y>=0&y<=100, 1,0)'
    [] 
    
    [eta2] # Initial Geometry condition setup for Order Parameter 2 (IMC Phase; Grain 1)
        variable = eta2
        type = FunctionIC
        function = 'if(y>100&y<=325 & x>50&x<=250, 1,0)'
    []

    [eta3] # Initial Geometry condition setup for Order Parameter 3 (IMC Phase; Grain 2)
        variable = eta3
        type = FunctionIC
        function = 'if(y>100&y<=275 & x>300&x<=500, 1,0)'
    []

    [eta4] # Initial Geometry condition setup for Order Parameter 4 (IMC Phase; Grain 3)
        variable = eta4
        type = FunctionIC
        function = 'if(y>100&y<=325 & x>550&x<=750, 1,0)'
    []

    [eta5] # Initial Geometry condition setup for Order Parameter 5 (IMC Phase; Grain 4)
        variable = eta5
        type = FunctionIC
        function = 'if(y>100&y<=300 & x>800&x<=1000, 1,0)'
    []
    
    [eta6] # Initial Geometry condition setup for Order Parameter 6 (LIQUID Phase)
        variable = eta6
        type = FunctionIC
        function = 'if(y>400&y<=500 & x>=50&x<=150, 1,0)'
    []

    [eta7]
        variable = eta7  # Initial Geometry condition setup for Order Parameter 7 (FCC Phase; it is set up as all remaining area that is left in computation domain except all the etas from 1 to 6)
        type = FunctionIC
        function = 'if(y>=0&y<=100, 0,
                    if(y>100&y<=325 & x>50&x<=250,  0,
                    if(y>100&y<=275 & x>300&x<=500, 0,
                    if(y>100&y<=325 & x>550&x<=750, 0,
                    if(y>100&y<=300 & x>800&x<=1000, 0,
                    if(y>400&y<=500 & x>=50&x<=150,  0, 1))))))'
    []

    [c] # Initial Global Composition in each of the four Phases
        variable = c
        type = FunctionIC
        function =  '0.07* if(y>=0&y<=100, 1, 0)
                    + 0.2* if(y>100&y<=325 & x>50&x<=250,    1, 
                            if(y>100&y<=275 & x>300&x<=500,  1, 
                            if(y>100&y<=325 & x>550&x<=750,  1, 
                            if(y>100&y<=300 & x>800&x<=1000, 1, 0))))
                    + 0.90* if(y>400&y<=500 & x>=50&x<=150,  1, 0)
                    + 0.70* if(y>=0&y<=100, 0,
                            if(y>100&y<=325 & x>50&x<=250,   0,
                            if(y>100&y<=275 & x>300&x<=500,  0,
                            if(y>100&y<=325 & x>550&x<=750,  0,
                            if(y>100&y<=300 & x>800&x<=1000, 0,
                            if(y>400&y<=500 & x>=50&x<=150,  0, 1))))))' 
    []

    [velocity_x] # Initial x-component flow velocity inside LIQUID Phase
        variable = vel_x
        type = FunctionIC
        function =  '1.0e7*if(y>400&y<=500 & x>=50&x<=150,  1, 0)'
    []

    [velocity_y] # Initial y-component flow velocity inside LIQUID Phase
        variable = vel_y
        type = FunctionIC
        function =  '1.0e7*if(y>400&y<=500 & x>=50&x<=150,  1, 0)'
    []


[]  

################ BOUNDARY CONDITION BLOCK ################
#########################################################
    
[BCs]

    [temp_fixed] # Boundary condtion for temperature which is fixed for bottom side to 300K to act as a sink (to reduce computational cost in simulating computational domain with large depth)
        type = ADDirichletBC
        variable = temp
        boundary = 'bottom'
        value = 300
    []

    [neumann1] # Boundary condition for top surface to avoide the periodic boundary condtion in y-direction
        type = NeumannBC
        boundary = 'top'
        variable = 'eta1'
        value = 0
    []

    [neumann2] # Boundary condition for top surface to avoide the periodic boundary condtion in y-direction
        type = NeumannBC
        boundary = 'bottom'
        variable = 'eta6'
        value = 0
    []

    [neumann3] # Boundary condition for top surface to avoide the periodic boundary condtion in y-direction
        type = NeumannBC
        boundary = 'bottom'
        variable = 'eta7'
        value = 0
    []

    [vel_x_bottom] # Boundary Condition for x-component of Fluid Flow inside Meltpool in 3 sides
        type = ADDirichletBC
        variable = vel_x
        boundary = 'left right bottom'
        value = 0
    []

    [vel_y_bottom] # Boundary Condition for x-component of Fluid Flow inside Meltpool in 3 sides
        type = ADDirichletBC
        variable = vel_y
        boundary = 'left right bottom'
        value = 0
    []


    # [Periodic] # If periodic Boundary condition is to be included in the simulation uncomment this block
    #     [all_exe_temp]
    #       auto_direction = 'x'
    #     #   variable = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 c c1 c2 c3 c4 w temp' 
    #     []
    # []
[]

################ AUXILIARY VARIABLE BLOCK ###############
#########################################################

[AuxVariables]

    [Energy] # Energy as an aux variable
        order = CONSTANT
        family = MONOMIAL
    []

    [bnds] # To visualize the boundaries of each phases in post processing of simulation
    []

    [gr_c] # To visualize the compositional boundaries of phases in post processing of simulation
        order = CONSTANT
        family = MONOMIAL
    []
[]

#################### FUNCTIONS BLOCK ####################
#########################################################


[Functions]

    [path_x] # Path of laser movement as a function of time and x-component of the domain
        type = ParsedFunction
        expression = 125+${scan_speed}*1e-4*t
    []

    [path_y] # Path of laser movement as a function of y-component of the domain
        type = ParsedFunction  
        expression = 500
    []

[]

#################### MATERIALS BLOCK ####################
#########################################################

[Materials]

    ## Anything associated with 1 is for HCP phase, 2 is for IMC Phase, 3 is for LIQUID Phase and 4 is for FCC phase in the material properties

    [scale] # scaling of the computational domain
        type = GenericConstantMaterial
        prop_names  = 'length_scale  energy_scale       time_scale      v_mol'      # m-->nm;  J-->eV;  s-->10mus   v_mol-->m^3/mol (Au->10.21e-6 m^3/mol)
        prop_values = '1.0e+9       6.24150943e+18      1.0e+7         10.21e-6'
    []

    [mu_values] # Viscosity values
        type = GenericConstantMaterial
        prop_names  = 'mu1  mu2   pseudo_mu3      mu4'
        prop_values = '1e4  1e4   0.8195498933    1e4' 
    []

    [mu_LIQ] # Viscosity of LIQUID Phase that is temperature dependent
        type = ParsedMaterial
        property_name = mu3
        material_property_names = 'pseudo_mu3'
        constant_names = 'factor_mu Q_mu'                   # factor for convergence
        constant_expressions = '1.139336879e-3 2669'
        coupled_variables = 'temp'
        expression = 'factor_mu*pseudo_mu3*exp(Q_mu/temp)'  # Source for the equation https://link.springer.com/article/10.1007/s10765-016-2104-7
    []

    [mu_NS] # Dynamic Viscosity of the Phases with interpolating function with length and time scaling
        type = ParsedMaterial
        property_name = mu_name
        material_property_names = 'length_scale time_scale energy_scale mu1 mu2 mu3 mu4 h1 h2 h3 h4 h5 h6 h7'
        expression = '(h1*mu1 + (h2+h3+h4+h5)*mu2 + h6*mu3 + h7*mu4 + 10*mu3*h6*(h1+h2+h3+h4+h5+h7)) / (length_scale*time_scale)'   # 10*mu3*h6*(h1+...h5+h7) is for mushy zone viscosity, using 10 times more viscous than LIQUID phase with interpolation function
    []    

    [conductivity] # Thermal conductivity of 4 phases
        type = GenericConstantMaterial
        prop_names  = 'k1       k2      pseudo_k3   pseudo_k4'
        prop_values = '16.07    91.585     100         338.91'      # https://doi.org/10.1016/j.jestch.2023.101413
    []

    [conductivity_LIQUID] # Function for themperature dependent thermal conductivity of LIQUID Phase
        type = DerivativeParsedMaterial
        property_name = k3
        material_property_names = 'pseudo_k3'
        constant_names =        'f_k3    alpha_k3'
        constant_expressions =  '1       0.027397'                  # Value of alpha_k3 = 0.027397 W/mK^2 reference: https://github.com/anilkunwar/temperature_dependent_material_properties/blob/main/getAuThermalConductivity.F90
        coupled_variables = 'temp'
        expression = 'f_k3*(pseudo_k3 + alpha_k3*temp)'
    []

    [conductivity_FCC]  # Function for themperature dependent thermal conductivity of FCC Phase
        type = DerivativeParsedMaterial
        property_name = k4
        material_property_names = 'pseudo_k4'
        constant_names =        'f_k4      alpha_k4'
        constant_expressions =  '1   -6.93e-2'                      # Value of alpha_k4 = -6.93e-2 W/mK reference: https://github.com/anilkunwar/temperature_dependent_material_properties/blob/main/getAuThermalConductivity.F90
        coupled_variables = 'temp'
        expression = 'f_k4*(pseudo_k4 + alpha_k4*temp)'
    []

    [thermal_conductivity_phases] # Thermal conductivity of all four phases with interpolating functions
        type = ParsedMaterial
        property_name = thermal_conductivity
        material_property_names = 'energy_scale length_scale time_scale k1 k2 k3 k4 h1 h2 h3 h4 h5 h6 h7'
        expression = '(h1*k1 + (h2+h3+h4+h5)*k2 + h6*k3 + h7*k4) * energy_scale/(length_scale*time_scale)'
    []

    [density_values] # Density of phases
        type = GenericConstantMaterial
        prop_names =    'rho1     rho2    pseudo_rho3     pseudo_rho4'
        prop_values =   '4504.7   8202.9  19325.28        19657.6'             # https://doi.org/10.1016/j.jestch.2023.101413
    []
    
    [density_LIQUID] # Function for themperature dependent density of LIQUID Phase
        type = DerivativeParsedMaterial
        property_name = rho3
        material_property_names = 'pseudo_rho3'
        constant_names =        'f_density3   alpha_rho3'
        constant_expressions =  '1              -1.44'                       # Value of alpha_rho3 = -1.44 reference: https://github.com/anilkunwar/temperature_dependent_material_properties/blob/main/getAuDensity.F90
        coupled_variables = 'temp'
        expression = 'f_density3*(pseudo_rho3 + alpha_rho3*temp)'
    []

    [density_FCC] # Function for themperature dependent density of FCC Phase
        type = DerivativeParsedMaterial
        property_name = rho4
        material_property_names = 'pseudo_rho4'
        constant_names = 'fos_density4 alpha_rho4'
        constant_expressions = '1 -1.2'                                     # Value of alpha_rho4 = -1.2 reference: https://github.com/anilkunwar/temperature_dependent_material_properties/blob/main/getAuDensity.F90
        coupled_variables = 'temp'
        expression = 'fos_density4*(pseudo_rho4 + alpha_rho4*temp)'
    []

    [density_phases] # Material property density of all four phases with interpolating functions
        type = ParsedMaterial
        property_name = density_name
        material_property_names = 'length_scale rho1 rho2 rho3 rho4 h1 h2 h3 h4 h5 h6 h7'
        expression = '(h1*rho1 + (h2+h3+h4+h5)*rho2 + h6*rho3 + h7*rho4) / length_scale^3'
    []

    [spec_heat_values] # Specific Heat of phases
        type = GenericConstantMaterial
        prop_names =  'sp1        sp2       pseudo_sp3  pseudo_sp4'
        prop_values = '528.049    428.77     158         132'                           # https://webbook.nist.gov/cgi/inchi?ID=C7440326&Mask=2&Type=JANAFS&Table=on (convert J/molK to J/kgK) &  https://doi.org/10.1016/S0364-5916(01)00026-8
    []

    [Sp_LIQUID] # Function for themperature dependent Specific heat of LIQUID Phase
        type = DerivativeParsedMaterial
        property_name = sp3
        material_property_names = 'pseudo_sp3'
        constant_names =       'f_sp3   alpha_sp3   beta_sp3'
        constant_expressions = '1       5.08e7      -0.0114'                           # Fitted the data from: .TDB file from https://doi.org/10.1016/S0364-5916(01)00026-8
        coupled_variables = 'temp'
        expression = 'f_sp3*(alpha_sp3*exp(beta_sp3*temp)+pseudo_sp3)'
    []

    [Sp_FCC] # Function for themperature dependent Specific heat of FCC Phase
        type = DerivativeParsedMaterial
        property_name = sp4
        material_property_names = 'pseudo_sp4'
        constant_names =        'f_sp4  alpha_sp4   beta_sp4'
        constant_expressions =  '1      2.5e-5      -0.011'                            # Fitted the data from: .TDB file from https://doi.org/10.1016/S0364-5916(01)00026-8
        coupled_variables = 'temp'
        expression = 'f_sp4*(alpha_sp4*temp^2 + beta_sp4*temp + pseudo_sp4)'
    []

    [sp_heat]  # Material property Specific Heat of all four phases with interpolating functions
        type = ParsedMaterial
        property_name = specific_heat
        material_property_names = 'energy_scale sp1 sp2 sp3 sp4 h1 h2 h3 h4 h5 h6 h7'
        expression = '(sp1*h1 + (h2+h3+h4+h5)*sp2 + h6*sp3 + h7*sp4) * energy_scale'
    []

    [absorptivity_value] # Thermal Absorptivity of Au 
        type = ParsedMaterial
        property_name = absorptivity
        material_property_names = 'length_scale'
        expression = '8.5e7/length_scale'                                   # https://doi.org/10.1016/j.jestch.2023.101413
    []

    [beam_radius] # Laser Beam Radius
        type = ParsedMaterial
        property_name = rG
        material_property_names = 'length_scale'
        expression = '120.0e-9*length_scale'
    []

    # [irradiance] # Irradiance (Power Density) of Laser in kW/cm^2
    #     type = ParsedMaterial
    #     property_name = Irradiance
    #     expression = '134.6' 
    # []

    [power_fun] # Calculation of Laser Power from Laser Irradiance with Unit converstion to J/s
        type = ParsedMaterial
        property_name = pow
        constant_names = 'f_irrad pi'                       # f_irrad for unit conversion of kW/cm^2
        constant_expressions = '5.747261834e-4 3.1415'
        material_property_names = 'energy_scale time_scale rG'
        expression = '${Irradiance}*(pi*rG*rG)*f_irrad*energy_scale/time_scale'  
    []

    [volumetric_heat] # Volumetric Heat 
        type = GaussianHS
        power = pow
        efficiency = 0.9    # Efficiency of laser heat absorption by the material taken as 90%
        Ca = 11.72646029    # Coefficient Constant Outside Exponential
        Cb = 3              # Coefficient Constant Inside Exponential
        rG = rG
        factor = 9.803921569e-7     # Factor for convergence
        alpha = absorptivity
        function_x= path_x
        function_y= path_y
    []  

    [Free_Energy_Phase_1] # For Ti rich HCP phase
        type = DerivativeParsedMaterial
        property_name = F1
        material_property_names = 'length_scale energy_scale v_mol'
        constant_names =       'F_HCP   A_HCP       B_HCP       C_HCP   D_HCP   c_eq_HCP    T_eq_HCP    m_HCP   n_HCP'
        constant_expressions = '1000    655.124     9225e-6     0       -6.95   0.1015      1100        2       2'
        coupled_variables = 'c1 temp'
        expression = 'F_HCP*(A_HCP*(c1-c_eq_HCP)^m_HCP + B_HCP*(temp-T_eq_HCP)^n_HCP + C_HCP*c1*temp + D_HCP)*energy_scale/(v_mol*length_scale^3)'
    []

    [Free_Energy_Phase_2] # Ti3Au rich (IMC) phase
        type = DerivativeParsedMaterial
        property_name = F2
        material_property_names = 'length_scale energy_scale v_mol'
        constant_names =        'F_IMC  A_IMC       B_IMC   C_IMC   D_IMC   c_eq_IMC    T_eq_IMC    m_IMC   n_IMC   '
        constant_expressions =  '1000   2328.75     2.2e-5  0       -12.6   0.261       859.11      2       2'
        coupled_variables = 'c2 temp'
        expression = 'F_IMC*(A_IMC*(c2-c_eq_IMC)^m_IMC + B_IMC*(temp-T_eq_IMC)^n_IMC + C_IMC*c2*temp + D_IMC)*energy_scale/(v_mol*length_scale^3)'
    []

    [Free_Energy_Phase_3]  # Au rich LIQUID phase
        type = DerivativeParsedMaterial
        property_name = F3
        material_property_names = 'length_scale energy_scale v_mol'
        constant_names =        'F_LIQ  A_LIQ       B_LIQ       C_LIQ   D_LIQ   c_eq_LIQ    T_eq_LIQ    m_LIQ   n_LIQ'
        constant_expressions =  '11.599 30013.79    0.00122     -0.75   -21     0.5272      1939.34     4       2'
        coupled_variables = 'c3 temp'
        expression = 'F_LIQ*(A_LIQ*(c3-c_eq_LIQ)^m_LIQ + B_LIQ*(temp-T_eq_LIQ)^n_LIQ + C_LIQ*c3*temp + D_LIQ)*energy_scale/(v_mol*length_scale^3)'
    []

    [Free_Energy_Phase_4] # Au rich FCC phase
        type = DerivativeParsedMaterial
        property_name = F4
        material_property_names = 'length_scale energy_scale v_mol'
        constant_names =        'F_FCC  A_FCC       B_FCC       C_FCC   D_FCC   c_eq_FCC    T_eq_FCC    m_FCC   n_FCC'
        constant_expressions =  '19.858 208000.75   0.001463    -0.75   -5.8    0.627       593.26      6       2'
        coupled_variables = 'c4 temp'
        expression = 'F_FCC*(A_FCC*(c4-c_eq_FCC)^m_FCC + B_FCC*(temp-T_eq_FCC)^n_FCC + C_FCC*c4*temp + D_FCC)*energy_scale/(v_mol*length_scale^3)'
    []

     #Switching Functions        ## Eq 10,11 of https://doi.org/10.1016/j.actamat.2010.10.038
    [h1]
        type = SwitchingFunctionMultiPhaseMaterial
        h_name = h1
        all_etas = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7'
        phase_etas = eta1
    []

    [h2]
        type = SwitchingFunctionMultiPhaseMaterial
        h_name = h2
        all_etas = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7'
        phase_etas = eta2
    []

    [h3]
        type = SwitchingFunctionMultiPhaseMaterial
        h_name = h3
        all_etas = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7'
        phase_etas = eta3
    []

    [h4]
        type = SwitchingFunctionMultiPhaseMaterial
        h_name = h4
        all_etas = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7'
        phase_etas = eta4
    []

    [h5]
        type = SwitchingFunctionMultiPhaseMaterial
        h_name = h5
        all_etas = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7'
        phase_etas = eta5
    []

    [h6]
        type = SwitchingFunctionMultiPhaseMaterial
        h_name = h6
        all_etas = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7'
        phase_etas = eta6
    []

    [h7]
        type = SwitchingFunctionMultiPhaseMaterial
        h_name = h7
        all_etas = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7'
        phase_etas = eta7
    []

     # Barrier functions for each phase
    [g1]
        type = BarrierFunctionMaterial
        g_order = SIMPLE
        eta = eta1
        function_name = g1
    []

    [g2]
        type = BarrierFunctionMaterial
        g_order = SIMPLE
        eta = eta2
        function_name = g2
    []

    [g3]
        type = BarrierFunctionMaterial
        g_order = SIMPLE
        eta = eta3
        function_name = g3
    []

    [g4]
        type = BarrierFunctionMaterial
        g_order = SIMPLE
        eta = eta4
        function_name = g4
    []

    [g5]
        type = BarrierFunctionMaterial
        g_order = SIMPLE
        eta = eta5
        function_name = g5
    []

    [g6]
        type = BarrierFunctionMaterial
        g_order = SIMPLE
        eta = eta6
        function_name = g6
    []

    [g7]
        type = BarrierFunctionMaterial
        g_order = SIMPLE
        eta = eta7
        function_name = g7
    []

    [constants] # Material Properties Constants
        type = GenericConstantMaterial
        # prop_names = 'sigma delta gamma M1 M2 M3 M4 R'
        prop_names  = 'sigma    delta       gamma   F_M1        F_M2        F_M3        F_M4        R       TAU' 
        prop_values = '0.4      10.0e-9     1.5     1.0e-21     1.0e-21     1.0e-15     1.0e-21     8.3145  0.1693' #sigma -> J/m^2; delta -> meter; M_si unit-> m^5/Js
    []
    
    [mu] # Model Parameter m
        type = ParsedMaterial
        property_name = mu
        material_property_names = 'sigma delta energy_scale length_scale'
        expression = '6*(sigma/delta)*(energy_scale/length_scale^3)' 
    []
    
    [kappa] # Model parameter k
        type = ParsedMaterial
        property_name = kappa
        material_property_names = 'sigma delta energy_scale length_scale'
        expression = '0.75*(sigma*delta)*(energy_scale/length_scale)' 
    []


     ## Note: All Mobilities should be in SI unit to avoid double scaling
    [Mobility_HCP]    # Temperature dependent Mobility of Phase HCP              # http://arfc.github.io/software/moltres/wiki/input_example/  # exponential expression for temp dep.
        type = DerivativeParsedMaterial
        property_name = M1
        material_property_names = 'F_M1 R'
        constant_names = 'fa_1 F_temp M_1 Q_1'
        constant_expressions = '0 1.0 3.72 1550'
        coupled_variables = 'temp'
        expression = 'fa_1+F_M1*F_temp*M_1*exp(-Q_1/(R*temp))'                    # Eqn 9 of https://www.sciencedirect.com/science/article/pii/S0026271417305449#s0010
    []

    [Mobility_IMC]    # Temperature dependent Mobility of Phase IMC
        type = DerivativeParsedMaterial
        property_name = M2
        material_property_names = 'F_M2 R'
        constant_names = 'fa_2 F_temp M_2 Q_2'
        constant_expressions = '0 1.0 6.28 1123'
        coupled_variables = 'temp'
        expression = 'fa_2+F_M2*F_temp*M_2*exp(-Q_2/(R*temp))'
    []

    [Mobility_LIQUID]    # Temperature dependent Mobility of Phase LIQUID
        type = DerivativeParsedMaterial
        property_name = M3
        material_property_names = 'F_M3 R'
        constant_names = 'fa_3 F_temp M_3 Q_3'
        constant_expressions = '0 1.0 2.908 3043'
        coupled_variables = 'temp'
        expression = 'fa_3+F_M3*F_temp*M_3*exp(-Q_3/(R*temp))'
    []

    [Mobility_FCC]    # Temperature dependent Mobility of Phase FCC
        type = DerivativeParsedMaterial
        property_name = M4
        material_property_names = 'F_M4 R'
        constant_names = 'fa_4 F_temp M_4 Q_4'
        constant_expressions = '0 1.0 4.88 2221'
        coupled_variables = 'temp'
        expression = 'fa_4+F_M4*F_temp*M_4*exp(-Q_4/(R*temp))'
    []

    [Mobility] # Mobilities of all Phases with interpolation function
        type = ParsedMaterial
        property_name = M
        material_property_names = 'length_scale energy_scale time_scale M1 M2 M3 M4 M_gb h1 h2 h3 h4 h5 h6 h7'
        expression = '(h1*M1 + (h2+h3+h4+h5)*M2 + h6*M3 + h7*M4 +(h2+h3+h4+h5)*(h2+h3+h4+h5)*M_gb) * (length_scale^5/(time_scale*energy_scale))'
    []

    [M_grain_bound] # Grain Boundary Mobilities for IMC Grain taken as 100 times that of IMC phase mobility
        type = ParsedMaterial
        material_property_names = 'M2' # M2-->IMC_grain
        property_name = M_gb
        expression = '100*M2'
    []

    [L1-2]
        type = ParsedMaterial
        property_name = L1_2 
        constant_names = 'factor_L'
        constant_expressions = '1.0'
        material_property_names = 'M1 M2 energy_scale time_scale length_scale mu kappa TAU'
        expression = 'factor_L*(0.5/TAU)*(M1+M2)*(4/3)*(mu/kappa)*(length_scale^3/(time_scale*energy_scale))'  ## 0.1693 ==> total values of dinominator of Eq 21 ( here we don't know ci,eq,i/j/k values) so 0.2 is taken as whole 
    []

    [L2-3]
        type = ParsedMaterial
        property_name = L2_3
        constant_names = 'factor_L'
        constant_expressions = '1.0'
        material_property_names = 'M2 M_gb energy_scale time_scale length_scale mu kappa h2 h3 h4 h5 TAU'
        expression = 'factor_L*((0.5/TAU)*(M2+M2)+M_gb*(h2*h3+h3*h4+h4*h5))*(4/3)*(mu/kappa)*(length_scale^3/(time_scale*energy_scale))' 
    []

    [L2-7]
        type = ParsedMaterial
        property_name = L2_7
        constant_names = 'factor_L'
        constant_expressions = '1.0'
        material_property_names = 'M2 M4 energy_scale time_scale length_scale mu kappa TAU'
        expression = 'factor_L*(0.5/TAU)*(M2+M4)*(4/3)*(mu/kappa)*(length_scale^3/(time_scale*energy_scale))' 
    []

    [L2-6]
        type = ParsedMaterial
        property_name = L2_6
        constant_names = 'factor_L'
        constant_expressions = '1.0'
        material_property_names = 'M2 M3 energy_scale time_scale length_scale mu kappa TAU'
        expression = 'factor_L*(1/TAU)*(M3)*(4/3)*(mu/kappa)*(length_scale^3/(time_scale*energy_scale))' ## Taking Mobility of Liquid phase as dictator
    []

    [L1-6]
        type = ParsedMaterial
        property_name = L1_6
        constant_names = 'factor_L'
        constant_expressions = '1.0'
        material_property_names = 'M1 M3 energy_scale time_scale length_scale mu kappa TAU'
        expression = 'factor_L*(1/TAU)*(M3)*(4/3)*(mu/kappa)*(length_scale^3/(time_scale*energy_scale))'  ## Taking Mobility of Liquid phase as dictator
    []

    [L1-7]
        type = ParsedMaterial
        property_name = L1_7
        constant_names = 'factor_L'
        constant_expressions = '1.0'
        material_property_names = 'M1 M4 energy_scale time_scale length_scale mu kappa TAU'
        expression = 'factor_L*(0.5/TAU)*(M1+M4)*(4/3)*(mu/kappa)*(length_scale^3/(time_scale*energy_scale))' 
    []

    [L6-7]
        type = ParsedMaterial
        property_name = L6_7
        constant_names = 'factor_L'
        constant_expressions = '1.0'
        material_property_names = 'M3 M4 energy_scale time_scale length_scale mu kappa TAU'
        expression = 'factor_L*(1/TAU)*(M3)*(4/3)*(mu/kappa)*(length_scale^3/(time_scale*energy_scale))'  ## Taking Mobility of Liquid phase as dictator
    []

    [Interface_Mobility]
        type = ParsedMaterial
        property_name = L
        coupled_variables = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7'
        material_property_names = 'L1_2 L2_3 L2_7 L1_6 L2_6 L1_7 L6_7 h1 h2 h3 h4 h5 h6 h7'
        expression = '(L1_2*(h1*h2+h1*h3+h1*h4+h1*h5)
                    +L2_3*(h2*h3+h3*h4+h4*h5)
                    +L2_7*(h2*h7+h3*h7+h4*h7+h5*h7)
                    +L2_6*(h2*h6+h3*h6+h4*h6+h5*h6) 
                    +L1_6*h1*h6
                    +L1_7*h1*h7 
                    +L6_7*h6*h7)'
    []
[]

#################### KERNELS BLOCK ####################
#########################################################

[Kernels]
  # Kernels for velocity and pressure
    [mass]
        type = INSMass
        variable = p
        rho_name = density_name
        u = vel_x
        v = vel_y
        pressure = p
    []

    [x_momentum_space]
        type = INSMomentumLaplaceForm
        variable = vel_x
        rho_name = density_name
        mu_name = mu_name
        u = vel_x
        v = vel_y
        pressure = p
        component = 0
    []

    [y_momentum_space]
        type = INSMomentumLaplaceForm
        variable = vel_y
        rho_name = density_name
        mu_name = mu_name
        u = vel_x
        v = vel_y
        pressure = p
        component = 1
    []

    [x_momentum_time]
        type = INSMomentumTimeDerivative
        variable = vel_x 
        rho_name = density_name
    []
    
    [y_momentum_time]
        type = INSMomentumTimeDerivative
        variable = vel_y
        rho_name = density_name
    []

    [conv_c1]
        type = CHConvection
        variable = c1
        u = vel_x
        v = vel_y
    []

    [conv_c2]
        type = CHConvection
        variable = c2
        u = vel_x
        v = vel_y
    []

    [conv_c3]
        type = CHConvection
        variable = c3
        u = vel_x
        v = vel_y
    []

    [conv_c4]
        type = CHConvection
        variable = c4
        u = vel_x
        v = vel_y
    []

    [heat_convection_x_y]
        type = CHConvection
        variable = temp
        u = vel_x
        v = vel_y
    []    

   # Kernels for Temperature 
    [time]
        type = HeatConductionTimeDerivative
        variable = temp
        density_name = density_name
    []

    [heat_conduct]
        type = HeatConduction
        variable = temp
        diffusion_coefficient = thermal_conductivity
    []

    [heat_source]
        type = ADMatHeatSource
        material_property = volumetric_heat
        variable = temp
    []

     # Phase concentration constraints
    [chempot12]
        type = KKSPhaseChemicalPotential
        variable = c1
        cb = c2
        fa_name = F1
        fb_name = F2
        args_a = temp
    []

    [chempot23]
        type = KKSPhaseChemicalPotential
        variable = c2
        cb       = c3
        fa_name  = F2
        fb_name  = F3
        args_a = temp
    []

    [chempot34]
        type = KKSPhaseChemicalPotential
        variable = c3
        cb       = c4
        fa_name  = F3
        fb_name  = F4
        args_a = temp
    []
    
    [chempot41]
        type = KKSPhaseChemicalPotential
        variable = c4
        cb       = c1
        fa_name  = F4   ## If kept F1 instead then residual of c4 is 0 that is error
        fb_name  = F1
        args_a = temp
    []

    [phaseconcentration]
        type = KKSMultiPhaseConcentration
        variable = c4
        cj = 'c1 c2 c2 c2 c2 c3 c4'
        hj_names = 'h1 h2 h3 h4 h5 h6 h7'
        etas = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7'
        c = c
    []

     ## Kernels for split Cahn-Hilliard type equation
       ## CHBulk known as KKSSplitCHCRes is here to replace SplitCHParsed
       ## because in KKS model , gradient energy term is not allowed in the C-H type equation [Tonks2018-ComputationalMaterialsScience,vol. 147, pp.353-362.]
       ## while SplitCHParsed kernel consists of the term k\nabla^2 c_i (thus making it unsuitable here), KKSSplitCHCRes fortunately avoids this term.
       ## Never use SplitCHParsed kernel with KKS model
       ## Because of the KKS condition 1 (equality of chemical potential between any two adjacent phases), one KKSSplitCHCRes kernel (either for c1, c2 or c3) is sufficient and there is no need to put three such kernels corresponding to c1, c2 and c3.

    [CHBulk] # Gives the residual for the concentration, dF/dc-mu
        type = KKSSplitCHCRes
        args_a = temp
        variable = c
        ca = c2       # Why only c2? coz of KKS condition equality of chem pot between phases. & only F2 is used
        fa_name = F2
        w = w
    []

    [dcdt] # Gives dc/dt
        type = CoupledTimeDerivative
        variable = w
        v = c
    []

    [ckernel] # Gives residual for chemical potential dc/dt+M\grad(mu)
        type = SplitCHWRes
        mob_name = M
        variable = w
        coupled_variables = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 temp'
    []

   # Kernels for Allen-Cahn equation for eta1
    [deta1dt]
        type = TimeDerivative
        variable = eta1
    []

    [ACBulkF1]
        type = KKSMultiACBulkF
        variable = eta1
        Fj_names = 'F1 F2 F2 F2 F2 F3 F4'
        hj_names = 'h1 h2 h3 h4 h5 h6 h7'
        gi_name = g1
        eta_i = eta1
        wi = 1.0
        coupled_variables = 'c1 c2 c3 c4 eta2 eta3 eta4 eta5 eta6 eta7 temp'
        mob_name = L
    []

    [ACBulkC1]
        type = KKSMultiACBulkC
        variable = eta1
        Fj_names = 'F1 F2 F2 F2 F2 F3 F4'
        hj_names = 'h1 h2 h3 h4 h5 h6 h7'
        cj_names = 'c1 c2 c2 c2 c2 c3 c4'
        eta_i = eta1
        coupled_variables = 'eta2 eta3 eta4 eta5 eta6 eta7 temp'
        mob_name = L
    []

    [ACInterface1]
        type = ACInterface
        variable = eta1
        kappa_name = kappa
        mob_name = L
    []

    [ACdfintdeta1]
        type = ACGrGrMulti
        variable = eta1
        v = 'eta2 eta3 eta4 eta5 eta6 eta7'
        gamma_names = 'gamma gamma gamma gamma gamma gamma'
        mob_name = L
        coupled_variables = 'eta2 eta3 eta4 eta5 eta6 eta7 temp'
    []

   # Kernels for Allen-Cahn equation for eta2
    [deta2dt]
        type = TimeDerivative
        variable = eta2
    []

    [ACBulkF2]
        type = KKSMultiACBulkF
        variable = eta2
        Fj_names = 'F1 F2 F2 F2 F2 F3 F4'
        hj_names = 'h1 h2 h3 h4 h5 h6 h7'
        gi_name = g2
        eta_i = eta2
        wi = 1.0
        coupled_variables = 'c1 c2 c3 c4 eta1 eta3 eta4 eta5 eta6 eta7 temp'
        mob_name = L
    []

    [ACBulkC2]
        type = KKSMultiACBulkC
        variable = eta2
        Fj_names = 'F1 F2 F2 F2 F2 F3 F4'
        hj_names = 'h1 h2 h3 h4 h5 h6 h7'
        cj_names = 'c1 c2 c2 c2 c2 c3 c4'
        eta_i = eta2
        coupled_variables = 'eta1 eta3 eta4 eta5 eta6 eta7 temp'
        mob_name = L
    []

    [ACInterface2]
        type = ACInterface
        variable = eta2
        kappa_name = kappa
        mob_name = L
    []

    [ACdfintdeta2]
        type = ACGrGrMulti
        variable = eta2
        v = 'eta1 eta3 eta4 eta5 eta6 eta7'
        gamma_names = 'gamma gamma gamma gamma gamma gamma'
        mob_name = L
        coupled_variables = 'eta1 eta3 eta4 eta5 eta6 eta7 temp'
    []

   # Kernels for Allen-Cahn equation for eta3
    [deta3dt]
        type = TimeDerivative
        variable = eta3
    []

    [ACBulkF3]
        type = KKSMultiACBulkF
        variable = eta3
        Fj_names = 'F1 F2 F2 F2 F2 F3 F4'
        hj_names = 'h1 h2 h3 h4 h5 h6 h7'
        gi_name = g3
        eta_i = eta3
        wi = 1.0
        coupled_variables = 'c1 c2 c3 c4 eta1 eta2 eta4 eta5 eta6 eta7 temp'
        mob_name = L
    []

    [ACBulkC3]
        type = KKSMultiACBulkC
        variable = eta3
        Fj_names = 'F1 F2 F2 F2 F2 F3 F4'
        hj_names = 'h1 h2 h3 h4 h5 h6 h7'
        cj_names = 'c1 c2 c2 c2 c2 c3 c4'
        eta_i = eta3
        coupled_variables = 'eta1 eta2 eta4 eta5 eta6 eta7 temp'
        mob_name = L
    []

    [ACInterface3]
        type = ACInterface
        variable = eta3
        kappa_name = kappa
        mob_name = L
    []

    [ACdfintdeta3]
        type = ACGrGrMulti
        variable = eta3
        v = 'eta1 eta2 eta4 eta5 eta6 eta7'
        gamma_names = 'gamma gamma gamma gamma gamma gamma'
        mob_name = L
        coupled_variables = 'eta1 eta2 eta4 eta5 eta6 eta7 temp'
    []

   # Kernels for Allen-Cahn equation for eta4
    [deta4dt]
        type = TimeDerivative
        variable = eta4
    []

    [ACBulkF4]
        type = KKSMultiACBulkF
        variable = eta4
        Fj_names = 'F1 F2 F2 F2 F2 F3 F4'
        hj_names = 'h1 h2 h3 h4 h5 h6 h7'
        gi_name = g4
        eta_i = eta4
        wi = 1.0
        coupled_variables = 'c1 c2 c3 c4 eta1 eta2 eta3 eta5 eta6 eta7 temp'
        mob_name = L
    []

    [ACBulkC4]
        type = KKSMultiACBulkC
        variable = eta4
        Fj_names = 'F1 F2 F2 F2 F2 F3 F4'
        hj_names = 'h1 h2 h3 h4 h5 h6 h7'
        cj_names = 'c1 c2 c2 c2 c2 c3 c4'
        eta_i = eta4
        coupled_variables = 'eta1 eta2 eta3 eta5 eta6 eta7 temp'
        mob_name = L
    []

    [ACInterface4]
        type = ACInterface
        variable = eta4
        kappa_name = kappa
        mob_name = L
    []

    [ACdfintdeta4]
        type = ACGrGrMulti
        variable = eta4
        v = 'eta1 eta2 eta3 eta5 eta6 eta7'
        gamma_names = 'gamma gamma gamma gamma gamma gamma'
        mob_name = L
        coupled_variables = 'eta1 eta2 eta3 eta5 eta6 eta7 temp'
    []    

   # Kernels for Allen-Cahn equation for eta5
    [deta5dt]
        type = TimeDerivative
        variable = eta5
    []

    [ACBulkF5]
        type = KKSMultiACBulkF
        variable = eta5
        Fj_names = 'F1 F2 F2 F2 F2 F3 F4'
        hj_names = 'h1 h2 h3 h4 h5 h6 h7'
        gi_name = g5
        eta_i = eta5
        wi = 1.0
        coupled_variables = 'c1 c2 c3 c4 eta1 eta2 eta3 eta4 eta6 eta7 temp'
        mob_name = L
    []

    [ACBulkC5]
        type = KKSMultiACBulkC
        variable = eta5
        Fj_names = 'F1 F2 F2 F2 F2 F3 F4'
        hj_names = 'h1 h2 h3 h4 h5 h6 h7'
        cj_names = 'c1 c2 c2 c2 c2 c3 c4'
        eta_i = eta5
        coupled_variables = 'eta1 eta2 eta3 eta4 eta6 eta7 temp'
        mob_name = L
    []

    [ACInterface5]
        type = ACInterface
        variable = eta5
        kappa_name = kappa
        mob_name = L
    []

    [ACdfintdeta5]
        type = ACGrGrMulti
        variable = eta5
        v = 'eta1 eta2 eta3 eta4 eta6 eta7'
        gamma_names = 'gamma gamma gamma gamma gamma gamma'
        mob_name = L
        coupled_variables = 'eta1 eta2 eta3 eta4 eta6 eta7 temp'
    []    

   # Kernels for Allen-Cahn equation for eta6
    [deta6dt]
        type = TimeDerivative
        variable = eta6
    []

    [ACBulkF6]
        type = KKSMultiACBulkF
        variable = eta6
        Fj_names = 'F1 F2 F2 F2 F2 F3 F4'
        hj_names = 'h1 h2 h3 h4 h5 h6 h7'
        gi_name = g6
        eta_i = eta6
        wi = 1.0
        coupled_variables = 'c1 c2 c3 c4 eta1 eta2 eta3 eta4 eta5 eta7 temp'
        mob_name = L
    []

    [ACBulkC6]
        type = KKSMultiACBulkC
        variable = eta6
        Fj_names = 'F1 F2 F2 F2 F2 F3 F4'
        hj_names = 'h1 h2 h3 h4 h5 h6 h7 '
        cj_names = 'c1 c2 c2 c2 c2 c3 c4'
        eta_i = eta6
        coupled_variables = 'eta1 eta2 eta3 eta4 eta5 eta7 temp'
        mob_name = L
    []

    [ACInterface6]
        type = ACInterface
        variable = eta6
        kappa_name = kappa
        mob_name = L
    []

    [ACdfintdeta6]
        type = ACGrGrMulti
        variable = eta6
        v = 'eta1 eta2 eta3 eta4 eta5 eta7'
        gamma_names = 'gamma gamma gamma gamma gamma gamma'
        mob_name = L
        coupled_variables = 'eta1 eta2 eta3 eta4 eta5 eta7 temp'
    []
    
   # Kernels for Allen-Cahn equation for eta7
    [deta7dt]
        type = TimeDerivative
        variable = eta7
    []

    [ACBulkF7]
        type = KKSMultiACBulkF
        variable = eta7
        Fj_names = 'F1 F2 F2 F2 F2 F3 F4'
        hj_names = 'h1 h2 h3 h4 h5 h6 h7'
        gi_name = g7
        eta_i = eta7
        wi = 1.0
        coupled_variables = 'c1 c2 c3 c4 eta1 eta2 eta3 eta4 eta5 eta6 temp'
        mob_name = L
    []

    [ACBulkC7]
        type = KKSMultiACBulkC
        variable = eta7
        Fj_names = 'F1 F2 F2 F2 F2 F3 F4'
        hj_names = 'h1 h2 h3 h4 h5 h6 h7'
        cj_names = 'c1 c2 c2 c2 c2 c3 c4'
        eta_i = eta7
        coupled_variables = 'eta1 eta2 eta3 eta4 eta5 eta6 temp'
        mob_name = L
    []

    [ACInterface7]
        type = ACInterface
        variable = eta7
        kappa_name = kappa
        mob_name = L
    []

    [ACdfintdeta7]
        type = ACGrGrMulti
        variable = eta7
        v = 'eta1 eta2 eta3 eta4 eta5 eta6'
        gamma_names = 'gamma gamma gamma gamma gamma gamma'
        mob_name = L
        coupled_variables = 'eta1 eta2 eta3 eta4 eta5 eta6 temp'
    []
    
[]


[AuxKernels]
    [Energy_total]
        type = KKSMultiFreeEnergy
        Fj_names = 'F1 F2 F2 F2 F2 F3 F4'
        hj_names = 'h1 h2 h3 h4 h5 h6 h7'
        gj_names = 'g1 g2 g3 g4 g5 g6 g7'
        variable = Energy
        w = 1
        interfacial_vars = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7'
        kappa_names = 'kappa kappa kappa kappa kappa kappa kappa'
    []

    [bnds]
        type = BndsCalcAux
        variable = bnds
        var_name_base = eta
        op_num = 7
        v = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7'
    []

    # [sumCdothsquare]
    #     type = SevenPhasesSumCdothsquare
    #     variable = gr_c
    #     var1=c
    #     h1_name = h1
    #     h2_name = h2
    #     h3_name = h3
    #     h4_name = h4
    #     h5_name = h5
    #     h6_name = h6
    #     h7_name = h7
    #   []

    [sumCdothsquare]
        type                = PhaseGlobalComposition
        variable            = gr_c
        global_composition  = c
        total_etas          = 2
        h_names             = 'h1 h2 h3 h4 h5 h6 h7'
    []
[]

[Executioner]
    type = Transient
    solve_type          = 'PJFNK'
    
    # https://mooseframework.inl.gov/modules/phase_field/Solving.html

    # petsc_options_iname = '-pc_type'
    # petsc_options_value = 'lu'

    petsc_options       = '-snes_converged_reason -ksp_converged_reason -options_left'
    petsc_options_iname = '-ksp_gmres_restart -pc_factor_shift_type -pc_factor_shift_amount -pc_type'
    petsc_options_value = '100 NONZERO 1e-15 ilu'
    
    # petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
    # petsc_options_value = 'asm      31      preonly     lu      2'

    # petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
    # petsc_options_value = 'hypre    boomeramg      31       0.7'

    # petsc_options = '-snes_converged_reason -ksp_converged_reason -options_left'
     
    l_max_its = 30
    nl_max_its = 50
    l_tol = 1e-04
    nl_rel_tol = 1e-08
    nl_abs_tol = 1e-09

    end_time = 8.0e+07
    dtmax = 2.5e4
    dt = 2.5e+4

    # [Adaptivity]
    #     initial_adaptivity = 1
    #     refine_fraction = 0.7
    #     coarsen_fraction = 0.1
    #     max_h_level = 1
    #     weight_names = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7'
    #     weight_values = '1 1 1 1 1 2 1'
    # []

    # [TimeStepper]
    #     type = IterationAdaptiveDT
    #     dt = 2.5e+04
    #     cutback_factor = 0.8
    #     growth_factor = 1.5 
    #     optimal_iterations = 7
    #     # num_steps = 55
    #     # end_time = 1.0e+10
    # []
[]

[Preconditioning]

    [SMP]
        type = SMP
        full = true
        petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_mat_solver_type'
        petsc_options_value = 'lu       NONZERO               strumpack'
    []

    active = 'full'

    [full]
        type = SMP
        full = true
    []

    [mydebug]
        type = FDP
        full = true
    []
[]

[Postprocessors]
    # [total_energy]
    #   type = ElementIntegralVariablePostprocessor
    #   variable = f_density
    #   execute_on = 'Initial TIMESTEP_END'
    # []
    
   # Area of Phases
    [area_h1]
        type = ElementIntegralMaterialProperty
        mat_prop = h1
        execute_on = 'Initial TIMESTEP_END'
    []
  
    [area_h2]
        type = ElementIntegralMaterialProperty
        mat_prop = h2
        execute_on = 'Initial TIMESTEP_END'
    []
  
    [area_h3]
        type = ElementIntegralMaterialProperty
        mat_prop = h3
        execute_on = 'Initial TIMESTEP_END'
    []
  
    [area_h4]
        type = ElementIntegralMaterialProperty
        mat_prop = h4
        execute_on = 'Initial TIMESTEP_END'
    []
  
    [area_h5]
        type = ElementIntegralMaterialProperty
        mat_prop = h5
        execute_on = 'Initial TIMESTEP_END'
    []
  
    [area_h6]
        type = ElementIntegralMaterialProperty
        mat_prop = h6
        execute_on = 'Initial TIMESTEP_END'
    []

    [area_h7]
        type = ElementIntegralMaterialProperty
        mat_prop = h7
        execute_on = 'Initial TIMESTEP_END'
    []

    [temp_max]
        type = ElementExtremeValue
        variable = temp
    []

    [temp_avg]
        type = ElementAverageValue
        variable = temp
    []

  []

[Outputs]
    exodus = true
    interval = 1
    file_base = exodus/7_NS
    csv = true
    [my_checkpoint]
        type = Checkpoint
        num_files = 2
        interval = 2
        file_base = exodus/7_NS
    []
[]

[Debug]
    show_var_residual_norms = true
[]

