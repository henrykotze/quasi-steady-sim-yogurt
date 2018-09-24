DPpol = [1055.2  73.946  -13.649]; //  DP = p0 + p1*Q + p2*Q^2

    Qmin_ref=       2.6 // [m^3/s]     Fan minimum reference flow rate
    Qmax_ref=      11.0; // [m^3/s]     Fan maximum reference flow rate
    rho_ref =       1.2; // [kg/m^3]    Fan reference air density
    rpm_ref =     900.0; // [rpm]       Fan reference speed
    rpm_max =    1300.0; // [rpm]       Fan maximum rpm
    rpm_min =     100.0; // [rpm]       Fan minimum rpm
    rpm_set =     700.0; // [rpm]       Fan setpoint
    n_fan   =       1;   // [ ]         Number of fans

    //---- Heat store -----------------------------------------

    L_b     =    10.0;   // [m]         Bed length 
    W_b     =    10.0;   // [m]         Bed width                 
    H_b     =     1.5;   // [m]         Bed height                
    eps_b   =     0.45;  // [ ]         Void ratio 
    n_b     =     3;     // [ ]         Number of lumped mass control volumes

    T_b0    =   326.15;  // [K]         Initial bed temperature at t=0  
    d_p     =     0.100; // [m]         Diameter of stones 
    psi     =     1.00;  // []          Spherisity
    rho_p   =  2640.0;   // [kg/m^3]    Density of stones
    c_pp    =   820.0;   // [J/(kg.K)]  Heat capacity of stones

    eta_r  =     1.2;    // [ ]         Rock store pressure loss factor

    //---- Electric Heater ------------------------------------

    q_hmax  = 50000.0;   // [W]         Electric heater max power        
    q_h0    =     0.0;   // [W]         Electric heater initial setting  

    KI      =     5.0;   // [W/K.s]     PI controller integral constant 
    KP      =  2000.0;   // [W/K]       PI controller proportional constant 
    T_set   =   328.15;  // [K]         PI controller temperature set point for incubator inlet 

    //---- Incubator -------------------------------------------

    D_y     =     0.075; // [m]         Yogurt canister diameter
    H_y     =     1.500; // [m]         Total yogurt cannister height
    c_py    =  3520.0;   // [J/(kg.K)]  Heat capacity of yogurt
    rho_y   =  1050.0;   // [kg/m^3]    Density of yogurt

    T_y0    =   293.15;  // [K]         Initial yogurt temperature at t=0 

    n_s     =     3      // [ ]         Number of stacks (control volumes)
    N_Ts    =    30;     // [ ]         Yogurt cannister rows - transeverse
    N_Ls    =    10;     // [ ]         Yogurt cannister rows - longitudinal
    S_Ts    =     0.095  // [m]         Pitch
    S_Ls    =     0.095  // [m]         Pitch

    F_s     =     0.98   // [ ]         Nusselt correction factor (linked to N_Li)
    f_s     =     0.40   // [ ]         Friction factor (linked to Re and S_L/D)
    X_s     =     1.00   // [ ]	        Friction correction (linked cannister layout)

    eta_s   =     1.2;    // [ ]        Yogurt stack store pressure loss factor

    //---- Environmental losses --------------------------------

    T_env   =   294.15;  // [K]          Environment external temperature
    A_env   =    100.0;  // [m^2]        Environment loss area  
    h_env   =     0.50;  // [W/(m^2.K)]  Environment heat transf coef (loss) 

    //---- System input data ----------------------------------

    T_0     =   293.15;  // [K]          Simulation initial temperature
    Dt      =    15.0;   // [s]          Simulation integration time step
    t_0     =     0.0;   // [s]          Initial time
    t_f     = 28800.0;   // [s]          Simulation final time 
