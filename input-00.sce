// Input function - 00

function [Dat] = InputDataSet()

    Dat = struct();

    //---- Fan input data -------------------------------------

    Dat.DPpol = [1055.2  73.946  -13.649]; //  DP = p0 + p1*Q + p2*Q^2

    Dat.Qmin_ref=       2.6; // [m^3/s]     Fan minimum reference flow rate
    Dat.Qmax_ref=      11.0; // [m^3/s]     Fan maximum reference flow rate
    Dat.rho_ref =       1.2; // [kg/m^3]    Fan reference air density
    Dat.rpm_ref =     900.0; // [rpm]       Fan reference speed
    Dat.rpm_max =    1300.0; // [rpm]       Fan maximum rpm
    Dat.rpm_min =     100.0; // [rpm]       Fan minimum rpm
    Dat.rpm_set =     700.0; // [rpm]       Fan setpoint
    Dat.n_fan   =       1;   // [ ]         Number of fans

    //---- Heat store -----------------------------------------

    Dat.L_b     =    10.0;   // [m]         Bed length 
    Dat.W_b     =    10.0;   // [m]         Bed width                 
    Dat.H_b     =     1.5;   // [m]         Bed height                
    Dat.eps_b   =     0.45;  // [ ]         Void ratio 
    Dat.n_b     =     3;     // [ ]         Number of lumped mass control volumes

    Dat.T_b0    =   326.15;  // [K]         Initial bed temperature at t=0  
    Dat.d_p     =     0.100; // [m]         Diameter of stones 
    Dat.psi     =     1.00;  // []          Spherisity
    Dat.rho_p   =  2640.0;   // [kg/m^3]    Density of stones
    Dat.c_pp    =   820.0;   // [J/(kg.K)]  Heat capacity of stones

    Dat.eta_r  =     1.2;    // [ ]         Rock store pressure loss factor

    //---- Electric Heater ------------------------------------

    Dat.q_hmax  = 50000.0;   // [W]         Electric heater max power        
    Dat.q_h0    =     0.0;   // [W]         Electric heater initial setting  

    Dat.KI      =     5.0;   // [W/K.s]     PI controller integral constant 
    Dat.KP      =  2000.0;   // [W/K]       PI controller proportional constant 
    Dat.T_set   =   328.15;  // [K]         PI controller temperature set point for incubator inlet 

    //---- Incubator -------------------------------------------

    Dat.D_y     =     0.075; // [m]         Yogurt canister diameter
    Dat.H_y     =     1.500; // [m]         Total yogurt cannister height
    Dat.c_py    =  3520.0;   // [J/(kg.K)]  Heat capacity of yogurt
    Dat.rho_y   =  1050.0;   // [kg/m^3]    Density of yogurt

    Dat.T_y0    =   293.15;  // [K]         Initial yogurt temperature at t=0 

    Dat.n_s     =     3      // [ ]         Number of stacks (control volumes)
    Dat.N_Ts    =    30;     // [ ]         Yogurt cannister rows - transeverse
    Dat.N_Ls    =    10;     // [ ]         Yogurt cannister rows - longitudinal
    Dat.S_Ts    =     0.095  // [m]         Pitch
    Dat.S_Ls    =     0.095  // [m]         Pitch

    Dat.F_s     =     0.98   // [ ]         Nusselt correction factor (linked to N_Li)
    Dat.f_s     =     0.40   // [ ]         Friction factor (linked to Re and S_L/D)
    Dat.X_s     =     1.00   // [ ]	        Friction correction (linked cannister layout)

    Dat.eta_s   =     1.2;    // [ ]        Yogurt stack store pressure loss factor

    //---- Environmental losses --------------------------------

    Dat.T_env   =   294.15;  // [K]          Environment external temperature
    Dat.A_env   =    100.0;  // [m^2]        Environment loss area  
    Dat.h_env   =     0.50;  // [W/(m^2.K)]  Environment heat transf coef (loss) 

    //---- System input data ----------------------------------

    Dat.T_0     =   293.15;  // [K]          Simulation initial temperature
    Dat.Dt      =    15.0;   // [s]          Simulation integration time step
    Dat.t_0     =     0.0;   // [s]          Initial time
    Dat.t_f     = 28800.0;   // [s]          Simulation final time 

endfunction
