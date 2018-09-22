//===========================================================================
//   TERMODYNAMIC FUNCTIONS (Kroger)
//        All Temperatures are in Kelvin
//===========================================================================

//--- Standard constants ----------------------------------------------------
R_air   = 287.08; // [J/(kg.K)]
P_atm   = 101325; // [Pa]
T_atm   = 293.15; // [K]
rho_atm = 1.204;  // [kg/m^3]


//--- Degrees Celcius to Kelvin ---------------------------------------------
//     TC      Temperature in degrees Celcius
//     returns Temperature in Kelvin
function TK = C2K(TC)
    TK = TC + 273.15;
endfunction


//--- Kelvin to degrees Celcius ---------------------------------------------
//     TC      Temperature in Kelvin
//     returns temperature in degrees Celcius
function TC = K2C(TK)
    TC = TK - 273.15;
endfunction


//--- Temperature limits for functions [K] ----------------------------------
function T =Tlim(T)
    T = max(T, 220.0)
    T = min(T, 380.0)
endfunction


//--- Air specific heat [J/(kg.K)] ------------------------------------------
function cp = cp_air(T)
    T = Tlim(T);
    cp =   ((-0.0000002705209*T + 0.0007083814) .* T - 0.3161783) .* T + 1045.356; 
endfunction


//--- Air thermal conductivity [W/(m.K)] ------------------------------------
function k = k_air(T)
    T = Tlim(T);
    k = ((1.250603E-11*T -0.00000004627937 ) .* T + 0.0001018087) .* T - 0.0004937787;
endfunction


//--- Air density [kg/m^3] --------------------------------------------------
function rho = rho_air(T)
    T = Tlim(T);
    rho = P_atm ./ (R_air * T);
endfunction


//--- Air dynamic viscosity [kg/(m.s)] --------------------------------------
function mu = mu_air(T)
    T = Tlim(T);
    mu =   ((8.15038E-15 * T  - 3.131956E-11) .* T + 0.00000006259793) .* T +  0.000002287973; 
endfunction


//--- Air kinematic viscosity [m^2/s] ---------------------------------------
function nu = nu_air(T)
    nu = mu_air(T) ./ rho_air(T);
endfunction


//--- Air Prandtl number [ ] ------------------------------------------------
function Pr = Pr_air(T)
   Pr =  mu_air(T) .* cp_air(T) ./ k_air(T);
endfunction
