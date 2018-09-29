//===========================================================================
// MAIN PROGRAM SCRIPTS
//
// Name:   Henry Kotzé
// Std nr: 19231865
//===========================================================================

clear()
clearglobal()

//--- Input and defs --------------------------------------------------------

thispath = get_absolute_file_path("main-script.sce")

exec(thispath+'func-thermo.sce', -1)   // Thermo funcs and consts
exec(thispath+'input-00.sce', -1)      // Input data (Will change)
//exec(thispath+'func-simul.sce', -1)    // Put all your functions here!!!!

//--- Global values ---------------------------------------------------------

global data                            // Input data
data = InputDataSet()                  // Declared in 'input-xx.sce'


// All variables are declared in input-00.sce and explained
X = [326.15 326.15-5 326.15-10,293.15+10 293.15+12 293.15+15,0]
DX = [0 0 0 0 0 0 0 ]
Tr = [data.T_b0-15 data.T_b0-10 data.T_b0-8 data.T_b0-4]
Ts = [data.T_y0 data.T_y0-5 data.T_y0-7 data.T_y0-10]
T = [308.16176   308.16176   316.15001   320.61251   320.2396   308.16176]

//---Determine the Fan Pressure Rise--------------------------------------
function [DP] = DP_Fan(Q,N)
      
    Qr = Q/data.n_fan;
    Q_star = Qr*data.rpm_ref/N;
    DP = data.DPpol(3)*Q_star^2+data.DPpol(2)*Q_star+data.DPpol(1);
    
    
endfunction


//----- Determine the Total Pressure loss in system----------------------------
function [DP]= DP_Sys(Tr,Ts,Q)
        Pd = 0
        Py = 0
        for j = 1:data.n_b,
            T_avg = 0.5*( Tr(j)+Tr(j+1) );                      // average temperature
            rho_avg = rho_air( T_avg) ;                         // density of air
            mu_avg = mu_air( T_avg );                           // dynamic viscosity of air
            G_avg = rho_avg*(Q/(data.W_b*data.H_b));                            // mass flow rate per unit area
            Re_avg = G_avg*data.d_p/(mu_avg);                        // Reynolds Num
            f_avg = 4.466*((Re_avg)^0.2)*((data.psi)^0.696)*((data.eps_b)^-2.945)*exp( 11.85*( log(data.psi) )^2 ); // friction factor
            Pd = -f_avg*G_avg^2*(data.L_b/data.n_b)/(rho_avg*data.d_p) + Pd;   //pressure drop
    end
    Pd = data.eta_r*Pd
    
    for i = 1:data.n_s,
        T_avg = 0.5*(Ts(i)+Ts(i+1));            // Average Temperature
        rho_avg = rho_air(T_avg);               // density of air
        mu_avg = mu_air(T_avg);                 // dynamic viscosity of air
        As = data.N_Ts*data.S_Ts*data.H_y;                     // Effective Area 
        mu_s = Q/As; 
        
        // Pressure drop 
        Py = Py + -data.N_Ls*0.3*(rho_avg*mu_s^2)/2;
    end
    Py = data.eta_s*Py;
    
    DP = Py+Pd;

endfunction

//------
function [X,T,Q,Tr,Ts]=Cycle_Temp(t,X,T,Q,Tr,Ts)
    
    
    I_i  = X(7);
    Tb = X(1:3);
    Ty = X(4:6);
    
    for k=10,
            T(2) = T(1);                // no temperature change accross fan
            Tr(1) = T(2);
            for j = 1:data.n_b,
                T_avg = 0.5*( Tr(j)+Tr(j+1) );                      // average temperature
                rho_avg = rho_air( T_avg) ;                         // density of air
                mu_avg = mu_air( T_avg );                           // dynamic viscosity of air
                cp_avg = cp_air( T_avg );                           // specific heat capacity of air
                k_avg = k_air( T_avg );                             // thermal conductivity of air
                G_avg = rho_avg*(Q/(data.W_b*data.H_b));                            // mass flow rate per unit area
                Re_avg = G_avg*data.d_p/(mu_avg);                        // Reynolds Num
                f_avg = 4.466*((Re_avg)^0.2)*((data.psi)^0.696)*((data.eps_b)^-2.945)*exp( 11.85*( log(data.psi) )^2 ); // friction factor
                Nu_avg = 0.437*((Re_avg)^0.75)*((data.psi)^3.35)*((data.eps_b)^-1.620)*exp( 29.03*( log(data.psi) )^2 ); // Nusselt Number
                h_avg = (k_avg*Nu_avg)/(data.d_p^2);                     // Convection coefficient 
                Tr(j+1) = Tb(j) - ( Tb(j)-Tr(j) )*exp( (-h_avg*( data.L_b/data.n_b ) )/( cp_avg*G_avg ) );

            end
        T(3) = Tr(data.n_b+1);
    
        // Heater
       T_avg = (T(3)+T(4))/2;                          // Average Temperature
        cp_avg = cp_air(T_avg);                         // Heat capacity air
        rho_avg = rho_air(T_avg);                       // density of air
        q_h = data.KP*(data.T_set-Ts(1))+data.KI*I_i                     // heat added
        T(4) = T(3) + q_h/(cp_avg*rho_avg*Q)           // temperature change
    
       // Environemental Losses
        T_avg = (T(4)+T(5))/2;                          // average temperature
        cp_avg = cp_air(T_avg);                         // Heat capacity air
        rho_avg = rho_air(T_avg);                       // density of air
        
        q = rho_avg*Q*cp_avg*(T(4)-T(5));
    
        T(5) = data.T_env + (T(4)-data.T_env)/(exp(data.h_env*data.A_env/(rho_avg*Q*cp_avg)));
        Ts(1) = T(5);
        
    // Yogurt Warming
        for i = 1:data.n_s,
            T_avg = 0.5*(Ts(i)+Ts(i+1));            // Average Temperature
            rho_avg = rho_air(T_avg);               // density of air
            mu_avg = mu_air(T_avg);                 // dynamic viscosity of air
            cp_avg = cp_air(T_avg);                 // Heat capacity air
            k_avg = k_air(T_avg);                   // thermal conductivity of air
            Pr_avg = Pr_air(T_avg);                 // Prantel Number of Air
            Pr_avg_e = Pr_air(Ts(i+1));             // Prantel Number of air for next control volume
            nu_avg = nu_air(T_avg);                 // Nusselt Number of air
            As = data.N_Ts*data.S_Ts*data.H_y;                     // Effective Area 
            mu_s = Q/As;                           
            Re_avg = mu_s*data.D_y/(nu_avg);                         // Reynolds Numbers
            
            if Re_avg > 10e3 || Re_avg < 2*10e5 then            // Determine correct Nusselt number equation based on Re  
                Nu_d =  0.72*data.F_s*Re_avg^0.63*Pr_avg^0.365*(Pr_avg/Pr_avg_e)^0.25;   // Nusselt Number
            elseif Re_avg > 2*10e5 || Re_avg < 2*10e6  then
                Nu_d = 0.033*data.F_s*Re_avg^0.8*Pr_avg^0.4*(Pr_avg/Pr_avg_e)^0.25;      // Nusselt Number
            end,
            m_y = (%pi)/4*rho_avg*(data.D_y)^2*data.H_y*data.N_Ls*data.N_Ts;
            h_avg = Nu_d*k_avg/(data.D_y);
            A_y = %pi*data.D_y*data.N_Ls*data.N_Ts;
            
            Ts(i+1) = Ty(i) - (Ty(i) - Ts(i))*exp( -1*(A_y*h_avg)/(rho_avg*Q*cp_avg) );
    
        end
    T(6) = Ts(data.n_s+1);
    if abs( (T(1)-T(6))/(T(6))) < 1e-8  then
        k =10;
    else
        T(1) = T(6);                // closed system
    end
    
    
    end
    
endfunction


//----- Determine the min and max volume flow rates  that the curve fit is valid for, using the fan laws
function [Q,T,Tr,Ts] = Workpoint(t,X,T,Tr,Ts,N)
    
    Qa = data.Qmin_ref/data.n_fan*(data.rpm_ref/N);
    Qb = data.Qmax_ref/data.n_fan*(data.rpm_ref/N);
    Qm = (Qa+Qb)/2;
    
    for i = 1:50
        pressure_pump = DP_Fan(Qm,N);
        pressure_sys = DP_Sys(Tr,Ts,Qm);
        if pressure_pump > abs(pressure_sys) then
            Qa = Qm;
        else
            Qb = Qm;
        end
        
        Qm = (Qa+Qb)/2;
        [X,T,Qm,Tr,Ts]=Cycle_Temp(t,X,T,Qm,Tr,Ts)
        
        if abs( (Qb-Qa )/(Qm) ) < 1e-5   then
            i = 50;
        elseif i == 50 then
            disp('Warning: Workpoint reached iteration limit')
        end
    end
    
    Q = Qm;
endfunction


function [DX,t,X,T,Tr,Ts,N] = Derivatives(t,X,T,Tr,Ts,N,DX)
    [Q,T,Tr,Ts] = Workpoint(t,X,T,Tr,Ts,N)
    for j = 1:data.n_b,
        T_avg = 0.5*( Tr(j)+Tr(j+1) );                      // average temperature
        rho_avg = rho_air( T_avg) ;                         // density of air
        cp_avg = cp_air( T_avg );                           // specific heat capacity of air
        G_avg = rho_avg*(Q/(data.W_b*data.H_b));                            // mass flow rate per unit area
        DX(j) = (cp_avg*G_avg*(Tr(j)-Tr(j+1)))/(2640*820*(1-Tr(j+1))*(data.L_b/data.n_b));
    end
    
    for i = 1:data.n_s,
        T_avg = 0.5*(Ts(i)+Ts(i+1));            // Average Temperature
        rho_avg = rho_air(T_avg);               // density of air
        mu_avg = mu_air(T_avg);                 // dynamic viscosity of air
        cp_avg = cp_air(T_avg);                 // Heat capacity air
        k_avg = k_air(T_avg);                   // thermal conductivity of air
        Pr_avg = Pr_air(T_avg);                 // Prantel Number of Air
        Pr_avg_e = Pr_air(Ts(i+1));             // Prantel Number of air for next control volume
        nu_avg = nu_air(T_avg);                 // Nusselt Number of air
        As = data.N_Ts*data.S_Ts*data.H_y;                     // Effective Area 

        mu_s = Q/As;                           
        Re_avg = mu_s*data.D_y/(nu_avg);                         // Reynolds Numbers
        
        if Re_avg > 10e3 || Re_avg < 2*10e5 then            // Determine correct Nusselt number equation based on Re  
            Nu_d =  0.72*data.F_s*Re_avg^0.63*Pr_avg^0.365*(Pr_avg/Pr_avg_e)^0.25;   // Nusselt Number
        elseif Re_avg > 2*10e5 || Re_avg < 2*10e6  then
            Nu_d = 0.033*data.F_s*Re_avg^0.8*Pr_avg^0.4*(Pr_avg/Pr_avg_e)^0.25;      // Nusselt Number
        end,
        m_y = (%pi)/4*rho_avg*(data.D_y)^2*data.H_y*data.N_Ls*data.N_Ts;
        h_avg = Nu_d*k_avg/(data.D_y);
        A_y = %pi*data.D_y*data.N_Ls*data.N_Ts;
        // change in yogurt stack temperature
        DX(3+i) = (rho_avg*Q)/(3520*m_y)*(Ts(i)-Ts(i+1))
    end
endfunction


//--- Simulation ------------------------------------------------------------
//fHndl = mopen(thispath+'output-00.txt','wt');       // Output file

//
// Enter your main program here
//



Q= 3.1734
N_set = 900
t = 1;


//[Q,T,Tr,Ts] = Workpoint(t,X,T,Tr,Ts,N_set)

//disp('workpoint Q')
//disp(Q)

//disp('Workpoing T')
//disp(T)


//[DP]= DP_Sys(Tr,Ts,Q);

//disp("sys DP")
//disp(DP)


//[DP] = DP_Fan(Q,N_set)
//disp('pump DP')
//disp(DP)

[DX,t,X,T,Tr,Ts,N] = Derivatives(t,X,T,Tr,Ts,N_set,DX)

disp(DX)


//mclose(fHndl)
//---------------------------------------------------------------------------

// please learn us something fundamentally new.... we have done similiar projects like these
// it was literally just copy paste.








