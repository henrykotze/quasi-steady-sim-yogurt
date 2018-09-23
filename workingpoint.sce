clear()
clearglobal()


//--- Input and defs --------------------------------------------------------

thispath = "C:\Users\Henry\Desktop\quasi-steady-sim-yogurt\";

exec(thispath+'func-thermo.sce', -1)   // Thermo funcs and consts
exec(thispath+'input_vars.sce', -1)      // Input data (Will change)
//exec(thispath+'func-simul.sce', -1)    // Put all your functions here!!!!


L_b
// Variables

// number of Rockbed control volumes
nb = 3;

// number of Yogurt stacks
ns = 3;

// Initial Values
T = [ 293.15 293.15 293.15 293.15 293.15 293.15];

// Rock bed air Temps:
Tr = [293.15 293.15 293.15 293.15]

//Yogurt Stack air temps:
Ts = [293.15 293.15 293.15 293.15];

//Yogurt mid points
Ty = [293.15 293.15 293.15]

// spmething
Tb = [293.15 293.15 293.15]

num_iteration = 10;

// Volume flow rate;
Qr =11;
Ab = W_b*H_b;
Nu_d = 2000;



for k=1:num_iteration
    T(2) = T(1);
    Tr(1) = T(2);
    
    for j = 1:n_b
        disp(j)
        T_avg = 0.5*(Tr(j)+Tr(j+1));
        rho_avg = rho_air(T_avg);
        mu_avg = mu_air(T_avg);
        cp_avg = cp_air(T_avg);
        k_avg = k_air(T_avg);
        G_avg = rho_avg*(Qr/Ab);
        Re_avg = G_avg*d_p/(mu_avg);                    // Reynolds Num
        f_avg = 4.466*(Re_avg)^0.75*(psi)^0.696*(eps_b)^-2.945*exp(11.85*(log(psi))^2);
        Nu_avg = 0.437*(Re_avg)^0.75*(psi)^3.35*(eps_b)^3.35*exp(29.03*((log(psi)^2)));
        h_avg = (k_avg*Nu_avg)/(d_p^2);
        Tr(j+1) = Tb(j) - (Tb(j)-Tr(j))*exp((h_avg*(L_b/n_b))/(-cp_avg*G_avg));
        //do somethinh
        //do something
    end
    T(3) = Tr(n_b+1);
    
    // Heater
    T_avg = (T(3)+T(4))/2;
    cp_avg = cp_air(T_avg);
    rho_avg = rho_air(T_avg);
    q_h = KP*(T_set-Ts(1))+KI*0
    T(4) = T(3) + q_h/(cp_avg*rho_avg*Qr)
    
    // Environemental Losses
    T_avg = (T(4)+T(5))/2;
    cp_avg = cp_air(T_avg);
    rho_avg = rho_air(T_avg);
    T(5) = (T(4) - T_env)/(exp((A_env*h_env)/(rho_avg*Qr*cp_avg)));
    Ts(1) = T(5);
    
    // Yogurt Warming
    for i = 1:n_s
        T_avg = 0.5*(Ts(j)+Ts(j+1));
        rho_avg = rho_air(T_avg);
        mu_avg = mu_air(T_avg);
        cp_avg = cp_air(T_avg);
        k_avg = k_air(T_avg);
        Pr_avg = Pr_air(T_avg);
        Pr_avg_e = Pr_air(Ts(i+1));
        nu_avg = nu_air(T_avg);
        As = N_Ts*S_Ts*H_y;
        Re_avg = mu_avg*D_y/(nu_avg);
        m_y = (%pi)/4*rho_avg*(D_y)^2*H_y*N_Ls*N_Ts;
        h_avg = Nu_d*k_avg/(D_y);
        A_y = %pi*D_y*N_Ls*N_Ts;
        
        Ts(i+1) = Ty(i) + (Ty(i) - Ts(i))*exp( -1*(A_y*h_avg)/(rho_avg*Qr*cp_avg));
    end
    
    T(6) = Ts(n_s+1);
end
disp('T6');
disp(T(6));
disp('T1');
disp(T(1));

disp('T')
disp(T);
disp('Ts')
disp(Ts);
disp('Tr')
disp(Tr);
disp('Ty')
disp(Ty)


