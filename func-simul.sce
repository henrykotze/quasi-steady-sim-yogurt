//===========================================================================
// SIMULATION FUNCTIONS
//
// Name:   HAR de Werker
// Std nr: 12345678
//===========================================================================

//--- Fan pressure increase ----------------------------------------------------
function DP = Fan_DP(T, Q, N) //.... You can add more or other params here

    global data cnst                     // Globals declared in main script

    rpm = max(N, data.rpm_min)           // Limit fan speed
    rpm = min(rpm, data.rpm_max)

	Tf = ....  
    rho = rho_air(Tf)                    // rho_air() declared in 'func-thermo.sci' 

    //
    //  DO THE REST eg. Check Q limits, etc.
    //
	// DPref = horner(cnst.Pfan, Qref)   // Use horner() function for polyfit
	// 

endfunction


//
//  OTHER FUNCTIONS
//
//  //--- Calculate .....
//  function [a,b] = afunc(T, Q, X, ....)
//     global data cnst  
//     //...............
//     // a = ...
//     // b = ...
//  endfunction
//
