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
exec(thispath+'func-simul.sce', -1)    // Put all your functions here!!!!

//--- Global values ---------------------------------------------------------

global data                            // Input data
data = InputDataSet()                  // Declared in 'input-xx.sce'

// Declare all your global variables and constants here
// For example:
//
//    global cnst                      // User constants 
//    cnst = struct()
//
//    cnst.Pfan   = poly(data.DPpol, 'Q', 'coeff')  // Fan polynomial 
//    cnst.err    = 1e-8               // Convergence limit 
//    cnst.nmax   = 50                 // Maximum number of iterations
//

//--- Simulation ------------------------------------------------------------
fHndl = mopen(thispath+'output-00.txt','wt');       // Output file

//
// Enter your main program here
//

mclose(fHndl)
//---------------------------------------------------------------------------

// please learn us something fundamentally new, rather than wasting our time
// a project.
