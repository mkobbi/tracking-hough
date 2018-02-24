#include <stdio.h>
#include <stdlib.h>
#include <tk.h>
#include "photo.hh"

/* extern int matherr();
int *tclDummyMathPtr = (int *) matherr;
*/


int Mon_AppInit(Tcl_Interp *interp);

int main(int argc, char *argv[])
{
	Tk_Main(argc, argv, (Tcl_AppInitProc *)Tcl_AppInit);
	return(0);
}

int Tcl_AppInit(Tcl_Interp *interp)
{
	if (Tcl_Init(interp) == TCL_ERROR) {
		return TCL_ERROR;
	}
	
   	if (Tk_Init(interp) == TCL_ERROR) {
		return TCL_ERROR;
    	}
	/* Initialisation du pointeur tampon image */
	Tcl_CreateCommand(interp, "Intinew", Inti_new, (ClientData)NULL,
		(Tcl_CmdDeleteProc *)NULL);
       init_retcc();
			
	Tcl_SetVar(interp, "tcl_rcFileName", "~/.tclshrc", TCL_GLOBAL_ONLY);
	return TCL_OK;
}
