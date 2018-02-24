#include <stdio.h>
#include <stdlib.h>
#include <tk.h>
#include "photo.hh"

extern "C" {
  int Inti_Init(Tcl_Interp *interp);
	   };

extern "C" {
  int Inti_SafeInit(Tcl_Interp *interp);
	   };


extern int Inti_Init(Tcl_Interp *interp) {
	Tcl_CreateCommand(interp, "Intinew", Inti_new, (ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);

	init_retcc();
	return TCL_OK;
}

extern int Inti_SafeInit(Tcl_Interp *interp) {
       Tcl_SetResult(interp,"Safe interpreter not supported",TCL_VOLATILE);
       //sprintf(interp->result,"Safe interpreter not supported");
       return TCL_ERROR;
}
