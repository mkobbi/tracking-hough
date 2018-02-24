#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tk.h>
#include "photo.hh"
#include "Image.hh"
#include "hough_line.hh"
#include "hough_circle.hh"
#include "hough_track.hh"
#include "contour.hh"

/* stockage des infos d'image */
Tk_PhotoImageBlock *blockPtr;
template <class T> struct infoimage{
  Image<int> img;
  Tk_PhotoHandle handle;
  
};

// stockage de l'espace des paramètres
double *HOUGH;
double vote_max;

// stockage de la R-Table (Générale)
deplacement_s **RTABLE;

// stockage des n meilleures positions
int *BestPos;

/* Conversion en une image Inti avec allocation eventuelle */

Image<int> tkimage2Inti(Tk_PhotoImageBlock *tkim){
       int i,j,n;
       int *pixrim;
       unsigned char *pixtk;
       int tw=tkim->width;
       int th=tkim->height;
       Image<int> p(tw,th);

       pixrim=p.PI();
       for (j=0;j<tkim->height;j++) {
          pixtk=tkim->pixelPtr + tkim->pitch*j;
          for (i=0;i<tkim->width;i++) {
	   //  if ((*(pixtk+tkim->offset[0]) != 255)&&(*(pixtk+tkim->offset[1]) != 255) 
// 		&&(*(pixtk+tkim->offset[2]) != 255))
              n = (int)((*(pixtk+tkim->offset[0]) + *(pixtk+tkim->offset[1]) + *(pixtk+tkim->offset[2]))/3);
	    // else 
// 	      n = 500+(*(pixtk+tkim->offset[0]) == 255)+2*(*(pixtk+tkim->offset[1]) == 255)+4*(*(pixtk+tkim->offset[2]) == 255);
              *pixrim = n;
              pixtk+=tkim->pixelSize;pixrim++;
          }
       }
       return p;
}

/* conversion en une image tk avec allocation eventuelle */

Tk_PhotoImageBlock *Inti2tkimage(Image<int>& pp,Tk_PhotoImageBlock *tkim) {
       unsigned char *tempotk;
       int *tempoplan;
       int i;
       if (tkim == NULL) {
              tkim = (Tk_PhotoImageBlock *)malloc(sizeof(Tk_PhotoImageBlock));
	      tkim->pixelPtr = NULL;
       }

       tkim->width=pp.PW();
       tkim->height=pp.PL();
       tkim->pixelSize=3;
       tkim->pitch=tkim->width*tkim->pixelSize;
       tkim->offset[0]=0;
       tkim->offset[1]=1;
       tkim->offset[2]=2;

       if (tkim->pixelPtr == NULL) {
	 tempotk= tkim->pixelPtr=(unsigned char *)malloc(tkim->pixelSize*tkim->width*tkim->height*sizeof(unsigned char));
       } else {
	 tempotk= tkim->pixelPtr=(unsigned char *)realloc(tkim->pixelPtr,tkim->pixelSize*tkim->width*tkim->height*sizeof(unsigned char));
       }
       tempoplan=pp.PI();
       for (i=0;i<tkim->width*tkim->height;i++) {
	 if (*tempoplan < 256) {
       	*tempotk = *(tempotk+tkim->offset[1]) = *(tempotk+tkim->offset[2]) = (unsigned char)(*tempoplan);
	 } else {
	   int c=(*tempoplan)-500;
	   *(tempotk+tkim->offset[0]) = 255*(c%2); 
	   *(tempotk+tkim->offset[1]) = 255*((c/2)%2);
	   *(tempotk+tkim->offset[2]) = 255*((c/4)%2);}
	tempoplan++;
	tempotk+=tkim->pixelSize;
       }
       return tkim;
}

void freetkim(Tk_PhotoImageBlock *tkim) {
       free(tkim->pixelPtr);
       free(tkim);
}



int Inti_cmd(ClientData clientData, Tcl_Interp *interp, int argc, char *argv[])
{
       infoimage<int>* cdata;
       cdata  = (infoimage<int>*) clientData;
       
       float x,y,t1,t2;
       double theta,sigma;
       int a,b,c,d,e,f,g,n;
       int chg=1;
       int colour,orig_x,orig_y,nb_ech;
       char* numb;
       char *message;
       float *BestCircles;
       int xbest,ybest;
       
       
       if (strcmp(argv[1],"change_bord") == 0) {
	 a = atoi(argv[2]);
	 ChangeBord(a);chg=0;}
       else
       if (strncmp(argv[1],"save",7) == 0) {
	 a = atoi(argv[2]);
	 cdata[a].img.Imagetopgm(argv[3]);chg=0;}
        else 
	  if (strcmp(argv[1],"Reset_Hough2d") == 0) {
	   a = atoi(argv[2]);b = atoi(argv[3]);
	   Reset_Hough2d(HOUGH,a,b);
	   chg = 0;}
	else
	  if (strcmp(argv[1],"MaJ_Hough2d") == 0) {
	    a = atoi(argv[2]);
	    MaJ_Hough2d(cdata[a].img,HOUGH,vote_max);
	 }
      else
	 if  (strcmp(argv[1],"Dessine_rectangle") == 0) {
	   a = atoi(argv[2]); b = atoi(argv[3]);
	   c = atoi(argv[4]); d = atoi(argv[5]); 
	   e = atoi(argv[6]); colour = atoi(argv[7]); 
	   Trace_rectangle(cdata[a].img,b,c,d,e,colour);
	 }
	else 
	 if (strcmp(argv[1],"Allocate_RTable") == 0) {
	   a = atoi(argv[2]);
	   Allocate_RTable(&RTABLE,a); 
	   chg = 0;
	 }
       else 
	 if (strcmp(argv[1],"Calcul_RTable_Proto_ord0") == 0) {
	   a = atoi(argv[2]);b = atoi(argv[3]);c = atoi(argv[4]); 
	   x = atof(argv[5]);d = atoi(argv[6]);
	   Create_RTable_ordre0(RTABLE,cdata[a].img,cdata[b].img,cdata[c].img,x,d);
	   chg = 23;
	 }
       else 
	 if (strcmp(argv[1],"Calcul_RTable_Proto_ord1") == 0) {
	   a = atoi(argv[2]);b = atoi(argv[3]);c = atoi(argv[4]); 
	   x = atof(argv[5]);d = atoi(argv[6]);
	   Create_RTable_ordre1(RTABLE,cdata[a].img,cdata[b].img,cdata[c].img,x,d);
	   chg = 23;
	 }
       else 
	 if (strcmp(argv[1],"Calcul_RTable_Proto_ord2") == 0) {
	   a = atoi(argv[2]);b = atoi(argv[3]);c = atoi(argv[4]); 
	   x = atof(argv[5]);d = atoi(argv[6]);
	   Create_RTable_ordre2(RTABLE,cdata[a].img,cdata[b].img,cdata[c].img,x,d);
	   chg = 23;
	 }
        else 
        if (strcmp(argv[1],"Allocate_Hough_2d") == 0) {
	   a = atoi(argv[2]);
           b = atoi(argv[3]);
	   printf("Allocation du tableau de Hough 2d\n");
	   printf("Taille du tableau : %d x %d\n",a,b);
	   HOUGH = (double *)calloc(a*b,sizeof(double));
	   chg = 0;
	 }
       else 
	 if (strcmp(argv[1],"Hough_General_ord0") == 0) {
	   a = atoi(argv[2]);
	   x = atof(argv[3]);b = atoi(argv[4]);
	   Hough_General_ordre0(RTABLE,HOUGH,cdata[a].img,x,b,&vote_max);
	   chg = 0;
	 }
       else 
	 if (strcmp(argv[1],"Hough_General_ord1") == 0) {
	   a = atoi(argv[2]);
	   x = atof(argv[3]);b = atoi(argv[4]);
	   Hough_General_ordre1(RTABLE,HOUGH,cdata[a].img,x,b,&vote_max);
	   chg = 0;
	 }
       else 
	 if (strcmp(argv[1],"Hough_General_ord2") == 0) {
	   a = atoi(argv[2]);
	   x = atof(argv[3]);b = atoi(argv[4]);
	   Hough_General_ordre2(RTABLE,HOUGH,cdata[a].img,x,b,&vote_max);
	   chg = 0;
	 }
         else 
	 if (strcmp(argv[1],"Allocate_BestPos") == 0) {
	   a = atoi(argv[2]);
           printf("Allocation du tableau de positions\n");
	   printf("Nombre de positions : %d\n",a);
	   BestPos = (int*)calloc(2*a + 1,sizeof(int)); 
	   chg = 0;
	 }
         else 
	 if (strcmp(argv[1],"Find_Best_THG") == 0) {
	   a = atoi(argv[2]);b = atoi(argv[3]);//Dimensions images
           c = atoi(argv[4]);//Nombre de positions
           d = atoi(argv[5]);e = atoi(argv[6]);//Dimensions proto
	   Find_Best_GHT(HOUGH,a,b,c,BestPos,d,e);    
	   numb = (char*)malloc(10*sizeof(char));
           for (n = 1;n <= 2*BestPos[0];n++) {
	     sprintf(numb,"%d",BestPos[n]);
	     Tcl_AppendElement(interp,numb);
	     }
	   chg = 0;}
         else 
	 if (strcmp(argv[1],"Update_Proto") == 0) {
            a = atoi(argv[2]);//Prototype
            b = atoi(argv[3]);//Index prototype
            c = atoi(argv[4]);//Poids prototype
            d = atoi(argv[5]);//Image courante
            e = atoi(argv[6]);f = atoi(argv[7]);//Position courante du prototype
            g = atoi(argv[8]);//Nombre de labels
            x = atof(argv[9]);
            Update_Proto_RTable(RTABLE,cdata[a].img,cdata[b].img,cdata[c].img,cdata[d].img,
                                BestPos,e,f,&xbest,&ybest,g,x);
            numb = (char*)malloc(10*sizeof(char));
            sprintf(numb,"%d",xbest);Tcl_AppendElement(interp,numb);
            sprintf(numb,"%d",ybest);Tcl_AppendElement(interp,numb);
            chg = 123;}
	 else {
                sprintf(message,"Commande %s inconnue !",argv[1]);
                Tcl_SetResult(interp,message,TCL_VOLATILE);
		return TCL_ERROR;
	}
       if (chg==1) {
       blockPtr=Inti2tkimage(cdata[a].img,blockPtr);
       Tk_PhotoPutBlock(cdata[a].handle, blockPtr, 0, 0, blockPtr->width, blockPtr->height);
       } else  if (chg==2) {
       blockPtr=Inti2tkimage(cdata[b].img,blockPtr);
       Tk_PhotoPutBlock(cdata[b].handle, blockPtr, 0, 0, blockPtr->width, blockPtr->height);
       } else  if (chg==23) {
	 blockPtr=Inti2tkimage(cdata[b].img,blockPtr);
	 Tk_PhotoPutBlock(cdata[b].handle, blockPtr, 0, 0, blockPtr->width, blockPtr->height);
       blockPtr=Inti2tkimage(cdata[c].img,blockPtr);
       Tk_PhotoPutBlock(cdata[c].handle, blockPtr, 0, 0, blockPtr->width, blockPtr->height);
       } else  if (chg==123) {
         blockPtr=Inti2tkimage(cdata[a].img,blockPtr);
         Tk_PhotoPutBlock(cdata[a].handle, blockPtr, 0, 0, blockPtr->width, blockPtr->height);
	 blockPtr=Inti2tkimage(cdata[b].img,blockPtr);
	 Tk_PhotoPutBlock(cdata[b].handle, blockPtr, 0, 0, blockPtr->width, blockPtr->height);
         blockPtr=Inti2tkimage(cdata[c].img,blockPtr);
         Tk_PhotoPutBlock(cdata[c].handle, blockPtr, 0, 0, blockPtr->width, blockPtr->height);
       }
        return TCL_OK;       
}

int img_new(ClientData clientData, Tcl_Interp *interp, int argc, char *argv[])
{
       infoimage<int> *cdata;
       Tk_PhotoImageBlock blockPtr;
       int x = atoi(argv[2]);
       cdata = (infoimage<int>*) clientData;
       char *error_message;
 
     if (argc != 3) {
                error_message = (char *)("Usage: Imgnew <image> <numero>");
                Tcl_SetResult(interp,error_message,TCL_VOLATILE);
		return TCL_ERROR;
	}
   cdata[x].handle=Tk_FindPhoto(interp,argv[1]);
   Tk_PhotoGetImage(cdata[x].handle, &blockPtr);
   cdata[x].img=tkimage2Inti(&blockPtr);
 return TCL_OK;
}


/* creation d'une image de type Inti */
/* Intinew <nomimage> <nomcommanderet> */

void init_retcc() {
blockPtr = NULL;
}

int Inti_new(ClientData clientData, Tcl_Interp *interp, int argc, char *argv[])
{
       infoimage<int> *cdata;
       int width,height;
       char *error_message;

       if (argc < 3) {
         error_message = (char *)("Usage: Intinew <nomcommande> <nb_images> <largeur(opt)> <hauteur(opt)>");
         Tcl_SetResult(interp,error_message,TCL_VOLATILE);
	 return TCL_ERROR;
       }
       int nbimages = atoi(argv[2]);
       if (argc > 4) {
	 width = atoi(argv[3]);
	 height = atoi(argv[4]);
       }
       cdata = new infoimage<int>[nbimages];
       Tcl_CreateCommand(interp, "Imgnew", img_new, (ClientData)cdata,(Tcl_CmdDeleteProc *)NULL);
       Tcl_CreateCommand(interp, argv[1], Inti_cmd, (ClientData)cdata,(Tcl_CmdDeleteProc *)NULL);
       return TCL_OK;
}


