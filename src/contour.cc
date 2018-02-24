#include "Image.hh"
#include "filtre.hh"
#include "morfo.hh"
#include "contour.hh"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Détection de contour, méthode de Marr et Hildreth :
// (1) Calcul du Laplacien de Gaussienne (LoG)
// (2) Calcul des passage par zéro du LoG
// (3) Seuillage par hystérésis (relativement au 
// module du gradient)
void Marr_Hildreth(Image<int>& p,float sigma,float th_haut,float th_bas) { 
  int h=p.PL();
  int w=p.PW();
  int i,j,val;
  int index=0;
  int *PIX = p.PI();
  Image<int> qq(w,h);
  int *QPIX = qq.PI();
  double *VAL,*DX,*DY,*GX,*GY,*GXX;
  double valmin,valmax;

  VAL = (double *)calloc(w*h,sizeof(double));
  GXX = (double *)calloc(w*h,sizeof(double));
  // Passage en double...
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      VAL[index] = (double) (PIX[index++]);
  // Calcul du gradient...
  
  DX = (double *)calloc(w*h,sizeof(double));
  DY = (double *)calloc(w*h,sizeof(double));
  Gauss_Rec(VAL,GXX,w,h,0,1,sigma,valmin,valmax);
  Gauss_Rec(GXX,DX,w,h,1,0,sigma,valmin,valmax);
  Gauss_Rec(VAL,GXX,w,h,0,0,sigma,valmin,valmax);
  Gauss_Rec(GXX,DY,w,h,1,1,sigma,valmin,valmax);
  // Calcul du module...
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++) {
      DY[index] = sqrt(DX[index]*DX[index] + DY[index]*DY[index]);
      index++;
    }
  // Calcul du laplacien...
  GX = (double *)calloc(w*h,sizeof(double));
  GY = (double *)calloc(w*h,sizeof(double));
  Gauss_Rec(VAL,GXX,w,h,0,1,sigma,valmin,valmax);
  Gauss_Rec(GXX,GX,w,h,2,0,sigma,valmin,valmax);
  Gauss_Rec(VAL,GXX,w,h,0,0,sigma,valmin,valmax);
  Gauss_Rec(GXX,GY,w,h,2,1,sigma,valmin,valmax);
  // Calcul du laplacien et détection des passages par zéro...
  // ...avec seuillage par hystérésis du gradient
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++) {
      GY[index] = GX[index] + GY[index];
      index++;
    }
  for (i=0;i<w-1;i++) {PIX[i]=0;QPIX[i]=0;}
  for (i=w;i<w*h;i+=w) {PIX[i]=0;QPIX[i]=0;}
  for (i=2*w-1;i<w*h;i+=w) {PIX[i]=0;QPIX[i]=0;}
  for (i=(h-1)*w;i<w*h;i++) {PIX[i]=0;QPIX[i]=0;}
  for (i=1;i<w-1;i++)
    for (j=1;j<h-1;j++) {
      index = j*w+i;
      if (((GY[index]*GY[index-1]<0)||
	   (GY[index]*GY[index-w]<0)||
	   (GY[index]*GY[index-1-w]<0))) {
	if (DY[index]>th_haut) {PIX[index]=255;QPIX[index]=255;}
	else if (DY[index]>th_bas) {PIX[index]=0;QPIX[index]=255;}
	else {PIX[index]=0;QPIX[index]=0;}
      } else {PIX[index]=0;QPIX[index]=0;}
    }
  Reconstruit(p,qq,4);
  p = !p;
  free(DX);
  free(DY);
  free(GX);
  free(GY);
  free(GXX);
  free(VAL); 
}

// Détection de contour, méthode de Canny
// Implantation de Jean-Baptiste Desmottes (ENSTA'09)
// (1) Calcul du module du gradient (échelle sigma)
// (2) Sélection des extrema locaux dans la direction du gradient
// (3) Seuillage par hystérésis

void Canny(Image<int>& p,float sigma,float th_haut,float th_bas) { 
   int h=p.PL();
   int w=p.PW();
   int i,j;
   int index=0;
   int *PIX = p.PI();
   Image<int> qq(w,h);
   int *QPIX = qq.PI();
   double *VAL,*DX,*DY,*TMP,*MGR;
   double valmin,valmax;
   double vAbsDY,vAbsDX;
   double a;
   
   double tgPIsur8 = 0.4142;
 
   VAL = (double *)calloc(w*h,sizeof(double));
   TMP = (double *)calloc(w*h,sizeof(double));
   MGR = (double *)calloc(w*h,sizeof(double));
   DX  = (double *)calloc(w*h,sizeof(double));
   DY  = (double *)calloc(w*h,sizeof(double));
   
   // Passage en double...
   index=0;
   for (i=0;i<w;i++)
     for (j=0;j<h;j++)
       VAL[index] = (double) (PIX[index++]);

   // Calcul du gradient...
   Gauss_Rec(VAL,TMP,w,h,0,1,sigma,valmin,valmax);
   Gauss_Rec(TMP,DX,w,h,1,0,sigma,valmin,valmax);
   Gauss_Rec(VAL,TMP,w,h,0,0,sigma,valmin,valmax);
   Gauss_Rec(TMP,DY,w,h,1,1,sigma,valmin,valmax);
   // Calcul du module
   valmin=255;valmax=0;
   index=0;
   for (i=0;i<w;i++)
     for (j=0;j<h;j++) 
     {
       MGR[index] = sqrt(DX[index]*DX[index] + DY[index]*DY[index]);
       if (MGR[index]>valmax)
         valmax=MGR[index];
       else if (MGR[index]<valmin)
         valmin=MGR[index];
       index++;
     }

   // Suppression des non-maxima
   //     on ignore les bords
   for (i=0;i<w-1;i++) TMP[i]=0;
   for (i=w;i<w*h;i+=w) TMP[i]=0;
   for (i=w-1;i<w*h;i+=w) TMP[i]=0;
   for (i=(h-1)*w;i<w*h;i++) TMP[i]=0;
   for (i=1;i<w-1;i++)
     for (j=1;j<h-1;j++)
     {
       index = j*w+i;
       vAbsDX = valabs(DX[index]) ;
       vAbsDY = valabs(DY[index]) ;
       
       if( vAbsDX<0.01 && vAbsDY<0.01 )
       {
         TMP[index]=0;
         continue;
       }

       if( vAbsDX < vAbsDY ) //orientation verticale
       {
         a=vAbsDX/vAbsDY;
         
         
         if( DX[index]*DY[index] < 0 )  // Nord-Est/Sud-Ouest
         {
           if( (a<tgPIsur8 && MGR[index]>MGR[index-w] && MGR[index]>MGR[index+w])
                || (a>=tgPIsur8 && MGR[index]>MGR[index-w+1] && MGR[index]>MGR[index+w-1]) )
             TMP[index]=MGR[index];
           else
             TMP[index]=0;
         }
         else // Nord-Ouest/Sud-Est
         {
           if( (a<tgPIsur8 && MGR[index]>MGR[index-w] && MGR[index]>MGR[index+w])
                || (a>=tgPIsur8 && MGR[index]>MGR[index-w-1] && MGR[index]>MGR[index+w+1]) )
             TMP[index]=MGR[index];
           else
             TMP[index]=0;
         }
       }
       else  //orientation horizontale
       {
         a=vAbsDY/vAbsDX;
         if( DX[index]*DY[index] < 0 )  // Nord-Est/Sud-Ouest
         {
           if( (a<tgPIsur8 && MGR[index]>MGR[index+1] && MGR[index]>MGR[index-1])
               || (a>=tgPIsur8 && MGR[index]>MGR[index-w+1] && MGR[index]>MGR[index+w-1]) )
             TMP[index]=MGR[index];
           else
             TMP[index]=0;
         }
         else // Nord-Ouest/Sud-Est
         {
           if( (a<tgPIsur8 && MGR[index]>MGR[index+1] && MGR[index]>MGR[index-1])
              || (a>=tgPIsur8 && MGR[index]>MGR[index-w-1] && MGR[index]>MGR[index+w+1]) )
             TMP[index]=MGR[index];
           else
             TMP[index]=0;
         }
       }
     }

   // Seuillage par hystérésis
   index=0;
   for (i=0;i<w;i++)
     for (j=0;j<h;j++) 
     {
       index = j*w+i;
       if (TMP[index]>th_haut) {PIX[index]=255; QPIX[index]=255;}
       else if (TMP[index]>th_bas) {PIX[index]=0; QPIX[index]=255;}
       else {PIX[index]=0; QPIX[index]=0;}
     }

   Reconstruit(p,qq,4);
   p = !p;

   free(DX);
   free(DY);
   free(TMP);
   free(MGR);
   free(VAL); 
}

// Méthode "repère local" (g,t) :
// (1) Calcul des dérivées à l'ordre 1,2 et 3
// dans la direction du gradient g
// (2) Sélection des passages par zéros de Igg
// avec la condition Iggg < 0
// (3) Seuil par hystérésis (relativement au 
// module du gradient)

void Contours_RepereLocal(Image<int>& p,float sigma,float th_haut,float th_bas) { 
  int h=p.PL();
  int w=p.PW();
  int i,j,val;
  int index=0;
  int *PIX = p.PI();
  Image<int> qq(w,h);
  int *QPIX = qq.PI();
  double *VAL,*DX,*DY,*DVV;
  double *TMP,*DXX,*DYY,*DXY;
  double *DXXX,*DXXY,*DXYY,*DYYY,*DVVV;
  double valmin,valmax;

  VAL = (double *)calloc(w*h,sizeof(double));
  TMP = (double *)calloc(w*h,sizeof(double));
  DVV = (double *)calloc(w*h,sizeof(double));
  // Passage en double...
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      VAL[index] = (double) (PIX[index++]);
  // Calcul des dérivées premières...
  DX = (double *)calloc(w*h,sizeof(double));
  DY = (double *)calloc(w*h,sizeof(double));
  Gauss_Rec(VAL,TMP,w,h,0,1,sigma,valmin,valmax);
  Gauss_Rec(TMP,DX,w,h,1,0,sigma,valmin,valmax);
  Gauss_Rec(VAL,TMP,w,h,0,0,sigma,valmin,valmax);
  Gauss_Rec(TMP,DY,w,h,1,1,sigma,valmin,valmax);
  // Calcul des dérivées secondes...
  DXX = (double *)calloc(w*h,sizeof(double));
  DYY = (double *)calloc(w*h,sizeof(double));
  DXY = (double *)calloc(w*h,sizeof(double));
  Gauss_Rec(VAL,TMP,w,h,0,1,sigma,valmin,valmax);
  Gauss_Rec(TMP,DXX,w,h,2,0,sigma,valmin,valmax);
  Gauss_Rec(VAL,TMP,w,h,0,0,sigma,valmin,valmax);
  Gauss_Rec(TMP,DYY,w,h,2,1,sigma,valmin,valmax);
  Gauss_Rec(VAL,TMP,w,h,0,1,sigma,valmin,valmax);
  Gauss_Rec(TMP,VAL,w,h,1,0,sigma,valmin,valmax);
  Gauss_Rec(VAL,DXY,w,h,1,1,sigma,valmin,valmax);
  // Calcul de la dérivée seconde dans la direction
  // du gradient Dvv = Dx^2 Dxx + 2 Dx Dy Dxy + Dyy Dy^2
  // + Calcul du module du gradient pour seuillage 
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++) {
      DVV[index] = DX[index]*DX[index]*DXX[index]
	+ 2*DX[index]*DY[index]*DXY[index]
	+ DY[index]*DY[index]*DYY[index];
      TMP[index] = sqrt(DX[index]*DX[index] + DY[index]*DY[index]);
      index++;
    }
  // Calcul des dérivées troisièmes...
  DXXX = (double *)calloc(w*h,sizeof(double));
  DYYY = (double *)calloc(w*h,sizeof(double));
  DXYY = (double *)calloc(w*h,sizeof(double));
  DXXY = (double *)calloc(w*h,sizeof(double));
  Gauss_Rec(DXX,DXXX,w,h,1,0,sigma,valmin,valmax);
  Gauss_Rec(DXX,DXXY,w,h,1,1,sigma,valmin,valmax);
  Gauss_Rec(DYY,DXYY,w,h,1,0,sigma,valmin,valmax);
  Gauss_Rec(DYY,DYYY,w,h,1,1,sigma,valmin,valmax);
  // Calcul de la dérivée troisieme dans la direction
  // du gradient Dvvv = Dx^3 Dxxx + 3 Dx^2 Dy Dxxy 
  // + 3 Dx Dy^2 Dxyy + Dyyy Dy^3 
  DVVV = (double *)calloc(w*h,sizeof(double));
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++) {
      DVVV[index] = DX[index]*DX[index]*DX[index]*DXXX[index]
	+ 3*DX[index]*DX[index]*DY[index]*DXXY[index]
	+ 3*DX[index]*DY[index]*DY[index]*DXYY[index]
	+ DY[index]*DY[index]*DY[index]*DYYY[index];
      index++;
    }

  for (i=0;i<w-1;i++) {PIX[i]=0;QPIX[i]=0;}
  for (i=w;i<w*h;i+=w) {PIX[i]=0;QPIX[i]=0;}
  for (i=2*w-1;i<w*h;i+=w) {PIX[i]=0;QPIX[i]=0;}
  for (i=(h-1)*w;i<w*h;i++) {PIX[i]=0;QPIX[i]=0;}
  for (i=1;i<w-1;i++)
    for (j=1;j<h-1;j++) {
      index = j*w+i;
      if (((DVV[index]*DVV[index-1]<0)||
	   (DVV[index]*DVV[index-w]<0)||
	   (DVV[index]*DVV[index-1-w]<0)) &&
	  ((DVVV[index]<0)&&(DVVV[index-1]<0)&&
	   (DVVV[index-w]<0)&&(DVVV[index-w-1]<0))){
	if (TMP[index]>th_haut) {PIX[index]=255;QPIX[index]=255;}
	else if (TMP[index]>th_bas) {PIX[index]=0;QPIX[index]=255;}
	else {PIX[index]=0;QPIX[index]=0;}
      } else {PIX[index]=0;QPIX[index]=0;}
    }
  Reconstruit(p,qq,4);
  p = !p;
  free(DX);
  free(DY);
  free(DXX);
  free(DYY);
  free(DXY);
  free(DVV);
  free(TMP);
  free(VAL); 
  free(DXXX);
  free(DXXY);
  free(DXYY);
  free(DYYY);
  free(DVVV);
}
