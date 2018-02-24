#include "Image.hh"
#include "histo.hh"
#include "morfo.hh"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void histogramme(Image<int> p,int* histo) {
  int index=0;
  int* PIX=p.PI();
  int L=p.PL();
  int W=p.PW();
  int i,j;
  
  for (i=0;i<257;i++) histo[i]=0;
  for (j=0;j<L;j++)
    for (i=0;i<W;i++) 
      histo[PIX[index++]]++;
  //histo[256] représente la valeur de l'effectif maximum
  for (i=0;i<256;i++) if (histo[256]<histo[i]) histo[256]=histo[i];    
}

void TransAffine(Image<int>& p,double a,double b){
  int index;
  int* PIX=p.PI();
  int L=p.PL();
  int W=p.PW();
  
  for (int j=0;j<L;j++)
    for (int i=0;i<W;i++) {
      index = j * W + i;
	  PIX[index] = (int)(a * PIX[index] + b);
	  if (PIX[index] > 255) PIX[index] = 255;
  	  if (PIX[index] < 0) PIX[index] = 0; 
	}
}

void Egalisation(Image<int>& p){
  int index=0;
  int* PIX=p.PI();
  int L=p.PL();
  int W=p.PW();
  int i,j,k;
  int *histo;
  int *histo_cumule;
  
  
  histo = (int*)calloc(257,sizeof(int));
  histo_cumule = (int*)calloc(257,sizeof(int));
  for (i=0;i<257;i++) {
	histo[i] = 0;
	histo_cumule[i] = 0;
  }
  for (j=0;j<L;j++)
    for (i=0;i<W;i++)
	  histo[PIX[index++]]++;
  for (i=0;i<256;i++) { 
	for (k=i;k<256;k++)
	  histo_cumule[k] +=  histo[i];
  }
  index=0;
  for (j=0;j<L;j++)
    for (i=0;i<W;i++) {
	  PIX[index] = 255 * histo_cumule[PIX[index]]/(L*W);
	  index++;
	}
  free(histo);
  free(histo_cumule);
}

void Histogramme_Cumule(Image<int> p,int* histo_cumule) {
  int index=0;
  int* PIX=p.PI();
  int L=p.PL();
  int W=p.PW();
  int i,j,k;
  int *histo;

  histo = (int*)calloc(257,sizeof(int));
  for (i=0;i<257;i++) {
	histo[i]=0;
	histo_cumule[i] = 0;
  }
  for (j=0;j<L;j++)
    for (i=0;i<W;i++) 
      histo[PIX[index++]]++;
  for (i=0;i<256;i++) { 
	for (k=i;k<256;k++)
	  histo_cumule[k] += histo[i];
  }
  //histo_cumule[256] représente la valeur de l'effectif maximum
 histo_cumule[256]=histo_cumule[255]; 
  
  free(histo);
}

void SeuilSimple(Image<int>& p,int seuil) {
  int index=0;
  int* PIX=p.PI();
  int L=p.PL();
  int W=p.PW();
  int i,j;

  for (j=0;j<L;j++)
    for (i=0;i<W;i++){ 
	  if (PIX[index] > seuil) PIX[index] = 255; 
	  else PIX[index] = 0 ;
	  index++;
	}
}

void SeuilDouble(Image<int>& p,int seuil1,int seuil2) {
  int index=0;
  int* PIX=p.PI();
  int L=p.PL();
  int W=p.PW();
  int i,j;

  for (j=0;j<L;j++) 
    for (i=0;i<W;i++){
      if ((PIX[index] > seuil1)&&(PIX[index] < seuil2)) PIX[index] = 255;
	  else PIX[index] = 0;
      index++;
    }
}


