#include "Image.hh"
#include "filtre.hh"
#include "hough_track.hh"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Conversion routine for displaying the GHT
void MaJ_Hough2d(Image<int>& p_transf,double *HOUGH,double vote_max) {
  int* PIX=p_transf.PI();
  int w = p_transf.PW();
  int h = p_transf.PL();
  int index;

  for (index = 0;index < w*h;index++) 
    PIX[index] = (int)(255.0*HOUGH[index]/vote_max);
}

// Reset routine for the GHT
void Reset_Hough2d(double *HOUGH,int width_img,int height_img) {
  int index;
  for (index = 0;index < width_img*height_img;index++) HOUGH[index]=0; 
} 

void Trace_rectangle(Image<int>& p,int x_centre,int y_centre,int largeur,int hauteur,int couleur) {
  int* PIX=p.PI();
  int w = p.PW();
  int h = p.PL();
  int i,index;
  int xmin,xmax,ymin,ymax;
  
  xmin = max(x_centre-largeur/2,0);
  xmax = min(x_centre+largeur/2,w-1);
  ymin = max(y_centre-hauteur/2,0);
  ymax = min(y_centre+hauteur/2,h-1);

  for (i=xmin;i<=xmax;i++) {
    PIX[ymin*w+i]=500+couleur;
    PIX[ymax*w+i]=500+couleur;
  }
  for (i=ymin;i<=ymax;i++) {
    PIX[i*w+xmin]=500+couleur;
    PIX[i*w+xmax]=500+couleur;
  }
}

void Allocate_RTable(deplacement_s*** RTABLE,int taille) {
  int i;
  
  if ((*RTABLE) != NULL) free(*RTABLE);
  (*RTABLE) = (deplacement_s**)calloc(taille,sizeof(deplacement_s*));
  for (i=0;i<taille;i++) {
    (*RTABLE)[i] = (deplacement_s*)malloc(sizeof(deplacement_s));
    (*RTABLE)[i]->next = NULL;
  }
  //printf("R-Table allouée, taille = %d\n",taille);
}

void Insert_RTable(deplacement_s** RTABLE,int index,int delta_x,int delta_y,float omega) {
  deplacement_s *p;

  p = (deplacement_s*)malloc(sizeof(deplacement_s));
  p->dx = delta_x;
  p->dy = delta_y;
  p->poids = omega;
  p->next = RTABLE[index]->next;
  RTABLE[index]->next = p;
}

// Creating the GHT Model:
// Filling the R-Table from the Prototype Image
// 1st function: R-Table of order 0 (pixel value)

void Create_RTable_ordre0(deplacement_s** RTABLE,Image<int>& p_proto,Image<int>& p_valeur,Image<int>& p_poids,float sigma,int nb_labels) {
  int* PIX=p_proto.PI();
  int w = p_proto.PW();
  int h = p_proto.PL();
  int* PIX_V=p_valeur.PI();
  int* PIX_W=p_poids.PI();
  double *VAL, *TEMP;
  int i,j,index;
  double valmin,valmax;
  int index_val,delta_x,delta_y;

  printf("Calculate R-Table (order 0), sigma = %f, nb_labels = %d\n",sigma,nb_labels);
  VAL = (double *)calloc(w*h,sizeof(double));
  TEMP = (double *)calloc(w*h,sizeof(double));
  // Conversion in double...
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      VAL[index] = (double) (PIX[index++]);
  // Gaussian smoothing
  Gauss_Rec(VAL,TEMP,w,h,0,0,sigma,valmin,valmax);
  Gauss_Rec(TEMP,VAL,w,h,0,1,sigma,valmin,valmax);
  index=0;
  for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      index_val = (int)((VAL[index]*(nb_labels-1))/(255.0));
      delta_x = w/2 - i;
      delta_y = h/2 - j;
      Insert_RTable(RTABLE,index_val,delta_x,delta_y,1.0);
      PIX_V[index] = (int)(((float)(index_val)*255.0)/nb_labels);
      PIX_W[index] = 1;	
      index++;
    }
  free(VAL);free(TEMP);
}

// 2nd function: R-Table of order 1 (gradient direction)

void Create_RTable_ordre1(deplacement_s** RTABLE,Image<int>& p_proto,Image<int>& p_direction,Image<int>& p_poids,float sigma,int nb_labels) {
  int* PIX=p_proto.PI();
  int w = p_proto.PW();
  int h = p_proto.PL();
  int* PIX_D=p_direction.PI();
  int* PIX_W=p_poids.PI();
  double *VAL,*GX,*GY,*TEMP;
  int i,j,index;
  double valmin,valmax;
  double argument,norme,norme_max;
  int index_arg,delta_x,delta_y;

  float PI = 3.1415926536;
    
  printf("Calculate R-Table (order 1), sigma = %f, nb_labels = %d\n",sigma,nb_labels);
  VAL = (double *)calloc(w*h,sizeof(double));
  GX = (double *)calloc(w*h,sizeof(double));
  GY = (double *)calloc(w*h,sizeof(double));
  TEMP = (double *)calloc(w*h,sizeof(double));
  // Conversion in double...
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      VAL[index] = (double) (PIX[index++]);
  // Calculate the gradient components
  Gauss_Rec(VAL,TEMP,w,h,0,1,sigma,valmin,valmax);
  Gauss_Rec(TEMP,GX,w,h,1,0,sigma,valmin,valmax);
  Gauss_Rec(VAL,TEMP,w,h,0,0,sigma,valmin,valmax);
  Gauss_Rec(TEMP,GY,w,h,1,1,sigma,valmin,valmax);
  index=0;
  for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      if ((GX[index]!=0)||(GY[index]!=0)) {
	argument = atan2(GY[index],GX[index]);
	index_arg = (int)(((argument+PI)*nb_labels)/(2*PI));
	TEMP[index] = sqrt(GX[index]*GX[index]+GY[index]*GY[index]);
	if (norme_max<TEMP[index]) norme_max = TEMP[index];
	delta_x = w/2 - i;
	delta_y = h/2 - j;
	Insert_RTable(RTABLE,index_arg,delta_x,delta_y,sigma*(float)TEMP[index]);
	PIX_D[index] = (int)(((float)(index_arg)*255.0)/nb_labels);
      } else {
	PIX_D[index] = 0;
	TEMP[index] = 0;
      }
      index++;
    }
  // normalising the weights for display (inverse video)
  for (index=0;index <w*h;index++) 
    PIX_W[index] = (int)(255.0 - (TEMP[index]*255.0)/norme_max);
  free(VAL);free(GX);free(GY);free(TEMP);
}

// 3rd function: R-Table of order 2 (isophote curvature)
void Create_RTable_ordre2(deplacement_s** RTABLE,Image<int>& p_proto,Image<int>& p_courbure,Image<int>& p_poids,float sigma,int nb_labels) {
  int* PIX=p_proto.PI();
  int w = p_proto.PW();
  int h = p_proto.PL();
  int* PIX_C=p_courbure.PI();
  int* PIX_W=p_poids.PI();
  double *VAL,*GX,*GY,*TEMP;
  double *GXX,*GYY,*GXY,*CURV;
  int i,j,index;
  double valmin,valmax,norme,poids_max;
  float omega;
  int index_curv,delta_x,delta_y;
    
  printf("Calculate R-Table (order 2), sigma = %f, nb_labels = %d\n",sigma,nb_labels);
  VAL = (double *)calloc(w*h,sizeof(double));
  GX = (double *)calloc(w*h,sizeof(double));
  GY = (double *)calloc(w*h,sizeof(double));
  TEMP = (double *)calloc(w*h,sizeof(double));
  GXX = (double *)calloc(w*h,sizeof(double));
  GXY = (double *)calloc(w*h,sizeof(double));
  GYY = (double *)calloc(w*h,sizeof(double));
  CURV = (double *)calloc(w*h,sizeof(double));
 
  // Conversion in double...
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      VAL[index] = (double) (PIX[index++]);
  // Calculate the gradient components
  Gauss_Rec(VAL,TEMP,w,h,0,1,sigma,valmin,valmax);
  Gauss_Rec(TEMP,GX,w,h,1,0,sigma,valmin,valmax);
  Gauss_Rec(VAL,TEMP,w,h,0,0,sigma,valmin,valmax);
  Gauss_Rec(TEMP,GY,w,h,1,1,sigma,valmin,valmax);
  // Calculate G_tt/G_g (i.e. isophote curvature)
  // Calculate G_xx
  Gauss_Rec(VAL,TEMP,w,h,0,1,sigma,valmin,valmax);
  Gauss_Rec(TEMP,GXX,w,h,2,0,sigma,valmin,valmax);
  // then G_yy
  Gauss_Rec(VAL,TEMP,w,h,0,0,sigma,valmin,valmax);
  Gauss_Rec(TEMP,GYY,w,h,2,1,sigma,valmin,valmax);
  // and then G_xy
  Gauss_Rec(VAL,TEMP,w,h,1,0,sigma,valmin,valmax);
  Gauss_Rec(TEMP,GXY,w,h,1,1,sigma,valmin,valmax);
  // Calculate G_tt/G_g, isophote curvature, saved in CURV, 
  // and the interest function, equal to the total curvature
  // (Frobenius norm of the Hessian matrix), saved in VAL
  poids_max = 0;
  for (index=0;index<w*h;index++) {
    norme = sqrt(GX[index]*GX[index] + GY[index]*GY[index]);
    if (norme > 0) {
      CURV[index] = - GY[index]*GY[index]*GXX[index]
	+ 2*GX[index]*GY[index]*GXY[index]
	- GX[index]*GX[index]*GYY[index];
      VAL[index] = sqrt(GXX[index]*GXX[index] + 2*GXY[index]*GXY[index] + GYY[index]*GYY[index]);
      if (poids_max<VAL[index]) poids_max = VAL[index];
      CURV[index] /= pow(norme,3);
      // The significant curvature is restricted within [-1,+1]
      if (CURV[index] > 1.0)  CURV[index] = 1.0;
      if (CURV[index] < -1.0)  CURV[index] = -1.0;
    } else CURV[index] = 0;
  }
  index=0;
  for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      index_curv = (int)((nb_labels-1)*(CURV[index]+1)/2);
      delta_x = w/2 - i;
      delta_y = h/2 - j;
      omega = (float)(VAL[index]);
      Insert_RTable(RTABLE,index_curv,delta_x,delta_y,sigma*sigma*omega);
      PIX_C[index] = (int)(128.0 + CURV[index]*127.0);
      PIX_W[index] = (int)(255.0 - (omega*255.0)/poids_max);
      index++;
    }
  free(VAL);free(GXX);free(GXY);free(GYY);free(TEMP);
  free(CURV);free(GX);free(GY);
}

// Calculate the GHT
// Order 0: value
void Hough_General_ordre0(deplacement_s** RTABLE,double *HOUGH,Image<int>& p_img,float sigma,int nb_labels,double *vote_max) {
  int* PIX=p_img.PI();
  int w = p_img.PW();
  int h = p_img.PL();
  double *VAL,*TEMP;
  int i,j,index;
  double valmin,valmax;
  int index_val,x_vote,y_vote,index_vote;
  deplacement_s* p;

  (*vote_max) = 0;
  printf("Detection (order 0) : sigma = %f, nb_labels = %d\n",sigma,nb_labels);
  VAL = (double *)calloc(w*h,sizeof(double));
  TEMP = (double *)calloc(w*h,sizeof(double));
  // Conversion in double...
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      VAL[index] = (double) (PIX[index++]);
  // Gaussian Smoothing
  Gauss_Rec(VAL,TEMP,w,h,0,0,sigma,valmin,valmax);
  Gauss_Rec(TEMP,VAL,w,h,0,1,sigma,valmin,valmax);
  index=0;
  for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
       index_val = (int)((VAL[index]*(nb_labels-1))/(255.0));
	//Looking up the R-Table at row "index_val"
	p = RTABLE[index_val];
	while (p->next != NULL) {
	  x_vote = i+p->dx;
	  y_vote = j+p->dy;
	  if ((x_vote>0)&&(x_vote<w)&&(y_vote>0)&&(y_vote<h)) {
	    index_vote = y_vote*w + x_vote;
	    HOUGH[index_vote] += (p->poids);
	    if ((*vote_max)<HOUGH[index_vote]) (*vote_max)=HOUGH[index_vote];
	  }
	  p = p->next;
      }
      index++;
    }
  free(VAL);free(TEMP);
}

// Calculate the GHT
// At order 1: gradient direction
void Hough_General_ordre1(deplacement_s** RTABLE,double *HOUGH,Image<int>& p_img,float sigma,int nb_labels,double *vote_max) {
  int* PIX=p_img.PI();
  int w = p_img.PW();
  int h = p_img.PL();
  double *VAL,*GX,*GY,*TEMP;
  int i,j,index;
  double valmin,valmax;
  double argument,norme;
  int index_arg,x_vote,y_vote,index_vote;
  deplacement_s* p;

  float PI = 3.1415926536;
  
  (*vote_max) = 0;
  printf("Detection (order 1): sigma = %f, nb_labels = %d\n",sigma,nb_labels);
  VAL = (double *)calloc(w*h,sizeof(double));
  GX = (double *)calloc(w*h,sizeof(double));
  GY = (double *)calloc(w*h,sizeof(double));
  TEMP = (double *)calloc(w*h,sizeof(double));
  // Conversion in double...
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      VAL[index] = (double) (PIX[index++]);
  // Calculate the gradient components
  Gauss_Rec(VAL,TEMP,w,h,0,1,sigma,valmin,valmax);
  Gauss_Rec(TEMP,GX,w,h,1,0,sigma,valmin,valmax);
  Gauss_Rec(VAL,TEMP,w,h,0,0,sigma,valmin,valmax);
  Gauss_Rec(TEMP,GY,w,h,1,1,sigma,valmin,valmax);
  index=0;
  for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      if ((GX[index]!=0)||(GY[index]!=0)) {
	argument = atan2(GY[index],GX[index]);
	index_arg = (int)(((argument+PI)*nb_labels)/(2*PI));
	norme = sqrt(GX[index]*GX[index]+GY[index]*GY[index]);
	//Looking up the R-Table at row "index_arg"
	p = RTABLE[index_arg];
	while (p->next != NULL) {
	  x_vote = i+p->dx;
	  y_vote = j+p->dy;
	  if ((x_vote>0)&&(x_vote<w)&&(y_vote>0)&&(y_vote<h)) {
	    index_vote = y_vote*w + x_vote;
	    // The vote can be weighted according to the current pixel
	    // or to the weight of the corresponding prototype, or a 
	    // combination of them, or even a constant weight...
	    //HOUGH[index_vote] += (p->poids)*norme;
	    //HOUGH[index_vote] += 1;
	    HOUGH[index_vote] += (p->poids);
	    if ((*vote_max)<HOUGH[index_vote]) (*vote_max)=HOUGH[index_vote];
	  }
	  p = p->next;
	}
      }
      index++;
    }
  free(VAL);free(GX);free(GY);free(TEMP);
}

// Calculate the GHT
// At order 2: isophote curvature
void Hough_General_ordre2(deplacement_s** RTABLE,double *HOUGH,Image<int>& p_img,float sigma,int nb_labels,double *vote_max) {
  int* PIX=p_img.PI();
  int w = p_img.PW();
  int h = p_img.PL();
  double *VAL,*GX,*GY,*TEMP,*CURV;
  double *GXX,*GXY,*GYY;
  int i,j,index;
  double valmin,valmax,norme;
  int index_curv,x_vote,y_vote,index_vote;
  deplacement_s* p;

  (*vote_max) = 0;
  printf("Detection (order 2): sigma = %f, nb_labels = %d\n",sigma,nb_labels);
  VAL = (double *)calloc(w*h,sizeof(double));
  GX = (double *)calloc(w*h,sizeof(double));
  GY = (double *)calloc(w*h,sizeof(double));
  TEMP = (double *)calloc(w*h,sizeof(double));
  GXX = (double *)calloc(w*h,sizeof(double));
  GXY = (double *)calloc(w*h,sizeof(double));
  GYY = (double *)calloc(w*h,sizeof(double));
  CURV = (double *)calloc(w*h,sizeof(double));
  // Conversion to double...
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      VAL[index] = (double) (PIX[index++]);
  // Calculate the gradient components
  Gauss_Rec(VAL,TEMP,w,h,0,1,sigma,valmin,valmax);
  Gauss_Rec(TEMP,GX,w,h,1,0,sigma,valmin,valmax);
  Gauss_Rec(VAL,TEMP,w,h,0,0,sigma,valmin,valmax);
  Gauss_Rec(TEMP,GY,w,h,1,1,sigma,valmin,valmax);
  // Calculate G_tt/Gv (i.e. isophote curvature)
  // Calculate G_xx
  Gauss_Rec(VAL,TEMP,w,h,0,1,sigma,valmin,valmax);
  Gauss_Rec(TEMP,GXX,w,h,2,0,sigma,valmin,valmax);
  // then G_yy
  Gauss_Rec(VAL,TEMP,w,h,0,0,sigma,valmin,valmax);
  Gauss_Rec(TEMP,GYY,w,h,2,1,sigma,valmin,valmax);
  // and then G_xy
  Gauss_Rec(VAL,TEMP,w,h,1,0,sigma,valmin,valmax);
  Gauss_Rec(TEMP,GXY,w,h,1,1,sigma,valmin,valmax);
  // Calculate G_tt/G_g, isophote curvature, saved in CURV, 
  // and the interest function, equal to the total curvature
  // (Frobenius norm of the Hessian matrix), saved in VAL
  for (index=0;index<w*h;index++) {
    norme = sqrt(GX[index]*GX[index] + GY[index]*GY[index]);
    if (norme > 0) {
      CURV[index] = -GY[index]*GY[index]*GXX[index]
	+ 2*GX[index]*GY[index]*GXY[index]
	- GX[index]*GX[index]*GYY[index];
      VAL[index] = sqrt(GXX[index]*GXX[index] + 2*GXY[index]*GXY[index] + GYY[index]*GYY[index]);
      CURV[index] /= pow(norme,3);
      // The significant curvature is restricted within [-1,+1]
      if (CURV[index] > 1.0)  CURV[index] = 1.0;
      if (CURV[index] < -1.0)  CURV[index] = -1.0;
    } else CURV[index] = 0;
  }
  // GHT itself!
  index=0;
  for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      index_curv = (int)((nb_labels-1)*(CURV[index]+1)/2);
	//Looking up the R-Table at row "index_curv"
      p = RTABLE[index_curv];
      while (p->next != NULL) {
	x_vote = i+p->dx;
	y_vote = j+p->dy;
	if ((x_vote>0)&&(x_vote<w)&&(y_vote>0)&&(y_vote<h)) {
	  index_vote = y_vote*w + x_vote;
	 // The vote can be weighted according to the current pixel
	    // or to the weight of the corresponding prototype, or a 
	    // combination of them, or even a constant weight...
	  //HOUGH[index_vote] += (p->poids)*VAL[index];
	  //HOUGH[index_vote] += 1;
	  HOUGH[index_vote] += (p->poids);
	  if ((*vote_max)<HOUGH[index_vote]) (*vote_max)=HOUGH[index_vote];
	}
	p = p->next;
      }
      index++;
    }
  free(VAL);free(GX);free(GY);free(TEMP);
  free(GXX);free(GXY);free(GYY);free(CURV);
}

// Searching the best detection in the THG
void Find_Best_GHT(double *HOUGH,int width_img,int height_img,int nbest,int *Best_Pos,int taille_x,int taille_y) {
  int i_best,j_best;
  int i,j,k,index;
  double valmax;
  
  int rayon_exclusion_dx = taille_x/4;
  int rayon_exclusion_dy = taille_y/4;
  
  Best_Pos[0] = nbest;

  for (k = 0; k < nbest; k++) {
    // Searching the maximal value in the Hough space
    index = 0;
    valmax = 0.0;
    for (j=0;j<height_img;j++)
      for (i=0;i<width_img;i++) {
        if (HOUGH[index]>valmax) {
	  valmax = HOUGH[index];
	  i_best = i;j_best = j;
        }
        index++;
      }
    // Resetting the values around the last max
    for (i=max(i_best-rayon_exclusion_dx,0);i<=min(i_best+rayon_exclusion_dx,width_img-1);i++)
      for (j=max(j_best-rayon_exclusion_dy,0);j<=min(j_best+rayon_exclusion_dy,height_img-1);j++)
        HOUGH[j*width_img + i] = 0;
    Best_Pos[2*k + 1] = i_best;
    Best_Pos[2*k + 2] = j_best;
  }
}

void Update_Proto_RTable(deplacement_s** RTABLE,Image<int>& p_proto,Image<int>& p_index,Image<int>& p_poids,
                         Image<int>& p_courante,int *BestPos,int old_position_x,int old_position_y,
                         int *new_position_x,int *new_position_y,int nb_labels,float sigma) {

  int Iw = p_courante.PW();
  int Ih = p_courante.PL();
  int Pw = p_proto.PW();
  int Ph = p_proto.PL();
  int *PIX_P = p_proto.PI();
  int dist_min = Iw*Iw + Ih*Ih;
  int diff_x,diff_y,sum_sq;
  int i, j, k;
  
  for (k = 0; k < BestPos[0]; k++) {//We search the "best" position that is closest from the previous one
      diff_x = old_position_x - BestPos[2*k+1];
      diff_y = old_position_y - BestPos[2*k+2];
      sum_sq = (diff_x*diff_x) + (diff_y*diff_y);
      if (sum_sq < dist_min) {
 	(*new_position_x) = BestPos[2*k+1];
        (*new_position_y) = BestPos[2*k+2];
        dist_min = sum_sq;
      }
  }

  int index = 0; // We simply replace the current proto by the new one
  for (j = 0;j < Ph; j++)
    for (i = 0;i < Pw; i++) 
       PIX_P[index++] = p_courante.X(*new_position_x - Pw/2 + i, 
                                     *new_position_y - Ph/2 + j);
  Create_RTable_ordre1(RTABLE,p_proto,p_index,p_poids,sigma,nb_labels);
}
