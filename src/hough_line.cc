#include "Image.hh"
#include "feature.hh"
#include "filtre.hh"
#include "morfo.hh"
#include "hough_line.hh"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Tracé d'une droite dans l'image
// (Algorithme de Bresenham)
void Trace_droite(Image<int>& p,float angle,int orig_x,int orig_y,int colour) {

  int* PIX=p.PI();
  int WIDTH = p.PW();
  int HEIGHT = p.PL();
  int delta_x,delta_y;
  float pente,pente_inverse;
  int index,index_new;
  
  pente = tan(angle);
  if (pente) pente_inverse = 1/pente;
  //printf("Pente = %f\n",pente);
  
  // Point de départ
  index = orig_y*WIDTH + orig_x;
  PIX[index]=500+colour;
  
  delta_x = 0; delta_y = 0;
  if (abs(pente)<1) { // cas 1, 1 pixel par colonne
    while ((orig_x+delta_x<WIDTH)&&(orig_y+delta_y<HEIGHT)&&(orig_y+delta_y>=0)) {
      delta_x += 1;
      index_new = index + delta_x + delta_y*WIDTH;
      PIX[index_new]=500+colour;
      if ((float)(abs(delta_y))/(float)(abs(delta_x)) < abs(pente))
	if (pente > 0) delta_y += 1; else delta_y -= 1;
      //printf("Delta_x = %d, Delta_y = %d\n",delta_x,delta_y);
    }
  } else { // cas 2, 1 pixel par ligne
    while ((orig_x+delta_x<WIDTH)&&(orig_y+delta_y<HEIGHT)&&(orig_y+delta_y>=0)) {
      if (pente_inverse > 0) delta_y += 1; else delta_y -= 1;
      index_new = index + delta_x + delta_y*WIDTH;
      PIX[index_new]=500+colour;
      if ((float)(abs(delta_x))/(float)(abs(delta_y)) < abs(pente_inverse))
	delta_x += 1;
      //printf("Delta_x = %d, Delta_y = %d\n",delta_x,delta_y);
    }
  }
}

// Calcul de la transformée de Hough directement à partir 
// du gradient : version 1 / Paramétrage Theta x Rho
// -PI < Theta < PI ;  0 < Rho < sqrt(w*w+h*h)
void Hough_gradient_Rho_Theta_old(Image<int>& p_img,Image<int>& p_transf,float sigma_init,int nb_ech) {
  int* PIX_T=p_transf.PI();
  int* PIX=p_img.PI();
  int w = p_img.PW();
  int h = p_img.PL();
  int T_theta = p_transf.PW();
  int T_rho = p_transf.PL();
  double *VAL,*TEMP,*GX,*GY,*GXGY;
  double *HOUGH;
  double valmin=10000;
  double valmax=-10000;
  int i,j,k,index;
  
  float PI = 3.1415926536;
  float sigma = sigma_init;
  double theta,rho,alpha,norme,cx,cy;
  double poids_rho,poids_theta,pos_theta,vote_total;
  int index_theta,index_rho,defined;

  VAL = (double *)calloc(w*h,sizeof(double));
  GX = (double *)calloc(w*h,sizeof(double));
  GY = (double *)calloc(w*h,sizeof(double));
  TEMP = (double *)calloc(w*h,sizeof(double));
  // Passage en double...
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      VAL[index] = (double) (PIX[index++]);
  // Allocation du tableau de Hough
  HOUGH = (double *)calloc(T_theta*T_rho,sizeof(double));
  
  for (k=0;k<nb_ech;k++) {
    // Calcul des composantes du gradient
    Gauss_Rec(VAL,TEMP,w,h,0,1,sigma,valmin,valmax);
    Gauss_Rec(TEMP,GX,w,h,1,0,sigma,valmin,valmax);
    Gauss_Rec(VAL,TEMP,w,h,0,0,sigma,valmin,valmax);
    Gauss_Rec(TEMP,GY,w,h,1,1,sigma,valmin,valmax);
    
    // Calcul de la transformée de Hough
    //printf("Coucou\n");
    index = 0;
    for (j=0;j<h;j++) {
      for (i=0;i<w;i++) {
	defined=0;
	// Calcul du gradient conjugué (isophote)
	if (GY[index]>0) {cx = GY[index];cy = -GX[index];}
	else {cx = -GY[index];cy = GX[index];}
	// Calcul de Rho : distance de l'origine à la droite
	// D d'équation cy*x - cx*y - i*cy + j*cx = 0
	norme = sqrt(cx*cx + cy*cy);
	if (norme) {
	  rho = fabs(cx*j - cy*i)/norme;
	  index_rho = (int)(trunc(rho));
	  poids_rho = 1 - rho + index_rho;
	  defined = 1;
	}
	// Calcul de l'angle de la tangente au contour
	if (cx) alpha = atanf(cy/cx); else alpha = PI/2;
	// Calcul de theta, l'angle entre Ox et la droite 
	// passant par O et perpendiculaire à la droite D
	if (alpha < 0) theta = PI/2 + alpha;
	else {
	  // 1er cas : D est en dessous de O
	  // -PI/2 < Theta < 0 ; cx*j - cy*i < 0
	  if (cx*j - cy*i < 0) theta = alpha - PI/2;
	  // 2eme cas : D est au dessus de O
	  // PI/2 < Theta < PI ; cx*j - cy*i > 0
	  else theta = alpha + PI/2;
	}
	pos_theta = ((theta + PI)*(T_theta-1))/(2*PI);
	index_theta = (int)(trunc(pos_theta));
	poids_theta = 1 - pos_theta + index_theta;
	// Le poids du vote est égal au carré du module du gradient
	// Le vote est interpolé sur les 4 cellules adjacentes
	if (defined) {
	  //vote_total = norme;//sans normalisation par l'échelle
	  vote_total = norme*sigma;//avec normalisation par l'échelle
	  
	  HOUGH[index_rho*T_theta + index_theta] += vote_total*poids_rho*poids_theta;
	  if (poids_rho < 1) 
	    HOUGH[(index_rho+1)*T_theta + index_theta] += vote_total*(1-poids_rho)*poids_theta;
	  if (poids_theta < 1) 
	    HOUGH[index_rho*T_theta + index_theta+1] += vote_total*poids_rho*(1-poids_theta);
	  if ((poids_rho < 1)&&(poids_theta<1)) 
	    HOUGH[(index_rho+1)*T_theta + index_theta+1] += vote_total*(1-poids_rho)*(1-poids_theta);
	}
	index++;
      }
    }
    sigma *= 2;
  }
  // Normalisation pour affichage
  valmax = 0;
  for (index=0;index<T_theta*T_rho;index++) 
  if (HOUGH[index]>valmax) valmax = HOUGH[index];
  for (index=0;index<T_theta*T_rho;index++) 
    PIX_T[index] = (int)((255.0*HOUGH[index])/valmax);
  free(VAL);
  free(GX);
  free(GY);
  free(HOUGH);
  free(TEMP);
}

// Calcul de la transformée de Hough directement à partir 
// du gradient : version 1 / Paramétrage Theta x Rho
// -PI < Theta < PI ;  0 < Rho < sqrt(w*w+h*h)
void Hough_gradient_Rho_Theta(Image<int>& p_img,Image<int>& p_transf,float sigma_init,int nb_ech) {
  int* PIX_T=p_transf.PI();
  int* PIX=p_img.PI();
  int w = p_img.PW();
  int h = p_img.PL();
  int T_theta = p_transf.PW();
  int T_rho = p_transf.PL();
  double *VAL,*TEMP,*GX,*GY;
  double *HOUGH;
  double valmin=10000;
  double valmax=-10000;
  int i,j,k,index;
  
  float PI = 3.1415926536;
  float sigma = sigma_init;
  double theta,rho,alpha,norme,dist_sign;
  double poids_rho,poids_theta,pos_theta,vote_total;
  int index_theta,index_rho,defined;

  VAL = (double *)calloc(w*h,sizeof(double));
  GX = (double *)calloc(w*h,sizeof(double));
  GY = (double *)calloc(w*h,sizeof(double));
  TEMP = (double *)calloc(w*h,sizeof(double));
  // Passage en double...
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      VAL[index] = (double) (PIX[index++]);
  // Allocation du tableau de Hough
  HOUGH = (double *)calloc(T_theta*T_rho,sizeof(double));
  
  for (k=0;k<nb_ech;k++) {
    // Calcul des composantes du gradient
    Gauss_Rec(VAL,TEMP,w,h,0,1,sigma,valmin,valmax);
    Gauss_Rec(TEMP,GX,w,h,1,0,sigma,valmin,valmax);
    Gauss_Rec(VAL,TEMP,w,h,0,0,sigma,valmin,valmax);
    Gauss_Rec(TEMP,GY,w,h,1,1,sigma,valmin,valmax);
    
    // Calcul de la transformée de Hough
    //printf("Coucou\n");
    index = 0;
    for (j=0;j<h;j++) {
      for (i=0;i<w;i++) {
	defined=0;
	// Calcul de Rho : distance de l'origine à la droite
	// D d'équation Ix*x + Iy*y - i*Ix - j*Iy = 0
	norme = sqrt(GX[index]*GX[index]+GY[index]*GY[index]);
	dist_sign = i*GX[index] + j*GY[index];
	if (norme) {
	  rho = fabs(dist_sign)/norme;
	  index_rho = (int)(trunc(rho));
	  poids_rho = 1 - rho + index_rho;
	  defined = 1;
	}
	// Calcul de l'argument du gradient
	if (GX[index]) 
	  alpha = atan(GY[index]/GX[index]);
	else 
	  alpha = PI/2;
	// Calcul de theta, l'angle entre Ox et la droite 
	// passant par O et perpendiculaire à la droite D
	if ((alpha < 0)&&(GY[index]*dist_sign > 0))
	  theta = PI + alpha;
	else 
	  theta = alpha;
	pos_theta = ((theta + PI)*(T_theta-1))/(2*PI);
	index_theta = (int)(trunc(pos_theta));
	poids_theta = 1 - pos_theta + index_theta;
	// Le poids du vote est égal au carré du module du gradient
	// Le vote est interpolé sur les 4 cellules adjacentes
	if (defined) {
	  //vote_total = norme;//sans normalisation par l'échelle
	  vote_total = norme*sigma;//avec normalisation par l'échelle
	  
	  HOUGH[index_rho*T_theta + index_theta] += vote_total*poids_rho*poids_theta;
	  if (poids_rho < 1) 
	    HOUGH[(index_rho+1)*T_theta + index_theta] += vote_total*(1-poids_rho)*poids_theta;
	  if (poids_theta < 1) 
	    HOUGH[index_rho*T_theta + index_theta+1] += vote_total*poids_rho*(1-poids_theta);
	  if ((poids_rho < 1)&&(poids_theta<1)) 
	    HOUGH[(index_rho+1)*T_theta + index_theta+1] += vote_total*(1-poids_rho)*(1-poids_theta);
	}
	index++;
      }
    }
    sigma *= 2;
  }
  // Normalisation pour affichage
  valmax = 0;
  for (index=0;index<T_theta*T_rho;index++) 
  if (HOUGH[index]>valmax) valmax = HOUGH[index];
  for (index=0;index<T_theta*T_rho;index++) 
    PIX_T[index] = (int)((255.0*HOUGH[index])/valmax);
  free(VAL);
  free(GX);
  free(GY);
  free(HOUGH);
  free(TEMP);
}


// Calcul de la transformée de Hough directement à partir 
// du gradient : version 2 / Paramétrage Alpha x Coordonnée à l'Origine
// -PI/2 < Alpha < PI/2 ;  0 < Coord_orig < w+h
void Hough_gradient_Alpha_Coord(Image<int>& p_img,Image<int>& p_transf,float sigma_init,int nb_ech) {
  int* PIX_T=p_transf.PI();
  int* PIX=p_img.PI();
  int w = p_img.PW();
  int h = p_img.PL();
  int T_angle = p_transf.PW();
  int T_coord = p_transf.PL();
  double *VAL,*TEMP,*GX,*GY,*GXGY;
  double *HOUGH;
  double valmin=10000;
  double valmax=-10000;
  int i,j,k,index;
  
  float PI = 3.1415926536;
  float sigma = sigma_init;
  double alpha,cx,cy,kx,ky;
  int index_angle,index_coord,defined;

  VAL = (double *)calloc(w*h,sizeof(double));
  GX = (double *)calloc(w*h,sizeof(double));
  GY = (double *)calloc(w*h,sizeof(double));
  TEMP = (double *)calloc(w*h,sizeof(double));
  // Passage en double...
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      VAL[index] = (double) (PIX[index++]);
  // Allocation du tableau de Hough
  HOUGH = (double *)calloc(T_angle*T_coord,sizeof(double));
  for (k=0;k<nb_ech;k++) {
    // Calcul des composantes du gradient
    Gauss_Rec(VAL,TEMP,w,h,0,1,sigma,valmin,valmax);
    Gauss_Rec(TEMP,GX,w,h,1,0,sigma,valmin,valmax);
    Gauss_Rec(VAL,TEMP,w,h,0,0,sigma,valmin,valmax);
    Gauss_Rec(TEMP,GY,w,h,1,1,sigma,valmin,valmax);
    // Calcul de la transformée de Hough
    //printf("Coucou\n");
    index = 0;
    for (j=0;j<h;j++) {
      for (i=0;i<w;i++) {
	defined=0;
	// Calcul du gradient conjugué (isophote)
	if (GY[index]>0) {cx = GY[index];cy = -GX[index];}
	else {cx = -GY[index];cy = GX[index];}
	// Calcul de l'angle de la tangente au contour
	if (cx) alpha = atanf(cy/cx); else alpha = PI/2;
	index_angle = (int)(((alpha + PI/2)*(T_angle-1))/PI);
	// 1er cas : Pente positive
	if (cy > 0) {
	  // Calcul de l'ordonnée à l'origine
	  if (cx) ky = (j*cx - i*cy)/cx; else ky = -1000;
	  // Vote de Hough dans le quadrant 4
	  if ((ky>=0)&&(ky<h)) {
	    index_coord = w + (int)(ky);
	    defined = 1;
	  } else {
	    // Calcul de l'abscisse à l'origine
	    kx = (i*cy - j*cx)/cy;
	    // Vote de Hough dans le quadrant 3
	    if ((kx>=0)&&(kx<w)) {
	      index_coord = (int)(w-kx);
	      defined = 1;
	    }
	  }
	} else {
	  // 2d cas : Pente négative
	  // Calcul de l'ordonnée à l'origine
	  if (cx) ky = (j*cx - i*cy)/cx; else ky = -1000;
	  // Vote de Hough dans le quadrant 2
	  if ((ky>=0)&&(ky<h)) {
	    index_coord = w + (int)(ky);
	    defined = 1;
	  } else {
	    // Calcul de l'abscisse finale
	    if (cy) kx = (cx*(h-1-j) + i*cy)/cy;
	    // Vote de Hough dans le quadrant 1
	    if ((cy)&&(kx>=0)&&(kx<w)) {
	      index_coord = (int)(w-kx);
	      defined = 1;
	    }
	  }
	}
	// Le poids du vote est égal au module du gradient
	if (defined)
	  HOUGH[index_coord*T_angle + index_angle] += sqrt(cx*cx + cy*cy);
	index++;
      }
    }
    sigma *= 2;
  }
  // Normalisation pour affichage
  valmax = 0;
  for (index=0;index<T_angle*T_coord;index++) 
  if (HOUGH[index]>valmax) valmax = HOUGH[index];
  for (index=0;index<T_angle*T_coord;index++) 
    PIX_T[index] = (int)((255.0*HOUGH[index])/valmax);
  free(VAL);
  free(GX);
  free(GY);
  free(HOUGH);
  free(TEMP);
}

// Trace dans l'espace des paramètres (theta,rho) l'image du point (x,y)
// i.e. la sinusoïde d'équation rho = x cos theta + y sin theta
void Trace_sinusoide(Image<int>& p_transf,int x,int y) {
  int* PIX_T=p_transf.PI();
  int T_theta = p_transf.PW();
  int T_rho = p_transf.PL();
  int i;
  double theta,rho;
  
  float PI = 3.1415926536;

  for (i = 0;i < T_theta;i++) {
    // Calcul de l'angle theta entre -PI et PI
    theta = (2*PI*i)/T_theta - PI;
    // Calcul de Rho
    rho = x*cos(theta) + y*sin(theta);
    if ((rho>=0)&&(rho<T_rho)) PIX_T[((int)(rho))*T_theta + i] = 255;
  }
}

// Trace dans l'espace des paramètres (alpha,coord) 
// l'image du point (x,y)
void Trace_courbeAC(Image<int>& p_img,Image<int>& p_transf,int x,int y) {
  int* PIX_T=p_transf.PI();
  int w = p_img.PW();
  int h = p_img.PL();
  int T_angle = p_transf.PW();
  int T_coord = p_transf.PL();
  int i,index_coord;
  double alpha,cx,cy,kx,ky;
  
  float PI = 3.1415926536;
  
  //printf("x= %d ; y = %d;\n",x,y);
  for (i = 0;i < T_angle;i++) {
    // Calcul de l'angle alpha entre -PI/2 et PI/2
    alpha = (PI*i)/(T_angle-1) - PI/2;
    cx = cos(alpha);
    cy = sin(alpha);
    // 1er cas : Pente positive
    if (cy > 0) {
      // Calcul de l'ordonnée à l'origine
      if (cx) ky = (y*cx - x*cy)/cx; else ky = -1000;
      if ((ky>=0)&&(ky<h)) index_coord = w + (int)(ky);
      else {
	// Calcul de l'abscisse à l'origine
	kx = (x*cy - y*cx)/cy;
	if ((kx>=0)&&(kx<w)) index_coord = (int)(w-kx);
      }
    } else {
      // 2d cas : Pente négative
      // Calcul de l'ordonnée à l'origine
      if (cx) ky = (y*cx - x*cy)/cx; else ky = -1000;
      if ((ky>=0)&&(ky<h)) index_coord = w + (int)(ky);
      else {
	// Calcul de l'abscisse finale
	if (cy) kx = (cx*(h-1-y) + x*cy)/cy;
	if ((cy)&&(kx>=0)&&(kx<w)) index_coord = (int)(w-kx);
      }
    }
    PIX_T[index_coord*T_angle+i]=255;
  }
}

void Hough_lines_contours(Image<int>& p_img,Image<int>& p_transf) {
  int* PIX=p_img.PI();
  int w = p_img.PW();
  int h = p_img.PL();
  int T_theta = p_transf.PW();
  int T_rho = p_transf.PL();
  int* PIX_T=p_transf.PI();
  float vote_max = 0;
  float *HOUGH;
  int index,i,j,index_hough,t;
  float theta,rho;
  
  float PI = 3.1415926536;
  // Allocation du tableau de Hough
  HOUGH = (float *)calloc(T_theta*T_rho,sizeof(float));
  // Procédure de vote
  index = 0;
  for (j = 0;j < h;j++) {
    for (i = 0;i < w;i++) {
      // Tous les points noirs votent
      if (PIX[index]==0) {
	for (t = 0;t < T_theta;t++) {
	  // Calcul de l'angle theta entre -PI et PI
	  theta = (2*PI*t)/T_theta - PI;
	  // Calcul de Rho
	  rho = i*cos(theta) + j*sin(theta);
	  if ((rho>=0)&&(rho<T_rho)) {
	    index_hough = ((int)(rho))*T_theta + t;
	    HOUGH[index_hough] += 1;
	    if (HOUGH[index_hough]>vote_max) vote_max = HOUGH[index_hough];
	  }
	}
      }
      index++;
    }
  }
  // Repassage en entier et normalisation
  for (index_hough=0;index_hough<T_theta*T_rho;index_hough++)
    PIX_T[index_hough] = (int)((HOUGH[index_hough]*255.0)/vote_max);
  free(HOUGH);
}

// Recherche des meilleurs droites dans l'espace de Hough
void Trace_Best_Hough_Lines_old(Image<int>& p_transf,int *xbest,int *ybest) {
  int T_theta = p_transf.PW();
  int T_rho = p_transf.PL();
  int* PIX_T=p_transf.PI();
  int i,j,i_best,j_best;
  int valmax,index;
  
  int rayon_exclusion = 10;
  
  // Recherche de la valeur max dans l'espace de Hough
  index = 0;
  valmax = 0;
  for (j=0;j<T_rho;j++)
    for (i=0;i<T_theta;i++) {
      if (PIX_T[index]>valmax) {
	valmax = PIX_T[index];
	i_best = i;j_best = j;
      }
      index++;
    }
  // Elimination des valeurs autour du dernier max
  for (i=max(i_best-rayon_exclusion,0);i<=min(i_best+rayon_exclusion,T_theta-1);i++)
    for (j=max(j_best-rayon_exclusion,0);j<=min(j_best+rayon_exclusion,T_rho-1);j++)
      PIX_T[j*T_theta + i] = 0;
  (*xbest) = i_best;
  (*ybest) = j_best;
}

// Recherche des meilleurs droites dans l'espace de Hough
void Trace_Best_Hough_Lines(Image<int>& p_transf,int *xbest,int *ybest) {
  int T_theta = p_transf.PW();
  int T_rho = p_transf.PL();
  int* PIX_T=p_transf.PI();
  int i,j,i_best,j_best;
  int valmax,index;
  
  int rayon_exclusion_theta = 15;
  int rayon_exclusion_rho = 12;

  // Préliminaire : suppression des non-maxima locaux
  index = 0;
  for (j=0;j<T_rho;j++)
    for (i=0;i<T_theta;i++) {
      if ((PIX_T[index]<p_transf.X(i-1,j-1))||
	  (PIX_T[index]<p_transf.X(i-1,j))||
	  (PIX_T[index]<p_transf.X(i-1,j+1))||
	  (PIX_T[index]<p_transf.X(i,j-1))||
	  (PIX_T[index]<p_transf.X(i,j+1))||
	  (PIX_T[index]<p_transf.X(i+1,j-1))||
	  (PIX_T[index]<p_transf.X(i+1,j))||
	  (PIX_T[index]<p_transf.X(i+1,j+1)))
	PIX_T[index]=0;
      index++;
    }
  // Recherche de la valeur max dans l'espace de Hough
  index = 0;
  valmax = 0;
  for (j=0;j<T_rho;j++)
    for (i=0;i<T_theta;i++) {
      if (PIX_T[index]>valmax) {
	valmax = PIX_T[index];
	i_best = i;j_best = j;
      }
      index++;
    }
  // Elimination des valeurs autour du dernier max
  for (i=max(i_best-rayon_exclusion_theta,0);i<=min(i_best+rayon_exclusion_theta,T_theta-1);i++)
    for (j=max(j_best-rayon_exclusion_rho,0);j<=min(j_best+rayon_exclusion_rho,T_rho-1);j++)
      PIX_T[j*T_theta + i] = 0;
  (*xbest) = i_best;
  (*ybest) = j_best;
}
