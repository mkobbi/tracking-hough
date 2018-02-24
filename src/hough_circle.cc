#include "Image.hh"
#include "feature.hh"
#include "filtre.hh"
#include "morfo.hh"
#include "hough_circle.hh"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Calcul de la transformée de Hough pour les cercles
// directement à partir du gradient et de la courbure
// Le tableau HOUGH est de taille rayon_max*w*h
void Hough_circles(double *HOUGH,Image<int> p,float sigma_init,int nb_echelles,float gamma_integr,int rayon_max,double *vote_max) {
  int* PIX=p.PI();
  int w = p.PW();
  int h = p.PL();
  double *VAL,*TEMP,*GX,*GY,*GXGY;
  double *GXX,*GYY,*GXY,*CURV;
  double valmin=10000;
  double valmax=-10000;
  double square_module;
  int i,j,k,index;
  double rayon,xc,yc,dx,dy;
  int index_rayon,index_x,index_y,index_hough;
  float sigma_diff = sigma_init;

  *vote_max = 0.0;
  printf("Sigma_init = %f Nb Echelles = %d Gamma_liss = %f Rayon_max = %d\n",sigma_init,nb_echelles,gamma_integr,rayon_max);
  VAL = (double *)calloc(w*h,sizeof(double));
  GX = (double *)calloc(w*h,sizeof(double));
  GY = (double *)calloc(w*h,sizeof(double));
  GXX = (double *)calloc(w*h,sizeof(double));
  GXY = (double *)calloc(w*h,sizeof(double));
  GYY = (double *)calloc(w*h,sizeof(double));
  CURV = (double *)calloc(w*h,sizeof(double));
  TEMP = (double *)calloc(w*h,sizeof(double));
  // Passage en double...
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      VAL[index] = (double) (PIX[index++]);
  for (k=0;k<nb_echelles;k++) {
    // Calcul des composantes du gradient
    Gauss_Rec(VAL,TEMP,w,h,0,1,sigma_diff,valmin,valmax);
    Gauss_Rec(TEMP,GX,w,h,1,0,sigma_diff,valmin,valmax);
    Gauss_Rec(VAL,TEMP,w,h,0,0,sigma_diff,valmin,valmax);
    Gauss_Rec(TEMP,GY,w,h,1,1,sigma_diff,valmin,valmax);
    // Calcul de G_tt (i.e. courbure de l'isophote)
    // Calcul de G_xx
    Gauss_Rec(VAL,TEMP,w,h,0,1,sigma_diff,valmin,valmax);
    Gauss_Rec(TEMP,GXX,w,h,2,0,sigma_diff,valmin,valmax);
    // puis G_yy
    Gauss_Rec(VAL,TEMP,w,h,0,0,sigma_diff,valmin,valmax);
    Gauss_Rec(TEMP,GYY,w,h,2,1,sigma_diff,valmin,valmax);
    // et enfin G_xy
    Gauss_Rec(VAL,TEMP,w,h,1,0,sigma_diff,valmin,valmax);
    Gauss_Rec(TEMP,GXY,w,h,1,1,sigma_diff,valmin,valmax);
    // Calcul de G_tt, courbure de l'isophote,
    // multiplié par la norme du gradient au cube
    for (index=0;index<w*h;index++) {
      CURV[index] = GY[index]*GY[index]*GXX[index]
	- 2*GX[index]*GY[index]*GXY[index]
	+ GX[index]*GX[index]*GYY[index];
    }
    // Début du calcul de la transformée de Hough circulaire
    index=0;
    for (j=0;j<h;j++)
      for (i=0;i<w;i++) {
	square_module = GX[index]*GX[index] + GY[index]*GY[index];
	if (CURV[index] != 0) {
	  dx = -GX[index]*square_module/CURV[index];
	  dy = -GY[index]*square_module/CURV[index];
	} else {
	  dx = 0;dy = 0;
	}
	// Calcul du rayon du cercle osculateur
	rayon = sqrt(dx*dx + dy*dy);
	index_rayon = (int)(rayon - 1);
	if ((index_rayon>=0)&&(index_rayon<rayon_max)) {
	  index_x = (int)(i+dx);
	  index_y = (int)(j+dy);
	  // Vote de Hough (poids = norme de Frobenius de la Hessienne
	  if ((index_x>=0)&&(index_y>=0)&&(index_x<w)&&(index_y<h)) {
	    index_hough = w*h*index_rayon + w*index_y + index_x;
	    // avec normalisation de l'échelle : multiplication par sigma^2)
	    HOUGH[index_hough] += sigma_diff*sigma_diff*sqrt(GXX[index]*GXX[index] + 2*GXY[index]*GXY[index] + GYY[index]*GYY[index]);
	    //HOUGH[index_hough] += sigma_diff*sqrt(square_module);
	    //HOUGH[index_hough] += 1;
	    // sans normalisation de l'échelle
	    //HOUGH[index_hough] += sqrt(GXX[index]*GXX[index] + 2*GXY[index]*GXY[index] + GYY[index]*GYY[index]);
	    if (HOUGH[index_hough]> *vote_max) *vote_max = HOUGH[index_hough];
	  }
	}
	index++;
      }
    sigma_diff *= 2;
  }
  // Filtrage spatial de l'espace de Hough
  if (gamma_integr < 6.0)
    FiltreExp_Hough_3d(HOUGH,gamma_integr,rayon_max,w,h,vote_max);
  free(VAL);free(GXX);free(GXY);free(GYY);free(TEMP);
  free(CURV);free(GX);free(GY);
}

// Calcul de la transformée de Hough pour les cercles
// à partir de l'image binaire de contours (points noirs à 0)
// Le tableau HOUGH est de taille rayon_max*w*h
void Hough_circles_contours(double *HOUGH,Image<int> p,int rayon_max,double *vote_max) {
  int* PIX=p.PI();
  int w = p.PW();
  int h = p.PL();
  int i,j,inew,jnew;
  int *DIST;
  double *POIDS;
  double distance;
  int index,index_tab,index_hough1,index_hough2;

  *vote_max = 0.0;

  // On commence à calculer la matrice des distances et des poids
  DIST = (int *)calloc((2*rayon_max+1)*(2*rayon_max+1),sizeof(int));
  POIDS = (double *)calloc((2*rayon_max+1)*(2*rayon_max+1),sizeof(double));
  index = 0;
  for (i=-rayon_max;i<=rayon_max;i++)
    for (j=-rayon_max;j<=+rayon_max;j++) {
      distance = sqrt(i*i+j*j);
      if (distance <= rayon_max) {
	DIST[index] = (int)(trunc(distance));
	POIDS[index] = 1 - (distance - DIST[index]);
      } else {
	DIST[index] = -1000;
      }
      index++;
    } 
  // Procédure de vote
  index = 0;
  for (j = 0;j < h;j++) {
    for (i = 0;i < w;i++) {
      // Tous les points noirs votent sur tous les points du cône
      // le poids des votes est interpolés selon les valeurs de rayon 
      if (PIX[index]==0) {
	index_tab = 0;
	for (inew=i-rayon_max;inew<=i+rayon_max;inew++)
	  for (jnew=j-rayon_max;jnew<=j+rayon_max;jnew++) {
	    if ((inew>=0)&&(jnew>=0)&&(inew<w)&&(jnew<h)
		&&(DIST[index_tab]>0)) {
	      index_hough1 = (DIST[index_tab]-1)*w*h + jnew*w + inew;
              //sans normalisation par le rayon
              //HOUGH[index_hough1] += POIDS[index_tab];
              //avec normalisation par le rayon
	      HOUGH[index_hough1] += POIDS[index_tab]/(1.0+DIST[index_tab]);
	      if (HOUGH[index_hough1]> *vote_max) *vote_max = HOUGH[index_hough1];
	      
	      if (POIDS[index_tab] < 1) {
		index_hough2 = DIST[index_tab]*w*h + jnew*w + inew;
                //sans normalisation par le rayon
		//HOUGH[index_hough2] += 1-POIDS[index_tab];
	        //avec normalisation par le rayon
                HOUGH[index_hough2] += (1-POIDS[index_tab])/(1.0+DIST[index_tab]);
		if (HOUGH[index_hough2]> *vote_max) *vote_max = HOUGH[index_hough2];
	      }
	    }
	    index_tab++;
	  }
      }
      index++;
    }
    //if ((j%10)==0) printf("Ligne n°%d\n",j);
  }
  free(POIDS);
  free(DIST);
}

// Routine de conversion en plans d'images de la transformée de Hough circulaire
 void MaJ_PlanHough(Image<int>& p_transf_plan,int numero_plan,double *HOUGH,double vote_max) {
  int* PIX=p_transf_plan.PI();
  int w = p_transf_plan.PW();
  int h = p_transf_plan.PL();
  int index,index_hough;
  
  index_hough = w*h*(numero_plan-1);
  for (index = 0;index < w*h;index++) 
    PIX[index] = (int)(255.0*HOUGH[index_hough++]/vote_max);
}

// Routine de remise à zéro de la transformée de Hough circulaire
void Reset_Hough(double *HOUGH,int rayon_max,int width_img,int height_img) {
  int index;
  for (index = 0;index < rayon_max*width_img*height_img;index++) HOUGH[index]=0; 
} 

void Trace_cercle(Image<int>& p,int x_centre,int y_centre,int rayon, int couleur) {
  int* PIX=p.PI();
  int WIDTH = p.PW();
  int HEIGHT = p.PL();
  Image<int> q(WIDTH,HEIGHT);
  int* QPIX=q.PI();
  int i,j,index,index_point;

  float *DIST;

  DIST = (float *)calloc((2*rayon+1)*(2*rayon+1),sizeof(float));
  index = 0;
  for (i=-rayon;i<=rayon;i++)
    for (j=-rayon;j<=+rayon;j++) 
      DIST[index++] = sqrt(i*i+j*j);
  index = 0;
  for (i=x_centre-rayon;i<=x_centre+rayon;i++)
    for (j=y_centre-rayon;j<=y_centre+rayon;j++) {
      if ((i>=0)&&(j>=0)&&(i<WIDTH)&&(j<HEIGHT)
	  &&(DIST[index]> rayon-0.5)&&(DIST[index]<=rayon+0.5)) {
	index_point = j*WIDTH + i;
	PIX[index_point]=500+couleur;
      }
      index++;
    }
  free(DIST);
}
 
// Filtrage exponentiel récursive de l'espace de Hough
// pour les cercles (rayon_max*w*h)

void FiltreExp_Hough_3d(double *HOUGH,double gamma,int rayon_max,int w,int h,double *vote_max) { 
  int i,j,k;
  int index;
  double alpha;

  alpha = 1 - exp(-gamma);
  // Passage en profondeur et causal
  index=w*h;
  for (k=1;k<rayon_max;k++)
    for (j=0;j<h;j++) 
      for (i=0;i<w;i++) {
	HOUGH[index] = alpha*HOUGH[index]+(1-alpha)*HOUGH[index-w*h];
	index++;
      }
  // Passage en profondeur et anticausal
  index -= w*h;
  for (k=1;k<rayon_max;k++)
    for (j=0;j<h;j++) 
      for (i=0;i<w;i++) {
	HOUGH[index] = alpha*HOUGH[index]+(1-alpha)*HOUGH[index+w*h];
	index--;
      }
  // Passage horizontal, causal
  for (k=0;k<rayon_max;k++)
    for (j=0;j<h;j++) {
      index++;
      for (i=1;i<w;i++) {
	HOUGH[index] = alpha*HOUGH[index]+(1-alpha)*HOUGH[index-1];
	index++;
      }
    }
  // Passage horizontal, anticausal
  for (k=0;k<rayon_max;k++)
    for (j=0;j<h;j++) {
      index--;
      for (i=1;i<w;i++) {
	HOUGH[index] = alpha*HOUGH[index]+(1-alpha)*HOUGH[index+1];
	index--;
      }
    }
  // Passage vertical, causal
  for (k=0;k<rayon_max;k++) {
    index += w;
    for (j=1;j<h;j++) 
      for (i=0;i<w;i++) {
	HOUGH[index] = alpha*HOUGH[index]+(1-alpha)*HOUGH[index-w];
	index++;
      }
  }
  // Passage vertical, anticausal, et calcul du max
  for (k=0;k<rayon_max;k++) {
    index -= w;
    for (j=1;j<h;j++) 
      for (i=0;i<w;i++) {
	HOUGH[index] = alpha*HOUGH[index]+(1-alpha)*HOUGH[index+w];
	if (HOUGH[index]>(*vote_max)) (*vote_max) = HOUGH[index];
	index--;
      }
  }
  
}

// Recherche des meilleurs cercles dans l'espace de Hough
void Trace_Best_Hough_Circles_old(double *HOUGH,int rayon_max,int w,int h,int *xbest,int *ybest,int *rhobest) {
  int i,j,k,i_best,j_best,k_best,index;
  double valmax;
  
  int rayon_excl_xy = 8;
  int rayon_excl_rho = 2;
  
  int rayon_minimum = 8;
  
  // Recherche de la valeur max dans l'espace de Hough
  index = rayon_minimum*w*h;
  valmax = 0;
  for (k=rayon_minimum;k<rayon_max;k++)
    for (j=0;j<h;j++)
      for (i=0;i<w;i++) {
	if (HOUGH[index]>valmax) {
	  valmax = HOUGH[index];
	  i_best = i;
	  j_best = j;
	  k_best = k;
	}
	index++;
      }
  // Elimination des valeurs autour du dernier max
  for (k=max(k_best-rayon_excl_rho,0);k<=min(k_best+rayon_excl_rho,rayon_max-1);k++)
    for (j=max(j_best-rayon_excl_xy,0);j<=min(j_best+rayon_excl_xy,h-1);j++)
      for (i=max(i_best-rayon_excl_xy,0);i<=min(i_best+rayon_excl_xy,w-1);i++)
	HOUGH[k*w*h + j*w + i] = 0;
  (*xbest) = i_best;
  (*ybest) = j_best;
  (*rhobest) = k_best;
}

// Recherche des meilleurs cercles dans l'espace de Hough
void Trace_Best_Hough_Circles(double *HOUGH,int rayon_max,int w,int h,int *xbest,int *ybest,int *rhobest) {
  int i,j,k,i_best,j_best,k_best,index;
  double valmax;
  
  int rayon_minimum = 3;

  int rayon_excl_xy = 3;
  int rayon_excl_rho = 8;

  // Préliminaire : suppression des non-maxima locaux
  for (k=1;k<rayon_max-1;k++)
    for (j=1;j<h-1;j++)
      for (i=1;i<w-1;i++) {
	index = k*w*h + j*w + i;
	if ((HOUGH[index]<HOUGH[index-w*h-1])||
	    (HOUGH[index]<HOUGH[index-w*h+1])||
	    (HOUGH[index]<HOUGH[index-w*h-w])||
	    (HOUGH[index]<HOUGH[index-w*h+w])||
	    (HOUGH[index]<HOUGH[index-w*h-1-w])||
	    (HOUGH[index]<HOUGH[index-w*h-1+w])||
	    (HOUGH[index]<HOUGH[index-w*h+1-w])||
	    (HOUGH[index]<HOUGH[index-w*h+1+w])||
	    (HOUGH[index]<HOUGH[index-w*h])||
	    (HOUGH[index]<HOUGH[index+w*h-1])||
	    (HOUGH[index]<HOUGH[index+w*h+1])||
	    (HOUGH[index]<HOUGH[index+w*h-w])||
	    (HOUGH[index]<HOUGH[index+w*h+w])||
	    (HOUGH[index]<HOUGH[index+w*h-1-w])||
	    (HOUGH[index]<HOUGH[index+w*h-1+w])||
	    (HOUGH[index]<HOUGH[index+w*h+1-w])||
	    (HOUGH[index]<HOUGH[index+w*h+1+w])||
	    (HOUGH[index]<HOUGH[index+w*h])||
	    (HOUGH[index]<HOUGH[index-1])||
	    (HOUGH[index]<HOUGH[index+1])||
	    (HOUGH[index]<HOUGH[index-w])||
	    (HOUGH[index]<HOUGH[index+w])||
	    (HOUGH[index]<HOUGH[index-1-w])||
	    (HOUGH[index]<HOUGH[index-1+w])||
	    (HOUGH[index]<HOUGH[index+1-w])||
	    (HOUGH[index]<HOUGH[index+1+w]))
	  HOUGH[index]=0;
      } 
  // Recherche de la valeur max dans l'espace de Hough
  index = rayon_minimum*w*h;
  valmax = 0;
  for (k=rayon_minimum;k<rayon_max;k++)
    for (j=0;j<h;j++)
      for (i=0;i<w;i++) {
	if (HOUGH[index]>valmax) {
	  valmax = HOUGH[index];
	  i_best = i;
	  j_best = j;
	  k_best = k;
	}
	index++;
      }
  // Elimination des valeurs autour du dernier max
  for (k=max(k_best-rayon_excl_rho,0);k<=min(k_best+rayon_excl_rho,rayon_max-1);k++)
    for (j=max(j_best-rayon_excl_xy,0);j<=min(j_best+rayon_excl_xy,h-1);j++)
      for (i=max(i_best-rayon_excl_xy,0);i<=min(i_best+rayon_excl_xy,w-1);i++)
	HOUGH[k*w*h + j*w + i] = 0;
  (*xbest) = i_best;
  (*ybest) = j_best;
  (*rhobest) = k_best;
}

// Calcul de la transformée de Hough en cercle, en un seul plan
// On cumule tous les rayons, on ne s'intéresse qu'au centre
void HCircleCentre(Image<int>& p_img,Image<int>& p_transf,float sigma_init,int nb_ech,int rayon_min) {
  int* PIX_T=p_transf.PI();
  int* PIX=p_img.PI();
  int w = p_img.PW();
  int h = p_img.PL();
  double *VAL,*TEMP,*GX,*GY,*GXGY;
  double *GXX,*GYY,*GXY,*CURV,*HOUGH;
  double valmin=10000;
  double valmax=-10000;
  double square_module;
  int i,j,k,index;
  double vote_max,dx,dy;
  double vote_total,poids_x,poids_y;
  int index_hough,index_x,index_y;
  float sigma_diff = sigma_init;
  
  printf("Sigma = %f Nb_Ech = %d Rayon_min = %d\n",sigma_init,nb_ech,rayon_min);
  vote_max = 0.0;
  VAL = (double *)calloc(w*h,sizeof(double));
  GX = (double *)calloc(w*h,sizeof(double));
  GY = (double *)calloc(w*h,sizeof(double));
  GXX = (double *)calloc(w*h,sizeof(double));
  GXY = (double *)calloc(w*h,sizeof(double));
  GYY = (double *)calloc(w*h,sizeof(double));
  CURV = (double *)calloc(w*h,sizeof(double));
  TEMP = (double *)calloc(w*h,sizeof(double));
  HOUGH = (double *)calloc(w*h,sizeof(double));
  // Passage en double...
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      VAL[index] = (double) (PIX[index++]);
  for (k=0;k<nb_ech;k++) {
    // Calcul des composantes du gradient
    Gauss_Rec(VAL,TEMP,w,h,0,1,sigma_diff,valmin,valmax);
    Gauss_Rec(TEMP,GX,w,h,1,0,sigma_diff,valmin,valmax);
    Gauss_Rec(VAL,TEMP,w,h,0,0,sigma_diff,valmin,valmax);
    Gauss_Rec(TEMP,GY,w,h,1,1,sigma_diff,valmin,valmax);
    // Calcul de G_tt (i.e. courbure de l'isophote)
    // Calcul de G_xx
    Gauss_Rec(VAL,TEMP,w,h,0,1,sigma_diff,valmin,valmax);
    Gauss_Rec(TEMP,GXX,w,h,2,0,sigma_diff,valmin,valmax);
    // puis G_yy
    Gauss_Rec(VAL,TEMP,w,h,0,0,sigma_diff,valmin,valmax);
    Gauss_Rec(TEMP,GYY,w,h,2,1,sigma_diff,valmin,valmax);
    // et enfin G_xy
    Gauss_Rec(VAL,TEMP,w,h,1,0,sigma_diff,valmin,valmax);
    Gauss_Rec(TEMP,GXY,w,h,1,1,sigma_diff,valmin,valmax);
    // Calcul de G_tt, courbure de l'isophote,
    // multiplié par la norme du gradient au cube
    for (index=0;index<w*h;index++) {
      CURV[index] = GY[index]*GY[index]*GXX[index]
	- 2*GX[index]*GY[index]*GXY[index]
	+ GX[index]*GX[index]*GYY[index];
    }
    // Début du calcul de la transformée de Hough partielle
    index=0;
    for (j=0;j<h;j++)
      for (i=0;i<w;i++) {
	square_module = GX[index]*GX[index] + GY[index]*GY[index];
	if (CURV[index] != 0) {
	  dx = -GX[index]*square_module/CURV[index];
	  dy = -GY[index]*square_module/CURV[index];
	  if ((i+dx > 0)&&(i+dx < w-1)&&
	      (j+dy > 0)&&(j+dy < h-1)&&
	      (dx*dx+dy*dy > rayon_min*rayon_min)) {
	    // Poids total du vote (norme de Frobenius de la Hessienne)
	    // normalisé par l'échelle
	    vote_total = sigma_diff*sigma_diff*sqrt(GXX[index]*GXX[index] + 2*GXY[index]*GXY[index] + GYY[index]*GYY[index]);
	    //vote_total = 1;
	    // Répartition du poids par interpolation bilinéaire 
	    // sur les 4 cases les plus proches
	    index_x = (int)(trunc(i+dx));
	    poids_x = 1 - (i+dx) + index_x;
	    index_y = (int)(trunc(j+dy));
	    poids_y = 1 - (j+dy) + index_y;
	    index_hough = index_y*w + index_x;
	    HOUGH[index_hough] += vote_total*poids_x*poids_y;
	    if (HOUGH[index_hough] > vote_max) vote_max = HOUGH[index_hough];
	    if (poids_x < 1) {
	      index_hough = (index_y+1)*w + index_x;
	      HOUGH[index_hough] += vote_total*(1-poids_x)*poids_y;
	      if (HOUGH[index_hough] > vote_max) vote_max = HOUGH[index_hough];
	    }
	    if (poids_y < 1) {
	      index_hough = index_y*w + index_x + 1;
	      HOUGH[index_hough] += vote_total*poids_x*(1-poids_y);
	      if (HOUGH[index_hough] > vote_max) vote_max = HOUGH[index_hough];
	    }
	    if ((poids_x < 1)&&(poids_y<1)) {
	      index_hough = (index_y+1)*w + index_x + 1;
	      HOUGH[index_hough] += vote_total*(1-poids_x)*(1-poids_y);
	      if (HOUGH[index_hough] > vote_max) vote_max = HOUGH[index_hough];
	    }
	  }
	}
	index++;
      }
    sigma_diff *= 2;
  }
  // Repassage en entier
  for (index=0;index<w*h;index++) 
    PIX_T[index] = (int)((HOUGH[index]*255.0)/vote_max);
  free(VAL);free(GXX);free(GXY);free(GYY);free(TEMP);
  free(CURV);free(GX);free(GY);free(HOUGH);

}

// Calcul des meilleurs cercles fondés sur la transformée de Hough en 2 temps
// on calcule d'abord la TH précedente (centre uniquement), puis on sélectionne
// les meilleurs centre et compte le nombre de points par rayon pour chaque centre
// La deuxième phase est calculée seulement à une échelle, la plus grande, mais
// pondéré selon la norme du gradient de la plus petite échelle...


void HCC_BestCircles(Image<int>& p_img,Image<int>& p_transf,float sigma_init,int nb_ech,int rayon_min,int rayon_max,int nb_cercles,float *BestCircles) {
  int* PIX_T=p_transf.PI();
  int* PIX=p_img.PI();
  int w = p_img.PW();
  int h = p_img.PL();
  double *VAL,*TEMP,*GX,*GY,*GXGY;
  double *GXX,*GYY,*GXY,*CURV,*HOUGH;
  double valmin=10000;
  double valmax=-10000;
  double square_module;
  int i,j,k,index;
  double vote_max;
  float dx,dy;
  double vote_total,poids_x,poids_y;
  int index_hough,index_x,index_y;
  float sigma_diff = sigma_init;
  Image<int> p_marqueur(w,h);
  int* PMARK = p_marqueur.PI();
  float **Best_coord;
  float **Best_histogramme;
  float rayon,best_score;
  int i_best,j_best,k_best,r_best;
  int rayon_exclusion_dx = 3;
  int rayon_exclusion_dy = 3;
  int rayon_exclusion_dr = 8;
  float sum_x,sum_y,sum_tot;
  
  printf("début procédure\n");
  printf("Paramètres : sigma = %f nb_ech = %d rayon_min = %d rayon_max = %d nb_cercles = %d\n",sigma_init,nb_ech,rayon_min,rayon_max,nb_cercles);
	 vote_max = 0.0;
  VAL = (double *)calloc(w*h,sizeof(double));
  GX = (double *)calloc(w*h,sizeof(double));
  GY = (double *)calloc(w*h,sizeof(double));
  GXX = (double *)calloc(w*h,sizeof(double));
  GXY = (double *)calloc(w*h,sizeof(double));
  GYY = (double *)calloc(w*h,sizeof(double));
  CURV = (double *)calloc(w*h,sizeof(double));
  TEMP = (double *)calloc(w*h,sizeof(double));
  HOUGH = (double *)calloc(w*h,sizeof(double));
  // Passage en double...
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      VAL[index] = (double) (PIX[index++]);
  for (k=0;k<nb_ech;k++) {
    // Calcul des composantes du gradient
    Gauss_Rec(VAL,TEMP,w,h,0,1,sigma_diff,valmin,valmax);
    Gauss_Rec(TEMP,GX,w,h,1,0,sigma_diff,valmin,valmax);
    Gauss_Rec(VAL,TEMP,w,h,0,0,sigma_diff,valmin,valmax);
    Gauss_Rec(TEMP,GY,w,h,1,1,sigma_diff,valmin,valmax);
    // Calcul de G_tt (i.e. courbure de l'isophote)
    // Calcul de G_xx
    Gauss_Rec(VAL,TEMP,w,h,0,1,sigma_diff,valmin,valmax);
    Gauss_Rec(TEMP,GXX,w,h,2,0,sigma_diff,valmin,valmax);
    // puis G_yy
    Gauss_Rec(VAL,TEMP,w,h,0,0,sigma_diff,valmin,valmax);
    Gauss_Rec(TEMP,GYY,w,h,2,1,sigma_diff,valmin,valmax);
    // et enfin G_xy
    Gauss_Rec(VAL,TEMP,w,h,1,0,sigma_diff,valmin,valmax);
    Gauss_Rec(TEMP,GXY,w,h,1,1,sigma_diff,valmin,valmax);
    // Calcul de G_tt, courbure de l'isophote,
    // multiplié par la norme du gradient au cube
    for (index=0;index<w*h;index++) {
      CURV[index] = GY[index]*GY[index]*GXX[index]
	- 2*GX[index]*GY[index]*GXY[index]
	+ GX[index]*GX[index]*GYY[index];
    }
    // Début du calcul de la transformée de Hough partielle
    index=0;
    for (j=0;j<h;j++)
      for (i=0;i<w;i++) {
	square_module = GX[index]*GX[index] + GY[index]*GY[index];
	if (CURV[index] != 0) {
	  dx = -GX[index]*square_module/CURV[index];
	  dy = -GY[index]*square_module/CURV[index];
	  if ((i+dx > 0)&&(i+dx < w-1)&&
	      (j+dy > 0)&&(j+dy < h-1)&&
	      (dx*dx+dy*dy > rayon_min*rayon_min)) {
	    // Poids total du vote (norme de Frobenius de la Hessienne)
	    // normalisé par l'échelle
	    vote_total = sigma_diff*sigma_diff*sqrt(GXX[index]*GXX[index] + 2*GXY[index]*GXY[index] + GYY[index]*GYY[index]);
	    //vote_total = 1;
	    // Répartition du poids par interpolation bilinéaire 
	    // sur les 4 cases les plus proches
	    index_x = (int)(trunc(i+dx));
	    poids_x = 1 - (i+dx) + index_x;
	    index_y = (int)(trunc(j+dy));
	    poids_y = 1 - (j+dy) + index_y;
	    index_hough = index_y*w + index_x;
	    HOUGH[index_hough] += vote_total*poids_x*poids_y;
	    if (HOUGH[index_hough] > vote_max) vote_max = HOUGH[index_hough];
	    if (poids_x < 1) {
	      index_hough = (index_y+1)*w + index_x;
	      HOUGH[index_hough] += vote_total*(1-poids_x)*poids_y;
	      if (HOUGH[index_hough] > vote_max) vote_max = HOUGH[index_hough];
	    }
	    if (poids_y < 1) {
	      index_hough = index_y*w + index_x + 1;
	      HOUGH[index_hough] += vote_total*poids_x*(1-poids_y);
	      if (HOUGH[index_hough] > vote_max) vote_max = HOUGH[index_hough];
	    }
	    if ((poids_x < 1)&&(poids_y<1)) {
	      index_hough = (index_y+1)*w + index_x + 1;
	      HOUGH[index_hough] += vote_total*(1-poids_x)*(1-poids_y);
	      if (HOUGH[index_hough] > vote_max) vote_max = HOUGH[index_hough];
	    }
	  }
	}
	index++;
      }
    sigma_diff *= 2;
  }
  // Lissage de la transformée pour améliorer la détection des maxima
  // Lissage réalisée à l'échelle minimum (sigma_init)
  Gauss_Rec(HOUGH,TEMP,w,h,0,1,sigma_init,valmin,valmax);
  Gauss_Rec(TEMP,HOUGH,w,h,0,0,sigma_init,valmin,valmax);
  
  // Repassage en entier pour affichage de la transformée
  for (index=0;index<w*h;index++) 
    PIX_T[index] = (int)((HOUGH[index]*255.0)/vote_max);
  printf("Fin 1ere phase\n");

  // Deuxième phase : calcul de N meilleurs centres
  Best_coord = (float**)calloc(nb_cercles,sizeof(float*));
  for (k=0;k<nb_cercles;k++) Best_coord[k] = (float*)calloc(2,sizeof(float));
  for (k=0;k<nb_cercles;k++) {
    // Recherche de la valeur max dans l'espace de Hough
    index = 0;
    vote_max = 0;
    for (j=0;j<h;j++)
      for (i=0;i<w;i++) {
	if (HOUGH[index]>vote_max) {
	  vote_max = HOUGH[index];
	  i_best = i;j_best = j;
	}
	index++;
      }
    // Elimination des valeurs autour du dernier max
    // En même temps on calcule le centre de gravité de la fonction de vote
    // dans son rayon d'exclusion, pour fournir une estimation sub-pixelique
    // du centre du cercle, et marquage des meilleurs centres sur p_marqueur
    sum_x = 0.0;
    sum_y = 0.0;
    sum_tot = 0.0;
    for (i=max(i_best-rayon_exclusion_dx,0);i<=min(i_best+rayon_exclusion_dx,w-1);i++)
      for (j=max(j_best-rayon_exclusion_dy,0);j<=min(j_best+rayon_exclusion_dy,h-1);j++) {
	index = j*w + i;
	sum_x += i*HOUGH[index];
	sum_y += j*HOUGH[index];
	sum_tot += HOUGH[index];
	HOUGH[index] = 0;
	PMARK[index] = 255;
      }
    Best_coord[k][0] = sum_x/sum_tot;
    Best_coord[k][1] = sum_y/sum_tot;
   
  }
  printf("Fin 2eme phase\n");
  for (k=0;k<nb_cercles;k++) printf("Meilleur cercle n.%d = (%f,%f)\n",k,Best_coord[k][0],Best_coord[k][1]);
  // Deuxième passe, on compte les points appartenant aux 
  // meilleurs cercles, pour chaque rayon
  // On réutilise GXX et GYY pour calculer le module du gradient
  // à la plus petite échelle (sigma_init)
  Gauss_Rec(VAL,TEMP,w,h,0,1,sigma_init,valmin,valmax);
  Gauss_Rec(TEMP,GXX,w,h,1,0,sigma_init,valmin,valmax);
  Gauss_Rec(VAL,TEMP,w,h,0,0,sigma_init,valmin,valmax);
  Gauss_Rec(TEMP,GYY,w,h,1,1,sigma_init,valmin,valmax);

  Best_histogramme = (float**)calloc(nb_cercles,sizeof(float*));
  for (k=0;k<nb_cercles;k++) Best_histogramme[k] = (float*)calloc(rayon_max+1,sizeof(float));
  index=0;
  for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      // calcul du module du gradient au carré à grande échelle
      square_module = GX[index]*GX[index] + GY[index]*GY[index];
      // calcul du module du gradient à petite échelle
      GXX[index] = sqrt(GXX[index]*GXX[index] + GYY[index]*GYY[index]);
      if (CURV[index] != 0) {
	dx = -GX[index]*square_module/CURV[index];
	dy = -GY[index]*square_module/CURV[index];
	rayon = sqrt(dx*dx + dy*dy);
	if ((rayon > rayon_min)&&(rayon < rayon_max)) {
	  //vote_total = 1.0;
	  //vote_total = sqrt(square_module);
	  vote_total = GXX[index];
	  for (k=0;k<nb_cercles;k++) {
	    if ((abs((i+dx)-Best_coord[k][0])<1)&&(abs((j+dy)-Best_coord[k][1])<1)) {
	      Best_histogramme[k][(int)(trunc(rayon))] += vote_total*(1 - rayon + trunc(rayon));
	      if (rayon != trunc(rayon))
		Best_histogramme[k][(int)(trunc(rayon) + 1)] += vote_total*(rayon - trunc(rayon));
	    }
	  }
	}
      }
      index++;
    }
  printf("Fin 3eme phase\n");
  // Dernière phase : on recalcule les meilleurs cercles à partir 
  // des histogrammes par centre et par rayon...
  
  for (k=0;k<nb_cercles;k++) {
    best_score = 0.0;
    for (i=0;i<nb_cercles;i++) {
      for (j=rayon_min;j<=rayon_max;j++) {
	// Le score est normalisé par le rayon pour ne pas 
	// trop désavantager les petits cercles
	// if (Best_histogramme[i][j] > best_score*j) {
// 	  best_score =  Best_histogramme[i][j]/((float)j);
	if (Best_histogramme[i][j] > best_score) {
	  best_score =  Best_histogramme[i][j];
	  r_best = j;
	  k_best = i;
	}
      }
    }
    // Elimination des valeurs autour du dernier max
    // En même temps on interpole la valeur du rayon autour de son rayon d'exclusion
    sum_x = 0.0;
    sum_tot = 0.0;
    for (j=max(r_best-rayon_exclusion_dr,0);j<=min(r_best+rayon_exclusion_dr,rayon_max);j++) {
      sum_x += j*Best_histogramme[k_best][j];
      sum_tot += Best_histogramme[k_best][j];
      Best_histogramme[k_best][j] = 0;
    }
    // On remplit le tableau de sortie
    BestCircles[3*k] = Best_coord[k_best][0];
    BestCircles[3*k+1] = Best_coord[k_best][1];
    BestCircles[3*k+2] = sum_x/sum_tot;
    printf("Meilleur cercle n.%d : (%f,%f), rayon : %f\n",k,BestCircles[3*k],BestCircles[3*k+1],BestCircles[3*k+2]);
  }
  //p_transf = p_marqueur;
  printf("Fin 4eme phase\n");
  free(VAL);free(GXX);free(GXY);free(GYY);free(TEMP);
  free(CURV);free(GX);free(GY);free(HOUGH);

}
