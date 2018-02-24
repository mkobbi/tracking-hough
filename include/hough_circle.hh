#ifndef _HOUGHCIRCLE_
#define _HOUGHCIRCLE_

//Transformée de Hough classique 1-to-many à partir des contours 
void Hough_circles_contours(double *HOUGH,Image<int> p,int rayon_max,double *vote_max);
//Transformée de Hough dense 1-to-1 à partir de la courbure et du gradient
void Hough_circles(double *HOUGH,Image<int> p,float sigma_init,int nb_echelles,float gamma_integr,int rayon_max,double *vote_max);
//Transformée de Hough dense 1-to-1 en un seul plan (seulement le centre) pour la THD 2-1
void HCircleCentre(Image<int>& p_img,Image<int>& p_transf,float sigma_init,int nb_ech,int rayon_min);
//Transformée de Hough dense en deux temps (THD 2-1), centre d'abord, 
//puis vote réduit pour les rayons
void HCC_BestCircles(Image<int>& p_img,Image<int>& p_transf,float sigma_init,int nb_ech,int rayon_min,int rayon_max,int nb_cercles,float *BestCircles);
//Routine de mise à jour d'un plan de la TH 3d entière à partir du tableau de double
void MaJ_PlanHough(Image<int>& p_transf_plan,int numero_plan,double *HOUGH,double vote_max);
//Routine de remise à zéro de la TH 3d
void Reset_Hough(double *HOUGH,int rayon_max,int width_img,int height_img);
//Filtrage exponentiel récursif pour la TH 3d
void FiltreExp_Hough_3d(double *HOUGH,double gamma,int rayon_max,int w,int h,double *vote_max);
//Routine pour le tracé d'un cercle
void Trace_cercle(Image<int>& p,int x_centre,int y_centre,int rayon, int couleur);
//Calcul du meilleur cercle et mise à jour de la TH 3d
void Trace_Best_Hough_Circles(double *HOUGH,int rayon_max,int w,int h,int *xbest,int *ybest,int *rhobest);

#endif
