#ifndef _CONTOUR_
#define _CONTOUR_

// Détecteur de contours de Marr et Hildreth
void Marr_Hildreth(Image<int>& p,float sigma,float th_haut,float th_bas);
// Détecteur de contours de Canny
void Canny(Image<int>& p,float sigma,float th_haut,float th_bas);
// Détecteur de contours dans le repère local (gradient,isophote)
void Contours_RepereLocal(Image<int>& p,float sigma,float th_haut,float th_bas);

#endif
