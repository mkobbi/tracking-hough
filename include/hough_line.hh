#ifndef _HOUGHLINE_
#define _HOUGHLINE_

//Transformée de Hough classique 1-to-many à partir du contour
void Hough_lines_contours(Image<int>& p_img,Image<int>& p_transf);
//Transformée de Hough dense 1-to-1 à partir du gradient
void Hough_gradient_Rho_Theta(Image<int>& p_img,Image<int>& p_transf,float sigma_init,int nb_ech);
void Hough_gradient_Alpha_Coord(Image<int>& p_img,Image<int>& p_transf,float sigma_init,int nb_ech);
//Routines de tracé des courbes image ou paramètre
void Trace_droite(Image<int>& p,float angle,int orig_x,int orig_y,int colour);
void Trace_sinusoide(Image<int>& p_transf,int x,int y);
void Trace_courbeAC(Image<int>& p_img,Image<int>& p_transf,int x,int y);
//Calcul de la meilleur droite et mise à jour de la TH
void Trace_Best_Hough_Lines(Image<int>& p_transf,int *xbest,int *ybest);

#endif
