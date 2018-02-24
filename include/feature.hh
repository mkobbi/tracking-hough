#ifndef _FEATURE_
#define _FEATURE_

// Détecteurs de points d'intérêt
void Morpho_interest(Image<int>& p,float sigma,int niveau_seuil,int display);
void Harris(Image<int>& p,float sigma1,float sigma2,int niveau_seuil,int display);

// Autour de la matrice Hessienne...
void Hessian_measures(Image<int>& p,float sigma,int type_measure);

// Routines de calcul des local jets
void Calcul_LJ(double *Gin, int ordre, float sigma, int w, int h, double &valmin, double &valmax, double *Gout,int normalize);
void Compute_LJ_Features(Image<int>& p_in,int ordre_derivation,int nb_echelles,float sigma_initial,double **FEAT_out);
float dist_LJ(double **FEATURES,int index1, int index2,int ordre_derivation, int nb_echelles);
float dist_vect_LJ(double *vecteur, double **FEATURES, int index,int ordre_derivation, int nb_echelles);
// Pour les LJ invariants en rotation
void Calcul_RI_LJ(double *Gin, int ordre, int produit_module, float sigma, int w, int h, double &valmin, double &valmax, double *Gout,int normalize);
// Calcul et affichage de différentielles vectorielles
void display_flow_vector(Image<int>& p,int index,int dx,int dy,int couleur);
void Diff_Flow_Measures(Image<int>& p,float sigma,int type_measure,int pas_affichage);
// Calcul du filtre de Frangi pour la détection des structures tubulaires
void Filter_Frangi_Complet(Image<int>& p,float sigma_min,float sigma_max,float pas_sigma,float alpha,float beta,int moymax);
void Filter_Frangi_1Echelle(double *IMG_IN,double *IMG_OUT,int w,int h,float sigma,float alpha,float beta);

#endif
