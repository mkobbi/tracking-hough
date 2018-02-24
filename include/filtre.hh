#ifndef _FILTRE_
#define _FILTRE_

float abs(float a);
double valabs(double a);
int max(int a,int b);
int min(int a,int b);
int module(int a,int b);

static void array_swap(int t[], int m, int n);
static int partition(int m, int n, int t[]);
void quicksort(int m, int n, int t[]);


// Petits traitements locaux

void FH(Image<int>& p,int c);
void FV(Image<int>& p,int c);

// Traitements à noyaux fixes classiques

void FiltreFixe(Image<int>& p,char *nom);

// Traitements à noyaux variables

int CalculCoeff(char *nom,float z,int precision,int x,int y);

void Filtre_Rang(Image<int>& p,int tableau[NMAX][NMAX],float rang);
void FiltreMasque(Image<int>& p,char *nom,int tableau[NMAX][NMAX],float z);

// Etiquetage en composantes connexes
void Comp_connexe(Image<int>& p);

// Filtre moyenneur 
void Mean(Image<int>& p,int rayon);

// Filtre exponentiel (implantation récursive)
void Exponential(Image<int>& p,float gamma);

// Calcul du gaussien (implantation récursive)
void Gauss_Rec(double *IN,double *OUT,int w,int h,int derivee,int direction,
	       float sigma,double &valmin,double &valmax);
void Gaussian(Image<int>& p,int derivee,int direction,float sigma);
// Calcul des dérivées de Gaussienne dans le repère local
void Display_RI_LJ(Image<int>& p,int ordre,int pdt_module,float sigma);
// Réhaussement de contraste multi-échelle
void Contrast_enhance(Image<int>& p,float sigma,float gain);
// Schémas de diffusion non-linéaire
void Diffusion_anisotrope(Image<int>& p,int type,float K,int nbiter);
void Mean_curvature(Image<int>& p,float sigma,int nbiter);
// NL-Mean calculé dans l'espace des Local jets
void NLM_LJ_generic (Image<int>& p,int searching_window, int ordre_derivation, int nb_echelles, float sigma_init, float facteur_bruit);
// Routines d'estimation de la variance du bruit
void Estim_bruit (Image<int>& p,float sigma,float proportion_pixels_norme);
float Fonction_Estim_bruit (Image<int>& p,float sigma,float proportion_pixels_norme);
void Estim_bruit_v2 (Image<int>& p);

// Filtre anisotrope généralisé
void FlatZone_Diffusion(Image<int>& p,int size_zone,int ordre_derivation,int nb_echelles,float sigma_init);

void Haar_Filters(Image<int>& p,int w_filt,int h_filt,int ord_filt);

#endif
