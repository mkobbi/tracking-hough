#ifndef _MORFO_
#define _MORFO_

// Traitements morphologiques élémentaires

void Morphobase(Image<int>& p,char *nom,int tableau[NMAX][NMAX],int cx);
void Appliq_morpho(Image<int>& p,int erodil,int sens,int tab[NMAX][NMAX],int xdeb,int ydeb,int xfin,int yfin);

// Traitements locaux dans le 8 ou 4-voisinage

void Morpho_48(Image<int>& p,int erodil,int conex);

// Reconstruction géodésique;
// Algo séquentiel, parcours sens video et inverse
void Reconstruit_seq(Image<int>& p,Image<int>& p_ref,int cx);
// Avec files d'attente 
void Reconstruit(Image<int>& p,Image<int>& p_ref,int cx);
void Max_regionaux(Image<int>& p,int cx);

// Quelques traitements binaires
void MaxPropag(Image<int>& p,int conex);
void Ero_ult(Image<int>& p,int conex);
// Calcul de fonctions distance
int PutDistance(Image<int>& p,int i,int j,char *nom,int sens);
void FctDistance(Image<int>& p,char *nom);

// Nivellements
void Nivellement(Image<int>& p,char *nom,int tab[NMAX][NMAX],float z,int taille,int cx);


// Filtres alternés séquentiels

void FAS(Image<int>& p,int taille,int sens);

void FGrain(Image<int>& p,int taille,int sens,int cx);


void LPE(Image<int>& p,int conex);

void Filtrage_dynamique(Image<int>& p,int hauteur,int conex);



#endif
