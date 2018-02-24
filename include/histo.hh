#ifndef _HISTO_
#define _HISTO_

void histogramme(Image<int> p,int* histo);
void TransAffine(Image<int>& p,double a,double b);
void Egalisation(Image<int>& p);
void Histogramme_Cumule(Image<int> p,int* histo_cumule);
void SeuilSimple(Image<int>& p,int seuil);
void SeuilDouble(Image<int>& p,int seuil1,int seuil2);
#endif
