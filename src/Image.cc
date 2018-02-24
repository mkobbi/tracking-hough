
#include <stdio.h>
#include <string.h>

#include "Image.hh"

#define ABS(a) (a)>0?(a):(-a)
#define MAX(a,b) (a)>(b)?(a):(b)
#define MIN(a,b) (a)<(b)?(a):(b) 


int ChangeBord(int val) {BORD=val;return 1;}


/**************************************************************/
/******************** Operateurs booleens *********************/
/**************************************************************/

BOOL NON(BOOL i) {return 255-i;} 

BOOL ET(BOOL i,BOOL j){
  if (i<j) return i; else return j;}

BOOL OU(BOOL i,BOOL j){
  if (i>j) return i; else return j;}
      
BOOL XOU(BOOL i,BOOL j) {if (i!=j) return(255); else return (0);}

/************************************************************/
/*     classe Info           -          Implementation      */
/************************************************************/
template <class T> PixLC<T>::PixLC() { NomOperation=NULL;Suivant=NULL;Precedent=NULL;pix=NULL;}
template <class T> PixLC<T>::~PixLC() {if (NomOperation) delete NomOperation;}
template <class T> void PixLC<T>::InsereAvant(T * pix)
{
  PixLC<T>* plc=new PixLC<T>();
  plc->pix=pix;
  PixLC<T>* prec=Precedent;
  Precedent=plc;
  plc->Suivant=this;
  plc->Precedent=prec;
  if (prec) prec->Suivant=plc;
}
// template <class T> void PixLC<T>::InsereApres(T * pix)
// {
//   PixLC<T>* plc=new PixLC<T>();
//   plc->pix=pix;
//   PixLC<T>* suiv=Suivant;
//   Suivant=plc;
//   plc->Precedent=this;
//   plc->Suivant=prec;
//   if (prec) prec->Precedent=plc;
// }

template <class T> Info<T>::Info(){ pix = NULL; PixListe=new PixLC<T>();S=0; }
template <class T> Info <T>::Info(int taille) {
            pix = new T[taille];
            S = taille;
	    PixListe=new PixLC<T>();
	    PixListe->pix=pix;
}
template <class T> Info<T>::Info(int t,T* im) {
            pix = new T[t];
            S = t;
	    // for (int i=0;i<S;i++) pix[i] = im[i];
            memcpy(pix,im,S*sizeof(T));
	    PixListe=new PixLC<T>();
	    PixListe->pix=pix;
}
template <class T> Info<T>::Info(const Info<T>& ii){ 
  S = ii.S;
  pix = ii.pix;
  PixListe=new PixLC<T>();
  PixListe->pix=pix;
}
template <class T> Info<T>::~Info() { 
  PixLC<T> *prec=PixListe->Precedent;
  for (PixLC<T> *pl=PixListe,*suiv=pl?pl->Suivant:NULL;pl;pl=suiv,suiv=pl?pl->Suivant:NULL)
    {
      delete [] pl->pix;
      delete pl;
    }
  for (PixLC<T> *pl=prec,*suiv=pl?pl->Precedent:NULL;pl;pl=suiv,suiv=pl?pl->Precedent:NULL)
    {
      delete [] pl->pix;
      delete pl;
    }
}
  
template <class T> Info<T>&  Info<T>::operator=(const Info<T>& ii){
   int i;
    if (pix != ii.pix){
      if (pix) delete[] pix;
      pix = new T[ii.S];
      S = ii.S;
      memcpy(pix,ii.pix,S*sizeof(T));
      //  for (i=0;i<S;i++) pix[i] = ii.pix[i];
        }
    return *this;} 
template <class T> T Info<T>::operator[](const int& i) const{
    return pix[i];}
template <class T> T* Info<T>::ref() const{
    return pix;}
template <class T> void Info<T>::Ecrit_info(const int& i,const T& z) {
    pix[i] = z; }

template <class T> void Info<T>::Empile(char * NomOperation) {
  if (!PixListe) return; //pas d'image
  for (PixLC<T> *pl=PixListe->Suivant,*suiv=pl?pl->Suivant:NULL;pl;pl=suiv,suiv=pl?pl->Suivant:NULL)
    {
      delete [] pl->pix;
       delete pl;
    }
  PixListe->Suivant=NULL;
  T * npix=new T[S];
  memcpy(npix,pix,sizeof(T)*S);
  PixListe->InsereAvant(npix);
  PixListe->NomOperation=strdup(NomOperation);
}
template <class T> bool Info<T>::DepileAvant() {
   if ( (!PixListe) || (!PixListe->Suivant)) return false;
   PixListe->pix=pix;
   PixListe=PixListe->Suivant;
   pix=PixListe->pix;
   return PixListe->Suivant;
 }

template <class T> bool Info<T>::DepileArriere() {
  if ((!PixListe) || (!PixListe->Precedent)) return false;
  PixListe->pix=pix;
  PixListe=PixListe->Precedent;
  pix=PixListe->pix;
  return PixListe->Precedent;
 }

template <class T> const char * Info<T>::DerniereOperation() {
  if (!PixListe) return NULL;
  return PixListe->NomOperation;
}
 
/************************************************************/
/* structure Vecteur        -     Implementation            */
/************************************************************/


  Vecteur::Vecteur() { x = y = 0;}
  Vecteur::Vecteur(int i, int j) { x = i; y = j;}
  Vecteur::Vecteur(const Vecteur& v) {x=v.x;y=v.y;}
  Vecteur::~Vecteur() {}

  Vecteur& Vecteur::operator=(const Vecteur& v) {
    x = v.x; y = v.y; return *this;}
  Vecteur Vecteur::operator+(const Vecteur& v) {
    return Vecteur(x + v.x,y + v.y);}
  Vecteur& Vecteur::operator+=(const Vecteur& v) {
    x += v.x; y += v.y; 
    return *this;}
  int Vecteur::operator*(const Vecteur& v) {
    return x*v.x + y*v.y;}
  Vecteur Vecteur::operator-(){
    return Vecteur(-x,-y); }
  Vecteur Vecteur::rot90() {
    return Vecteur(-y,x);}
  Vecteur Vecteur::rot270() {
    return Vecteur(y,-x);}
  int Vecteur::operator==(const Vecteur& t)
    {
      return ((x == t.x) && (y == t.y));
    }


/************************************************************/
/* classe Noyau          -        Implementation            */
/************************************************************/
 
 
  Noyau::Noyau() {}  
  Noyau::~Noyau() {}
 
Noyau& Noyau::operator= (const Noyau& k)
{
  debx=k.debx;deby=k.deby;finx=k.finx;finy=k.finy;
  int larg=(finx-debx);
  int haut=(finy-deby);
  if (val != k.val)
    memcpy(val,k.val,larg*haut*sizeof(int));
  return *this; }
 
int Noyau::KX0(){ return debx;}
 
int Noyau::KX1(){ return finx;}
 
int Noyau::KY0(){ return deby;}     

int Noyau::KY1(){ return finy;}
 
Noyau& Noyau::Noyau_efface ()
   {
   for (int i=debx;i<=finx;i++)
     for (int j=deby;j<=finy;j++)
       val[i][j]=0;
   return *this;
   }

Noyau& Noyau::Noyau_setmasque(int tab[NMAX][NMAX])
{
int first=1;
 for (int i=0;i<NMAX;i++)
     for (int j=0;j<NMAX;j++) {
          if (tab[i][j]) {
            if (first) {
                first=0;
                debx=i-NMAX/2;deby=j-NMAX/2;
                finx=i-NMAX/2;finy=j-NMAX/2;
                } else {
                finx=i-NMAX/2;finy=j-NMAX/2;
                }
               }
             val[i-NMAX/2][j-NMAX/2]=tab[i][j];
             }       
return *this;
}

Noyau Noyau::Noyau_transpose ()
{
  Noyau q;
  q.debx=-finx;q.finx=-debx;
  q.deby=-finy;q.finy=-deby;
  for (int i=q.debx;i<=q.finx;i++)
    for (int j=q.deby;j<=q.finy;j++)
      q.val[i][j]=val[-i][-j];
   return q;
   }   
         
Noyau Noyau::operator+ (const Noyau& ker)
   {
   Noyau q;
   for (int i=debx;i<=finx;i++)
     for (int j=deby;j<=finy;j++)
       q.val[i][j] = ker.val[i][j]+val[i][j];
   return q;
    }
 
Noyau Noyau::operator* (const int x)
  {
    Noyau q;
    for (int i=debx;i<=finx;i++)
      for (int j=deby;j<=finy;j++)
        q.val[i][j] = x*val[i][j];
      return q;
    }       
        
/************************************************************/
/* classe FIFO           -        Implementation            */
/************************************************************/

FIFO::FIFO() {}

FIFO::FIFO(int taille)
{listx = new int[taille];
 listy = new int[taille];
 first = 0; last = 0;
 size = taille;}
 
FIFO::~FIFO()
{if (listx) delete[] listx;
 if (listy) delete[] listy;}

FIFO& FIFO::operator=(const FIFO& ff) {
  //size = ff.size;
  if (ff.last-ff.first > size) {
    printf("Opération irrégulière sur FIFO\n");
    return *this;
  } else {
    last = ff.last;
    first = ff.first;
    memcpy(listx,ff.listx,size*sizeof(int));
    memcpy(listy,ff.listy,size*sizeof(int));
    return *this;
  }
}

FIFO& FIFO::FIFO_delete () {
  first = 0; last = 0;
  return *this;
}

FIFO& FIFO::FIFO_alloc(int taille)
{listx = new int[taille];
 listy = new int[taille];
 first = 0; last = 0;
 size = taille;
 return *this;
}

void FIFO::FIFO_free()
{
  if (listx) delete[] listx;
  if (listy) delete[] listy;
}

FIFO& FIFO::FIFO_add (int x,int y) {
listx[last]=x;listy[last]=y;
last++;
if (last == size) this->FIFO_reinit();
return *this;
}

FIFO& FIFO::FIFO_reinit () {
for (int i=first;i<=last;i++) {
  listx[i-first]=listx[i];listy[i-first]=listy[i];}
last -=first;
first = 0;
return *this;
}

void FIFO::FIFO_get(int& x,int& y) {
x = listx[first];
y = listy[first++];}

int FIFO::FIFO_empty() const {
return (first==last);}

/************************************************************/
/* classe Image           -        Implementation            */
/************************************************************/


template <> Image<int>::Image() {}
template <> Image<int>::Image(const int& i,const int& j) : I(i*j)
{ W = i; L = j; }
template <> Image<int>::Image(const Image<int>& pp)
    {W = pp.W; 
     L = pp.L;
     I = pp.I;
     }
template <> Image<int>::~Image() {}
  
template <> Image<int>& Image<int>::operator= (const Image<int>& p)
    {W=p.W; L=p.L; I= p.I;
    return *this; }

template <> int* Image<int>::PI(){ return I.ref();}

template <> int Image<int>::PW(){ return W;}

template <> int Image<int>::PL(){ return L;}

template <> int Image<int>::X(int i,int j) const {
  switch(BORD) {
  case 0 : {
    if ((i<0)||(i>=W)||(j<0)||(j>=L)) return 0;
    else return I[j * W + i];
    break;
  }
  case 1 : {
    if (i < 0) i = 0;
    else if (i >= W) i = W-1;
    if (j < 0) j = 0;
    else if (j >= L) j = L-1;
    return I[j * W + i];
    break;
  }
  case 2 : {
    if (i < 0) i = W+i;
    else if (i >= W) i = i-W;
    if (j < 0) j = L+j;
    else if (j >= L) j = j-L;
    return I[j * W + i];
    break;
  }                        
}
}  

template <> Image<int> Image<int>::Image_translat(const Vecteur& v) const {

  Image<int> q(W,L);
  for (int j=0;j<q.L;j++)
    for (int i=0;i < q.W; i++)
      q.I.Ecrit_info(j * W + i,X(i-v.x,j-v.y));
  return q; 
}

template <> Image<int> Image<int>::Nord() const
{ return Image_translat(Vecteur(0,1));}

template <> Image<int> Image<int>::Sud() const
{ return Image_translat(Vecteur(0,-1));}

template <> Image<int> Image<int>::Est() const
{ return Image_translat(Vecteur(-1,0));}

template <> Image<int> Image<int>::Ouest() const
{ return Image_translat(Vecteur(1,0));}

template <> Image<int> Image<int>::NE() const
{ return Image_translat(Vecteur(-1,1));}

template <> Image<int> Image<int>::SE() const
{ return Image_translat(Vecteur(-1,-1));}

template <> Image<int> Image<int>::NW() const
{ return Image_translat(Vecteur(1,1));}

template <> Image<int> Image<int>::SW() const
{ return Image_translat(Vecteur(1,-1));}

template <> int* Image<int>::Image_hist ()
{
  return 0;
        } 

template <> Image<int>& Image<int>::Image_efface () 
   {
   for (int i=0;i<W*L;i++) I.Ecrit_info(i,0);
   return *this;
   }

template <> Image<int> Image<int>::operator+ (const Image<int>& p)
{
  Image<int> q(W,L);
  for (int j=0;j<L;j++)
    for (int i=0;i<W;i++) 
      q.I.Ecrit_info(j * W + i,p.X(i,j)+X(i,j)); 
  return q;
}

template <> Image<BOOL> Image<BOOL>::operator& (const Image<BOOL>& p) const
   {
     int index;
   Image<BOOL> q(W,L);
  for (int j=0;j<L;j++)
    for (int i=0;i<W;i++){
      index = j * W + i;
      q.I.Ecrit_info(index,ET(I[index],p.I[index]));
    } 
  return q;
   }

template <> Image<BOOL> Image<BOOL>::operator| (const Image<BOOL>& p) const
{
  int index;
  Image<BOOL> q(W,L);
  for (int j=0;j<L;j++)
    for (int i=0;i<W;i++) { 
      index = j * W + i;
      q.I.Ecrit_info(index,OU(I[index],p.I[index]));
    } 
  return q;
}

template <> Image<BOOL> Image<BOOL>::operator! () const
{
  int index;
  Image<BOOL> q(W,L);
  for (int j=0;j<L;j++)
    for (int i=0;i<W;i++) { 
      index = j * W + i;
      q.I.Ecrit_info(index,NON(I[index]));
    }
  return q;
}

template <> Image<BOOL> Image<BOOL>::operator^ (const Image<BOOL>& p) const
{
  int index;
  Image<BOOL> q(W,L);
  for (int j=0;j<L;j++)
    for (int i=0;i<W;i++) { 
      index = j * W + i;
      q.I.Ecrit_info(index,XOU(I[index],p.I[index]));
    }
  return q;
}

template <> Image<BOOL> Image<BOOL>::operator== (const Image<BOOL>& p) const
{
  int index;
  Image<BOOL> q(W,L);
  for (int j=0;j<L;j++)
    for (int i=0;i<W;i++) { 
      index = j * W + i;
      q.I.Ecrit_info(index,NON(XOU(I[index],p.I[index])));
    }
  return q;
}

template <> Image<BOOL> Image<BOOL>::operator- (const Image<BOOL>& p) const
{
  int index;
  Image<BOOL> q(W,L);
  for (int j=0;j<L;j++)
    for (int i=0;i<W;i++) {
      index = j * W + i;
      q.I.Ecrit_info(index,ET(I[index],NON(p.I[index])));
    }
  return q;
}

template <> Image<BOOL> Image<BOOL>::operator<= (const Image<BOOL>& p) const
{
  int index;
  Image<BOOL> q(W,L);
  for (int j=0;j<L;j++)
    for (int i=0;i<W;i++) {
      index = j * W + i;
      q.I.Ecrit_info(index,OU(NON(I[index]),p.I[index]));
    }     
  return q;
}

template <> Image<BOOL> Image<BOOL>::operator>= (const Image<BOOL>& p) const
{
  int index;
  Image<BOOL> q(W,L);
  for (int j=0;j<L;j++)
    for (int i=0;i<W;i++) {
      index = j * W + i;
      q.I.Ecrit_info(index,OU(I[index],NON(p.I[index]))); 
    } 
  return q;
}

template <> Image<BOOL> Image<int>::Image_seuil (const int s)
{
  int index=0;
  Image<BOOL> q(W,L);
  for (int j=0;j<L;j++)
    for (int i=0;i<W;i++) {
      q.I.Ecrit_info(index,255*(I[index] >= s));
      index++;
    }
  return q;
}

template <> Image<BOOL> Image<int>::Image_levelsets (const int step)
{
  int index=0;
  int k,l,val;
  int reste_min,reste_max;
  Image<BOOL> q(W,L);
  for (int j=0;j<L;j++)
    for (int i=0;i<W;i++) {
      reste_min=step;
      reste_max=0;
      for (k=-1;k<=1;k++)
	for (l=-1;l<=1;l++) {
	  val = (X(i-k,j+l)%step);
	  if (val<reste_min) reste_min = val;
	  else 
	    if (val>reste_max) reste_max = val;
	}
      q.I.Ecrit_info(index,255*((reste_min > step/2)||(reste_max < step/2)));
      index++;
    }
  return q;
}

template <> Image<int> Image<int>::Image_diff(const Image<int>& p)
{
  int index;
  Image<int> q(W,L);
  for (int j=0;j<L;j++)
    for (int i=0;i<W;i++) {
      index = j * W + i;
      q.I.Ecrit_info(index,I[index]-p.I[index]);
    }  
  return q;
}

template <> Image<int> Image<int>::Image_diffsignee(const Image<int>& p)
{
  int index=0;
  int valeur;
  Image<int> q(W,L);
  for (int j=0;j<L;j++)
    for (int i=0;i<W;i++) {
      valeur = 128+I[index]-p.I[index];
      if (valeur<0) valeur = 0;
      else if (valeur>255) valeur = 255;
      q.I.Ecrit_info(index,valeur);
      index++;
    }  
  return q;
}

template <> Image<int> Image<int>::Image_superpose(const Image<BOOL>& p,const int nb,const int coul)
{
  int index;
  Image<int> q(W,L);
  for (int j=0;j<L;j++)
    for (int i=0;i<W;i++) {
      index = j * W + i;
      if (p.I[index] == nb)
	q.I.Ecrit_info(index,500 + coul);
      else q.I.Ecrit_info(index,I[index]);
    }  
  return q;
}

template <> int Image<int>::Imagetopgm (char *nom) const
{
	FILE *fp;

        if ((fp=fopen(nom,"w"))==NULL) {
                perror("fopen\n");
                return -1; 
        };
        fprintf(fp,"P5\n");
        fprintf(fp,"# CREATOR : God\n");
        fprintf(fp,"%d %d\n%d\n",W,L,255);
        for (int i=0;i<W*L;i++) fprintf(fp,"%c",I[i]);
        fclose(fp);
	return 0;
}

template <> int Image<int>::Image_compte_pixel () const
  {
    int index;
    int compt = 0;
    for (int j=0;j<L;j++)
      for (int i=0;i<W;i++) {
	index = j * W + i;
	if (I[index] != 0) compt++;
      }
    return compt;
  } 

template <> void Image<int>::Empile(char * NomOperation) {I.Empile(NomOperation);}
template <> bool Image<int>::DepileAvant() {return I.DepileAvant();}
template <> bool Image<int>::DepileArriere() {return I.DepileArriere();}
template <> const char * Image<int>::DerniereOperation() {return I.DerniereOperation();}
