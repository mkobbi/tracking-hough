#ifndef _IMAGE_
#define _IMAGE_

#define NMAX 9


/************************************************************/
/* definitions de types                                     */
/************************************************************/

template <class T> class Info;
struct Vecteur;       
class Noyau; 
class FIFO;   
template <class T> class Image;     
typedef int BOOL;

/************************************************************/
/* variables globales et leurs procédures                   */
/************************************************************/
 
static int BORD;   // 0 : bords nuls; 1 : bords geodesiques;
                   // 2 : bords wrappes (image torique);         
int ChangeBord(int val); 

/************************************************************/
/* classe Info                                              */
/************************************************************/
template <class T> class PixLC {
public:
  PixLC(); 
  ~PixLC(); 
  void InsereApres(T* pix);
  void InsereAvant(T * pix);
private:
  T * pix;
  char * NomOperation;
  PixLC<T> * Suivant, *Precedent;
  friend class Info<T>;
  
};
template <class T> class Info {
private :
  T* pix;
  PixLC<T> * PixListe; //element courant de la liste
   int S;
  Info<T>* Precedent;
  Info<T> * Suivant;
public :
  Info<T>();
  Info<T>(int taille);
  Info<T>(int t,T* im);
  Info<T>(const Info<T>& ii);
  ~Info();
  
  Info<T>& operator=(const Info<T>& ii);
  void Info_decale(const int& i);
  T operator[](const int& i) const;
  T index(const int& i) const;
  T* ref() const;
  void Ecrit_info(const int& i,const T& z);
  void Empile(char * NomOperation);
  bool DepileAvant();
  bool DepileArriere();
  const char * DerniereOperation();
};

/************************************************************/
/* structure Vecteur                                        */
/************************************************************/

struct Vecteur {

  int x;
  int y;

  Vecteur();
  Vecteur(int i, int j);
  Vecteur(const Vecteur& v);
  ~Vecteur();

  Vecteur& operator=(const Vecteur& v);
  Vecteur operator+(const Vecteur& v);
  Vecteur& operator+=(const Vecteur& v);
  int operator*(const Vecteur& v);
  Vecteur operator-();
  Vecteur rot90();
  Vecteur rot270();
  int operator==(const Vecteur& t);
};

/************************************************************/
/* classe Noyau                                             */
/************************************************************/
 
class Noyau {

public :
  int debx,finx;                /*   Position en x    */
  int deby,finy;                /*   Position en y    */
  int val[NMAX][NMAX];          /*      Valeurs       */

  Noyau();
  ~Noyau();
 
 int KX0();
 int KX1();
 int KY0();
 int KY1();
 
Noyau& operator= (const Noyau& k);
Noyau Noyau_transpose ();
Noyau& Noyau_efface ();
Noyau& Noyau_setmasque(int tab[NMAX][NMAX]);
       
// Somme et multiplication par un scalaire
Noyau operator+ (const Noyau& p);
Noyau operator* (const int x);
 
};
  
/************************************************************/
/* classe FIFO                                              */
/************************************************************/
 
class FIFO {

private :
  int first;               /*   Début de liste    */
  int last;                /*    Fin de liste     */
  int size;                /*    Taille admise    */
  int *listx;              /*     Liste des x     */
  int *listy;              /*     Liste des y     */
   
public :
  FIFO();
  FIFO(int taille);
  ~FIFO();

  FIFO& operator=(const FIFO& ff); 
  FIFO& FIFO_delete ();
  FIFO& FIFO_alloc(int taille);
  void FIFO_free();
  FIFO& FIFO_add (int x,int y);
  FIFO& FIFO_reinit (); 
  void FIFO_get(int& x,int& y);
  int FIFO_empty() const;

}; 
 
/************************************************************/
/* classe Image                                              */
/************************************************************/

template<class T> class Image {

private :
  int W,L;          /*     Dimensions Statiques      */
  Info<T> I;        /*      Information totale       */
public :
  Image();
  Image(const int& i,const int& j);
  Image(const Image<int>& pp);
  ~Image();

int* PI();
 int PW();
 int PL();
int X (int i,int j) const; 
Image<int>& operator= (const Image<int>& p);
Image<int> Image_translat(const Vecteur& v) const;
Image<int> Nord() const;
Image<int> Sud() const;
Image<int> Ouest() const;
Image<int> Est() const;
Image<int> NE() const;
Image<int> SE() const;
Image<int> NW() const;
Image<int> SW() const;
Image<int>& Image_efface ();
int* Image_hist ();
/*                     Operations booleennes                      */ 
Image<BOOL> operator& (const Image<BOOL>& p) const;
Image<BOOL> operator| (const Image<BOOL>& p) const;
Image<BOOL> operator! () const;
Image<BOOL> operator^ (const Image<BOOL>& p) const;
Image<BOOL> operator- (const Image<BOOL>& p) const;
Image<BOOL> operator== (const Image<BOOL>& p) const; 
Image<BOOL> operator<= (const Image<BOOL>& p) const;
Image<BOOL> operator>= (const Image<BOOL>& p) const;

//    Operations de base sur les niveaux de gris

Image<BOOL> Image_seuil (const int s);
Image<BOOL> Image_levelsets (const int step);
Image<int> operator+ (const Image<int>& p);
Image<int> Image_superpose (const Image<BOOL>& p,const int nb,const int coul);
Image<int> Image_diff(const Image<int>& p);
  Image<int> Image_diffsignee(const Image<int>& p);

//    Ecriture sous differents formats
int Imagetopgm (char *nom) const;  

int Image_compte_pixel () const;

  //gestion de la pile des modifications (operations undo et redo)
  void Empile(char * NomOperation);  /*empile l'image actuelle dans la pile arriere, efface la pile avant*/ 
  bool DepileArriere(); /*undo:remplace l'image actuelle par la premiere de la pile ariere, empile l'image actuelle dans la pile avant*/
  bool DepileAvant(); /*redo:remplace l'image actuelle par la premiere de la pile avant, empile l'image actuelle dans la pile arriere*/
  const char * DerniereOperation(); /*renvoit le nom de la derniere operation effectuee*/
};

#endif

