#include "Image.hh"
#include "filtre.hh"
#include "morfo.hh"
#include "feature.hh"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

float abs(float a) {if (a<0) return -a; else return a;}
double valabs (double x) {if (x<0) return -x; else return x;}
double sqr(float a) {return a*a;}
int ABS(int a) {if (a<0) return -a; else return a;}
int max(int a,int b) {if (a>b) return a; else return b;}
int min(int a,int b) {if (a<b) return a; else return b;}
int module(int a,int b) {return (int)sqrt(a*a+b*b);}

void FH(Image<int>& p,int c)
{
  int index=0;
  Image<int> q(p);
  int* PIX=p.PI();
  int L=p.PL();
  int W=p.PW();
  
for (int j=0;j<L;j++)
  for (int i=0;i<W;i++) 
    PIX[index++]=q.X(i-1,j)+c*q.X(i,j)+q.X(i+1,j);
}

void FV(Image<int>& p,int c)
{
  int index=0;
  Image<int> q(p);
  int* PIX=p.PI();
  int L=p.PL();
  int W=p.PW();
  
for (int j=0;j<L;j++)
  for (int i=0;i<W;i++)
    PIX[index++]=q.X(i,j-1)+c*q.X(i,j)+q.X(i,j+1);
}

static void array_swap(int t[], int m, int n) {
  int temp = t[m];
  t[m] = t[n];
  t[n] = temp;
}

static int partition(int m, int n, int t[]) {
  int v = t[m];                 /* valeur pivot */
  int i = m-1;
  int j = n+1;         /* indice final du pivot */
  while (1) {
    do {
      j--;
    } while (t[j] > v);
    do {
      i++;
    } while (t[i] < v);
    if (i<j) {
      array_swap(t, i, j);
    } else {
      return j;
    }
  }
}

void quicksort(int m, int n, int t[]) {
  if (m<n) {
    int p = partition(m,n,t);
    quicksort(m,p,t);
    quicksort(p+1,n,t);
  }
}

void FiltreFixe(Image<int>& p,char *nom)
{
  int index=0;
  Image<int> q(p);
  int* PIX=p.PI();
  int L=p.PL();
  int W=p.PW();
  
  
if (strcmp(nom,"laplace4")==0) {
  for (int j=0;j<L;j++)
    for (int i=0;i<W;i++) {
      PIX[index]=128+q.X(i-1,j)+q.X(i,j-1)+q.X(i+1,j)+q.X(i,j+1)-4*q.X(i,j);
      if (PIX[index] > 255) PIX[index] = 255;
      else if (PIX[index] < 0) PIX[index] = 0;
      index++;
    }
} else 
if (strcmp(nom,"laplace8")==0) {
  for (int j=0;j<L;j++)
    for (int i=0;i<W;i++) {
      PIX[index]=128+q.X(i-1,j-1)+q.X(i,j-1)+q.X(i+1,j-1)
	+q.X(i-1,j)+q.X(i+1,j)+q.X(i-1,j+1)
	+q.X(i,j+1)+q.X(i+1,j+1)-8*q.X(i,j);
      if (PIX[index] > 255) PIX[index] = 255;
      else if (PIX[index] < 0) PIX[index] = 0;
      index++;
    }
} else 
if (strncmp(nom,"gradh",5)==0) {
  if (strcmp(nom,"gradhsobel")==0) FV(q,2); else FV(q,1);
  for (int j=0;j<L;j++)
    for (int i=0;i<W;i++) {
      PIX[index]=128+q.X(i+1,j)-q.X(i-1,j);
      if (PIX[index] > 255) PIX[index] = 255;
      if (PIX[index] < 0) PIX[index] = 0;
      index++;
    }
} else 
if (strncmp(nom,"gradv",5)==0) {
  if (strcmp(nom,"gradvsobel")==0) FH(q,2); else FH(q,1);
  for (int j=0;j<L;j++)
    for (int i=0;i<W;i++) {
      PIX[index]=128+q.X(i,j+1)-q.X(i,j-1);
      if (PIX[index] > 255) PIX[index] = 255;
      if (PIX[index] < 0) PIX[index] = 0; 
      index++;
    }
} else 
  if ((strcmp(nom,"sobel")==0) || (strcmp(nom,"prewitt")==0)){
    Image<int> r(p);
    if (strcmp(nom,"sobel")==0) {FH(q,2);FV(r,2);}   
    else {FH(q,1);FV(r,1);}
    for (int j=0;j<L;j++)
      for (int i=0;i<W;i++)
	PIX[index++] = min(255,module(q.X(i,j+1)-q.X(i,j-1),r.X(i+1,j)-r.X(i-1,j)));
  } else
   if ((strncmp(nom,"arg",3)==0)){
     Image<int> r(p);
     int norme,gx,gy;
     double PI=3.141592;
     if (strcmp(nom,"argsobel")==0) 
       {FH(q,2);FV(r,2);} else {FH(q,1);FV(r,1);}
     index=0;
     for (int j=0;j<L;j++)
       for (int i=0;i<W;i++) {
	 gy = q.X(i,j+1)-q.X(i,j-1);
	 gx = r.X(i+1,j)-r.X(i-1,j);
	 norme = module(gx,gy);
	 if (norme)
	   PIX[index] = (int)(((asin((float)gy/(float)norme)+PI/2)*255)/PI); 
	 else PIX[index] = 0;
	 index++;
       }
   } else
if (strcmp(nom,"nagao")==0) {
  int ecart_min,num;
  int a,b;
  int liste[9];
  for (int j=0;j<L;j++)
   for (int i=0;i<W;i++) {
     // Fenetre centrale et initialisation
     num=0;
     for (a=-1;a<=1;a++)
       for (b=-1;b<=1;b++)
         liste[num++] = q.X(i+a,j+b);
     quicksort(0,8,liste);
     ecart_min = liste[8]-liste[0];
     PIX[index] = liste[5];
     // Fenetres N,S,E,W
       num=0;
       liste[num++] = q.X(i,j);
       for (a=-1;a<=1;a++) 
          for (b=-2;b<=-1;b++)
            liste[num++] = q.X(i+a,j+b);
       quicksort(0,6,liste);
       if (liste[6]-liste[0] < ecart_min) {ecart_min = liste[6]-liste[0];
                                           PIX[index] = liste[4];}
       num=0;
       liste[num++] = q.X(i,j);
       for (a=-1;a<=1;a++) 
          for (b=1;b<=2;b++)
            liste[num++] = q.X(i+a,j+b);
       quicksort(0,6,liste);
       if (liste[6]-liste[0] < ecart_min) {ecart_min = liste[6]-liste[0];
                                           PIX[index] = liste[4];}
       num=0;
       liste[num++] = q.X(i,j);
       for (a=-2;a<=-1;a++) 
          for (b=-1;b<=1;b++)
            liste[num++] = q.X(i+a,j+b);
       quicksort(0,6,liste);
       if (liste[6]-liste[0] < ecart_min) {ecart_min = liste[6]-liste[0];
                                           PIX[index] = liste[4];}
        num=0;
       liste[num++] = q.X(i,j);
       for (a=1;a<=2;a++) 
          for (b=-1;b<=1;b++)
            liste[num++] = q.X(i+a,j+b);
       quicksort(0,6,liste);
       if (liste[6]-liste[0] < ecart_min) {ecart_min = liste[6]-liste[0];
                                           PIX[index] = liste[4];}
       // Fenetres NW,SE,NE,SW
       num=0;
       for (a=-2;a<=0;a++) 
          for (b=-2;b<=0;b++)
            if (!(abs(a-b)==2)) liste[num++] = q.X(i+a,j+b);
       quicksort(0,6,liste);
       if (liste[6]-liste[0] < ecart_min) {ecart_min = liste[6]-liste[0];
                                           PIX[index] = liste[4];}
        num=0;
       for (a=0;a<=2;a++) 
          for (b=0;b<=2;b++)
            if (!(abs(a-b)==2)) liste[num++] = q.X(i+a,j+b);
       quicksort(0,6,liste);
       if (liste[6]-liste[0] < ecart_min) {ecart_min = liste[6]-liste[0];
                                           PIX[index] = liste[4];}
         num=0;
       for (a=0;a<=2;a++) 
          for (b=-2;b<=0;b++)
            if (!(abs(a-b)==2)) liste[num++] = q.X(i+a,j+b);
       quicksort(0,6,liste);
       if (liste[6]-liste[0] < ecart_min) {ecart_min = liste[6]-liste[0];
                                           PIX[index] = liste[4];}
         num=0;
       for (a=-2;a<=0;a++) 
          for (b=0;b<=2;b++)
            if (!(abs(a-b)==2)) liste[num++] = q.X(i+a,j+b);
       quicksort(0,6,liste);
       if (liste[6]-liste[0] < ecart_min) {ecart_min = liste[6]-liste[0];
                                           PIX[index] = liste[4];} 
       index++;
   }  
    
  } 
}

void Filtre_Rang(Image<int>& p,int tableau[NMAX][NMAX],float rang) {

int x1=NMAX,x2=0,y1=NMAX,y2=0;
int longueur=0;
Image<int> q(p);

  for (int i=0;i<NMAX;i++)
     for (int j=0;j<NMAX;j++)
          if (tableau[i][j]) {
                longueur++;
                if (i<x1) x1=i;if (j<y1) y1=j;
                if (i>x2) x2=i;if (j>y2) y2=j;
                }
                
  int* liste=new int[longueur];
  int index,valeur,num;
  int* PIX=p.PI();
  int L=p.PL();
  int W=p.PW();

for (int j=0;j<L;j++)
    for (int i=0;i<W;i++) {
      index = j * W + i;
      num=0;
      for (int a = x1;a <= x2;a++)
        for (int b = y1;b <= y2;b++) {
            if (tableau[a][b]) {
                    valeur = q.X(i-(a-NMAX/2),j-(b-NMAX/2));
                    liste[num++]=valeur;
                     }
             }
      quicksort(0,longueur-1,liste);
      PIX[index] = liste[(int)(rang*(longueur-1))];
            }
  
}

int CalculCoeff(char *nom,float z,int precision,int x,int y) {

  if (strcmp(nom,"exponentiel")==0)
    return (int)((pow(10,precision))*((z*z)/4)*exp(-z*(abs(x)+abs(y))));
  else 
    if (strcmp(nom,"gauss")==0) {
      float PI = 3.141582635;
      return (int)((pow(10,precision)/(2*PI*z*z))*exp(-(x*x + y*y)/(2*z*z)));}
    else
      if (strcmp(nom,"moyenne")==0) return 1;
}

void FiltreMasque(Image<int>& p,char *nom,int tableau[NMAX][NMAX],float z) {

 int x1=NMAX,x2=0,y1=NMAX,y2=0;
 int longueur=0;
 int a,b;

  for (int i=0;i<NMAX;i++)
     for (int j=0;j<NMAX;j++)
          if (tableau[i][j]) {
            tableau[i][j]=CalculCoeff(nom,z,3,i-NMAX/2,j-NMAX/2);
	    longueur+=tableau[i][j];
	    if (i<x1) x1=i;if (j<y1) y1=j;
	    if (i>x2) x2=i;if (j>y2) y2=j;
	  }
                
  int index,valeur;
  Image<int> q(p);
  int* PIX=p.PI();
  int L=p.PL();
  int W=p.PW();

for (int j=0;j<L;j++)
    for (int i=0;i<W;i++) {
      index = j * W + i;
      valeur=0;
      for (int a = x1;a <= x2;a++)
        for (int b = y1;b <= y2;b++)
            if (tableau[a][b]) valeur += tableau[a][b]*q.X(i-(a-NMAX/2),j-(b-NMAX/2));
      if (longueur)  PIX[index] = valeur/longueur;
            }
    }

void Comp_connexe(Image<int>& p) { 

  int h=p.PL();
  int w=p.PW();
  Image<int> Etiq(p);
  int i,j,x,y,a,b;
  int num = 0;
  int minx,miny,maxx,maxy;
  int* LAB=Etiq.PI();
  int* PIX=p.PI();
  int xmin[100];
  int xmax[100];
  int ymin[100];
  int ymax[100];
  FIFO f(w*h);
  
  for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      if (Etiq.X(i,j) == 255) {
	LAB[j*w+i] = ++num;
	xmin[num] = i;
	xmax[num] = i;
	ymin[num] = j;
	ymax[num] = j;
	f.FIFO_add(i,j);
	while (!f.FIFO_empty()) {
	  f.FIFO_get(x,y);
	  for (int a=-1;a<=1;a++)
	    for (int b=-1;b<=1;b++) {
	      if ((a!=0)||(b!=0)) {
		if ((Etiq.X(x+a,y+b)==255)&&(x+a>=0)&&(y+b>=0)&&(x+a<w)&&(y+b<h)) {
		  f.FIFO_add(x+a,y+b);
		  LAB[(y+b)*w+x+a]=num;
		  if (x+a<xmin[num]) xmin[num]= x+a; else if (x+a>xmax[num]) xmax[num]= x+a;
		  if (y+b<ymin[num]) ymin[num]= y+b; else if (y+b>ymax[num]) ymax[num]= y+b;
		}
	      }
	    }
	}
      }
    }
   for (j=0;j<h;j++)
    for (i=0;i<w;i++) 
      if (p.X(i,j)==255) PIX[j*w+i] = (LAB[j*w+i]*255)/num;
   printf("%d composantes connexes\n",num);
   for (i=1;i<=num;i++) {
     for (j=xmin[i];j<=xmax[i];j++) {PIX[ymin[i]*w+j]=503;PIX[ymax[i]*w+j]=503;}
     for (j=ymin[i];j<=ymax[i];j++) {PIX[j*w+xmin[i]]=503;PIX[j*w+xmax[i]]=503;}
   }       
}

void Mean(Image<int>& p,int rayon) { 
  int h=p.PL();
  int w=p.PW();
  int i,j,k,sum;
  int index=0;
  float area;
  int *PIX = p.PI();
  Image<int> q(p);
  int *QPIX = q.PI();

  // Passage horizontal
  index=0;
  for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      for (k = 1;k <= rayon;k++)
	QPIX[index] += p.X(i-k,j)+p.X(i+k,j);
      index++;
    }
      // Passage vertical
  index=0;
  area = (float)(2*rayon+1)*(2*rayon+1);
  for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      for (k = 1;k <= rayon;k++)
	PIX[index] += q.X(i,j-k)+q.X(i,j+k);
      PIX[index] = (int)round(((float)(PIX[index]))/area);
      index++;
    }
}

// Implantation récursive du noyau exponentiel

void Exponential(Image<int>& p,float gamma) { 
  int h=p.PL();
  int w=p.PW();
  int i,j;
  int index=0;
  float area;
  int *PIX = p.PI();
  double *VAL;
  double alpha;

  VAL = (double *)calloc(w*h,sizeof(double));
  alpha = 1 - exp(-gamma);
  // Passage horizontal, causal et conversion en double
  index=0;
  for (j=0;j<h;j++) {
    VAL[index++] = p.X(0,j);
    for (i=1;i<w;i++) {
      VAL[index] = alpha*p.X(i,j)+(1-alpha)*VAL[index-1];
      index++;
    }
  }
  // Passage horizontal, anticausal
  for (j=0;j<h;j++) {
    index = j*w + w-2;
    for (i=w-2;i>=0;i--) {
      VAL[index] = alpha*VAL[index]+(1-alpha)*VAL[index+1];
      index--;
    }
  }
  // Passage vertical, causal
  for (i=0;i<w;i++) {
    index = i;
    for (j=1;j<h;j++) {
      index +=w;
      VAL[index] = alpha*VAL[index]+(1-alpha)*VAL[index-w];
    }
  }
  // Passage vertical, anticausal et conversion en entier
  for (i=0;i<w;i++) {
    index = (h-1)*w + i;
    for (j=h-2;j>=0;j--) {
      index -=w;
      VAL[index] = alpha*VAL[index]+(1-alpha)*VAL[index+w];
      PIX[index] = (int)(VAL[index]);
    }
  }
  free(VAL);
}

// Implantation récursive du noyau gaussien et de ses dérivées
// Cf. "Recursive implementation of the Gaussian filter"
// Ian T. Young and Lucas J. van Vliet 
// in Signal Processing 44(95) 139--151

// Parametres :
// * dérivée = 0 --> Convolution par filtre gaussien G
// direction: 0 = selon y; 1 = selon x; 2 = x puis y;
// * dérivée = 1 --> Dérivées premières de G
// direction: 0 = conv. par Gy; 1 = conv. par Gx; 2 = module du gradient;
// * dérivée = 2 --> Dérivées secondes de G
// direction: 0 = conv. par Gyy; 1 = conv. par Gxx; 2 = conv. par Gxy; 3 = laplacien;



void Gauss_Rec(double *IN,double *OUT,int w,int h,int derivee,int direction,
	       float sigma,double &valmin,double &valmax) { 
  int i,j,val;
  int index=0;
  double q,b0,b1,b2,b3,bb;
  double *GX,*GY;

  // Calcul des coefficients du filtrage récursif
  if (sigma<2.5) q = 3.97156 - 4.14554*sqrt(1-0.26891*sigma);
  else q = 0.98711*sigma - 0.9633;
  
  b0 = 1.57825 + 2.44413*q + 1.4281*q*q + 0.422205*q*q*q;
  b1 = 2.44413*q + 2.85619*q*q + 1.26661*q*q*q;
  b2 = -1.4281*q*q - 1.26661*q*q*q;
  b3 = 0.422205*q*q*q;
  bb = 1 - (b1+b2+b3)/b0;
  valmin = 1000;
  valmax = -1000;
  
  if (derivee==0) {
    if (direction==0) {
      for (j=0;j<h;j++) {
	OUT[j*w] = IN[j*w];
	OUT[j*w+1] = bb*IN[j*w+1]+((b1+b2+b3)*OUT[j*w])/b0;
	OUT[j*w+2] = bb*IN[j*w+2]+(b1*OUT[j*w+1]+(b2+b3)*OUT[j*w])/b0; 
	for (i=3;i<w;i++) {
	  index=j*w+i;
	  OUT[index] = bb*IN[index]+(b1*OUT[index-1]+b2*OUT[index-2]+b3*OUT[index-3])/b0;
	}
      }
      for (j=h-1;j>=0;j--) {
	OUT[(j+1)*w-2] = bb*OUT[(j+1)*w-2]+((b1+b2+b3)*OUT[(j+1)*w-1])/b0;
	OUT[(j+1)*w-3] = bb*OUT[(j+1)*w-3]+(b1*OUT[(j+1)*w-2]+(b2+b3)*OUT[(j+1)*w-1])/b0;
	for (i=w-4;i>=0;i--) {
	  index=j*w+i;
	  OUT[index] = bb*OUT[index]+(b1*OUT[index+1]+b2*OUT[index+2]+b3*OUT[index+3])/b0;
	}
      }
    } else {
      for (i=0;i<w;i++) {
	OUT[i] = IN[i];
	OUT[w+i] = bb*IN[w+i]+((b1+b2+b3)*OUT[i])/b0;
	OUT[2*w+i] = bb*IN[2*w+i]+(b1*OUT[w+i]+(b2+b3)*OUT[i])/b0;
	for (j=3;j<h;j++) {
	  index=j*w+i;
	  OUT[index] = bb*IN[index]+(b1*OUT[index-w]+b2*OUT[index-2*w]+b3*OUT[index-3*w])/b0;
	}
      }
      for (i=w-1;i>=0;i--) {
	OUT[w*(h-2)+i] = bb*OUT[w*(h-2)+i]+((b1+b2+b3)*OUT[w*(h-1)+i])/b0;
	OUT[w*(h-3)+i] = bb*OUT[w*(h-3)+i]+(b1*OUT[w*(h-2)+i]+(b2+b3)*OUT[w*(h-1)+i])/b0;
	for (j=h-4;j>=0;j--) {
	  index=j*w+i;
	  OUT[index] = bb*OUT[index]+(b1*OUT[index+w]+b2*OUT[index+2*w]+b3*OUT[index+3*w])/b0;
	}
      }
    }
  } else if (derivee==1) {
    if (direction==0) {
      GX = (double *)calloc(w*h,sizeof(double));
      for (j=0;j<h;j++) {
	GX[j*w] = 0;
	GX[j*w+1] = 0.5*bb*(IN[j*w+2]-IN[j*w]);
	GX[j*w+2] = 0.5*bb*(IN[j*w+3]-IN[j*w+1])+(b1*GX[j*w+1])/b0; 
	for (i=3;i<w-1;i++) {
	  index=j*w+i;
	  GX[index] = 0.5*bb*(IN[index+1]-IN[index-1])+(b1*GX[index-1]+b2*GX[index-2]+b3*GX[index-3])/b0;
	}
	GX[(j+1)*w-1] = (b1*GX[(j+1)*w-2]+b2*GX[(j+1)*w-3]+b3*GX[(j+1)*w-4])/b0;
      }
      for (j=h-1;j>=0;j--) {
	OUT[(j+1)*w-1] = GX[(j+1)*w-1];
	OUT[(j+1)*w-2] = bb*GX[(j+1)*w-2]+((b1+b2+b3)*OUT[(j+1)*w-1])/b0;
	OUT[(j+1)*w-3] = bb*GX[(j+1)*w-3]+(b1*OUT[(j+1)*w-2]+(b2+b3)*OUT[(j+1)*w-1])/b0;
	for (i=w-4;i>=0;i--) {
	  index=j*w+i;
	  OUT[index] = bb*GX[index]+(b1*OUT[index+1]+b2*OUT[index+2]+b3*OUT[index+3])/b0;
	  if (OUT[index]>valmax) valmax=OUT[index]; else if (OUT[index]<valmin) valmin=OUT[index];
	}
      }
      free(GX);
    } else {
      GY = (double *)calloc(w*h,sizeof(double));
      for (i=0;i<w;i++) {
	GY[i] = 0;
	GY[w+i] = 0.5*bb*(IN[2*w+i]-IN[i]);
	GY[2*w+i] = 0.5*bb*(IN[3*w+i]-IN[w+i])+(b1*GY[w+i])/b0;
	for (j=3;j<h-1;j++) {
	  index=j*w+i;
	  GY[index] = 0.5*bb*(IN[index+w]-IN[index-w])+(b1*GY[index-w]+b2*GY[index-2*w]+b3*GY[index-3*w])/b0;
	}
	GY[(h-1)*w+i] = (b1*GY[index-w]+b2*GY[index-2*w]+b3*GY[index-3*w])/b0;
      }
      for (i=w-1;i>=0;i--) {
	OUT[w*(h-1)+i] = GY[w*(h-1)+i];
	OUT[w*(h-2)+i] = bb*GY[w*(h-2)+i]+((b1+b2+b3)*OUT[w*(h-1)+i])/b0;
	OUT[w*(h-3)+i] = bb*GY[w*(h-3)+i]+(b1*OUT[w*(h-2)+i]+(b2+b3)*OUT[w*(h-1)+i])/b0;
	for (j=h-4;j>=0;j--) {
	  index=j*w+i;
	  OUT[index] = bb*GY[index]+(b1*OUT[index+w]+b2*OUT[index+2*w]+b3*OUT[index+3*w])/b0;
	  if (OUT[index]>valmax) valmax=OUT[index]; else if (OUT[index]<valmin) valmin=OUT[index];
	}
      }
      free(GY);
    }
  } else if (derivee==2) {
    if (direction==0) {
      GX = (double *)calloc(w*h,sizeof(double));
      for (j=0;j<h;j++) {
	GX[j*w] = 0;
	GX[j*w+1] = bb*(IN[j*w+1]-IN[j*w]);
	GX[j*w+2] = bb*(IN[j*w+2]-IN[j*w+1])+(b1*GX[j*w+1])/b0; 
	for (i=3;i<w;i++) {
	  index=j*w+i;
	  GX[index] = bb*(IN[index]-IN[index-1])+(b1*GX[index-1]+b2*GX[index-2]+b3*GX[index-3])/b0;
	}
      }
      for (j=h-1;j>=0;j--) {
	OUT[(j+1)*w-1] = 0;
	OUT[(j+1)*w-2] = bb*(GX[(j+1)*w-1]-GX[(j+1)*w-2]);
	OUT[(j+1)*w-3] = bb*(GX[(j+1)*w-2]-GX[(j+1)*w-3])+(b1*OUT[(j+1)*w-2])/b0;
	for (i=w-4;i>=0;i--) {
	  index=j*w+i;
	  OUT[index] = bb*(GX[index+1]-GX[index])+(b1*OUT[index+1]+b2*OUT[index+2]+b3*OUT[index+3])/b0;
	  if (OUT[index]>valmax) valmax=OUT[index]; else if (OUT[index]<valmin) valmin=OUT[index];
	}
      }
      free(GX);
    } else {
      GY = (double *)calloc(w*h,sizeof(double));
      for (i=0;i<w;i++) {
	GY[i] = 0;
	GY[w+i] = bb*(IN[w+i]-IN[i]);
	GY[2*w+i] = bb*(IN[2*w+i]-IN[w+i])+(b1*GY[w+i])/b0;
	for (j=3;j<h-1;j++) {
	  index=j*w+i;
	  GY[index] = bb*(IN[index]-IN[index-w])+(b1*GY[index-w]+b2*GY[index-2*w]+b3*GY[index-3*w])/b0;
	}
      }
      for (i=w-1;i>=0;i--) {
	OUT[w*(h-1)+i] = 0;
	OUT[w*(h-2)+i] = bb*(GY[w*(h-1)+i]-GY[w*(h-2)+i]);
	OUT[w*(h-3)+i] = bb*(GY[w*(h-2)+i]-GY[w*(h-3)+i])+(b1*OUT[w*(h-2)+i])/b0;
	for (j=h-4;j>=0;j--) {
	  index=j*w+i;
	  OUT[index] = bb*(GY[index+w]-GY[index])+(b1*OUT[index+w]+b2*OUT[index+2*w]+b3*OUT[index+3*w])/b0;
	  if (OUT[index]>valmax) valmax=OUT[index]; else if (OUT[index]<valmin) valmin=OUT[index];
	}
      }
      free(GY);
    }
  }
}

void Gaussian(Image<int>& p,int derivee,int direction,float sigma) { 
  int h=p.PL();
  int w=p.PW();
  int i,j,val;
  int index=0;
  int *PIX = p.PI();
  double *VAL,*GX,*GY,*GXX,*GYY;
  double valmin=10000;
  double valmax=-10000;

  VAL = (double *)calloc(w*h,sizeof(double));
  GX = (double *)calloc(w*h,sizeof(double));
  GY = (double *)calloc(w*h,sizeof(double));
  // Passage en double...
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      VAL[index] = (double) (PIX[index++]);
  
  if (derivee==0) {
    if (direction==0) Gauss_Rec(VAL,GX,w,h,0,0,sigma,valmin,valmax);
    else if (direction==1) Gauss_Rec(VAL,GX,w,h,0,1,sigma,valmin,valmax);
    else if (direction==2) {
      Gauss_Rec(VAL,GY,w,h,0,0,sigma,valmin,valmax);
      Gauss_Rec(GY,GX,w,h,0,1,sigma,valmin,valmax);
    }
    // Repassage en entier...
    index=0;
    for (i=0;i<w;i++)
      for (j=0;j<h;j++)
	PIX[index] = (int) (GX[index++]);
  } else if (derivee==1) {
    GXX = (double *)calloc(w*h,sizeof(double));
    GYY = (double *)calloc(w*h,sizeof(double));
    if (direction==0) {
      Gauss_Rec(VAL,GXX,w,h,0,1,sigma,valmin,valmax);
      Gauss_Rec(GXX,GX,w,h,1,0,sigma,valmin,valmax);
    } else if (direction==1) {
      Gauss_Rec(VAL,GXX,w,h,0,0,sigma,valmin,valmax);
      Gauss_Rec(GXX,GX,w,h,1,1,sigma,valmin,valmax);
    } else if (direction==2) {
      Gauss_Rec(VAL,GXX,w,h,0,1,sigma,valmin,valmax);
      Gauss_Rec(GXX,GX,w,h,1,0,sigma,valmin,valmax);
      Gauss_Rec(VAL,GYY,w,h,0,0,sigma,valmin,valmax);
      Gauss_Rec(GYY,GY,w,h,1,1,sigma,valmin,valmax);
      valmin=255;valmax=0;
      index=0;
      for (i=0;i<w;i++)
	for (j=0;j<h;j++) {
	  GX[index] = sqrt(GX[index]*GX[index] + GY[index]*GY[index]);
	  if (GX[index]>valmax) valmax=GX[index]; else if (GX[index]<valmin) valmin=GX[index];
	  index++;
	}
    }
    // Repassage en entier avec équilibrage des niveaux de gris...
    index=0;
    if (direction != 2) {
      for (i=0;i<w;i++)
	for (j=0;j<h;j++) {
	  if (GX[index]>0) 
	    PIX[index] = 128 + (int) ((GX[index]*127)/valmax);
	  else 
	    PIX[index] = 128 - (int) ((GX[index]*127)/valmin);
	  index++;
	}
    } else {
      for (i=0;i<w;i++)
	for (j=0;j<h;j++) 
	  PIX[index] = (int) ((GX[index++]-valmin)*255/(valmax-valmin));
    }
    free(GXX);
    free(GYY);
  } else if (derivee==2) {
    GXX = (double *)calloc(w*h,sizeof(double));
    GYY = (double *)calloc(w*h,sizeof(double));
    if (direction==0) {
      Gauss_Rec(VAL,GXX,w,h,0,1,sigma,valmin,valmax);
      Gauss_Rec(GXX,GX,w,h,2,0,sigma,valmin,valmax);
    } else if (direction==1) {
      Gauss_Rec(VAL,GXX,w,h,0,0,sigma,valmin,valmax);
      Gauss_Rec(GXX,GX,w,h,2,1,sigma,valmin,valmax);
    } else if (direction==2) {
      Gauss_Rec(VAL,GXX,w,h,0,1,sigma,valmin,valmax);
      Gauss_Rec(GXX,GY,w,h,1,0,sigma,valmin,valmax);
      Gauss_Rec(GY,GYY,w,h,0,0,sigma,valmin,valmax);
      Gauss_Rec(GYY,GX,w,h,1,1,sigma,valmin,valmax);
    } else if (direction==3) {
      Gauss_Rec(VAL,GXX,w,h,0,1,sigma,valmin,valmax);
      Gauss_Rec(GXX,GX,w,h,2,0,sigma,valmin,valmax);
      Gauss_Rec(VAL,GYY,w,h,0,0,sigma,valmin,valmax);
      Gauss_Rec(GYY,GY,w,h,2,1,sigma,valmin,valmax);
      valmin=10000;valmax=-10000;
      index=0;
      for (i=0;i<w;i++)
	for (j=0;j<h;j++) {
	  GX[index] = GX[index] + GY[index];
	  if (GX[index]>valmax) valmax=GX[index]; else if (GX[index]<valmin) valmin=GX[index];
	  index++;
	}
    }
    // Repassage en entier avec équilibrage des niveaux de gris...
    index=0;
    for (i=0;i<w;i++)
      for (j=0;j<h;j++) {
	if (GX[index]>0) 
	  PIX[index] = 128 + (int) ((GX[index]*127)/valmax);
	else 
	  PIX[index] = 128 - (int) ((GX[index]*127)/valmin);
	index++;
      }
    free(GXX);
    free(GYY);
  }
}

void Display_RI_LJ(Image<int>& p,int ordre,int pdt_module,float sigma) { 
  int h=p.PL();
  int w=p.PW();
  int index;
  int *PIX = p.PI();
  double *Vin,*Vout;
  double valmin=10000;
  double valmax=-10000;

  Vin = (double *)calloc(w*h,sizeof(double));
  Vout = (double *)calloc(w*h,sizeof(double));
  // Passage en double...
  for (index=0;index<w*h;index++)
    Vin[index] = (double) (PIX[index]);
  // Calcul du LJ en double
  int normalisation = 2;//Normalisation scale-space
  Calcul_RI_LJ(Vin,ordre,pdt_module,sigma,w,h,valmin,valmax,Vout,normalisation);
  // repassage en entier 
  if (ordre == 0) {// le domaine est inchangé
    for (index=0;index<w*h;index++) 
      PIX[index] = (int)(Vout[index]);
  } else if (ordre == 1) {// module : non signé
    for (index=0;index<w*h;index++) {
      PIX[index] = (int) ((Vout[index]*255.0)/valmax);
    }  
  } else {//grandeurs d'ordre 2 : signées
    for (index=0;index<w*h;index++) {
      if (Vout[index]>0) 
	PIX[index] = 128 + (int) ((Vout[index]*127)/valmax);
      else 
	PIX[index] = 128 - (int) ((Vout[index]*127)/valmin);
    }
  }
  free(Vin);free(Vout);
}


void Contrast_enhance(Image<int>& p,float sigma,float gain) { 
  int h=p.PL();
  int w=p.PW();
  int i,j,val;
  int index=0;
  int *PIX = p.PI();
  Image<int> qq(w,h);
  int *QPIX = qq.PI();
  double *VAL,*GX,*GY,*GXX;
  double valmin,valmax;

  VAL = (double *)calloc(w*h,sizeof(double));
  GXX = (double *)calloc(w*h,sizeof(double));
  // Passage en double...
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      VAL[index] = (double) (PIX[index++]);
  // Calcul du laplacien...
  GX = (double *)calloc(w*h,sizeof(double));
  GY = (double *)calloc(w*h,sizeof(double));
  Gauss_Rec(VAL,GXX,w,h,0,1,sigma,valmin,valmax);
  Gauss_Rec(GXX,GX,w,h,2,0,sigma,valmin,valmax);
  Gauss_Rec(VAL,GXX,w,h,0,0,sigma,valmin,valmax);
  Gauss_Rec(GXX,GY,w,h,2,1,sigma,valmin,valmax);
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++) {
      PIX[index] = (int) (PIX[index] -gain*(GX[index] + GY[index]));
      if (PIX[index]<0) PIX[index]=0; else 
	 if (PIX[index]>255) PIX[index]=255;
      index++;
    }
  free(GX);
  free(GY);
  free(GXX);
  free(VAL); 
}

void Diffusion_anisotrope(Image<int>& p,int type,float K,int nbiter) { 
  int h=p.PL();
  int w=p.PW();
  int n,i,j,index;
  double *GN,*GS,*GE,*GW,*PV,*PH,*P;
  int *PIX = p.PI();
  float lambda = 0.1;
  
  
  GN = (double *)calloc(w*h,sizeof(double));
  GS = (double *)calloc(w*h,sizeof(double));
  GE = (double *)calloc(w*h,sizeof(double));
  GW = (double *)calloc(w*h,sizeof(double));
  PV = (double *)calloc(w*h,sizeof(double));
  PH = (double *)calloc(w*h,sizeof(double));
  P = (double *)calloc(w*h,sizeof(double));

  // Passage en double...
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      P[index] = (double) (PIX[index++]);
  
  for (n=1;n<=nbiter;n++) {
    // Calcul des gradients (Sobel) 
    index = 0;
    for (j=0;j<h;j++)
      for (i=0;i<w;i++) {
	if ((i>0)&&(i<w-1))
	  PH[index] = P[index-1]+2*P[index]+P[index+1];
	else PH[index] = P[index];
	index++;
      }
    index=0;
    for (j=0;j<h;j++)
      for (i=0;i<w;i++) {
	if ((j>0)&&(j<h-1))
	  PV[index] = P[index-w]+2*P[index]+P[index+w];
	else PV[index] = P[index];
	index++;
      }
    index=0;
    for (j=0;j<h;j++)
      for (i=0;i<w;i++) {
	if (i>0) GW[index] = PV[index-1]-PV[index]; else GW[index] = 0;
	if (i<w-1) GE[index] = PV[index+1]-PV[index]; else GE[index] = 0;
	if (j>0) GN[index] = PH[index-w]-PH[index]; else GN[index] = 0;
	if (j<h-1) GS[index] = PH[index+w]-PH[index]; else GS[index] = 0;
	index++;
      }
    // Diffusion
    index=0;
    if (type ==0) {
      for (j=0;j<h;j++)
	for (i=0;i<w;i++) { 
	  P[index] += lambda*(exp(-sqr(valabs(GW[index])/K))*GW[index]
			      + exp(-sqr(valabs(GE[index])/K))*GE[index]
			      + exp(-sqr(valabs(GN[index])/K))*GN[index]
			      + exp(-sqr(valabs(GS[index])/K))*GS[index]);
	  index++;
	}
    } else {
      for (j=0;j<h;j++)
	for (i=0;i<w;i++) {
	  P[index] += lambda*((1/(1+sqr(valabs(GW[index])/K)))*GW[index]
			      + (1/(1+sqr(valabs(GE[index])/K)))*GE[index]
			      + (1/(1+sqr(valabs(GN[index])/K)))*GN[index]
			      + (1/(1+sqr(valabs(GS[index])/K)))*GS[index]);
	 index++;
	}
    }
  }
  // repassage en entier...
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++) {
      PIX[index] = (int) (P[index]);
      if (PIX[index]<0) PIX[index]=0;
      if (PIX[index]>255) PIX[index]=255;
      index++;
    }
  free(GN);free(GS);free(GE);free(GW);free(PV);free(PH);free(P);
}

void Mean_curvature(Image<int>& p,float sigma,int nbiter) { 
  int h=p.PL();
  int w=p.PW();
  int i,j,k,val;
  int index=0;
  int *PIX = p.PI();
  double *VAL,*GX,*GY,*TEMP,*TEMP1,*GXX,*GYY,*GXY;
  double valmin=10000;
  double valmax=-10000;
  double seuil;
  double alpha = 1.0;
  
  VAL = (double *)calloc(w*h,sizeof(double));
  TEMP = (double *)calloc(w*h,sizeof(double));
  TEMP1 = (double *)calloc(w*h,sizeof(double));
  GX = (double *)calloc(w*h,sizeof(double));
  GY = (double *)calloc(w*h,sizeof(double));
  GXX = (double *)calloc(w*h,sizeof(double));
  GYY = (double *)calloc(w*h,sizeof(double));
  GXY = (double *)calloc(w*h,sizeof(double));
  // Passage en double...
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h;j++)
      VAL[index] = (double) (PIX[index++]);
  for (k=1;k<=nbiter;k++) {
    Gauss_Rec(VAL,TEMP,w,h,0,1,sigma,valmin,valmax);
    Gauss_Rec(TEMP,GX,w,h,1,0,sigma,valmin,valmax);
    Gauss_Rec(TEMP,GXX,w,h,2,0,sigma,valmin,valmax);
    Gauss_Rec(VAL,TEMP,w,h,0,0,sigma,valmin,valmax);
    Gauss_Rec(TEMP,GY,w,h,1,1,sigma,valmin,valmax);
    Gauss_Rec(TEMP,GYY,w,h,2,1,sigma,valmin,valmax);
    Gauss_Rec(VAL,TEMP,w,h,0,1,sigma,valmin,valmax);
    Gauss_Rec(TEMP,TEMP1,w,h,1,0,sigma,valmin,valmax);
    Gauss_Rec(TEMP1,TEMP,w,h,0,0,sigma,valmin,valmax);
    Gauss_Rec(TEMP,GXY,w,h,1,1,sigma,valmin,valmax);
    index=0;
    for (i=0;i<w;i++)
      for (j=0;j<h;j++) {
	VAL[index] += alpha*(GXX[index]*GY[index]*GY[index]-2*GXY[index]*GX[index]*GY[index]+GYY[index]*GX[index]*GX[index])/(GX[index]*GX[index]+GY[index]*GY[index]);
	index++;
      }
  }
  index=0;
  for (i=0;i<w;i++)
    for (j=0;j<h-1;j++) {
      PIX[index]=(int)(VAL[index]);
      if (PIX[index]<0) PIX[index]=0;
      if (PIX[index]>255) PIX[index]=255;
      index++;
    }
  free(VAL);
  free(GX);
  free(GY);
  free(GXY);
  free(GXX);
  free(GYY);
  free(TEMP);
  free(TEMP1);
}

void Estim_bruit (Image<int>& p,float sigma,float proportion_pixels_norme) {

  int index=0;
  int i,j;
  Image<int> gpx(p);
  Image<int> mask(p);
  int* PIX=p.PI();
  int* GNORM=gpx.PI();
  int* PMSK=mask.PI();
  int HEIGHT=p.PL();
  int WIDTH=p.PW();
  float sum=0;
  float sumsq=0;
  float moy,ectype,diff;
  int effectif_estimateur;
  int valeur_theorique;

  // Calcul de la norme du gradient a l'echelle sigma
  Gaussian(gpx,1,2,sigma);
  // Calcul de la moyenne de la norme du gradient
  for (j=0;j<HEIGHT;j++)
    for (i=0;i<WIDTH;i++) {
      sum += GNORM[index];
      index++;
    }
  moy = sum/(HEIGHT*WIDTH);
  printf("Moyenne du module du gradient : %f\n",moy);
  // Calcul de la valeur theorique y tq F(y) = proportion_pixels_norme
  // avec F fonction de repartition de la loi de Rayleigh de moyenne "moy"
  valeur_theorique = (int)(sqrt((-4*moy*moy*logf(1-proportion_pixels_norme))/3.1415926535));
  // Affichage : je mets en rouge tous les points de
  // contraste hors norme, et prepare le masque de calcul
  // de l'estimateur du bruit
  index = 0;
  sum = 0;
  for (j=0;j<HEIGHT;j++)
    for (i=0;i<WIDTH;i++) {
      if (GNORM[index] > valeur_theorique) {
	PIX[index] = 501;
	PMSK[index] = 0;
      } else {
	sum++;
	PMSK[index] = 1;
      }
      index++;		       
    }
  // Verification de la conformite a la loi de Rayleigh
  printf("Nombre de pixels dans la norme : %f\n",(sum*100.0)/(WIDTH*HEIGHT));
  printf("(Prévision loi de Rayleigh : %f)\n",100.0*proportion_pixels_norme);
  // Estimation du bruit 
  index = 0;
  effectif_estimateur = 0;
  sumsq = 0;
  for (j=0;j<HEIGHT;j++)
    for (i=0;i<WIDTH;i++) {
      if (mask.X(i,j)*mask.X(i-1,j)*mask.X(i,j-1)*mask.X(i+1,j)*mask.X(i,j+1)) {
	effectif_estimateur += 1;
	diff = PIX[index] - 0.25*(p.X(i-1,j)+p.X(i,j-1)+p.X(i+1,j)+p.X(i,j+1));
	sumsq += diff*diff;
      }
      index++;
    }
  printf("Ecart-type du bruit estime : %f\n",sqrt(sumsq/effectif_estimateur));
      
}

// Le meme, en fonction...

float Fonction_Estim_bruit (Image<int>& p,float sigma,float proportion_pixels_norme) {

  int index=0;
  int i,j;
  Image<int> gpx(p);
  Image<int> mask(p);
  int* PIX=p.PI();
  int* GNORM=gpx.PI();
  int* PMSK=mask.PI();
  int HEIGHT=p.PL();
  int WIDTH=p.PW();
  float sum=0;
  float sumsq=0;
  float moy,ectype,diff;
  int effectif_estimateur;
  int valeur_theorique;
  float ec_bruit;

  // Calcul de la norme du gradient a l'echelle sigma
  Gaussian(gpx,1,2,sigma);
  // Calcul de la moyenne de la norme du gradient
  for (j=0;j<HEIGHT;j++)
    for (i=0;i<WIDTH;i++) {
      sum += GNORM[index];
      index++;
    }
  moy = sum/(HEIGHT*WIDTH);
  printf("Moyenne du module du gradient : %f\n",moy);
  // Calcul de la valeur theorique y tq F(y) = proportion_pixels_norme
  // avec F fonction de repartition de la loi de Rayleigh de moyenne "moy"
  valeur_theorique = (int)(sqrt((-4*moy*moy*logf(1-proportion_pixels_norme))/3.1415926535));
  // Masque de calcul de l'estimateur du bruit
  index = 0;
  sum = 0;
  for (j=0;j<HEIGHT;j++)
    for (i=0;i<WIDTH;i++) {
      if (GNORM[index] > valeur_theorique) {
	PMSK[index] = 0;
      } else {
	sum++;
	PMSK[index] = 1;
      }
      index++;		       
    }
  // Verification de la conformite a la loi de Rayleigh
  printf("Nombre de pixels dans la norme : %f\n",(sum*100.0)/(WIDTH*HEIGHT));
  printf("(Prévision loi de Rayleigh : %f)\n",100.0*proportion_pixels_norme);
  // Estimation du bruit 
  index = 0;
  effectif_estimateur = 0;
  sumsq = 0;
  for (j=0;j<HEIGHT;j++)
    for (i=0;i<WIDTH;i++) {
      if (mask.X(i,j)*mask.X(i-1,j)*mask.X(i,j-1)*mask.X(i+1,j)*mask.X(i,j+1)) {
	effectif_estimateur += 1;
	diff = PIX[index] - 0.25*(p.X(i-1,j)+p.X(i,j-1)+p.X(i+1,j)+p.X(i,j+1));
	sumsq += diff*diff;
      }
      index++;
    }
  ec_bruit = sqrt(sumsq/effectif_estimateur);
  printf("Ecart-type du bruit estime : %f\n",ec_bruit);
  return(ec_bruit);
      
}

// En utilisant la technique de Immerkaer (CVIU'96, vol64-2, pp300-302)
void Estim_bruit_v2 (Image<int>& p) {

  int index=0;
  int i,j;
  Image<int> Msk(p);
  int* PIX=p.PI();
  int* MASK=Msk.PI();
  int HEIGHT=p.PL();
  int WIDTH=p.PW();
  float sum=0;
  float sumsq=0;
  float moy,ectype,diff;

  // Calcul du masque pour supprimer la structure
  for (j=0;j<HEIGHT;j++)
    for (i=0;i<WIDTH;i++) {
      MASK[index++] = 4*p.X(i,j) -2*p.X(i-1,j) -2*p.X(i,j-1) -2*p.X(i+1,j) -2*p.X(i,j+1) + p.X(i-1,j-1) + p.X(i-1,j+1) + p.X(i+1,j-1) + p.X(i+1,j+1);
    }
  // Affichage : on affiche le masque sur 0--255
  index = 0;
  for (j=0;j<HEIGHT;j++)
    for (i=0;i<WIDTH;i++) {
      sumsq += MASK[index]*MASK[index];
      PIX[index] = 128 + MASK[index];
      if (PIX[index]>255) PIX[index]=255;
      else if (PIX[index]<0) PIX[index]=0;
      index++;		       
    }
  // Estimation du bruit 
  
  printf("Ecart-type du bruit estime : %f\n",sqrt(sumsq/(36*HEIGHT*WIDTH)));
      
}


void NLM_LJ_generic (Image<int>& p,int searching_window,int ordre_derivation, int nb_echelles, float sigma_init, float facteur_bruit) {

  int index=0;
  Image<int> copiep(p);
  int* COPIX=copiep.PI(); 
  int* PIX=p.PI();
  int HEIGHT=p.PL();
  int WIDTH=p.PW();

  int SEARCH_WINDOW = (int)floor(searching_window/2) ;

  int i0, j0, ic, jc, i, j, k;
  int mini,minj,maxi,maxj;
  int index2;

  float wmax,dist,weight,average,sweight;

  int nb_derivees = (ordre_derivation + 1)*(ordre_derivation + 2)/2;

  double **FEAT;
  
  // Allocation du vecteur de features
  FEAT = (double **)calloc(nb_echelles*nb_derivees,sizeof(double*));

  for (i = 0;i < nb_echelles*nb_derivees; i++) 
    FEAT[i] = (double *)calloc(WIDTH*HEIGHT,sizeof(double));

  Compute_LJ_Features(p,ordre_derivation,nb_echelles,sigma_init,FEAT);
  
  // Avec estimation automatique du bruit
  float h2 = facteur_bruit*Fonction_Estim_bruit(p,2.0,0.7);
  h2 = h2*h2;
  
  //on parcourt tous les points de l'image dans BitMapAg et on écrit le résultat dans BitMap
  index=0;
  for (j0=0;j0<HEIGHT ; j0++)
    {
      //printf("j0 = %d\n",j0);
      for (i0=0;i0<WIDTH ;i0++ )
    	{
	  wmax = 0; //poids max dans le voisinage
	  average = 0; //moyenne pondérée
	  sweight = 0; //poids total

	  //on parcours les points dans la fenetre de comparaison
	  mini=max(0,i0-SEARCH_WINDOW);
	  maxi=min(i0+SEARCH_WINDOW+1, WIDTH);
	  minj=max(0,j0-SEARCH_WINDOW);
	  maxj=min(j0+SEARCH_WINDOW+1, HEIGHT);

	  for (ic=mini;ic<maxi ;ic++ )
	    {
	      for (jc=minj;jc<maxj ;jc++ )
		{
		  if ((ic!=0) || (jc!=0))
		    {
		      //on calcule la distance entre (i0,j0) et (ic,jc)
		      index2 = jc*WIDTH+ic;
		      dist = dist_LJ(FEAT,index,index2,ordre_derivation,nb_echelles);
		      //on calcule le poids associé
		      weight = exp(- dist / h2);
		      if (weight > wmax) wmax = weight ;
		      average += weight * copiep.X(ic,jc);
		      sweight += weight;
		    }
		}
	    }
	  
	  //on traite à part le cas de y = x
	  average += wmax * PIX[index];
	  sweight += wmax;

	  PIX[index++]=(unsigned char)(average / sweight);
    	}
    }
  for (i = 0;i < nb_echelles*nb_derivees; i++) free(FEAT[i]);
  free(FEAT);
}


void FlatZone_Diffusion(Image<int>& p,int size_zone,int ordre_derivation,int nb_echelles,float sigma_init) {
  int index,i,j;
  int* PIX=p.PI();
  int HEIGHT=p.PL();
  int WIDTH=p.PW();
  
  float *DistTab_N,*DistTab_E,*DistTab_W,*DistTab_S;
  
  int nb_derivees = (ordre_derivation + 1)*(ordre_derivation + 2)/2;
  
  double **FEAT;
  
  // Allocation du vecteur de features
  FEAT = (double **)calloc(nb_echelles*nb_derivees,sizeof(double*));
  for (i = 0;i < nb_echelles*nb_derivees; i++) 
    FEAT[i] = (double *)calloc(WIDTH*HEIGHT,sizeof(double));
  // Calcul de la projection dans l'espace des features
  Compute_LJ_Features(p,ordre_derivation,nb_echelles,sigma_init,FEAT);
  

  // Pré-calcul des 4 tableaux de distance (carte de potentiels)
  DistTab_N = (float*)calloc(WIDTH*HEIGHT,sizeof(float));
  DistTab_E = (float*)calloc(WIDTH*HEIGHT,sizeof(float));
  DistTab_W = (float*)calloc(WIDTH*HEIGHT,sizeof(float));
  DistTab_S = (float*)calloc(WIDTH*HEIGHT,sizeof(float));
  for (index = 0;index < WIDTH*HEIGHT;index++) {
    if (index >= WIDTH) DistTab_N[index]=dist_LJ(FEAT,index,index-WIDTH,ordre_derivation, nb_echelles); else DistTab_N[index]=255;
    if ((index% WIDTH)!=0) DistTab_E[index]=dist_LJ(FEAT,index,index-1,ordre_derivation, nb_echelles); else DistTab_E[index]=255;
    if ((index% WIDTH)!=WIDTH-1) DistTab_W[index]=dist_LJ(FEAT,index,index+1,ordre_derivation, nb_echelles); else DistTab_W[index]=255;
    if (index < WIDTH*(HEIGHT-1)) DistTab_S[index]=dist_LJ(FEAT,index,index+WIDTH,ordre_derivation, nb_echelles); else DistTab_S[index]=255;
  }
  // On peut désallouer la mémoire des features
  for (i = 0;i < nb_echelles*nb_derivees; i++) free(FEAT[i]);
  free(FEAT);
  
  // Affichage des tableaux de distance (debug)
  float valmax = 0.0;
  for (index = 0;index < WIDTH*HEIGHT;index++)
    if (DistTab_W[index]>valmax) valmax = DistTab_W[index];
  for (index = 0;index < WIDTH*HEIGHT;index++)
    PIX[index] = (int)((DistTab_W[index]*255.0)/valmax);


  // désallocation des cartes de potentiels
  free(DistTab_N);free(DistTab_E);free(DistTab_W);free(DistTab_S);
}


void Haar_Filters(Image<int>& p,int w_filt,int h_filt,int ord_filt) {
  int index,i,j;
  int* PIX=p.PI();
  int HEIGHT=p.PL();
  int WIDTH=p.PW();
  
  float *SV,*S2;
  float *OUT;
  float valmax;

  int x1,x2,y1,y2;
  int xm1,ym1,xm2,ym2;
  
  SV = (float*)calloc(WIDTH*HEIGHT,sizeof(float));
  S2 = (float*)calloc(WIDTH*HEIGHT,sizeof(float));
  OUT = (float*)calloc(WIDTH*HEIGHT,sizeof(float));
  
  // Calcul de l'image intégrale
  SV[0] = S2[0] = (float)(PIX[0]);
  for (index=1;index<WIDTH;index++) {
    SV[index] = (float)(PIX[index]);
    S2[index] = S2[index-1] + (float)(PIX[index]);
  }
  for (index=WIDTH;index<WIDTH*HEIGHT;index++) {
    SV[index] = SV[index-WIDTH] + (float)(PIX[index]);
    if ((index % WIDTH)==0) S2[index] = SV[index];
    else S2[index] = S2[index-1] + SV[index];
  }
  // Calcul du filtre en fonction de l'ordre
  if (ord_filt == 0) {//Filtre intégral (moyenne)
    index = 0;
    valmax = 0;
    for (j=0;j<HEIGHT;j++)
      for (i=0;i<WIDTH;i++) {
	if (i<WIDTH-1-w_filt/2) x1 = i+w_filt/2; else x1 = WIDTH-1;
	if (j<HEIGHT-1-h_filt/2) y1 = j+h_filt/2; else y1 = HEIGHT-1;
	if (i>w_filt/2) x2 = i-w_filt/2; else x2 = 0;
	if (j>h_filt/2) y2 = j-h_filt/2; else y2 = 0;
	OUT[index] = S2[y1*WIDTH+x1]-S2[y1*WIDTH+x2]
	  -S2[y2*WIDTH+x1]+S2[y2*WIDTH+x2];
	if (valmax < OUT[index]) valmax = OUT[index];
	index++;
      }
  } else if (ord_filt == 2) {//Filtre dérivateur vertical
    index = 0;
    valmax = 0;
    for (j=0;j<HEIGHT;j++)
      for (i=0;i<WIDTH;i++) {
	if (i<WIDTH-1-w_filt/2) x1 = i+w_filt/2; else x1 = WIDTH-1;
	if (j<HEIGHT-1-h_filt/2) y1 = j+h_filt/2; else y1 = HEIGHT-1;
	if (i>w_filt/2) x2 = i-w_filt/2; else x2 = 0;
	if (j>h_filt/2) y2 = j-h_filt/2; else y2 = 0;
	OUT[index] = S2[y1*WIDTH+x1]-S2[y1*WIDTH+x2]
	  -2*S2[j*WIDTH+x1]+2*S2[j*WIDTH+x2]
	  +S2[y2*WIDTH+x1]-S2[y2*WIDTH+x2];
	if (valmax < abs(OUT[index])) valmax = abs(OUT[index]);
	index++;
      }
  } else if (ord_filt == 1) {//Filtre dérivateur horizontal
    index = 0;
    valmax = 0;
    for (j=0;j<HEIGHT;j++)
      for (i=0;i<WIDTH;i++) {
	if (i<WIDTH-1-w_filt/2) x1 = i+w_filt/2; else x1 = WIDTH-1;
	if (j<HEIGHT-1-h_filt/2) y1 = j+h_filt/2; else y1 = HEIGHT-1;
	if (i>w_filt/2) x2 = i-w_filt/2; else x2 = 0;
	if (j>h_filt/2) y2 = j-h_filt/2; else y2 = 0;
	OUT[index] = S2[y1*WIDTH+x1]-S2[y2*WIDTH+x1]
	  -2*S2[y1*WIDTH+i]+2*S2[y2*WIDTH+i]
	  +S2[y1*WIDTH+x2]-S2[y2*WIDTH+x2];
	if (valmax < abs(OUT[index])) valmax = abs(OUT[index]);
	index++;
      }
  } else if (ord_filt == 4) {//Filtre dérivée seconde verticale
    index = 0;
    valmax = 0;
    for (j=0;j<HEIGHT;j++)
      for (i=0;i<WIDTH;i++) {
	if (i<WIDTH-1-w_filt/2) x1 = i+w_filt/2; else x1 = WIDTH-1;
	if (j<HEIGHT-1-h_filt/2) y1 = j+h_filt/2; else y1 = HEIGHT-1;
	if (j<HEIGHT-1-h_filt/6) ym1 = j+h_filt/6; else ym1 = HEIGHT-1;
	if (i>w_filt/2) x2 = i-w_filt/2; else x2 = 0;
	if (j>h_filt/2) y2 = j-h_filt/2; else y2 = 0;
	if (j>h_filt/6) ym2 = j-h_filt/6; else ym2 = 0;
	OUT[index] = S2[y1*WIDTH+x1]-S2[y1*WIDTH+x2]-3*S2[ym1*WIDTH+x1]+3*S2[ym1*WIDTH+x2]
	  +3*S2[ym2*WIDTH+x1]-3*S2[ym2*WIDTH+x2]-S2[y2*WIDTH+x1]+S2[y2*WIDTH+x2];
	if (valmax < abs(OUT[index])) valmax = abs(OUT[index]);
	index++;
      }
  } else if (ord_filt == 3) {//Filtre dérivée seconde horizontale
    index = 0;
    valmax = 0;
    for (j=0;j<HEIGHT;j++)
      for (i=0;i<WIDTH;i++) {
	if (i<WIDTH-1-w_filt/2) x1 = i+w_filt/2; else x1 = WIDTH-1;
	if (j<HEIGHT-1-h_filt/2) y1 = j+h_filt/2; else y1 = HEIGHT-1;
	if (i<WIDTH-1-w_filt/6) xm1 = i+w_filt/6; else xm1 = WIDTH-1;
	if (i>w_filt/2) x2 = i-w_filt/2; else x2 = 0;
	if (j>h_filt/2) y2 = j-h_filt/2; else y2 = 0;
	if (i>w_filt/6) xm2 = i-w_filt/6; else xm2 = 0;
	OUT[index] = S2[y1*WIDTH+x1]-S2[y1*WIDTH+x2]-3*S2[y1*WIDTH+xm1]+3*S2[y1*WIDTH+xm2]
	  +3*S2[y2*WIDTH+xm1]-3*S2[y2*WIDTH+xm2]-S2[y2*WIDTH+x1]+S2[y2*WIDTH+x2];
	if (valmax < abs(OUT[index])) valmax = abs(OUT[index]);
	index++;
      }
  } else if (ord_filt == 5) {//Dérivée seconde diagonale
    index = 0;
    valmax = 0;
    for (j=0;j<HEIGHT;j++)
      for (i=0;i<WIDTH;i++) {
	if (i<WIDTH-1-w_filt/2) x1 = i+w_filt/2; else x1 = WIDTH-1;
	if (j<HEIGHT-1-h_filt/2) y1 = j+h_filt/2; else y1 = HEIGHT-1;
	if (i>w_filt/2) x2 = i-w_filt/2; else x2 = 0;
	if (j>h_filt/2) y2 = j-h_filt/2; else y2 = 0;
	OUT[index] = S2[y1*WIDTH+x1]-2*S2[j*WIDTH+x1]-2*S2[y1*WIDTH+i]
	  +4*S2[j*WIDTH+i]+S2[y1*WIDTH+x2]-2*S2[j*WIDTH+x2]
	  +S2[y2*WIDTH+x1]-2*S2[y2*WIDTH+i]+S2[y2*WIDTH+x2];	  
	if (valmax < abs(OUT[index])) valmax = abs(OUT[index]);
	index++;
      }
  }
  // Repassage en entier avec normalisation
  if (ord_filt == 0) {
  for (index=WIDTH;index<WIDTH*HEIGHT;index++) 
    PIX[index] = (int)((OUT[index]*255.0)/valmax);
  } else {
    for (index=WIDTH;index<WIDTH*HEIGHT;index++) 
      PIX[index] = (int)(128.0+(OUT[index]*127.0)/valmax);
  }  
    
  }
