#include "Image.hh"
#include "filtre.hh"
#include "morfo.hh"
#include <string.h>
#include <stdio.h>
#include <math.h>

void Morphobase(Image<int>& p,char *nom,int tableau[NMAX][NMAX],int cx) {

int x1=NMAX,x2=0,y1=NMAX,y2=0;

 for (int i=0;i<NMAX;i++)
     for (int j=0;j<NMAX;j++)
          if (tableau[i][j]) {
                if (i<x1) x1=i;if (j<y1) y1=j;
                if (i>x2) x2=i;if (j>y2) y2=j;
                }

if (strcmp(nom,"erode")==0) Appliq_morpho(p,0,-1,tableau,x1,y1,x2,y2);
else
if (strcmp(nom,"dilate")==0) Appliq_morpho(p,1,-1,tableau,x1,y1,x2,y2);
else
if (strcmp(nom,"ouvre")==0) {Appliq_morpho(p,0,-1,tableau,x1,y1,x2,y2);
                             Appliq_morpho(p,1,1,tableau,x1,y1,x2,y2);}
else
if (strcmp(nom,"ferme")==0) {Appliq_morpho(p,1,-1,tableau,x1,y1,x2,y2);
                             Appliq_morpho(p,0,1,tableau,x1,y1,x2,y2);}
else
if (strcmp(nom,"grad")==0) {Image<int> q(p);
                            Appliq_morpho(q,0,-1,tableau,x1,y1,x2,y2);
                            Appliq_morpho(p,1,-1,tableau,x1,y1,x2,y2);
                            p = p.Image_diff(q);}
else
if (strcmp(nom,"lapl")==0) {Image<int> q(p);
                            Image<int> r(p);
                            Appliq_morpho(q,0,-1,tableau,x1,y1,x2,y2);
			    q = p.Image_diff(q);
                            Appliq_morpho(r,1,-1,tableau,x1,y1,x2,y2);
                            p = r.Image_diff(p);
			    p = p.Image_diffsignee(q);
                           }
else
if (strcmp(nom,"amincit")==0) {Image<int> q(p);
                               Appliq_morpho(q,0,-1,tableau,x1,y1,x2,y2);
			       p = p.Image_diff(q);}
else
if (strcmp(nom,"epaissit")==0) {Image<int> q(p);
                                Appliq_morpho(q,0,-1,tableau,x1,y1,x2,y2);
				p = p|q;}
else
if (strcmp(nom,"tophat")==0) {Image<int> q(p);
                              Appliq_morpho(q,0,-1,tableau,x1,y1,x2,y2); 
                              Appliq_morpho(q,1,1,tableau,x1,y1,x2,y2);
                              p = p.Image_diff(q);}  
else
if (strcmp(nom,"hattop")==0) {Image<int> q(p);
                              Appliq_morpho(p,1,-1,tableau,x1,y1,x2,y2);
                              Appliq_morpho(p,0,1,tableau,x1,y1,x2,y2);
                              p = p.Image_diff(q);} 
else
if (strcmp(nom,"opreco")==0) {Image<int> q(p);
                              Appliq_morpho(p,0,-1,tableau,x1,y1,x2,y2);
                              Appliq_morpho(p,1,1,tableau,x1,y1,x2,y2);
                              Reconstruit(p,q,cx);} 
else
if (strcmp(nom,"cloreco")==0) {Image<int> q(p);
                              Appliq_morpho(p,1,-1,tableau,x1,y1,x2,y2);
                              Appliq_morpho(p,0,1,tableau,x1,y1,x2,y2);
			      p = !p; q = !q;
                              Reconstruit(p,q,cx);
			      p = !p;
                              }   
else
if (strcmp(nom,"contrast")==0) {Image<int> q(p);
                                Image<int> r(p);
				int* PIX=p.PI();
				int* QPIX=q.PI();
				int* RPIX=r.PI();
				int L=p.PL();
				int W=p.PW();
				int index=0;
				Appliq_morpho(q,0,-1,tableau,x1,y1,x2,y2);
				Appliq_morpho(r,1,-1,tableau,x1,y1,x2,y2);
				for (int j=0;j<L;j++)
				  for (int i=0;i<W;i++) {
				    if ((PIX[index]-QPIX[index])<(RPIX[index]-PIX[index])) PIX[index]=QPIX[index]; else PIX[index]=RPIX[index];
				    index++;
				  }           
}
}

void Appliq_morpho(Image<int>& p,int erodil,int sens,int tab[NMAX][NMAX],int xdeb,int ydeb,int xfin,int yfin) {

  int index,valeur;
  Image<int> q(p);
  int* PIX=p.PI();
  int L=p.PL();
  int W=p.PW();

for (int j=0;j<L;j++)
    for (int i=0;i<W;i++) {
      index = j * W + i;
      for (int a = xdeb;a <= xfin;a++)
        for (int b = ydeb;b <= yfin;b++) {
            if (tab[a][b]) {
              valeur = q.X(i-sens*(a-NMAX/2),j-sens*(b-NMAX/2));
              if (tab[a][b]==2) valeur = 255-valeur;
              if ((a==xdeb)&&(b==ydeb)) PIX[index]=valeur;
                else {
                 if (erodil == 0)
                  PIX[index]=min(PIX[index],valeur);
                  else  PIX[index]=max(PIX[index],valeur); }
              }
            }
       }  
     }

void Morpho_48(Image<int>& p,int erodil,int conex) {

  int index,valeur;
  Image<int> q(p);
  int* PIX=p.PI();
  int L=p.PL();
  int W=p.PW();

  if (conex == 4) {
    for (int j=0;j<L;j++)
      for (int i=0;i<W;i++) {
	index = j * W + i;
                 if (erodil == 0)
                 PIX[index]=min(PIX[index],min(q.X(i-1,j),min(q.X(i,j-1),min(q.X(i+1,j),q.X(i,j+1)))));
		 else  PIX[index]=max(PIX[index],max(q.X(i-1,j),max(q.X(i,j-1),max(q.X(i+1,j),q.X(i,j+1)))));}
  } else {
    for (int a=0;a<=1;a++)
      for (int b=-1;b<=1;b+=2) {
	for (int j=0;j<L;j++)
	  for (int i=0;i<W;i++) {
	    index = j * W + i;
	    if (erodil == 0)
	      PIX[index]=min(PIX[index],q.X(i-a*b,j-(1-a)*b));
	    else PIX[index]=max(PIX[index],q.X(i-a*b,j-(1-a)*b));
	  }
	q=p;
      }
  }
}
      

void Reconstruit_seq(Image<int>& p,Image<int>& p_ref,int cx){
 
int w = p.PW();
int h = p.PL();
int* PIX=p.PI();
Image<int> p_tamp(w,h);
int pas_fini = 1;
int Nb,index;
 
while (pas_fini) {
p_tamp = p; 
for (int j=0;j<h;j++)
  for (int i=0;i<w;i++){
    index = j * w + i;
    if (cx == 4) {PIX[index]=min(p_ref.X(i,j),max(PIX[index],max(p.X(i-1,j),p.X(i,j-1))));}
    else {PIX[index]=min(p_ref.X(i,j),max(PIX[index],max(p.X(i-1,j-1),max(p.X(i,j-1),p.X(i-1,j)))));}
  }
for (int j=h-1;j>=0;j--)
  for (int i=w-1;i>=0;i--){
    index = j * w + i;
    if (cx == 4) {PIX[index]=min(p_ref.X(i,j),max(PIX[index],max(p.X(i+1,j),p.X(i,j+1))));}
    else {PIX[index]=min(p_ref.X(i,j),max(PIX[index],max(p.X(i+1,j+1),max(p.X(i,j+1),p.X(i+1,j)))));}
  }
p_tamp = p.Image_diff(p_tamp);    
Nb = p_tamp.Image_compte_pixel();
pas_fini = (Nb != 0);
     }
  }      
  
void Reconstruit(Image<int>& p,Image<int>& p_ref,int cx){
 
int w = p.PW();
int h = p.PL();
int* PIX=p.PI();
int index,newi,newj,newindex,initial,x,y;
FIFO f(w*h);

Max_regionaux(p,cx);

if (cx == 4) { 
 for (int j=0;j<h;j++)
   for (int i=0;i<w;i++){
     index = j * w + i;
     if ((PIX[index]) && !(p.X(i-1,j) && p.X(i,j-1) && p.X(i+1,j) && p.X(i,j+1))) f.FIFO_add(i,j);
     }
 while (!f.FIFO_empty()) {
      f.FIFO_get(x,y);
      for (int a=0;a<=1;a++)
        for (int b=-1;b<=1;b+=2) {
              newi = x-a*b; newj = y-(1-a)*b;
              if ((p.X(newi,newj) < p.X(x,y)) && (p_ref.X(newi,newj) != p.X(newi,newj)) && (newi>=0) && (newj>=0) && (newi<w) && (newj<h)) {
                   newindex = newj * w + newi;
                   PIX[newindex] = min(p.X(x,y),p_ref.X(newi,newj));
                   f.FIFO_add(newi,newj);}
               }
       }
 }
else {
for (int j=0;j<h;j++)
  for (int i=0;i<w;i++){
    index = j * w + i;
    initial=1;
    for (int a=-1;a<=1;a++)
      for (int b=-1;b<=1;b++)
        initial = initial && p.X(i-a,j-b);
    if (PIX[index] && !initial) f.FIFO_add(i,j);
    }
while (!f.FIFO_empty()) {
           f.FIFO_get(x,y);
           for (int a=-1;a<=1;a++)
            for (int b=-1;b<=1;b++) {
              if ((a!=0)||(b!=0)) {
                 newi = x-a; newj = y-b;
                 if ((p.X(newi,newj) < p.X(x,y)) && (p_ref.X(newi,newj) != p.X(newi,newj)) && (newi>=0) && (newj>=0) && (newi<w) && (newj<h)) {
                       newindex = newj * w + newi;
                       PIX[newindex] = min(p.X(x,y),p_ref.X(newi,newj));
                       f.FIFO_add(newi,newj);}
	          }
	     }
	}   
 }   
}
	                           
void FAS(Image<int>& p,int taille,int sens){
 
  for (int i = 1;i <= taille;i++)
    for (int opclo = 0;opclo <= 1;opclo++)
      for (int erdil = 0;erdil <= 1;erdil++) {
	int k=(sens!=(erdil!=opclo));
	for (int j=1;j <= i;j++)
	  if (j%2) Morpho_48(p,k,4); else Morpho_48(p,k,8);
      }
}  

void FGrain(Image<int>& p,int taille,int sens,int cx){

Image<int> q(p);

  for (int i = 1;i <= taille;i++)
    for (int opclo = 0;opclo <= 1;opclo++) {
      for (int erdil = 0;erdil <= 1;erdil++) {
	int k=(sens!=(erdil!=opclo));
	for (int j=1;j <= i;j++)
	  if (j%2) Morpho_48(p,k,4); else Morpho_48(p,k,8);
      }
      if (opclo) {Reconstruit(p,q,cx);} 
      else {p = !p;q = !q;Reconstruit(p,q,cx);p = !p;}
    }
}  

void Max_regionaux(Image<int>& p,int cx){
 
int w = p.PW();
int h = p.PL();
int* PIX=p.PI();
Image<int> q(p);
int index,newi,newj,indexnew,x,y,maxloc,newindex;
FIFO f(w*h);
 
if (cx == 4) { 
for (int j=0;j<h;j++)
  for (int i=0;i<w;i++){
    index = j * w + i;
    if ((PIX[index] < max(q.X(i-1,j),max(q.X(i,j-1),max(q.X(i+1,j),q.X(i,j+1))))) && (q.X(i,j))) {PIX[index]=0;f.FIFO_add(i,j);}
    while (!f.FIFO_empty()) {
                            f.FIFO_get(x,y);
                            for (int a=0;a<=1;a++)
                              for (int b=-1;b<=1;b+=2) {
                                 newi = x-a*b; newj = y-(1-a)*b;
                                 if ((q.X(newi,newj) <= q.X(x,y)) && (newi>=0) && (newj>=0) && (newi<w) && (newj<h)) {
                                  newindex = newj * w + newi;
                                  if (PIX[newindex]) {
                                 PIX[newindex] = 0;f.FIFO_add(newi,newj);}
                                 }
	                          }
	                     }
     }
   }
else {
  for (int j=0;j<h;j++)
  for (int i=0;i<w;i++){
    index = j * w + i;
    maxloc=1;
    for (int a=-1;a<=1;a++)
      for (int b=-1;b<=1;b++)
        maxloc = maxloc && (PIX[index] >= q.X(i-a,j-b));
    if ((!maxloc) && (q.X(i,j))) {PIX[index]=0;f.FIFO_add(i,j);}
    while (!f.FIFO_empty()) {
                            f.FIFO_get(x,y);
                            for (int a=-1;a<=1;a++)
                              for (int b=-1;b<=1;b++) {
                                if ((a!=0)||(b!=0)) {
                                   newi = x-a; newj = y-b;
                                   if ((q.X(newi,newj) <= q.X(x,y)) && (newi>=0) && (newj>=0) && (newi<w) && (newj<h)) {
                                   newindex = newj * w + newi;
                                   if (PIX[newindex]) {
                                            PIX[newindex] = 0;
                                            f.FIFO_add(newi,newj);}
                                  }
	                        }
	                      }
	                   }   
     }
   }   
 }    	   


void Nivellement(Image<int>& p,char *nom,int tab[NMAX][NMAX],float z,int taille,int cx) {

Image<int> p_ref1(p);
Image<int> p_ref2(p);

 if (strcmp(nom,"fas")==0) FAS(p,taille,0);
 else if (strcmp(nom,"median")==0) Filtre_Rang(p,tab,0.5);
 else FiltreMasque(p,nom,tab,z);

 p_ref1 = p_ref1 | p; 
 Reconstruit(p,p_ref1,cx);
 p = !p;
 p_ref2 = !p_ref2;
 p_ref2 = p_ref2 | p; 
 Reconstruit(p,p_ref2,12-cx);
 p = !p;
}

void MaxPropag(Image<int>& p,int conex){

int w = p.PW();
int h = p.PL();
Image<int> q(p);
Image<int> qtemp(p);
Image<int> pres(p);
int* Qpix=qtemp.PI();
int* Ppix=pres.PI();
int index,x,y,interior;
FIFO f(w*h);
FIFO g(w*h);
 
 pres = p;
 for (int j=0;j<h;j++)
   for (int i=0;i<w;i++){
     index = j*w + i;
     if (p.X(i,j)) {
       interior = p.X(i-1,j)&&p.X(i+1,j)&&p.X(i,j-1)&&p.X(i,j+1);
       if (conex == 8) interior = interior &&
			 p.X(i-1,j-1)&&p.X(i-1,j+1)&&p.X(i+1,j-1)&&p.X(i+1,j+1);
       if (!interior) {f.FIFO_add(i,j); Qpix[index]=0;Ppix[index]=0;}
     }
   }

 while (!f.FIFO_empty()) { 
   q = qtemp;
   while (!f.FIFO_empty()) {
     f.FIFO_get(x,y);
     if (qtemp.X(x-1,y)) {g.FIFO_add(x-1,y);Qpix[y*w+x-1]=0;Ppix[y*w+x-1]=0;}
     if (qtemp.X(x+1,y)) {g.FIFO_add(x+1,y);Qpix[y*w+x+1]=0;Ppix[y*w+x+1]=0;}
     if (qtemp.X(x,y-1)) {g.FIFO_add(x,y-1);Qpix[(y-1)*w+x]=0;Ppix[(y-1)*w+x]=0;}
     if (qtemp.X(x,y+1)) {g.FIFO_add(x,y+1);Qpix[(y+1)*w+x]=0;Ppix[(y+1)*w+x]=0;}
     if (conex == 8) {
       if (qtemp.X(x-1,y-1)) {g.FIFO_add(x-1,y-1);Qpix[(y-1)*w+x-1]=0;Ppix[(y-1)*w+x-1]=0;}
       if (qtemp.X(x+1,y-1)) {g.FIFO_add(x+1,y-1);Qpix[(y-1)*w+x+1]=0;Ppix[(y-1)*w+x+1]=0;}
       if (qtemp.X(x-1,y+1)) {g.FIFO_add(x-1,y+1);Qpix[(y+1)*w+x-1]=0;Ppix[(y+1)*w+x-1]=0;}
       if (qtemp.X(x+1,y+1)) {g.FIFO_add(x+1,y+1);Qpix[(y+1)*w+x+1]=0;Ppix[(y+1)*w+x+1]=0;}
     }
     interior = q.X(x-1,y)||q.X(x+1,y)||q.X(x,y-1)||q.X(x,y+1);
     if (conex != 4) interior = interior ||
		       q.X(x-1,y-1)||q.X(x+1,y-1)||q.X(x-1,y+1)||q.X(x+1,y+1);
     if (!interior) Ppix[y*w+x]=255;
   }
   f = g;
   g.FIFO_delete();
 }
 p = pres;
}

int PutDistance(Image<int>& p,int i,int j,char *nom,int sens) {
  if (sens == 1) {
    if (strcmp(nom,"d4")==0) return min(p.X(i-1,j),p.X(i,j-1))+1;
    if (strcmp(nom,"d8")==0) return min(p.X(i-1,j-1),min(p.X(i,j-1),min(p.X(i+1,j-1),p.X(i-1,j))))+1;
    if (strcmp(nom,"ch3-4")==0) return min(p.X(i-1,j-1)+4,min(p.X(i,j-1)+3,min(p.X(i+1,j-1)+4,p.X(i-1,j)+3)));
    if (strcmp(nom,"ch5-7-11")==0)
      return min(min(min(min(min(p.X(i-1,j),p.X(i,j-1))+5,
                             min(p.X(i-1,j-1),p.X(i+1,j-1))+7),
                         min(p.X(i-2,j),p.X(i,j-2))+10),
                     min(p.X(i-2,j-1),min(p.X(i-1,j-2),min(p.X(i+1,j-2),p.X(i+2,j-1))))+11),
                 min(p.X(i-2,j-2),p.X(i+2,j-2))+14);
    if (strcmp(nom,"ch14-20-31-44")==0) 
      return min(min(min(min(min(min(min(min(min(p.X(i-1,j),p.X(i,j-1))+14,
					     min(p.X(i-1,j-1),p.X(i+1,j-1))+20),
					 min(p.X(i-2,j),p.X(i,j-2))+28),
				     min(p.X(i-2,j-1),min(p.X(i-1,j-2),min(p.X(i+1,j-2),p.X(i+2,j-1))))+31),
				 min(p.X(i-2,j-2),p.X(i+2,j-2))+40),
			     min(p.X(i,j-3),p.X(i-3,j))+42),
			 min(p.X(i-3,j-1),min(p.X(i-1,j-3),min(p.X(i+1,j-3),p.X(i+3,j-1))))+44),
		     min(p.X(i-3,j-2),min(p.X(i-2,j-3),min(p.X(i+2,j-3),p.X(i+3,j-2))))+51),
		 min(p.X(i-3,j-3),p.X(i+3,j-3))+60);
  } else {
    if (strcmp(nom,"d4")==0) return min(p.X(i,j),(min(p.X(i+1,j),p.X(i,j+1))+1));
    if (strcmp(nom,"d8")==0) return min(p.X(i,j),min(p.X(i+1,j+1),min(p.X(i,j+1),min(p.X(i-1,j+1),p.X(i+1,j))))+1);
    if (strcmp(nom,"ch3-4")==0) return min(p.X(i,j),min(p.X(i+1,j+1)+4,min(p.X(i,j+1)+3,min(p.X(i-1,j+1)+4,p.X(i+1,j)+3))));
    if (strcmp(nom,"ch5-7-11")==0)
      return min(p.X(i,j),min(min(min(min(min(p.X(i+1,j),p.X(i,j+1))+5,
                                          min(p.X(i+1,j+1),p.X(i-1,j+1))+7),
                                      min(p.X(i+2,j),p.X(i,j+2))+10),
                                  min(p.X(i+2,j+1),min(p.X(i+1,j+2),min(p.X(i-1,j+2),p.X(i-2,j+1))))+11),
                              min(p.X(i+2,j+2),p.X(i-2,j+2))+14));
    if (strcmp(nom,"ch14-20-31-44")==0) 
      return min(p.X(i,j),min(min(min(min(min(min(min(min(min(p.X(i+1,j),p.X(i,j+1))+14,
							  min(p.X(i+1,j+1),p.X(i-1,j+1))+20),
						      min(p.X(i+2,j),p.X(i,j+2))+28),
						  min(p.X(i+2,j+1),min(p.X(i+1,j+2),min(p.X(i-1,j+2),p.X(i-2,j+1))))+31),
					      min(p.X(i+2,j+2),p.X(i-2,j+2))+40),
					  min(p.X(i,j+3),p.X(i+3,j))+42),
				      min(p.X(i+3,j+1),min(p.X(i+1,j+3),min(p.X(i-1,j+3),p.X(i-3,j+1))))+44),
				  min(p.X(i+3,j+2),min(p.X(i+2,j+3),min(p.X(i-2,j+3),p.X(i-3,j+2))))+51),
			      min(p.X(i+3,j+3),p.X(i-3,j+3))+60));
  }
}
 
void FctDistance(Image<int>& p,char *nom) {
  int h=p.PL();
  int w=p.PW();
  int i,j;
  int index = 0;
  int *PIX = p.PI();
  int valmax = 0;
 
  for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      if (*PIX > 0) *PIX = PutDistance(p,i,j,nom,1);
      PIX++;
        }             
  for (j=h-1;j>=0;j--)  
    for (i=w-1;i>=0;i--) { 
      PIX--;
      if (*PIX > 0) { 
        *PIX =  PutDistance(p,i,j,nom,-1); 
        if (*PIX > valmax) valmax = *PIX; 
         }
       }  
  // Affichage : réglage de la dynamique
  for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      if (*PIX > 0) *PIX = (*PIX * 255)/valmax;
      PIX++;
       }
} 

void Ero_ult(Image<int>& p,int conex) {
   int h=p.PL();
   int w=p.PW();
   int i,j;
   int *PIX = p.PI();
   int index = 0;
  char nom_dist[8];
   
   if (conex == 4) sprintf(nom_dist,"d4");    
        else sprintf(nom_dist,"d8");
   FctDistance(p,nom_dist);
   Max_regionaux(p,conex);
   for (j=0;j<h;j++)
     for (i=0;i<w;i++) {
       if (PIX[index] > 0) PIX[index] = 255;
       index++;
     }
 }

// LPE version 1 : LPE fine : frontière arbitraire
void LPE_version1(Image<int>& p,int conex) {

  int w = p.PW();
  int h = p.PL();
  Image<int> p_minima(p);
  Image<int> p_label(w,h);
  Image<int> p_done(w,h);
  Image<int> p_lpe(w,h);
  int* MIN=p_minima.PI();
  int* LAB=p_label.PI();
  int* DONE=p_done.PI();
  int* LPE=p_lpe.PI();
  int* PIX=p.PI();
  int index,newi,newj,x,y,minloc,newindex,indice,num;
  FIFO f(w*h);
  int i,j,a,b;
 
  printf("Connexité : %d\n",conex);
// (1) Calcul des minima régionaux
if (conex == 4) {
  index=0; 
  for (j=0;j<h;j++)
    for (i=0;i<w;i++){
      if (MIN[index] > min(p.X(i-1,j),min(p.X(i,j-1),min(p.X(i+1,j),p.X(i,j+1))))) {MIN[index]=506;f.FIFO_add(i,j);}
      while (!f.FIFO_empty()) {
	f.FIFO_get(x,y);
	for (a=0;a<=1;a++) {
	  for (b=-1;b<=1;b+=2) {
	    newi = x-a*b; newj = y-(1-a)*b;
	    if ((p.X(newi,newj) >= p.X(x,y)) && (newi>=0) && (newj>=0) && (newi<w) && (newj<h)) {
	      newindex = newj * w + newi;
	      if (MIN[newindex]!=506) {
		MIN[newindex] = 506;f.FIFO_add(newi,newj);}
	    }
	  }
	}
      }
      index++;
    }
} else {
  index=0;
  for (j=0;j<h;j++)
    for (i=0;i<w;i++){
      minloc=1;
      for (a=-1;a<=1;a++)
	for (b=-1;b<=1;b++)
	  minloc = minloc && (MIN[index] <= p.X(i-a,j-b));
      if (!minloc) {MIN[index]=506;f.FIFO_add(i,j);}
      while (!f.FIFO_empty()) {
	f.FIFO_get(x,y);
	for (a=-1;a<=1;a++)
	  for (b=-1;b<=1;b++) {
	    if ((a!=0)||(b!=0)) {
	      newi = x-a; newj = y-b;
	      if ((p.X(newi,newj) >= p.X(x,y)) && (newi>=0) && (newj>=0) && (newi<w) && (newj<h)) {
		newindex = newj * w + newi;
		if (MIN[newindex]!=506) {
		  MIN[newindex] = 506;
		  f.FIFO_add(newi,newj);}
	      }
	    }
	  }
      } 
      index++;
    }
}

// (2) Etiquettage des minima

 p_label = p_minima;
 p_done = p_minima;
 num = 0;
 if (conex == 4) {
  for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      if (p_done.X(i,j) != 506) {
	indice = 1;
	LAB[j*w+i] = ++num;
	DONE[j*w+i] = 506;
	f.FIFO_add(i,j);
	while (!f.FIFO_empty()) {
	  f.FIFO_get(x,y);
	  for (a=0;a<=1;a++) {
	    for (b=-1;b<=1;b+=2) {
	      newi = x-a*b; newj = y-(1-a)*b;
	      if ((p_done.X(newi,newj)!=506)&&(newi>=0)&&(newj>=0)&&(newi<w)&&(newj<h)) {
		  f.FIFO_add(newi,newj);
		  indice++;
		  index = newj*w+newi;
		  LAB[index]=num;
		  DONE[index]=506;
		}
	      }
	    }
	}
	printf("minimum n°%d : %d points\n",num,indice);
      }
    }
  } else {
    for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      if (p_done.X(i,j) != 506) {
	indice = 1;
	LAB[j*w+i] = ++num;
	DONE[j*w+i] = 506;
	f.FIFO_add(i,j);
	while (!f.FIFO_empty()) {
	  f.FIFO_get(x,y);
	  for (a=-1;a<=1;a++)
	    for (b=-1;b<=1;b++) {
	      if ((a!=0)||(b!=0)) {
		if ((p_done.X(x+a,y+b)!=506)&&(x+a>=0)&&(y+b>=0)&&(x+a<w)&&(y+b<h)) {
		  f.FIFO_add(x+a,y+b);
		  indice++;
		  index = (y+b)*w+x+a;
		  LAB[index]=num;
		  DONE[index]=506;
		}
	      }
	    }
	}
	printf("minimum n°%d : %d points\n",num,indice);
      }
    }
 }
   // Nombre de labels
  printf("%d labels\n",num);  
  
  // (3) Innondation à partir des minima

  num = 1;
 p_lpe = p_minima; 
if (conex == 4) {
  for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      if (p_label.X(i,j) == num) {
	f.FIFO_add(i,j);
	DONE[j*w+i]=507+num;
	while (!f.FIFO_empty()) {
	  f.FIFO_get(x,y);
	  for (a=0;a<=1;a++) {
	    for (b=-1;b<=1;b+=2) {
	      newi = x-a*b; newj = y-(1-a)*b;
	      if ((p_done.X(newi,newj)!=507+num)&&
		  (newi>=0)&&(newj>=0)&&(newi<w)&&(newj<h)&&
		  (p.X(newi,newj)>=p.X(x,y))) {
		indice = newj*w+newi;
		f.FIFO_add(newi,newj);
		MIN[indice] = p_minima.X(x,y);
		if (LAB[indice] == 506)
		  LAB[indice] = num;
		LPE[indice]= MIN[indice];
		DONE[indice] = 507+num;
	      }
	    }
	  }
	}
	num++;
      }
    }
} else {
   for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      if (p_label.X(i,j) == num) {
	f.FIFO_add(i,j);
	DONE[j*w+i]=507+num;
	while (!f.FIFO_empty()) {
	  f.FIFO_get(x,y);
	  for (a=-1;a<=1;a++) {
	    for (b=-1;b<=1;b++) {
	      if ((a!=0)||(b!=0)) {
		if ((p_done.X(x+a,y+b)!=507+num)&&
		    (x+a>=0)&&(y+b>=0)&&(x+a<w)&&(y+b<h)&&
		    (p.X(x+a,y+b)>=p.X(x,y))) {
		  indice = (y+b)*w+x+a;
		  f.FIFO_add(x+a,y+b);
		  MIN[indice] = p_minima.X(x,y);
		  if (LAB[indice] == 506)
		    LAB[indice] = num;
		  LPE[indice]= MIN[indice];
		  DONE[indice] = 507+num;
		}
	      }
	    }
	  }
	}
	num++;
      }
    }
}
  // Affichage des frontières
  indice=0;
 //  for (j=0;j<h;j++)
//     for (i=0;i<w;i++) {
//       LPE[indice] = 0;
//       if (LAB[indice]<p_label.X(i-1,j))
// 	LPE[indice] = min(PIX[indice]-MIN[indice],PIX[indice]-p_minima.X(i-1,j));
//       if (LAB[indice]<p_label.X(i+1,j))
// 	LPE[indice] = max(LPE[indice],min(PIX[indice]-MIN[indice],PIX[indice]-p_minima.X(i+1,j)));
//       if (LAB[indice]<p_label.X(i,j-1))
// 	LPE[indice] = max(LPE[indice],min(PIX[indice]-MIN[indice],PIX[indice]-p_minima.X(i,j-1)));
//       if (LAB[indice]<p_label.X(i,j+1)) 
// 	LPE[indice] = max(LPE[indice],min(PIX[indice]-MIN[indice],PIX[indice]-p_minima.X(i,j+1)));
//       indice++;
//     }

if (conex == 4) {
 for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      if ((LAB[indice]<p_label.X(i-1,j)) ||
	  (LAB[indice]<p_label.X(i+1,j)) ||
	  (LAB[indice]<p_label.X(i,j-1)) ||
	  (LAB[indice]<p_label.X(i,j+1)))
	LPE[indice] = 503;
      indice++;
    }
} else {
  for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      if ((LAB[indice]<p_label.X(i-1,j)) ||
	  (LAB[indice]<p_label.X(i+1,j)) ||
	  (LAB[indice]<p_label.X(i,j-1)) ||
	  (LAB[indice]<p_label.X(i,j+1)) ||
	  (LAB[indice]<p_label.X(i-1,j-1)) ||
	  (LAB[indice]<p_label.X(i+1,j-1)) ||
	  (LAB[indice]<p_label.X(i-1,j+1)) ||
	  (LAB[indice]<p_label.X(i+1,j+1)))
	LPE[indice] = 503;
      indice++;
    }
}
  
  p = p_lpe;
}

// LPE version 2 : LPE épaisse : pas de choix arbitraire 
// de frontière en cas de plateaux 
void LPE_version2(Image<int>& p,int conex) {

  int w = p.PW();
  int h = p.PL();
  Image<int> p_minima(p);
  Image<int> p_label(w,h);
  Image<int> p_done(w,h);
  Image<int> p_lpe(w,h);
  int* MIN=p_minima.PI();
  int* LAB=p_label.PI();
  int* DONE=p_done.PI();
  int* LPE=p_lpe.PI();
  int* PIX=p.PI();
  int index,newi,newj,x,y,minloc,newindex,indice,num;
  FIFO f(w*h);
  int i,j,a,b;
 
  printf("Connexité : %d\n",conex);
// (1) Calcul des minima régionaux
if (conex == 4) {
  index=0; 
  for (j=0;j<h;j++)
    for (i=0;i<w;i++){
      if (MIN[index] > min(p.X(i-1,j),min(p.X(i,j-1),min(p.X(i+1,j),p.X(i,j+1))))) {MIN[index]=506;f.FIFO_add(i,j);}
      while (!f.FIFO_empty()) {
	f.FIFO_get(x,y);
	for (a=0;a<=1;a++) {
	  for (b=-1;b<=1;b+=2) {
	    newi = x-a*b; newj = y-(1-a)*b;
	    if ((p.X(newi,newj) >= p.X(x,y)) && (newi>=0) && (newj>=0) && (newi<w) && (newj<h)) {
	      newindex = newj * w + newi;
	      if (MIN[newindex]!=506) {
		MIN[newindex] = 506;f.FIFO_add(newi,newj);}
	    }
	  }
	}
      }
      index++;
    }
} else {
  index=0;
  for (j=0;j<h;j++)
    for (i=0;i<w;i++){
      minloc=1;
      for (a=-1;a<=1;a++)
	for (b=-1;b<=1;b++)
	  minloc = minloc && (MIN[index] <= p.X(i-a,j-b));
      if (!minloc) {MIN[index]=506;f.FIFO_add(i,j);}
      while (!f.FIFO_empty()) {
	f.FIFO_get(x,y);
	for (a=-1;a<=1;a++)
	  for (b=-1;b<=1;b++) {
	    if ((a!=0)||(b!=0)) {
	      newi = x-a; newj = y-b;
	      if ((p.X(newi,newj) >= p.X(x,y)) && (newi>=0) && (newj>=0) && (newi<w) && (newj<h)) {
		newindex = newj * w + newi;
		if (MIN[newindex]!=506) {
		  MIN[newindex] = 506;
		  f.FIFO_add(newi,newj);}
	      }
	    }
	  }
      } 
      index++;
    }
}

// (2) Etiquettage des minima

 p_label = p_minima;
 p_done = p_minima;
 num = 0;
 if (conex == 4) {
  for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      if (p_done.X(i,j) != 506) {
	indice = 1;
	LAB[j*w+i] = ++num;
	DONE[j*w+i] = 506;
	f.FIFO_add(i,j);
	while (!f.FIFO_empty()) {
	  f.FIFO_get(x,y);
	  for (a=0;a<=1;a++) {
	    for (b=-1;b<=1;b+=2) {
	      newi = x-a*b; newj = y-(1-a)*b;
	      if ((p_done.X(newi,newj)!=506)&&(newi>=0)&&(newj>=0)&&(newi<w)&&(newj<h)) {
		  f.FIFO_add(newi,newj);
		  indice++;
		  index = newj*w+newi;
		  LAB[index]=num;
		  DONE[index]=506;
		}
	      }
	    }
	}
	printf("minimum n°%d : %d points\n",num,indice);
      }
    }
  } else {
    for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      if (p_done.X(i,j) != 506) {
	indice = 1;
	LAB[j*w+i] = ++num;
	DONE[j*w+i] = 506;
	f.FIFO_add(i,j);
	while (!f.FIFO_empty()) {
	  f.FIFO_get(x,y);
	  for (a=-1;a<=1;a++)
	    for (b=-1;b<=1;b++) {
	      if ((a!=0)||(b!=0)) {
		if ((p_done.X(x+a,y+b)!=506)&&(x+a>=0)&&(y+b>=0)&&(x+a<w)&&(y+b<h)) {
		  f.FIFO_add(x+a,y+b);
		  indice++;
		  index = (y+b)*w+x+a;
		  LAB[index]=num;
		  DONE[index]=506;
		}
	      }
	    }
	}
	printf("minimum n°%d : %d points\n",num,indice);
      }
    }
 }
   // Nombre de labels
  printf("%d labels\n",num);  
  
  // (3) Innondation à partir des minima

  num = 1;
 p_lpe = p_minima; 
if (conex == 4) {
  for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      if (p_label.X(i,j) == num) {
	f.FIFO_add(i,j);
	DONE[j*w+i]=507+num;
	while (!f.FIFO_empty()) {
	  f.FIFO_get(x,y);
	  for (a=0;a<=1;a++) {
	    for (b=-1;b<=1;b+=2) {
	      newi = x-a*b; newj = y-(1-a)*b;
	      if ((p_done.X(newi,newj)!=507+num)&&
		  (newi>=0)&&(newj>=0)&&(newi<w)&&(newj<h)&&
		  (p.X(newi,newj)>=p.X(x,y))) {
		indice = newj*w+newi;
		f.FIFO_add(newi,newj);
		MIN[indice] = p_minima.X(x,y);
		if (LAB[indice] == 506)
		  LAB[indice] = num;
		if (DONE[indice] > 507) LPE[indice]= 503; else LPE[indice]= MIN[indice];
		DONE[indice] = 507+num;
	      }
	    }
	  }
	}
	num++;
      }
    }
} else {
   for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      if (p_label.X(i,j) == num) {
	f.FIFO_add(i,j);
	DONE[j*w+i]=507+num;
	while (!f.FIFO_empty()) {
	  f.FIFO_get(x,y);
	  for (a=-1;a<=1;a++) {
	    for (b=-1;b<=1;b++) {
	      if ((a!=0)||(b!=0)) {
		if ((p_done.X(x+a,y+b)!=507+num)&&
		    (x+a>=0)&&(y+b>=0)&&(x+a<w)&&(y+b<h)&&
		    (p.X(x+a,y+b)>=p.X(x,y))) {
		  indice = (y+b)*w+x+a;
		  f.FIFO_add(x+a,y+b);
		  MIN[indice] = p_minima.X(x,y);
		  if (LAB[indice] == 506)
		    LAB[indice] = num;
		  if (DONE[indice] > 507) LPE[indice]= 503; else LPE[indice]= MIN[indice];
		  DONE[indice] = 507+num;
		}
	      }
	    }
	  }
	}
	num++;
      }
    }
}
  p = p_lpe;
}

// LPE version 3 : LPE fine : frontière =
// frontière intérieure de la LPE épaisse
 
void LPE_version3(Image<int>& p,int conex) {

  int w = p.PW();
  int h = p.PL();
  Image<int> p_minima(p);
  Image<int> p_label(w,h);
  Image<int> p_done(w,h);
  Image<int> p_lpe(w,h);
  int* MIN=p_minima.PI();
  int* LAB=p_label.PI();
  int* DONE=p_done.PI();
  int* LPE=p_lpe.PI();
  int* PIX=p.PI();
  int index,newi,newj,x,y,minloc,newindex,indice,num;
  FIFO f(w*h);
  int i,j,a,b;
 
  printf("Connexité : %d\n",conex);
// (1) Calcul des minima régionaux
if (conex == 4) {
  index=0; 
  for (j=0;j<h;j++)
    for (i=0;i<w;i++){
      if (MIN[index] > min(p.X(i-1,j),min(p.X(i,j-1),min(p.X(i+1,j),p.X(i,j+1))))) {MIN[index]=506;f.FIFO_add(i,j);}
      while (!f.FIFO_empty()) {
	f.FIFO_get(x,y);
	for (a=0;a<=1;a++) {
	  for (b=-1;b<=1;b+=2) {
	    newi = x-a*b; newj = y-(1-a)*b;
	    if ((p.X(newi,newj) >= p.X(x,y)) && (newi>=0) && (newj>=0) && (newi<w) && (newj<h)) {
	      newindex = newj * w + newi;
	      if (MIN[newindex]!=506) {
		MIN[newindex] = 506;f.FIFO_add(newi,newj);}
	    }
	  }
	}
      }
      index++;
    }
} else {
  index=0;
  for (j=0;j<h;j++)
    for (i=0;i<w;i++){
      minloc=1;
      for (a=-1;a<=1;a++)
	for (b=-1;b<=1;b++)
	  minloc = minloc && (MIN[index] <= p.X(i-a,j-b));
      if (!minloc) {MIN[index]=506;f.FIFO_add(i,j);}
      while (!f.FIFO_empty()) {
	f.FIFO_get(x,y);
	for (a=-1;a<=1;a++)
	  for (b=-1;b<=1;b++) {
	    if ((a!=0)||(b!=0)) {
	      newi = x-a; newj = y-b;
	      if ((p.X(newi,newj) >= p.X(x,y)) && (newi>=0) && (newj>=0) && (newi<w) && (newj<h)) {
		newindex = newj * w + newi;
		if (MIN[newindex]!=506) {
		  MIN[newindex] = 506;
		  f.FIFO_add(newi,newj);}
	      }
	    }
	  }
      } 
      index++;
    }
}

// (2) Etiquettage des minima

 p_label = p_minima;
 p_done = p_minima;
 num = 0;
 if (conex == 4) {
  for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      if (p_done.X(i,j) != 506) {
	indice = 1;
	LAB[j*w+i] = ++num;
	DONE[j*w+i] = 506;
	f.FIFO_add(i,j);
	while (!f.FIFO_empty()) {
	  f.FIFO_get(x,y);
	  for (a=0;a<=1;a++) {
	    for (b=-1;b<=1;b+=2) {
	      newi = x-a*b; newj = y-(1-a)*b;
	      if ((p_done.X(newi,newj)!=506)&&(newi>=0)&&(newj>=0)&&(newi<w)&&(newj<h)) {
		  f.FIFO_add(newi,newj);
		  indice++;
		  index = newj*w+newi;
		  LAB[index]=num;
		  DONE[index]=506;
		}
	      }
	    }
	}
	printf("minimum n°%d : %d points\n",num,indice);
      }
    }
  } else {
    for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      if (p_done.X(i,j) != 506) {
	indice = 1;
	LAB[j*w+i] = ++num;
	DONE[j*w+i] = 506;
	f.FIFO_add(i,j);
	while (!f.FIFO_empty()) {
	  f.FIFO_get(x,y);
	  for (a=-1;a<=1;a++)
	    for (b=-1;b<=1;b++) {
	      if ((a!=0)||(b!=0)) {
		if ((p_done.X(x+a,y+b)!=506)&&(x+a>=0)&&(y+b>=0)&&(x+a<w)&&(y+b<h)) {
		  f.FIFO_add(x+a,y+b);
		  indice++;
		  index = (y+b)*w+x+a;
		  LAB[index]=num;
		  DONE[index]=506;
		}
	      }
	    }
	}
	printf("minimum n°%d : %d points\n",num,indice);
      }
    }
 }
   // Nombre de labels
  printf("%d labels\n",num);  
  
  // (3) Innondation à partir des minima

  num = 1;
 p_lpe = p_minima; 
if (conex == 4) {
  for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      if (p_label.X(i,j) == num) {
	f.FIFO_add(i,j);
	DONE[j*w+i]=507+num;
	while (!f.FIFO_empty()) {
	  f.FIFO_get(x,y);
	  for (a=0;a<=1;a++) {
	    for (b=-1;b<=1;b+=2) {
	      newi = x-a*b; newj = y-(1-a)*b;
	      if ((p_done.X(newi,newj)!=507+num)&&
		  (newi>=0)&&(newj>=0)&&(newi<w)&&(newj<h)&&
		  (p.X(newi,newj)>=p.X(x,y))) {
		indice = newj*w+newi;
		f.FIFO_add(newi,newj);
		MIN[indice] = p_minima.X(x,y);
		if (LAB[indice] == 506)
		  LAB[indice] = num;
		if (DONE[indice] > 507) LPE[indice]= 503; else LPE[indice]= MIN[indice];
		DONE[indice] = 507+num;
	      }
	    }
	  }
	}
	num++;
      }
    }
} else {
   for (j=0;j<h;j++)
    for (i=0;i<w;i++) {
      if (p_label.X(i,j) == num) {
	f.FIFO_add(i,j);
	DONE[j*w+i]=507+num;
	while (!f.FIFO_empty()) {
	  f.FIFO_get(x,y);
	  for (a=-1;a<=1;a++) {
	    for (b=-1;b<=1;b++) {
	      if ((a!=0)||(b!=0)) {
		if ((p_done.X(x+a,y+b)!=507+num)&&
		    (x+a>=0)&&(y+b>=0)&&(x+a<w)&&(y+b<h)&&
		    (p.X(x+a,y+b)>=p.X(x,y))) {
		  indice = (y+b)*w+x+a;
		  f.FIFO_add(x+a,y+b);
		  MIN[indice] = p_minima.X(x,y);
		  if (LAB[indice] == 506)
		    LAB[indice] = num;
		  if (DONE[indice] > 507) LPE[indice]= 503; else LPE[indice]= MIN[indice];
		  DONE[indice] = 507+num;
		}
	      }
	    }
	  }
	}
	num++;
      }
    }
}
  // Calcul des frontières
 p_done = p_lpe; 
 indice=0;
 if (conex == 4) {
   for (j=0;j<h;j++)
     for (i=0;i<w;i++) {
       if ((p_done.X(i-1,j)==503)&&(p_done.X(i+1,j)==503)&&(p_done.X(i,j-1)==503)&&(p_done.X(i,j+1)==503))
	 LPE[indice] = MIN[indice];
       indice++;
     }
 } else {
   for (j=0;j<h;j++)
     for (i=0;i<w;i++) {
       if ((p_done.X(i-1,j)==503)&&(p_done.X(i+1,j)==503)&&
	   (p_done.X(i,j-1)==503)&&(p_done.X(i,j+1)==503)&&
	   (p_done.X(i-1,j-1)==503)&&(p_done.X(i-1,j+1)==503)&&
	   (p_done.X(i+1,j-1)==503)&&(p_done.X(i+1,j+1)==503))
	 LPE[indice] = MIN[indice];
       indice++;
     }
 }
  
  p = p_lpe;
}

void LPE(Image<int>& p,int conex) {
  LPE_version1(p,conex);
}

void Filtrage_dynamique(Image<int>& p,int hauteur,int conex) {
   int h=p.PL();
   int w=p.PW();
   int i,j;
   Image<int> q(p);
   int *QPIX = q.PI();
   int index = 0;
   
   
   for (j=0;j<h;j++)
     for (i=0;i<w;i++) 
       QPIX[index++] += hauteur;
   Reconstruit(p,q,conex);
}

