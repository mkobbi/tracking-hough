#ifndef _HOUGHTRACK_
#define _HOUGHTRACK_

//Update routine of the integer GHT from the double GHT (for display)
void MaJ_Hough2d(Image<int>& p_transf,double *HOUGH,double vote_max);
//Reset routine of the GHT
void Reset_Hough2d(double *HOUGH,int width_img,int height_img);
//Drawing rectangle routine in the image space
void Trace_rectangle(Image<int>& p,int x_centre,int y_centre,int largeur,int hauteur,int couleur);
//Data structure for a voting instruction of the GHT (= a box of the R-Table)
typedef struct deplacement
{
  int dx;
  int dy;
  float poids;
  struct deplacement *next;
} deplacement_s;
//Allocating the R-table (size: nunber of distinct indexes)
void Allocate_RTable(deplacement_s*** RTABLE,int taille);
//Insertiing a vote instruction (displacement) within the R-Table 
void Insert_RTable(deplacement_s** RTABLE,int index,int delta_x,int delta_y,float omega);
//Creating a R-Table from a unique prototype image (p_proto) 
void Create_RTable_ordre0(deplacement_s** RTABLE,Image<int>& p_proto,Image<int>& p_valeur,Image<int>& p_poids,float sigma,int nb_labels);
void Create_RTable_ordre1(deplacement_s** RTABLE,Image<int>& p_proto,Image<int>& p_direction,Image<int>& p_poids,float sigma,int nb_labels);
void Create_RTable_ordre2(deplacement_s** RTABLE,Image<int>& p_proto,Image<int>& p_courbure,Image<int>& p_poids,float sigma,int nb_labels);
//Generalised Hough Transforms
void Hough_General_ordre1(deplacement_s** RTABLE,double *HOUGH,Image<int>& p_img,float sigma,int nb_labels,double *vote_max);
void Hough_General_ordre2(deplacement_s** RTABLE,double *HOUGH,Image<int>& p_img,float sigma,int nb_labels,double *vote_max);
void Hough_General_ordre0(deplacement_s** RTABLE,double *HOUGH,Image<int>& p_img,float sigma,int nb_labels,double *vote_max);
//Calculate the "nbest" best detections from the GHT "HOUGH"
void Find_Best_GHT(double *HOUGH,int width_img,int height_img,int nbest,int *Best_Pos,int taille_x,int taille_y);
//Update the model according to the best detections (BestPos) and the previous position (old_position)
void Update_Proto_RTable(deplacement_s** RTABLE,Image<int>& p_proto,Image<int>& p_index,Image<int>& p_poids,
                         Image<int>& p_courante,int *BestPos,int old_position_x,int old_position_y,
                         int *new_position_x,int *new_position_y,int nb_labels,float sigma);

#endif
