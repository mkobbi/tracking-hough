#include "Image.hh"
#include "filtre.hh"
#include "morfo.hh"
#include "feature.hh"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

inline int puissance(int n, int p) { if (!p) return 1; else return n * puissance(n, p - 1); }

inline double valabs(double x) { if (x < 0) return -x; else return x; }

void Hessian_measures(Image<int> &p, float sigma, int type_measure) {
    int h = p.PL();
    int w = p.PW();
    int i, j, val;
    int index = 0;
    int *PIX = p.PI();
    double *VAL, *TEMP, *GXX, *GYY, *GXY;
    double valmin = 10000;
    double valmax = -10000;
    double sum_diag, diff_diag, delta;

    VAL = (double *) calloc(w * h, sizeof(double));
    TEMP = (double *) calloc(w * h, sizeof(double));
    GXX = (double *) calloc(w * h, sizeof(double));
    GYY = (double *) calloc(w * h, sizeof(double));
    GXY = (double *) calloc(w * h, sizeof(double));
    // Passage en double...
    for (index = 0; index < w * h; index++)
        VAL[index] = (double) (PIX[index]);
    // Calcul des composantes de la matrice Hessienne
    Gauss_Rec(VAL, TEMP, w, h, 0, 1, sigma, valmin, valmax);
    Gauss_Rec(TEMP, GXX, w, h, 2, 0, sigma, valmin, valmax);
    Gauss_Rec(VAL, TEMP, w, h, 0, 0, sigma, valmin, valmax);
    Gauss_Rec(TEMP, GYY, w, h, 2, 1, sigma, valmin, valmax);
    Gauss_Rec(VAL, TEMP, w, h, 1, 0, sigma, valmin, valmax);
    Gauss_Rec(TEMP, GXY, w, h, 1, 1, sigma, valmin, valmax);
    valmin = 10000;
    valmax = -10000;
    if (type_measure == 0) {// Trace (Laplacien)
        for (index = 0; index < w * h; index++) {
            VAL[index] = GXX[index] + GYY[index];
            if (VAL[index] > valmax) valmax = VAL[index]; else if (VAL[index] < valmin) valmin = VAL[index];
        }
    } else if (type_measure == 1) {// Déterminant (Fonction d'intérêt)
        for (index = 0; index < w * h; index++) {
            VAL[index] = GXX[index] * GYY[index] - GXY[index] * GXY[index];
            if (VAL[index] > valmax) valmax = VAL[index]; else if (VAL[index] < valmin) valmin = VAL[index];
        }
    } else if (type_measure == 2) {// Norme de Frobénius ("Unflatness")
        for (index = 0; index < w * h; index++) {
            VAL[index] = sqrt(GXX[index] * GXX[index] + GYY[index] * GYY[index] + 2 * GXY[index] * GXY[index]);
            if (VAL[index] > valmax) valmax = VAL[index];
        }
    } else if (type_measure == 3) {// Plus grande valeur propre (Courbure principale)
        for (index = 0; index < w * h; index++) {
            sum_diag = GXX[index] + GYY[index];
            diff_diag = GXX[index] - GYY[index];
            delta = sqrt(4 * GXY[index] * GXY[index] + diff_diag * diff_diag);
            VAL[index] = 0.5 * (sum_diag + delta);
            if (VAL[index] > valmax) valmax = VAL[index]; else if (VAL[index] < valmin) valmin = VAL[index];
        }
    } else if (type_measure == 4) {// Plus petite valeur propre (Courbure secondaire)
        for (index = 0; index < w * h; index++) {
            sum_diag = GXX[index] + GYY[index];
            diff_diag = GXX[index] - GYY[index];
            delta = sqrt(4 * GXY[index] * GXY[index] + diff_diag * diff_diag);
            VAL[index] = 0.5 * (sum_diag - delta);
            if (VAL[index] > valmax) valmax = VAL[index]; else if (VAL[index] < valmin) valmin = VAL[index];
        }
    }
    // Repassage en entier avec équilibrage des niveaux de gris...
    if (type_measure == 2) {// Valeurs positives
        for (index = 0; index < w * h; index++)
            PIX[index] = (int) ((VAL[index] * 255.0) / valmax);
    } else {// Valeurs signées
        for (index = 0; index < w * h; index++) {
            if (VAL[index] > 0)
                PIX[index] = 128 + (int) ((VAL[index] * 127.0) / valmax);
            else
                PIX[index] = 128 - (int) ((VAL[index] * 127.0) / valmin);
        }
    }
    free(GXY);
    free(GXX);
    free(GYY);
    free(TEMP);
}

void Harris(Image<int> &p, float sigma1, float sigma2, int niveau_seuil, int display) {
    int h = p.PL();
    int w = p.PW();
    int i, j, val;
    int index = 0;
    int *PIX = p.PI();
    double *VAL, *TEMP, *GX, *GY, *GXGY;
    double valmin = 10000;
    double valmax = -10000;
    double seuil;

    double alpha = 0.06;

    VAL = (double *) calloc(w * h, sizeof(double));
    GX = (double *) calloc(w * h, sizeof(double));
    GY = (double *) calloc(w * h, sizeof(double));
    GXGY = (double *) calloc(w * h, sizeof(double));
    TEMP = (double *) calloc(w * h, sizeof(double));
    // Passage en double...
    index = 0;
    for (i = 0; i < w; i++)
        for (j = 0; j < h; j++)
            VAL[index] = (double) (PIX[index++]);
    // Calcul des dérivées spatiales
    Gauss_Rec(VAL, TEMP, w, h, 0, 1, sigma1, valmin, valmax);
    Gauss_Rec(TEMP, GX, w, h, 1, 0, sigma1, valmin, valmax);
    Gauss_Rec(VAL, TEMP, w, h, 0, 0, sigma1, valmin, valmax);
    Gauss_Rec(TEMP, GY, w, h, 1, 1, sigma1, valmin, valmax);
    // Calcul des termes de la matrice d'autocorrelation
    index = 0;
    for (i = 0; i < w; i++)
        for (j = 0; j < h; j++) {
            GXGY[index] = GX[index] * GY[index];
            GX[index] = GX[index] * GX[index];
            GY[index] = GY[index] * GY[index];
            index++;
        }
    // Calcul des moyennes locales des dérivées
    Gauss_Rec(GX, TEMP, w, h, 0, 1, sigma2, valmin, valmax);
    Gauss_Rec(TEMP, GX, w, h, 0, 0, sigma2, valmin, valmax);
    Gauss_Rec(GY, TEMP, w, h, 0, 1, sigma2, valmin, valmax);
    Gauss_Rec(TEMP, GY, w, h, 0, 0, sigma2, valmin, valmax);
    Gauss_Rec(GXGY, TEMP, w, h, 0, 1, sigma2, valmin, valmax);
    Gauss_Rec(TEMP, GXGY, w, h, 0, 0, sigma2, valmin, valmax);
    // Calcul de la fonction d'interet (det(A)-alpha*trace(A)*trace(A))
    valmin = 10000;
    valmax = -10000;
    index = 0;
    for (i = 0; i < w; i++)
        for (j = 0; j < h; j++) {
            VAL[index] = GX[index] * GY[index] - GXGY[index] * GXGY[index]
                         - alpha * (GX[index] + GY[index]) * (GX[index] + GY[index]);
            if (VAL[index] > valmax) valmax = VAL[index]; else if (VAL[index] < valmin) valmin = VAL[index];
            index++;
        }
    if (display == 0) {
        // Repassage en entier avec équilibrage des niveaux de gris...
        index = 0;
        for (i = 0; i < w; i++)
            for (j = 0; j < h; j++) {
                if (VAL[index] > 0)
                    PIX[index] = 128 + (int) ((VAL[index] * 127) / valmax);
                else
                    PIX[index] = 128 - (int) ((VAL[index] * 127) / valmin);
                index++;
            }
    } else {
        // Calcul des maxima locaux et seuillage
        index = 0;
        seuil = 0.5 * niveau_seuil * valmax / 100;
        if (display == 2) { Gaussian(p, 0, 2, sigma1); }
        for (i = 1; i < w - 1; i++)
            for (j = 1; j < h - 1; j++) {
                index = j * w + i;
                if (VAL[index] > seuil) {
                    if ((VAL[index] >= VAL[index + 1]) &&
                        (VAL[index] >= VAL[index - 1]) &&
                        (VAL[index] >= VAL[index + w]) &&
                        (VAL[index] >= VAL[index + 1]) &&
                        (VAL[index] >= VAL[index - 1 - w]) &&
                        (VAL[index] >= VAL[index - 1 + w]) &&
                        (VAL[index] >= VAL[index + 1 - w]) &&
                        (VAL[index] >= VAL[index + 1 + w])) {
                        PIX[index] = 503;
                        PIX[index - w + 1] = 503;
                        PIX[index + w + 1] = 503;
                        PIX[index - w - 1] = 503;
                        PIX[index + w - 1] = 503;
                    }
                }
            }
    }
    free(GX);
    free(GXGY);
    free(GY);
    free(TEMP);
    free(VAL);
}

void Morpho_interest(Image<int> &p, float sigma, int niveau_seuil, int display) {
    int h = p.PL();
    int w = p.PW();
    int i, j, val;
    int index = 0;
    int *PIX = p.PI();
    double *VAL, *GX, *GY, *TEMP, *TEMP1, *GXX, *GYY, *GXY;
    double valmin = 10000;
    double valmax = -10000;
    double seuil;

    VAL = (double *) calloc(w * h, sizeof(double));
    TEMP = (double *) calloc(w * h, sizeof(double));
    TEMP1 = (double *) calloc(w * h, sizeof(double));
    GX = (double *) calloc(w * h, sizeof(double));
    GY = (double *) calloc(w * h, sizeof(double));
    GXX = (double *) calloc(w * h, sizeof(double));
    GYY = (double *) calloc(w * h, sizeof(double));
    GXY = (double *) calloc(w * h, sizeof(double));
    // Passage en double...
    index = 0;
    for (i = 0; i < w; i++)
        for (j = 0; j < h; j++)
            VAL[index] = (double) (PIX[index++]);

    Gauss_Rec(VAL, TEMP, w, h, 0, 1, sigma, valmin, valmax);
    Gauss_Rec(TEMP, GX, w, h, 1, 0, sigma, valmin, valmax);
    Gauss_Rec(TEMP, GXX, w, h, 2, 0, sigma, valmin, valmax);
    Gauss_Rec(VAL, TEMP, w, h, 0, 0, sigma, valmin, valmax);
    Gauss_Rec(TEMP, GY, w, h, 1, 1, sigma, valmin, valmax);
    Gauss_Rec(TEMP, GYY, w, h, 2, 1, sigma, valmin, valmax);
    Gauss_Rec(VAL, TEMP, w, h, 0, 1, sigma, valmin, valmax);
    Gauss_Rec(TEMP, TEMP1, w, h, 1, 0, sigma, valmin, valmax);
    Gauss_Rec(TEMP1, TEMP, w, h, 0, 0, sigma, valmin, valmax);
    Gauss_Rec(TEMP, GXY, w, h, 1, 1, sigma, valmin, valmax);
    valmin = 10000;
    valmax = -10000;
    index = 0;
    for (i = 0; i < w; i++)
        for (j = 0; j < h; j++) {
            VAL[index] = (GXX[index] * GY[index] * GY[index] - 2 * GXY[index] * GX[index] * GY[index] +
                          GYY[index] * GX[index] * GX[index])
                         / (GX[index] * GX[index] + GY[index] * GY[index]);
            if (VAL[index] > valmax) valmax = VAL[index]; else if (VAL[index] < valmin) valmin = VAL[index];
            index++;
        }
    if (display == 0) {
        // Repassage en entier avec équilibrage des niveaux de gris...
        index = 0;
        for (i = 0; i < w; i++)
            for (j = 0; j < h; j++) {
                if (VAL[index] > 0)
                    PIX[index] = 128 + (int) ((VAL[index] * 127) / valmax);
                else
                    PIX[index] = 128 - (int) ((VAL[index] * 127) / valmin);
                index++;
            }
    } else {
        // Calcul des maxima locaux et seuillage
        index = 0;
        seuil = (valmax - valmin) * niveau_seuil / 100;
        if (display == 2) { Gaussian(p, 0, 2, sigma); }
        for (i = 0; i < w; i++)
            for (j = 0; j < h - 1; j++) {
                VAL[index] = valabs(VAL[index]);
                index++;
            }
        for (i = 1; i < w - 1; i++)
            for (j = 1; j < h - 1; j++) {
                index = j * w + i;
                if (VAL[index] > seuil) {
                    if ((VAL[index] >= VAL[index + 1]) &&
                        (VAL[index] >= VAL[index - 1]) &&
                        (VAL[index] >= VAL[index + w]) &&
                        (VAL[index] >= VAL[index + 1]) &&
                        (VAL[index] >= VAL[index - 1 - w]) &&
                        (VAL[index] >= VAL[index - 1 + w]) &&
                        (VAL[index] >= VAL[index + 1 - w]) &&
                        (VAL[index] >= VAL[index + 1 + w])) {
                        PIX[index] = 503;
                        PIX[index - w + 1] = 503;
                        PIX[index + w + 1] = 503;
                        PIX[index - w - 1] = 503;
                        PIX[index + w - 1] = 503;
                    }
                }
            }
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

void Calcul_LJ(double *Gin, int ordre, float sigma, int w, int h, double &valmin, double &valmax, double *Gout,
               int normalize) {

    double *TEMP, *TEMP1;
    int i, j;
    int index = 0;

    TEMP = (double *) calloc(w * h, sizeof(double));

    switch (ordre) {
        case 0 : {
            // Calcul de la gaussienne
            Gauss_Rec(Gin, TEMP, w, h, 0, 0, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gout, w, h, 0, 1, sigma, valmin, valmax);
            break;
        }
        case 1 : {
            // Calcul de G_x
            Gauss_Rec(Gin, TEMP, w, h, 0, 1, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gout, w, h, 1, 0, sigma, valmin, valmax);
            break;
        }
        case 2 : {
            // Calcul de G_y
            Gauss_Rec(Gin, TEMP, w, h, 0, 0, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gout, w, h, 1, 1, sigma, valmin, valmax);
            break;
        }
        case 3 : {
            // Calcul de G_xx
            Gauss_Rec(Gin, TEMP, w, h, 0, 1, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gout, w, h, 2, 0, sigma, valmin, valmax);
            break;
        }
        case 4 : {
            // Calcul de G_yy
            Gauss_Rec(Gin, TEMP, w, h, 0, 0, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gout, w, h, 2, 1, sigma, valmin, valmax);
            break;
        }
        case 5 : {
            TEMP1 = (double *) calloc(w * h, sizeof(double));
            // Calcul de G_xy
            Gauss_Rec(Gin, TEMP, w, h, 0, 1, sigma, valmin, valmax);
            Gauss_Rec(TEMP, TEMP1, w, h, 1, 0, sigma, valmin, valmax);
            Gauss_Rec(TEMP1, Gout, w, h, 1, 1, sigma, valmin, valmax);
            free(TEMP1);
            break;
        }
    }
    // Normalisation
    if (normalize == 1) {
        printf("sigma = %f - derivee n° %d : Min = %f, Max = %f\n", sigma, ordre, valmin, valmax);
        // Basee histogramme
        for (i = 0; i < w; i++)
            for (j = 0; j < h; j++) {
                Gout[index] = (int) ((Gout[index] - valmin) * 255.0 / (valmax - valmin));
                index++;
            }
    } else if (normalize == 2) {
        // Basee scale-space
        float scale_factor;
        if (ordre > 2) scale_factor = sigma * sigma;
        else if (ordre > 0) scale_factor = sigma;
        else scale_factor = 1.0;
        for (i = 0; i < w; i++)
            for (j = 0; j < h; j++)
                Gout[index++] *= scale_factor;
    }
    free(TEMP);
}

void
Compute_LJ_Features(Image<int> &p_in, int ordre_derivation, int nb_echelles, float sigma_initial, double **FEAT_out) {

    double *G;
    int i, j, k;
    int index = 0;
    int *PIX = p_in.PI();
    int HEIGHT = p_in.PL();
    int WIDTH = p_in.PW();
    double valmin, valmax;
    int nb_derivees = (ordre_derivation + 1) * (ordre_derivation + 2) / 2;

    G = (double *) calloc(WIDTH * HEIGHT, sizeof(double));

    //  Passage en double...
    for (i = 0; i < WIDTH; i++)
        for (j = 0; j < HEIGHT; j++)
            G[index] = (double) (PIX[index++]);

    float sigma = sigma_initial;
    k = 0;
    for (i = 0; i < nb_echelles; i++) {
        for (j = 0; j < nb_derivees; j++)
            Calcul_LJ(G, j, sigma, WIDTH, HEIGHT, valmin, valmax, FEAT_out[k++], 2);
        sigma = 2 * sigma;
    }

    printf("Dimension des local jets : %d\n", nb_echelles * nb_derivees);
    free(G);
}


float dist_LJ(double **FEATURES, int index1, int index2, int ordre_derivation, int nb_echelles) {
    float dist, diff;
    int i;
    float sum_sq = 0;
    int nb_derivees = (ordre_derivation + 1) * (ordre_derivation + 2) / 2;

    for (i = 0; i < nb_echelles; i++) {
        diff = FEATURES[nb_derivees * i][index1] - FEATURES[nb_derivees * i][index2];
        sum_sq += 6 * diff * diff;
        if (ordre_derivation > 0) {
            diff = FEATURES[nb_derivees * i + 1][index1] - FEATURES[nb_derivees * i + 1][index2];
            sum_sq += 3 * diff * diff;
            diff = FEATURES[nb_derivees * i + 2][index1] - FEATURES[nb_derivees * i + 2][index2];
            sum_sq += 3 * diff * diff;
            if (ordre_derivation > 1) {
                diff = FEATURES[nb_derivees * i + 3][index1] - FEATURES[nb_derivees * i + 3][index2];
                sum_sq += 2 * diff * diff;
                diff = FEATURES[nb_derivees * i + 4][index1] - FEATURES[nb_derivees * i + 4][index2];
                sum_sq += 2 * diff * diff;
                diff = FEATURES[nb_derivees * i + 5][index1] - FEATURES[nb_derivees * i + 5][index2];
                sum_sq += 2 * diff * diff;
            }
        }
    }

    dist = sum_sq / (6 * (ordre_derivation + 1) * nb_echelles);
    return dist;
}


float dist_vect_LJ(double *vecteur, double **FEATURES, int index, int ordre_derivation, int nb_echelles) {
    float dist, diff;
    int i;
    float sum_sq = 0;
    int nb_derivees = (ordre_derivation + 1) * (ordre_derivation + 2) / 2;

    for (i = 0; i < nb_echelles; i++) {
        diff = vecteur[nb_derivees * i] - FEATURES[nb_derivees * i][index];
        sum_sq += 6 * diff * diff;
        if (ordre_derivation > 0) {
            diff = vecteur[nb_derivees * i + 1] - FEATURES[nb_derivees * i + 1][index];
            sum_sq += 3 * diff * diff;
            diff = vecteur[nb_derivees * i + 2] - FEATURES[nb_derivees * i + 2][index];
            sum_sq += 3 * diff * diff;
            if (ordre_derivation > 1) {
                diff = vecteur[nb_derivees * i + 3] - FEATURES[nb_derivees * i + 3][index];
                sum_sq += 2 * diff * diff;
                diff = vecteur[nb_derivees * i + 4] - FEATURES[nb_derivees * i + 4][index];
                sum_sq += 2 * diff * diff;
                diff = vecteur[nb_derivees * i + 5] - FEATURES[nb_derivees * i + 5][index];
                sum_sq += 2 * diff * diff;
            }
        }
    }

    dist = sum_sq / (6 * (ordre_derivation + 1) * nb_echelles);
    return dist;
}

// Calcul des Composantes dans le Repère Local (g,t)
// (Pour Invariance en Rotation)
void Calcul_RI_LJ(double *Gin, int ordre, int produit_module, float sigma, int w, int h, double &valmin, double &valmax,
                  double *Gout, int normalize) {

    double *TEMP, *Gx, *Gy, *Gxx, *Gxy, *Gyy, *TEMP2;
    int i, j;
    int index = 0;
    double square_module;

    TEMP = (double *) calloc(w * h, sizeof(double));

    switch (ordre) {
        case 0 : {
            // Calcul de la gaussienne
            Gauss_Rec(Gin, TEMP, w, h, 0, 0, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gout, w, h, 0, 1, sigma, valmin, valmax);
            break;
        }
        case 1 : {
            Gx = (double *) calloc(w * h, sizeof(double));
            Gy = (double *) calloc(w * h, sizeof(double));
            // Calcul de G_g (i.e. module du gradient)
            // On calcule G_x
            Gauss_Rec(Gin, TEMP, w, h, 0, 1, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gx, w, h, 1, 0, sigma, valmin, valmax);
            // puis G_y
            Gauss_Rec(Gin, TEMP, w, h, 0, 0, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gy, w, h, 1, 1, sigma, valmin, valmax);
            // puis le module euclidien
            valmax = 0;
            for (index = 0; index < w * h; index++) {
                Gout[index] = sqrt(Gx[index] * Gx[index] + Gy[index] * Gy[index]);
                if (Gout[index] > valmax) valmax = Gout[index];
            }
            //printf("Valeur max du gradient = %f\n",valmax);
            free(Gx);
            free(Gy);
            break;
        }
        case 2 : {
            Gx = (double *) calloc(w * h, sizeof(double));
            Gy = (double *) calloc(w * h, sizeof(double));
            Gxx = (double *) calloc(w * h, sizeof(double));
            Gxy = (double *) calloc(w * h, sizeof(double));
            Gyy = (double *) calloc(w * h, sizeof(double));
            TEMP2 = (double *) calloc(w * h, sizeof(double));
            // Calcul de G_gg (i.e. courbure du gradient)
            // On calcule G_x
            Gauss_Rec(Gin, TEMP, w, h, 0, 1, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gx, w, h, 1, 0, sigma, valmin, valmax);
            // puis G_y
            Gauss_Rec(Gin, TEMP, w, h, 0, 0, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gy, w, h, 1, 1, sigma, valmin, valmax);
            // puis G_xx
            Gauss_Rec(Gin, TEMP, w, h, 0, 1, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gxx, w, h, 2, 0, sigma, valmin, valmax);
            // puis G_yy
            Gauss_Rec(Gin, TEMP, w, h, 0, 0, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gyy, w, h, 2, 1, sigma, valmin, valmax);
            // et enfin G_xy
            Gauss_Rec(Gin, TEMP, w, h, 0, 1, sigma, valmin, valmax);
            Gauss_Rec(TEMP, TEMP2, w, h, 1, 0, sigma, valmin, valmax);
            Gauss_Rec(TEMP2, TEMP, w, h, 0, 0, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gxy, w, h, 1, 1, sigma, valmin, valmax);
            // Calcul de G_gg
            valmin = 1000;
            valmax = -1000;
            for (index = 0; index < w * h; index++) {
                Gout[index] = Gx[index] * Gx[index] * Gxx[index]
                              + 2 * Gx[index] * Gy[index] * Gxy[index]
                              + Gy[index] * Gy[index] * Gyy[index];
                if (produit_module == 0) {
                    square_module = Gx[index] * Gx[index] + Gy[index] * Gy[index];
                    if (square_module > 0)
                        Gout[index] /= square_module;
                }
                if (Gout[index] > valmax) valmax = Gout[index]; else if (Gout[index] < valmin) valmin = Gout[index];
            }
            free(Gx);
            free(Gy);
            free(Gxx);
            free(Gxy);
            free(Gyy);
            break;
        }
        case 3 : {
            Gx = (double *) calloc(w * h, sizeof(double));
            Gy = (double *) calloc(w * h, sizeof(double));
            Gxx = (double *) calloc(w * h, sizeof(double));
            Gxy = (double *) calloc(w * h, sizeof(double));
            Gyy = (double *) calloc(w * h, sizeof(double));
            TEMP2 = (double *) calloc(w * h, sizeof(double));
            // Calcul de G_tt (i.e. courbure de l'isophote)
            // On calcule G_x
            Gauss_Rec(Gin, TEMP, w, h, 0, 1, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gx, w, h, 1, 0, sigma, valmin, valmax);
            // puis G_y
            Gauss_Rec(Gin, TEMP, w, h, 0, 0, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gy, w, h, 1, 1, sigma, valmin, valmax);
            // puis G_xx
            Gauss_Rec(Gin, TEMP, w, h, 0, 1, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gxx, w, h, 2, 0, sigma, valmin, valmax);
            // puis G_yy
            Gauss_Rec(Gin, TEMP, w, h, 0, 0, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gyy, w, h, 2, 1, sigma, valmin, valmax);
            // et enfin G_xy
            Gauss_Rec(Gin, TEMP, w, h, 0, 1, sigma, valmin, valmax);
            Gauss_Rec(TEMP, TEMP2, w, h, 1, 0, sigma, valmin, valmax);
            Gauss_Rec(TEMP2, TEMP, w, h, 0, 0, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gxy, w, h, 1, 1, sigma, valmin, valmax);
            // Calcul de G_gg
            valmin = 1000;
            valmax = -1000;
            for (index = 0; index < w * h; index++) {
                Gout[index] = Gy[index] * Gy[index] * Gxx[index]
                              - 2 * Gx[index] * Gy[index] * Gxy[index]
                              + Gx[index] * Gx[index] * Gyy[index];
                if (produit_module == 0) {
                    square_module = Gx[index] * Gx[index] + Gy[index] * Gy[index];
                    if (square_module > 0)
                        Gout[index] /= square_module;
                }
                if (Gout[index] > valmax) valmax = Gout[index]; else if (Gout[index] < valmin) valmin = Gout[index];
            }
            free(Gx);
            free(Gy);
            free(Gxx);
            free(Gxy);
            free(Gyy);
            break;
        }
        case 4 : {
            Gx = (double *) calloc(w * h, sizeof(double));
            Gy = (double *) calloc(w * h, sizeof(double));
            Gxx = (double *) calloc(w * h, sizeof(double));
            Gxy = (double *) calloc(w * h, sizeof(double));
            Gyy = (double *) calloc(w * h, sizeof(double));
            TEMP2 = (double *) calloc(w * h, sizeof(double));
            // Calcul de G_gt
            // On calcule G_x
            Gauss_Rec(Gin, TEMP, w, h, 0, 1, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gx, w, h, 1, 0, sigma, valmin, valmax);
            // puis G_y
            Gauss_Rec(Gin, TEMP, w, h, 0, 0, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gy, w, h, 1, 1, sigma, valmin, valmax);
            // puis G_xx
            Gauss_Rec(Gin, TEMP, w, h, 0, 1, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gxx, w, h, 2, 0, sigma, valmin, valmax);
            // puis G_yy
            Gauss_Rec(Gin, TEMP, w, h, 0, 0, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gyy, w, h, 2, 1, sigma, valmin, valmax);
            // et enfin G_xy
            Gauss_Rec(Gin, TEMP, w, h, 0, 1, sigma, valmin, valmax);
            Gauss_Rec(TEMP, TEMP2, w, h, 1, 0, sigma, valmin, valmax);
            Gauss_Rec(TEMP2, TEMP, w, h, 0, 0, sigma, valmin, valmax);
            Gauss_Rec(TEMP, Gxy, w, h, 1, 1, sigma, valmin, valmax);
            // Calcul de G_gg
            valmin = 1000;
            valmax = -1000;
            for (index = 0; index < w * h; index++) {
                Gout[index] = Gx[index] * Gy[index] * (Gyy[index] - Gxx[index])
                              + Gxy[index] * (Gx[index] * Gx[index] - Gy[index] * Gy[index]);
                if (produit_module == 0) {
                    square_module = Gx[index] * Gx[index] + Gy[index] * Gy[index];
                    if (square_module > 0)
                        Gout[index] /= square_module;
                }
                if (Gout[index] > valmax) valmax = Gout[index]; else if (Gout[index] < valmin) valmin = Gout[index];
            }
            free(Gx);
            free(Gy);
            free(Gxx);
            free(Gxy);
            free(Gyy);
            break;
        }
    }
    // Normalisation
    if (normalize == 1) {
        printf("sigma = %f - derivee RI n° %d : Min = %f, Max = %f\n", sigma, ordre, valmin, valmax);
        // Basee histogramme
        for (index = 0; index < w * h; index++)
            Gout[index] = (int) ((Gout[index] - valmin) * 255.0 / (valmax - valmin));
    } else if (normalize == 2) {
        // Basee scale-space
        float scale_factor;
        if (ordre > 1) scale_factor = sigma * sigma;
        else if (ordre > 0) scale_factor = sigma;
        else scale_factor = 1.0;
        for (index = 0; index < w * h; index++)
            Gout[index] *= scale_factor;
        valmin *= scale_factor;
        valmax *= scale_factor;
    }
    free(TEMP);
}

// Routine pour tracer un vecteur en couleur
// à partir d'un pixel donné (index)
// (Algorithme de Bresenham)
void display_flow_vector(Image<int> &p, int index, int dx, int dy, int couleur) {

    int *PIX = p.PI();
    int WIDTH = p.PW();
    int HEIGHT = p.PL();
    int delta_x, delta_y;
    float pente, pente_inverse;
    int index_new;

    if (dx) pente = (float) (abs(dy)) / (float) (abs(dx));
    if (dy) pente_inverse = (float) (abs(dx)) / (float) (abs(dy));
    // On marque une croix au départ du vecteur
    //PIX[index]=0;  PIX[index-1]=0;  PIX[index+1]=0; PIX[index-w]=0;  PIX[index+w]=0;
    // Ou bien on ne marque que le départ du vecteur
    PIX[index] = 0;

    delta_x = 0;
    delta_y = 0;
    //printf("Dx = %d, Dy = %d\n",dx,dy);
    if (abs(dx) > abs(dy)) { // cas 1, 1 pixel par colonne
        while ((abs(delta_x) < abs(dx)) || (abs(delta_y) < abs(dy))) {
            delta_x += dx / abs(dx);
            index_new = index + delta_x + delta_y * WIDTH;
            if ((index_new > 0) && (index_new < WIDTH * HEIGHT))
                PIX[index_new] = 500 + couleur;
            if ((float) (abs(delta_y)) / (float) (abs(delta_x)) < pente)
                delta_y += dy / abs(dy);
            //printf("Delta_x = %d, Delta_y = %d\n",delta_x,delta_y);
        }
    } else { // cas 2, 1 pixel par ligne
        while ((abs(delta_x) < abs(dx)) || (abs(delta_y) < abs(dy))) {
            delta_y += dy / abs(dy);
            index_new = index + delta_x + delta_y * WIDTH;
            if ((index_new > 0) && (index_new < WIDTH * HEIGHT))
                PIX[index_new] = 500 + couleur;
            if ((float) (abs(delta_x)) / (float) (abs(delta_y)) < pente_inverse)
                delta_x += dx / abs(dx);
            //printf("Delta_x = %d, Delta_y = %d\n",delta_x,delta_y);
        }
    }
}

void Diff_Flow_Measures(Image<int> &p, float sigma, int type_measure, int pas_affichage) {
    int h = p.PL();
    int w = p.PW();
    int i, j, val;
    int index = 0;
    int *PIX = p.PI();
    int dx, dy;
    double *VAL, *TEMP, *TEMP1, *GX, *GY, *GXX, *GYY, *GXY;
    double valmin = 10000;
    double valmax = -10000;
    double coord_max, norm_grad, norm_vect;
    double sum_diag, diff_diag, delta;

    VAL = (double *) calloc(w * h, sizeof(double));
    TEMP = (double *) calloc(w * h, sizeof(double));
    // Passage en double...
    for (index = 0; index < w * h; index++)
        VAL[index] = (double) (PIX[index]);
    valmin = 10000;
    valmax = -10000;
    if (type_measure == 0) {// Gradient complet (avec sens et module)
        GX = (double *) calloc(w * h, sizeof(double));
        GY = (double *) calloc(w * h, sizeof(double));
        // Calcul de G_x
        Gauss_Rec(VAL, TEMP, w, h, 0, 1, sigma, valmin, valmax);
        Gauss_Rec(TEMP, GX, w, h, 1, 0, sigma, valmin, valmax);
        coord_max = max(valmax, valabs(valmin));
        // puis G_y
        valmin = 10000;
        valmax = -10000;
        Gauss_Rec(VAL, TEMP, w, h, 0, 0, sigma, valmin, valmax);
        Gauss_Rec(TEMP, GY, w, h, 1, 1, sigma, valmin, valmax);
        // On fixe la valeur max des coordonnées à 20 pixels
        coord_max = max(coord_max, max(valmax, valabs(valmin)));
        for (j = 0; j < h; j++)
            for (i = 0; i < w; i++)
                if ((i % pas_affichage == pas_affichage - 1) && (j % pas_affichage == pas_affichage - 1)) {
                    index = j * w + i;
                    dx = (int) ((GX[index] * 20.0) / coord_max);
                    dy = (int) ((GY[index] * 20.0) / coord_max);
                    display_flow_vector(p, index, dx, dy, 2);
                }
        free(GX);
        free(GY);
    } else if (type_measure == 1) {// Direction du gradient
        GX = (double *) calloc(w * h, sizeof(double));
        GY = (double *) calloc(w * h, sizeof(double));
        // Calcul de G_x
        Gauss_Rec(VAL, TEMP, w, h, 0, 1, sigma, valmin, valmax);
        Gauss_Rec(TEMP, GX, w, h, 1, 0, sigma, valmin, valmax);
        // puis G_y
        Gauss_Rec(VAL, TEMP, w, h, 0, 0, sigma, valmin, valmax);
        Gauss_Rec(TEMP, GY, w, h, 1, 1, sigma, valmin, valmax);
        for (j = 0; j < h; j++)
            for (i = 0; i < w; i++)
                if ((i % pas_affichage == pas_affichage - 1) && (j % pas_affichage == pas_affichage - 1)) {
                    index = j * w + i;
                    // Calcul du gradient normalisé
                    norm_grad = sqrt(GX[index] * GX[index] + GY[index] * GY[index]);
                    // seuillage du gradient pour affichage
                    if (norm_grad > 0.1) {
                        GX[index] /= norm_grad;
                        GY[index] /= norm_grad;
                        dx = (int) (GX[index] * 6.0);
                        dy = (int) (GY[index] * 6.0);
                    } else dx = dy = 0;
                    display_flow_vector(p, index, dx, dy, 2);
                    display_flow_vector(p, index, -dx, -dy, 2);
                }
        free(GX);
        free(GY);
    } else if (type_measure == 2) {// Direction de l'isophote
        GX = (double *) calloc(w * h, sizeof(double));
        GY = (double *) calloc(w * h, sizeof(double));
        // Calcul de G_x
        Gauss_Rec(VAL, TEMP, w, h, 0, 1, sigma, valmin, valmax);
        Gauss_Rec(TEMP, GX, w, h, 1, 0, sigma, valmin, valmax);
        // puis G_y
        Gauss_Rec(VAL, TEMP, w, h, 0, 0, sigma, valmin, valmax);
        Gauss_Rec(TEMP, GY, w, h, 1, 1, sigma, valmin, valmax);
        for (j = 0; j < h; j++)
            for (i = 0; i < w; i++)
                if ((i % pas_affichage == pas_affichage - 1) && (j % pas_affichage == pas_affichage - 1)) {
                    index = j * w + i;
                    // Calcul du gradient normalisé
                    norm_grad = sqrt(GX[index] * GX[index] + GY[index] * GY[index]);
                    // seuillage du gradient pour affichage
                    if (norm_grad > 0.1) {
                        GX[index] /= norm_grad;
                        GY[index] /= norm_grad;
                        dx = (int) (-GY[index] * 6.0);
                        dy = (int) (GX[index] * 6.0);
                    } else dx = dy = 0;
                    display_flow_vector(p, index, dx, dy, 2);
                    display_flow_vector(p, index, -dx, -dy, 2);
                }
        free(GX);
        free(GY);
    } else if (type_measure == 3) {// Premier vecteur propre du Hessien
        // module normalisé (direction)
        GXX = (double *) calloc(w * h, sizeof(double));
        GYY = (double *) calloc(w * h, sizeof(double));
        GXY = (double *) calloc(w * h, sizeof(double));
        TEMP1 = (double *) calloc(w * h, sizeof(double));
        // Calcul des composantes de la matrice Hessienne
        Gauss_Rec(VAL, TEMP, w, h, 0, 1, sigma, valmin, valmax);
        Gauss_Rec(TEMP, GXX, w, h, 2, 0, sigma, valmin, valmax);
        Gauss_Rec(VAL, TEMP, w, h, 0, 0, sigma, valmin, valmax);
        Gauss_Rec(TEMP, GYY, w, h, 2, 1, sigma, valmin, valmax);
        Gauss_Rec(VAL, TEMP, w, h, 1, 0, sigma, valmin, valmax);
        Gauss_Rec(TEMP, GXY, w, h, 1, 1, sigma, valmin, valmax);
        for (j = 0; j < h; j++)
            for (i = 0; i < w; i++)
                if ((i % pas_affichage == pas_affichage - 1) && (j % pas_affichage == pas_affichage - 1)) {
                    index = j * w + i;
                    diff_diag = GXX[index] - GYY[index];
                    delta = sqrt(4 * GXY[index] * GXY[index] + diff_diag * diff_diag);
                    norm_vect = delta - diff_diag;
                    norm_vect *= norm_vect;
                    norm_vect += 4 * GXY[index] * GXY[index];
                    norm_vect = sqrt(norm_vect);
                    if (norm_vect > 0.1) {
                        TEMP[index] = (2 * GXY[index]) / norm_vect;
                        TEMP1[index] = (delta - diff_diag) / norm_vect;
                        dx = (int) (TEMP[index] * 6.0);
                        dy = (int) (TEMP1[index] * 6.0);
                    } else dx = dy = 0;
                    display_flow_vector(p, index, dx, dy, 2);
                    display_flow_vector(p, index, -dx, -dy, 2);
                }
        free(GXX);
        free(GYY);
        free(GXY);
        free(TEMP1);
    } else if (type_measure == 4) {// Second vecteur propre du Hessien
        // module normalisé (direction)
        GXX = (double *) calloc(w * h, sizeof(double));
        GYY = (double *) calloc(w * h, sizeof(double));
        GXY = (double *) calloc(w * h, sizeof(double));
        TEMP1 = (double *) calloc(w * h, sizeof(double));
        // Calcul des composantes de la matrice Hessienne
        Gauss_Rec(VAL, TEMP, w, h, 0, 1, sigma, valmin, valmax);
        Gauss_Rec(TEMP, GXX, w, h, 2, 0, sigma, valmin, valmax);
        Gauss_Rec(VAL, TEMP, w, h, 0, 0, sigma, valmin, valmax);
        Gauss_Rec(TEMP, GYY, w, h, 2, 1, sigma, valmin, valmax);
        Gauss_Rec(VAL, TEMP, w, h, 1, 0, sigma, valmin, valmax);
        Gauss_Rec(TEMP, GXY, w, h, 1, 1, sigma, valmin, valmax);
        for (j = 0; j < h; j++)
            for (i = 0; i < w; i++)
                if ((i % pas_affichage == pas_affichage - 1) && (j % pas_affichage == pas_affichage - 1)) {
                    index = j * w + i;
                    diff_diag = GXX[index] - GYY[index];
                    delta = sqrt(4 * GXY[index] * GXY[index] + diff_diag * diff_diag);
                    norm_vect = -delta - diff_diag;
                    norm_vect *= norm_vect;
                    norm_vect += 4 * GXY[index] * GXY[index];
                    norm_vect = sqrt(norm_vect);
                    if (norm_vect > 0.1) {
                        TEMP[index] = (2 * GXY[index]) / norm_vect;
                        TEMP1[index] = (-delta - diff_diag) / norm_vect;
                        dx = (int) (TEMP[index] * 6.0);
                        dy = (int) (TEMP1[index] * 6.0);
                    } else dx = dy = 0;
                    display_flow_vector(p, index, dx, dy, 2);
                    display_flow_vector(p, index, -dx, -dy, 2);
                }
        free(GXX);
        free(GYY);
        free(GXY);
        free(TEMP1);
    } else if (type_measure == 5) {// Direction de courbure principale
        // i.e. direction du vecteur propre correspondant à la
        // valeur propre de plus grande valeur absolue
        GXX = (double *) calloc(w * h, sizeof(double));
        GYY = (double *) calloc(w * h, sizeof(double));
        GXY = (double *) calloc(w * h, sizeof(double));
        TEMP1 = (double *) calloc(w * h, sizeof(double));
        // Calcul des composantes de la matrice Hessienne
        Gauss_Rec(VAL, TEMP, w, h, 0, 1, sigma, valmin, valmax);
        Gauss_Rec(TEMP, GXX, w, h, 2, 0, sigma, valmin, valmax);
        Gauss_Rec(VAL, TEMP, w, h, 0, 0, sigma, valmin, valmax);
        Gauss_Rec(TEMP, GYY, w, h, 2, 1, sigma, valmin, valmax);
        Gauss_Rec(VAL, TEMP, w, h, 1, 0, sigma, valmin, valmax);
        Gauss_Rec(TEMP, GXY, w, h, 1, 1, sigma, valmin, valmax);
        for (j = 0; j < h; j++)
            for (i = 0; i < w; i++)
                if ((i % pas_affichage == pas_affichage - 1) && (j % pas_affichage == pas_affichage - 1)) {
                    index = j * w + i;
                    diff_diag = GXX[index] - GYY[index];
                    sum_diag = GXX[index] + GYY[index];
                    delta = sqrt(4 * GXY[index] * GXY[index] + diff_diag * diff_diag);
                    if (valabs(sum_diag + delta) > valabs(sum_diag - delta)) {
                        norm_vect = delta - diff_diag;
                        norm_vect *= norm_vect;
                        norm_vect += 4 * GXY[index] * GXY[index];
                        norm_vect = sqrt(norm_vect);
                        if (norm_vect > 0.1) {
                            TEMP1[index] = (delta - diff_diag) / norm_vect;
                        } else TEMP1[index] = 0;
                    } else {
                        norm_vect = -delta - diff_diag;
                        norm_vect *= norm_vect;
                        norm_vect += 4 * GXY[index] * GXY[index];
                        norm_vect = sqrt(norm_vect);
                        if (norm_vect > 0.1) {
                            TEMP1[index] = (-delta - diff_diag) / norm_vect;
                        } else TEMP1[index] = 0;
                    }
                    if (norm_vect > 0.1) {
                        TEMP[index] = (2 * GXY[index]) / norm_vect;
                    } else TEMP[index] = 0;
                    dx = (int) (TEMP[index] * 6.0);
                    dy = (int) (TEMP1[index] * 6.0);
                    display_flow_vector(p, index, dx, dy, 2);
                    display_flow_vector(p, index, -dx, -dy, 2);
                }
        free(GXX);
        free(GYY);
        free(GXY);
        free(TEMP1);
    }
    free(VAL);
    free(TEMP);
}


void Filter_Frangi_Complet(Image<int> &p, float sigma_min, float sigma_max, float pas_sigma, float alpha, float beta,
                           int moymax) {
    int h = p.PL();
    int w = p.PW();
    int i, j, val;
    int index = 0;
    int *PIX = p.PI();
    Image<int> q(p);
    float sigma;
    int nombre_echelle;

    double *ORIGINAL, *FRANGI_1E, *FRANGI_COMPLET;

    ORIGINAL = (double *) calloc(w * h, sizeof(double));
    FRANGI_1E = (double *) calloc(w * h, sizeof(double));
    FRANGI_COMPLET = (double *) calloc(w * h, sizeof(double));

    // Passage en double...
    for (index = 0; index < w * h; index++)
        ORIGINAL[index] = (double) (PIX[index]);
    nombre_echelle = 1;
    Filter_Frangi_1Echelle(ORIGINAL, FRANGI_COMPLET, w, h, sigma_min, alpha, beta);
    printf("Frangi, echelle = %f\n", sigma_min);

    if (moymax == 1) {

        for (sigma = sigma_min + pas_sigma; sigma <= sigma_max; sigma += pas_sigma) {
            printf("Frangi, echelle = %f, accumulation par maximum\n", sigma);
            Filter_Frangi_1Echelle(ORIGINAL, FRANGI_1E, w, h, sigma, alpha, beta);
            for (index = 0; index < w * h; index++)
                if (FRANGI_1E[index] > FRANGI_COMPLET[index]) FRANGI_COMPLET[index] = FRANGI_1E[index];
        }
    } else {
        for (sigma = sigma_min + pas_sigma; sigma <= sigma_max; sigma += pas_sigma) {
            printf("Frangi, echelle = %f, accumulation par moyenne\n", sigma);
            Filter_Frangi_1Echelle(ORIGINAL, FRANGI_1E, w, h, sigma, alpha, beta);
            nombre_echelle += 1;
            for (index = 0; index < w * h; index++)
                FRANGI_COMPLET[index] += FRANGI_1E[index];
        }
        for (index = 0; index < w * h; index++)
            FRANGI_COMPLET[index] /= nombre_echelle;
    }
    // Repassage en entier...
    for (index = 0; index < w * h; index++)
        PIX[index] = (int) (255.0 * FRANGI_COMPLET[index]);

    free(ORIGINAL);
    free(FRANGI_1E);
    free(FRANGI_COMPLET);
}

void Filter_Frangi_1Echelle(double *IMG_IN, double *IMG_OUT, int w, int h, float sigma, float alpha, float beta) {
    int i, j, val;
    int index = 0;
    double *TEMP, *GXX, *GYY, *GXY;
    double valmin = 10000;
    double valmax = -10000;
    double lambda1, lambda2;
    double sum_diag, diff_diag, delta;

    TEMP = (double *) calloc(w * h, sizeof(double));
    GXX = (double *) calloc(w * h, sizeof(double));
    GYY = (double *) calloc(w * h, sizeof(double));
    GXY = (double *) calloc(w * h, sizeof(double));
    // Calcul des composantes de la matrice Hessienne
    Gauss_Rec(IMG_IN, TEMP, w, h, 0, 1, sigma, valmin, valmax);
    Gauss_Rec(TEMP, GXX, w, h, 2, 0, sigma, valmin, valmax);
    Gauss_Rec(IMG_IN, TEMP, w, h, 0, 0, sigma, valmin, valmax);
    Gauss_Rec(TEMP, GYY, w, h, 2, 1, sigma, valmin, valmax);
    Gauss_Rec(IMG_IN, TEMP, w, h, 1, 0, sigma, valmin, valmax);
    Gauss_Rec(TEMP, GXY, w, h, 1, 1, sigma, valmin, valmax);

    // Calcul des valeurs propres et de la réponse de Frangi (1 échelle)
    for (index = 0; index < w * h; index++) {
        sum_diag = GXX[index] + GYY[index];
        diff_diag = GXX[index] - GYY[index];
        delta = sqrt(4 * GXY[index] * GXY[index] + diff_diag * diff_diag);
        lambda2 = 0.5 * (sum_diag + delta);
        lambda1 = 0.5 * (sum_diag - delta);
        if (lambda2 < 0) IMG_OUT[index] = 0;
        else {
            lambda2 *= lambda2;
            lambda1 *= lambda1;
            beta *= beta;
            alpha *= alpha;
            IMG_OUT[index] = exp(-lambda1 / (2.0 * beta * lambda2)) * (1 - exp(-(lambda1 + lambda2) / (2.0 * alpha)));
            //if ((index%1000)==0) printf("val = %f ",IMG_OUT[index]);
        }
    }
    free(GXY);
    free(GXX);
    free(GYY);
    free(TEMP);
}
