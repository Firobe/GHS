// -*- coding: utf-8 -*-
/*

   Graph Generator                                   © Cyril Gavoille


     Use: gengraph [-options] graph [parameters]


     A free command-line program to generate graphs in many formats:
     plain text, .dot, .pdf, .fig, .svg … Formats like .pdf or .fig
     are produced thanks to GraphViz and allow visualization. Easy to
     install, there is a single .c file to compile. There is a on-line
     manual available in the command (for the moment the manual is in
     French only and included at the end of the source).

                                -----

     Petit programme libre en ligne de commande permettant de générer
     des graphes dans plein de formats différents: texte, .dot, .pdf,
     .fig, .svg … Certains formats (comme .pdf ou .jpg), produits
     grâce à GraphViz, permettent de visualiser les graphes.  Très
     simple à installer, il n'y a qu'un seul fichier .c à compiler.
     Pour le reste il faut voir l'aide en ligne qui est aussi incluse
     à la fin du source.


   Comment l'installer / le compiler ?

     (MacOs)  gcc gengraph.c -o gengraph
     (Linux)  gcc gengraph.c -lm -o gengraph

*/

#define DEBUG(I) //if(1){ I }
#define _GNU_SOURCE   /* pour asprintf() */

#include <stdio.h>    /* pour printf(), sprintf() ... */
#include <stdlib.h>   /* pour system(), strtod(), RAND_MAX ... */
#include <string.h>   /* pour strcomp() ... */
#include <unistd.h>   /* pour getpid() ... */
#include <math.h>     /* pour sqrt(), cos(), acos() ... */
#include <float.h>    /* pour DBL_MAX, DBL_DIG ... */
#include <limits.h>   /* pour INT_MAX, LONG_MAX ... */
#include <time.h>     /* pour time() */
#include <sys/time.h> /* pour gettimeofday() */
#include <assert.h>   /* pour assert() */

/* graphe ou famille de graphes */

typedef struct _graph {
  int id;    // numéro du graphe, utilisé pour les familles de graphes 
  int n;     // nb de sommets, <0 si non défini
  int m;     // nb d'arêtes, <0 si non défini
  int sort;  // vrai ssi les listes d'adjacences sont triées
  int *d;    // d[u]=degré du sommet u
  int **L;   // L[u][i]=i-ème voisin du sommet u, i=0...d[u]-1
  int sym;   // vrai ssi les listes d'adjacence du graphe sont symétriques
  double **W;// W[u][i]=poids du i-ème voisin du sommet u, i=0...d[u]-1
  double *xpos,*ypos; // tableau de positions des sommets (graphe géométrique)
  int f;     // nombre de graphes de la famille, =0 si graphe normal
  struct _graph **G; // G[i]=i-ème graphe de la famille, i=0..f-1

  // les champs suivants sont utilisés pour communiquer des valeurs
  // ou des résultats à travers les appels de fonctions

  int int1;  // paramètre: entier
  int* pint1;// paramètre: tableau d'entiers
} graph;

/* fonction de test sur un graphe */

typedef int test(graph*);

/* chemin simple d'un graphe G */

typedef struct{
  int n;  /* nombre de sommets du chemin */
  int *P; /* liste ordonnée des sommets du chemin */
  int *V; /* V[u]=i si le sommet u de G est le i-ème (i dans [0,G->n[)
	     sommet du chemin (V[P[i]]=i), V[u]=-1 si u n'est dans le
	     chemin  */
  int **aux; /* tableau auxiliaire utilisé (et géré) par NextPath() */
  int nG;    /* nombre de sommets du graphe G, sert pour free_path() */
} path;

/* code pour le type "list" */

enum{
  T_NODE,  // item est un sommet u
  T_EDGE,  // item est un sommet v connecté au précédent de la liste par une arête: u-v
  T_ARC,   // item est un sommet v connecté au précédent de la liste par un arc: u->v
  T_ID,    // item est l'identifiant d'un nouveau graphe d'une famille
  T_NB,    // item est le nombre de graphes de la famille, item toujours en tête de liste
  T_OPENE, // item est le premier sommet v d'un groupe arête: u-(v ...
  T_OPENA, // item est le premier sommet v d'un groupe arc: u->(v ...
  /* token non encore géré */
  T_UNIV,  // u-* = item est un sommet universel u-(0 ... n)
  T_UNOUT, // u->* item est un sommet universel u->(0 ... n-1)
  T_UNIN,  // *->u item est un sommet universel u<-(0 ... n-1)
};

/* code pour les fonctions de hashage */

enum{
  H_PRIME,
  H_SHUFFLE,
  H_MIX,
  H_MOD,
};


/* code pour les normes et graphes géométriques */

enum{
  NORM_FAIL,  // norme indéterminée
  NORM_L1,    // norme L1
  NORM_L2,    // norme L2
  NORM_LMAX,  // norme Lmax
  NORM_LMIN,  // norme Lmin
  NORM_POLY,  // norme polygonale
  NORM_HYPER, // norme hyperbolique
};


typedef struct _list { int item,type; struct _list *next; } list; /* liste chaînée */
typedef struct{ double x,y; } point; /* un point du plan */
typedef struct{ int u,v; } edge; /* une arête */
typedef struct{ int x,y; long z; } triplet; /* un triplet d'entiers (avec lon long) */
typedef struct{ unsigned int r:8,g:8,b:8; } color;
/* une couleur est un int accessible par les champs r,g,b de 8 bits */


/* constantes pour la ligne de commande */

#define PARMAX      16
/* nb maximum de paramètres pour un graphe ayant un nombre de
   paramètres fixés (non variable). Le plus grand nombre de paramètres
   (non variable) est pour le graphe de Chvatal. Il est de 12. */

#define DIMAX       64 /* dimension maximum pour un graphe (hypercube, pancake, ...) */
#define CMDMAX    1024 /* nb maximum de caractères sur la ligne de commande */
#define NAMEMAX    256 /* nb maximum de caractères pour un nom de sommet */
#define PARAMSIZE 1024 /* taille mémoire des buffers FPARAM et CPARAM (en octets) */
#define SURFACEMAX 256 /* nombre maximum de coutures d'une surface */

/* constantes pour le format dot et -visu */

double VSIZESTD=0.05; /* taille standard des sommets */
double VSIZEXY=0.12;  /* taille des sommets pour les graphes géométriques */
double LEN=1.0; /* longueur des arêtes avec dot neato */
#define URL_vis_js1 "vis.min.js"
#define URL_vis_js2 "https://cdnjs.cloudflare.com/ajax/libs/vis/4.20.1/vis.min.js"
#define GRAPHPDF g.pdf /* nom du graphe de sortie par défaut */
#define VSIZEK 5.5 /* rapport entre taille max et min des sommets (-vsize) */
#define C32 32
/* Coefficient par lequel sont multipliées les coordonnées des points
   dans le cas des graphes géométriques et le format dot. Plus la
   valeur est élevée, plus les sommets paraissent petits et les arêtes
   fines et longues */

/* constantes diverses */

#define STR_DISTRIB "▪" /* caractère pour affichage des distributions */
#define LEN_DISTRIB (60.0) /* longueur max pour l'affichage des distributions */
#define LEN_BARRE 80 /* longueur d'une barre de séparation */
#define CHRONOMAX 5 /* nombre maximum de chronomètres */
#define RS_NI_FP "name-independent, fixed-port model"
#define RS_L_FP  "labeled, fixed-port model"


/* codes pour les formats de sorties possibles */

enum{
  F_standard,
  F_list,
  F_matrix,
  F_smatrix,
  F_dot,
  F_userdot,
  F_html,
  F_xy,
  F_no,
};


/* codes pour les algorithmes via l'option -check. La variable CHECK
   vaut l'une de ces constantes. Par défaut, CHECK=CHECK_OFF(=0). Si
   CHECK>0, alors le graphe sera stocké en mémoire. Si CHECK>1, alors
   en plus l'algorithme spécifique sera appliqué. */

enum{
  CHECK_OFF=0, // valeur par défaut
  CHECK_ON,    // graphe chargé en mémoire
  CHECK_BFS,
  CHECK_DFS,
  CHECK_NCC,
  CHECK_BELLMAN,
  CHECK_STRETCH,
  CHECK_DEGENERATE,
  CHECK_GCOLOR,
  CHECK_DEG,
  CHECK_ISO,
  CHECK_SUB,
  CHECK_ISUB,
  CHECK_MINOR,
  CHECK_TWDEG,
  CHECK_TW,
  CHECK_PS1,
  CHECK_PS1b,
  CHECK_PS1c,
  CHECK_PS1x,
  CHECK_GIRTH,
  CHECK_RADIUS,
  CHECK_DIAMETER,
  CHECK_VOLM,
  CHECK_PATHS,
  CHECK_MAINCC,
  CHECK_KCOLOR,
  CHECK_KCOLORSAT,
  CHECK_KINDEPSAT,
  CHECK_INFO,
  CHECK_SIMPLIFY,
  CHECK_RS_CLUSTER,
  CHECK_RS_DCR,
  CHECK_RS_AGMNT,
  CHECK_RS_TZRPLG,
  CHECK_RS_BC,
  CHECK_RS_HDLBR,
};


/* codes pour la variable XYtype indiquant le type de
   distribution/génération des points (xpos,ypos) des graphes
   géométriques */

enum{
  XY_FILE,    // points lus à partir d'un fichier
  XY_UNIF,    // points selon une loi uniforme
  XY_PLAW,    // points selon une loi en puissance
  XY_PERM,    // points selon une permutation π, (i,π(i))
  XY_MESH,    // points sur une grille Xmesh x Ymesh
  XY_CIRCLE,  // points aléatoires sur un cercle de rayon 1
  XY_CYCLE,   // points réguliers sur un cercle de rayon 1
  XY_RPOLY,   // points aléatoire dans un polygone convexe régulier
  XY_DISK,    // points dans le disque unité triés selon l'angle
  XY_HYPER,   // points selon une loi en puissance dans le disque unité
  XY_CONVEX,  // points dans le disque en position convexe
  XY_CONVEX2, // variante de XY_CONVEX
  XY_USER,    // points définis par la fonction d'adjacence
};


/* fonctions inline */

/* Attention! utiliser fmin() pour des doubles, car
   imin(double,double) ne provoque par de "warning" à la
   compilation */

static inline int imin(const int i,const int j){ return (i<j)?i:j; }
static inline int imax(const int i,const int j){ return (i>j)?i:j; }

/* renvoie l'entier ceil(x/y) lorsque x,y sont entiers et y<>0. Marche
   même pour les entiers négatifs. */
static inline int iceil(const int x,const int y){ return (x/y) + (((x<0)^(y>0))&&(x%y)); }

/* macros & pseudo-fonction utiles */

#define EQUAL(s)   (strcmp(ARGV[i],s)==0)
#define PREFIX(s)  (strncmp(ARGV[i],s,strlen(s))==0)
#define SUFFIX(s)  (strncmp(ARGV[i]+imax(strlen(ARGV[i])-strlen(s),0),s,strlen(s))==0)
#define PLURIEL(n) (((n)>1)?"s":"")
#define SWAP(x,y,z)  (z)=(x),(x)=(y),(y)=(z) /* échange x et y, via z */
#define XORSWAP(x,y) ((x)^=(y),(y)^=(x),(x)^=(y)) /* comme SWAP() mais sans z */
#define VIDE(s)    *s='\0'    /* vide le tableau de caractères s */ 
#define ESTVIDE(s) (*s=='\0') /* vrai ssi s est vide */
#define NONVIDE(s) (*s!='\0') /* vrai ssi s est non vide */
#define MEM(mem,pos,type) (*(type*)(mem+pos)) /* écrit/lit la mémoire mem[pos] */
#define RAND01  ((double)random()/(double)RAND_MAX) /* réel aléatoire dans [0,1[ */

/* appelle qsort(T,...) avec la bonne taille d'élément (taille de type T[]) */
#define QSORT(T,n,f) qsort(T,n,sizeof(*(T)),f)

/* permet l'expansion d'une macro. Ex: scanf("%"xstr(DMAX)"s",buffer); */
#define xstr(s) str(s)
#define str(s) #s

/* Alloue à la variable T un tableau de n valeurs du type de T[] */
/* Si n=0, alors ne fait rien, mais free(T) ne fera pas d'erreur */
#define _ALLOC(T,n) if(((T)=malloc((n)*sizeof(*(T))))==NULL) Erreur(3)
#define ALLOC(T,n) do{ _ALLOC(T,n); }while(0)
/*
  #define _ALLOC(T,n) do{				\
  if(!((T)=malloc((n)*sizeof(*(T))))) Erreur(3);	\
  printf("malloc("#T ","#n ")=%p\n",T);			\
  }while(0)
  #define free(p) printf("free("#p ")\n"),free(p)
*/

/* Comme ALLOC mais initialise en plus le tableau avec le term z */
#define ALLOCZ(T,n,z)						\
  do{								\
    _ALLOC(T,n);						\
    int _i;for(_i=0;_i<(n);_i++) (T)[_i]=(z);			\
  }while(0)

/* Un realloc() qui évite de faire un free() non souhaité à un
   pointeur si n=0, et donc qui évite les erreurs "double free". */
#define REALLOC(P,n) P=realloc(P,imax(n,1)*sizeof(*P))

/* Alloue à la variable T une matrice de n x s valeurs, c'est-à-dire
   un tableau de n tableaux de s valeurs du type de T[][]. */
#define ALLOCMAT(T,n,s)			\
  do{					\
    int _k=0,_n=n;			\
    _ALLOC(T,_n);			\
    for(;_k<_n;_k++) _ALLOC(T[_k],s);	\
  }while(0)

#define ALLOCMATZ(T,n,s,z)		\
  do{					\
    int _k=0,_n=n;			\
    _ALLOC(T,_n);			\
    for(;_k<_n;_k++) ALLOCZ(T[_k],s,z);	\
  }while(0)

/* Les variantes NALLOC permettent, en plus de ALLOC, de déclarer le
   pointeur T sur des valeurs de type t.

   Ex: int *T;      ->   NALLOC(int,T,n);
       ALLOC(T,n);

   Attention de ne pas mettre ces macros dans un bloc { ... }, car
   sinon la portée de T se limitera à ce block. Les variantes NALLOC
   ne peuvent remplacer des ALLOC qui sont dans les configurations
   suivantes:

   - la déclaration et le malloc() ne sont pas dans le même bloc, comme:
     int *T;
     if(...)
       ALLOC(T,...);

   - la déclaration est après l'étiquette d'un goto ou d'un case, comme:
     label:
       int *T;
       ALLOC(T,...);
*/
#define NALLOC(t,T,n) t *T;ALLOC(T,n);
#define NALLOCZ(t,T,n,z) t *T;ALLOCZ(T,n,z);
#define NALLOCMAT(t,T,n,s) t **T;ALLOCMAT(T,n,s);
#define NALLOCMATZ(t,T,n,s,z) t **T;ALLOCMATZ(T,n,s,z);

/* Libère un tableau de pointeurs. On ne peut pas faire de fonction,
   le type de X n'étant pas pré-déterminé on ne peut écrire le
   prototype de la fonction. */
#define FREE2(X,n)			\
  do{					\
    if(X){				\
      int _j=0;				\
      for(;_j<(n);_j++) free(X[_j]);	\
      free(X);				\
    }					\
  }while(0)

/* Comme FREE2() mais pour un tableau de dimension 3.
   Ex: FREE3(X,n,m[_i]); */
#define FREE3(X,n,m)			\
  do{					\
    if(X){				\
      int _i=0;				\
      for(;_i<(n);_i++) FREE2(X[_i],m);	\
      free(X);				\
    }					\
  }while(0)

/* ajoute une arête uv au graphe G en fin de liste de u et de v, en
   supposant qu'il y a assez de place. Attention ! seul le champs G->d
   est mis à jour. Les champs comme G->m, G->sort peuvent ne pas être
   cohérents */
#define ADD_EDGE(G,u,v)	   	\
  do{				\
    G->L[u][G->d[u]++]=v;	\
    G->L[v][G->d[v]++]=u;	\
  }while(0)

/* comme ADD_EDGE, mais sans ajouter u à la liste de v */
#define ADD_ARC(G,u,v)	   	\
  do{				\
    G->L[u][G->d[u]++]=v;	\
  }while(0)

/* macros de débuggage */

/* affiche une ligne vide */
#define PRINTN printf("\n")

/* affiche une variable v de type int */
#define PRINT(v)					\
  do{							\
    printf(#v " = %i\n",v);				\
  }while(0)

/* affiche une variable v de type double */
#define PRINTD(v)					\
  do{							\
    printf(#v " = %g\n",v);				\
  }while(0)

/* affiche une variable v de type pointeur */
#define PRINTP(v)					\
  do{							\
    printf(#v " = %p\n",v);				\
  }while(0)

/* affiche un tableau d'entiers */
#define PRINTT(T,n)					\
  do{							\
    int _i;						\
    printf(#T " =");					\
    if((T)==NULL) printf(" NULL"); else			\
      for(_i=0;_i<(n);_i++) printf(" %i",(T)[_i]);	\
    printf("\n");					\
  }while(0)

/* affiche une liste chaînée */
#define PRINTLIST(L)		        	\
  do{						\
    list *_L=L;                                 \
    printf(#L " = {");				\
    while(_L!=NULL){				\
      printf(" (%i,%i)",_L->item,_L->type);     \
      _L=_L->next;                              \
      } printf(" }\n");				\
  }while(0)

// affiche un aperçu d'une liste de n valeurs
// ex: APERCU(_i*_i,10,3,2) affichera
// "[ 0 1 2 ... 64 81 ]"
#define APERCU(F,n,K1,K2)						\
  do{									\
    int _i;								\
    printf("[ ");							\
    for(_i=0;_i<n;_i++)							\
      if((_i<K1)||(_i>n-K2-1)||(K1==n-K2-1))				\
	printf("%i ",F);						\
      else{ if((_i>=K1)&&(_i<n-K2)){ printf("... "); _i=n-K2-1; }}	\
    printf("]\n");							\
  }while(0)

#define PAUSE scanf("%*c") /* appuyer sur [RETURN] pour continuer */
#define STRTOI(s) ((int)strtol(s,NULL,10)) 
#define STRTOL(s) ((long)strtol(s,NULL,10)) 
#define STRTOD(s) strtod(s,NULL)
#define FAIL_ROUTING -1000000 // constante int suffisament négative pour provoquer un échec de routage

/* variables globales */

int NF;          /* nb de sommets final du graphes (donc après suppression) */
int *V;          /* étiquette des sommets, en principe V[i]=i */
int *VF;         /* VF[j] est l'indice i du j-ème sommet non supprimé */
int *INC;        /* INC[i]=deg(i). Si =0, alors i est un sommet isolé */
int ARGC;        /* variable globale pour argc */
char **ARGV;     /* variable globale pour argv */
char PARAM_PAL[64]; /* mot définissant le dégradé de couleur pour -vcolor pal */
void* CPARAM=NULL; /* liste de paramètres (pointeur tout type, en octets) pour -check */
void* FPARAM=NULL; /* liste de paramètres (pointeur tout type, en octets) pour -filter */
int CVALUE;      /* sert pour la valeur dans -filter */
int PVALUE;      /* =1 ssi on affiche la valeur du paramètre dans -filter */
test *FTEST;     /* pour l'option -filter */
double DELE=0;   /* proba de supprimer une arêtes */
double DELV=0;   /* proba de supprimer un sommet */
double REDIRECT=0; /* proba de rediriger une arête */
int VERTEX0=-1;  /* voisinage du sommet à afficher pour -format vertex */
int SHIFT=0;     /* début de la numérotation des sommets */
int DIRECTED=0;  /* vrai ssi graphe orienté for(i=0 ... for(j=...) ) */
int LOOP=0;      /* vrai ssi ont n'affiche pas les boucles */
int PERMUTE=0;   /* vrai ssi -permute */
int NOT=0;       /* vrai ssi -not */
int POS=0;       /* vrai ssi -pos */
int FAST=0;      /* vrai ssi -fast */
int CHECK=0;     /* vrai ssi option -check */
int VARIANT=0;   /* variante de l'option -variant */
int STAR=0;      /* paramètre de -star */
int APEX=0;      /* paramètre de -apex */
int LABEL=0;     /* vrai ssi affiche les labels des sommets (-label) */
int NORM=NORM_L2;/* norme pour les graphes géométriques: L2 par défaut */
int NORM_poly=3; /* nombre de cotés de la norme polygonale */
int FORMAT=F_standard; /* type de la sortie, par défaut le format standard */
int HEADER=0;    /* par défaut pas de préambule, sinon HEADER=1 */
int WIDTH=12;    /* nb maximum d'arêtes ou de sommets isolés affichés par ligne */
unsigned SEED;   /* graîne du générateur aléatoire */
char *DOTFILTER="neato"; /* nom du filtre "dot" par défaut */
char *DOTSCALE=NULL; /* facteur d'échelle des points et arêtes du graphe (dot) */
char *CAPTION=NULL; /* légende du graphe (dot) */

char *FILEXY;    /* nom de fichier des coordonnées */
double *XPOS=NULL,*YPOS=NULL; /* tableaux (X,Y) pour les graphes géométriques */
double BOXX=-1,BOXY; /* pour le redimensionement: <0 signifie pas de redim. */
double XMIN=0,YMIN=0,XMAX=1,YMAX=1; /* Bounding Box par défaut */
int ROUND=DBL_DIG; /* arrondi de XPOS/YPOS à 10^-ROUND près */
/* DBL_DIG (=15) est considérée comme valeur impossible pour ROUND et signifie aucun arrondi */
double XYnoiser=-1,XYnoisep; /* paramètres pour -xy noise: <0 signifie pas de "noise" */
int XYtype=XY_UNIF; /* type de génération des points, uniforme par défaut */
int XYunique=0; /* =1 ssi on élimine les points doubles */
int XYgrid=0; /* <>0 si on affiche une grille grisée */
int XYpoly; /* nombre de cotés du polygone régulier pour -xy polygon */
double XYratio=1; /* ratio pour les options -xy polar, circle, convex */
int XYzero=0; /* =1 ssi il faut ajouter le point (0,0) en rouge pour le format dot */
int XYborder=0; /* =1 ssi il faut ajouter un bord pour le format dot */
int Xmesh=0,Ymesh=0; /* dimension de la grille pour l'option -xy mesh */
double XYvsize=1; /* facteur d'échelle pour la taille des sommets dans F_dot */
int XYseedk;  /* nombre de graînes pour génération des points */
double XYpower; /* exposant pour la loi puissance de la génération des points */
double *XSEED=NULL,*YSEED=NULL; /* tableaux de doubles pour les graînes */
int XYsurfacesize=0; /* nombre de coutures de la surface pour les graphes géométriques */
int XYsurface[SURFACEMAX]; /* signature de la surface */
int LOADC=0; /* vrai ssi graphe "loadc file" */
graph *GF=NULL;     /* graphe pour l'option -check */
graph *FAMILY=NULL; /* graphe pour l'option -filter */
int HASH=H_PRIME; /* fonction de hashage par défaut */
int VSIZE=0; /* code pour la taille des sommets */
int VCOLOR=0; /* code pour la couleur les sommets */
char COLORCHAR[]="redjykugfocatbhsqvmpinzxlw";
color COLORBASE[]={ /* HTML Color Names */
  /* palette des couleurs de bases, l'ordre doit être celui de COLORCHAR */
  {255,  0,  0}, // r=red
  {210,105, 30}, // e=chocolate
  {255,140,  0}, // d=darkorange
  {255,165,  0}, // j=orange
  {255,255,  0}, // y=yellow
  {240,230,140}, // k=khaki
  {154,205, 50}, // u=yellowgreen
  {  0,255,  0}, // g=green (lime)
  { 34,139, 34}, // f=forestgreen
  {128,128,  0}, // o=olive
  {  0,255,255}, // c=cyan
  {127,255,212}, // a=aquamarine
  {  0,128,128}, // t=teal
  {  0,  0,255}, // b=blue
  {255,105,180}, // h=hotpink
  {250,128,114}, // s=salmon
  {255,192,203}, // q=pink
  {238,130,238}, // v=violet
  {255,  0,255}, // m=magenta
  {128,  0,128}, // p=purple
  { 75,  0,130}, // i=indigo
  {  0,  0,128}, // n=navy
  {  0,  0,  0}, // z=black
  {128,128,128}, // x=gray
  {230,230,250}, // l=lavender
  {255,255,255}  // w=white
};
color *PALETTE=COLORBASE; /* palette par défaut */
#define COLORNB    ((int)sizeof(COLORCHAR)-1) /* nb de couleurs de la palette de base */
int NPAL=COLORNB; /* taille de la palette par défaut */
struct{
  int mode,u,v,dist;
  double stretch;
} SCENARIO; /* scenario pour l'option "-check routing" */


/* structure utilisée pour l'évaluation d'un graphe et aussi Out() */
/* tous les calculs, notamment les entrées/sorties, effectués par les
   fonctions d'adjacence ne devraient utiliser que les éléments de
   cette structure */

struct _query {
  int code; // type de requête, voir enum{ QUERY_... }
  int i,j; // sommet i et j
  double *xpos,*ypos; // coodonnées des points pour les graphes géométriques
  char *sparam; // paramètre chaînes de caractères, libérer par free_query()
  
  // Les paramètres du graphes, param[] et dparam[], sont des tableaux
  // de taille PARMAX sauf pour les graphes ayant un nombre de
  // paramètres variables (bdrg, ...). Dans ce cas ils sont
  // réalloués. Ils seront libérés dans le free_query() final.

  int *param; // paramètres entiers du graphes, alloué/libéré par new/free_query()
  double *dparam; // paramètres réels du graphes, alloué/libéré par new/free_query()

  /* utilisée pour une valeur de retour */

  int n; // nombre de sommet du graphe
  int a; // est-ce que i est voisin de j ?
  char name[NAMEMAX+20]; // nom du sommet i, avec une marge de 20 pour
			 // pouvoir écrire au moins un "item" (entier,
			 // double, caractère, ...)
  int **rep; // pour la représentation implicite du graphe
  double **drep; // comme rep mais avec des doubles

  /* pas encore utilisés */

  int deg; // degré du sommet i
  int *list; // liste des voisins du sommet i
  graph *G; // graphe complet ou partiel
  int error; // code d'erreur, 0 si tout est ok, >0 sinon

  /* divers */

  int *wrap; // tableau annexe (cf. grid, rpartite, permutation,...)
  struct _query* query;
  int (*adj)(struct _query* const);
};


typedef struct _query query;
typedef int adjacence(query* const);


/* codes de requête/erreur pour les fonctions d'adjacence */

enum{
  QUERY_INIT,  // résultat dans ->n, initialise la fonction
  QUERY_END,   // termine la fonction
  QUERY_ADJ,   // résultat dans ->a
  QUERY_NAME,  // résultat dans ->name
  QUERY_DEG,   // résutlat dans ->deg
  QUERY_LIST,  // résultat dans ->list
  QUERY_GRAPH, // résultat dans ->G
  QUERY_ISOL,  // pour l'affichage de sommets isolés dans Out()
  QUERY_DOT,   // pour le dessin des arêtes (A FINIR)
  QUERY_LNAME, // pour la liste des noms des sommets (A FINIR)
};


struct{ /* pour définir une arête avec Out(i,j) sous le format F_userdot */
  adjacence *adj; // nom de la fonction d'adjacence
  int i,j;        // dernière arête calculée par adj(i,j)
  void *ptr;      // informations permettant d'afficher le dessin de l'arête i-j 
} USERDOT;


/***********************************

         ROUTINES EN VRAC

***********************************/


void Erreur(const int erreur){ /* affiche l'erreur et termine avec exit() */
  char *s;
  switch(erreur){
  case  1: s="option -xy non reconnue."; break;
  case  2: s="option non reconnue."; break;
  case  3: s="espace mémoire insuffisant."; break;
  case  4: s="nombre trop grand de paramètres."; break;
  case  5: s="format de sortie inconnu."; break;
  case  6: s="paramètre incorrect."; break;
  case  7: s="ouverture du fichier impossible."; break;
  case  8: s="tableau de coordonnées inexistant."; break;
  case  9: s="option -vcolor non reconnue."; break;
  case 10: s="nom de graphe inconnu ou erreur dans les paramètres."; break;
  case 11: s="le graphe doit être connexe."; break;
  case 12: s="option -check non reconnue."; break;
  case 13: s="format de famille de graphes invalide."; break;
  case 14: s="option -filter non reconnue."; break;
  case 15: s="graphe(s) non trouvé(s)."; break;
  case 16: s="plage de valeurs incorrecte."; break;
  case 17: s="nom ou identifiant de sommets trop grand."; break;
  case 18: s="opérateurs -star et -apex incompatibles."; break;
  case 19: s="nom de fichier trop long."; break;
  case 20: s="nombre de couleurs dans la palette trop important."; break;
  case 21: s="code inconnue dans la fonction SortInt()."; break;
  case 22: s="sommet de valeur négative dans le format standard."; break;
  case 23: s="code inconnu dans la fonction routing_test()."; break;
  case 24: s="option -visu incompatible avec -format no."; break;
  case 25: s="la variante -fast n'est pas implémentée pour ce graphe."; break;
  case 26: s="numéro de chronomètre incorrect."; break;
  case 27: s="plusieurs options -check sur la ligne de commande."; break;
  case 28: s="mauvais format du fichier d'entrée."; break;
  case 29: s="loadc et -not/-permute sont incompatibles."; break;
  case 30: s="loadc devrait être suivi d'une option -check."; break;
  case 31: s="le graphe doit être simple et non-orienté, essayez -check info."; break;
  case 32: s="problème dans la génération du graphe."; break;
  case 33: s="dépassement arithmétique."; break;
  case 34: s="la séquence de degré n'est pas graphique."; break;
  case 35: s="-caption ne devrait contenir qu'une occurence du format %XXX."; break;
  case 36: s="le graphe doit comporter au moins une arête."; break;
  case 37: s="sommet n'appartenant pas au graphe."; break;
  case 38: s="la probabilité devrait être dans l'intervalle [0,1]."; break;
  case 39: s="schéma de routage non reconnu."; break;
  case 40: s="fonction de hachage non reconnue."; break;
  case 41: s="scenario non reconnu."; break;
  case 42: s="nombre de coté du polygone incorrect."; break;
  case 43: s="option -norm non reconnue."; break;
  case 44: s="option définie seulement pour les graphes géométriques."; break;
  case 45: s="signature incorrecte dans l'option -xy surface."; break;
  case 46: s="l'entier doit être inférieur à INT_MAX."; break;
  case 47: s="algorithme non implémenté."; break;
  case 48: s="graphe vide ou dfs()/bfs() à partir d'un sommet non valide."; break;
  case 49: s="option -dot non reconnue."; break;
  default: s="code d'erreur inconnue."; /* ne devrait jamais arriver */
  }
  fprintf(stderr,"Erreur : %s\n",s);
  exit(EXIT_FAILURE);
}


void Permute(int *T,const int n){
/*
  Permute aléatoirement les n premiers éléments de T.
*/
  int i,j,k;
  for(i=0;i<n;i++){
    j=i+(random()%(n-i));
    SWAP(T[i],T[j],k);
  }
}


void name_vector(char* const S,
		 const int *R,const int n,
		 const char *sep,const char *par,const int d,const char *f){
/*
  Ecrit dans la variable dans S le nom représenté par les n premiers
  entiers (lettres) du tableau R. Les lettres (chiffres) sont séparées
  par la chaîne "sep" (qui peut être vide). Le mot est parenthésé par
  la chaîne "par" qui peut soit être vide soit avoir exactement deux
  caractères: par[0]=parenthèse gauche, par[1]=parenthèse droite. Si
  d>0, les lettres sont écrites dans le sens croissant des indices
  (vers la droite), sinon c'est dans le sens des indices
  décroissant. La chaîne f indique le format d'affichage (avec printf)
  de l'entier R[i]. En général, il s'agit de "%i".

  Ex: R={3,6,1}, n=3, f="%i".
  d=1, sep="," et par="{}" alors S="{3,6,1}"
  d=1, sep="-" et par=""   alors S="3-6-1"
  d=1, sep=""  et par=""   alors S="361"
  d=0, sep=""  et par=""   alors S="163"
*/
  int i,b,c;
  int p=!ESTVIDE(par); /* p=0 ou 1=nombre de caractère écrit dans S */
  
  if(d>0) c=1,i=0,b=n;
  else c=-1,i=n-1,b=-1;

  /* parcoure R dans un sens ou l'autre */
  VIDE(S); /* vide la chaîne */
  while(i!=b){
    p+=sprintf(S+p,f,R[i]);
    if(p>NAMEMAX) Erreur(17);
    i += c;
    if(i!=b) p+=sprintf(S+p,"%s",sep);
  }

  /* met les parenthèses ? */
  if(ESTVIDE(par)) return;
  S[0]=par[0];
  S[p]=par[1];
  S[p+1]='\0';
}


void name_base(char* const S,
	       int u,const int b,const int n,
	       const char *sep,const char *par,const int d){
/*
  Comme name_vector(...,n,sep,par,d,"%i") sauf que S est l'écriture de
  l'entier u en base b.

  Ex: R=361, b=10, n=3 et d=1.
  sep="," et par="{}" alors S="{3,6,1}"
  sep="-" et par=""   alors S="3-6-1"
  sep=""  et par=""   alors S="361"
*/
  if(b<2) return;
  int R[NAMEMAX],i;
  for(i=0;i<n;i++) R[i]=u%b, u /= b;
  name_vector(S,R,n,sep,par,d,"%i");
}


int LoadXY(const char *file){
/*
  Remplit les tableaux XPOS et YPOS à partir d'un fichier (file), et
  renvoie le nombre de points lues. Le programme peut être amené à
  lire un point de trop (et des coordonnées vides) si le "n" du
  fichier est > aux nombres de points qu'il contient réellement et que
  le fichier contient un retour de ligne finale.
*/
  FILE *f=stdin;
  int n,i;

  if(strcmp(file,"-")) f=fopen(file,"r"); /* ouvre en lecture */
  if(f==NULL) Erreur(7);
  i=fscanf(f,"%i",&n); /* lit n */
  if((n<0)||(i<0)) n=0;
  ALLOC(XPOS,n);
  ALLOC(YPOS,n);
  for(i=0;(i<n)&&(!feof(f));i++){
    if(fscanf(f,"//%*[^\n]\n")>0) continue;
    fscanf(f,"%lf %lf",XPOS+i,YPOS+i);
  }
  fclose(f);
  return i; /* i=minimum entre n et le nb de valeurs lues */
}


list *new_list(void){
/*
  Crée une nouvelle liste (qui est renvoyée) comprenant une célulle
  (sentinelle) avec le champs ->next=NULL.
*/
  NALLOC(list,L,1);
  L->next=NULL;
  return L;
}


static inline list *Insert(list *p,int v,int t){
/*
  Écrit (v,t) dans l'élément de la liste chaînée pointée par p qui ne
  doit pas être NULL, crée un nouvel élément qui est chaîné à la suite
  de p. On renvoit le nouvel élément crée.
*/
  p->item=v;
  p->type=t;
  return p->next=new_list(); /* nouvel élément vide */
}


graph *new_graph(const int n){
/*
  Renvoie un objet de type "graph". Les champs sont initialisés à
  leurs valeurs par défaut. Si n>0, alors les tableaux ->d et ->L de
  taille n sont alloués, mais pas les n tableaux ->L[u]. Les tableaux
  ->W, ->xpos et ->ypos ne sont pas alloués.
*/

  NALLOC(graph,G,1);
  G->n=0;
  G->m=-1;
  G->sort=0;
  G->id=-1;
  G->d=NULL;
  G->L=NULL;
  G->W=NULL;
  G->xpos=NULL;
  G->ypos=NULL;
  G->f=0;
  G->sym=1;
  G->G=NULL;

  G->pint1=NULL;
  G->int1=-1;

  if(n>0){
    G->n=n;
    ALLOC(G->d,n);
    ALLOC(G->L,n);
  }

  return G;
}


void free_graph(graph *G)
/*
  Libère G et tous ses tableaux.  Dans le cas d'une famille G, chaque
  graphe est aussi libéré (de manière récursive).
*/
{
  if(G==NULL) return;

  /* Remarque: ce n'est pas grave de faire free() sur un ptr NULL */

  int i;
  free(G->d);
  free(G->pint1);
  free(G->xpos);
  free(G->ypos);
  FREE2(G->L,G->n);
  FREE2(G->W,G->n);
  for(i=0;i<G->f;i++) free_graph(G->G[i]);
  free(G->G);
  free(G);
}


long SizeOfGraph(const graph *G){
/*
  Donne la taille mémoire (nombre d'octets) utilisé par le graphe G
  (ou une famille de graphe).
*/
  int u;
  const int n=G->n;
  long t=sizeof(*G) + 2*n*sizeof(int); // taille de G->L et G->d
  if(G->W) t += n*sizeof(*(G->W)); // ajoute G->W
  if(G->xpos) t += n*sizeof(*(G->xpos)); // ajoute G->W
  if(G->ypos) t += n*sizeof(*(G->ypos)); // ajoute G->W
  for(u=0;u<n;u++) t += G->d[u]*sizeof(int);
  if(G->f==0) return t; // cas d'un graphe simple

  for(u=0;u<G->f;u++) t += SizeOfGraph(G->G[u]);
  return t;
}


path *new_path(graph *G,int *L,int k){
/*
  Créer un chemin d'un graphe G, défini par une liste L de k sommets
  de G. Attention, le champs P du chemin renvoyé est utilisé en
  interne. Il ne faut pas détruire P après cet appel. P sera libéré
  par free_path(). Si L=NULL, alors le champs P de taille k est
  alloué, et le champs n=0. C'est une façon de créer un chemin vide
  d'au plus k sommets. Le champs V, de taille G->n, est initialisé
  suivant L (si L<>NULL), ou bien à -1 (si L=NULL).
*/
  if(G==NULL) return NULL;

  NALLOC(path,Q,1); // Q=nouveau chemin qui sera renvoyé
  ALLOCZ(Q->V,G->n,-1);
  Q->aux=NULL;
  Q->nG=G->n;

  if(L){
    int i;
    for(i=0;i<k;i++) Q->V[L[i]]=i;
    Q->P=L;
    Q->n=k;
  }
  else{
    ALLOC(Q->P,k);
    Q->n=0;
  }

  return Q;
}


void free_path(path *P){
  if(P==NULL) return;
  free(P->P);
  free(P->V);
  FREE2(P->aux,P->nG);
  free(P);
}


query *new_query(void){
  NALLOC(query,Q,1);

  Q->error=0; // par défaut tout est ok
  Q->i=Q->j=-1;
  Q->n=-1;
  Q->a=-1;
  Q->deg=-1;

  Q->rep=NULL;
  Q->list=NULL;
  Q->xpos=Q->ypos=NULL;
  Q->sparam=NULL;
  Q->wrap=NULL;
  Q->adj=NULL;
  Q->query=NULL;
  Q->G=NULL;
  VIDE(Q->name);

  ALLOC(Q->param,PARMAX);
  ALLOC(Q->dparam,PARMAX);

  return Q;
}


void free_query(query *Q){
  if(Q==NULL) return;
  free(Q->param);
  free(Q->dparam);
  free(Q->sparam);
  free_query(Q->query);
  free_graph(Q->G);
  free(Q);
  return;
}


/* structure de données pour une forêt enracinée */
typedef struct{
  int n;       /* nombre de sommets */
  int nroot;   /* nombre de racines, c'est-à-dire d'arbres de la forêt */
  int *lroot;  /* lroot[i]=i-ème racine de la forêt, i=0..nroot-1 */
  int *height; /* height[i]=hauteur du i-ème arbre, i=0..nroot-1 */
  int *root;   /* root[u]=racine de u (root[u]=u si u racine) */
  int *parent; /* parent[u]=parent de u, -1 si u racine */
  int *nchild; /* nchild[u]=nombre de fils de u */
  int **child; /* child[u][i]=i-ème fils de u, i=0..nchild[u]-1 */
  int *dfs;    /* dfs[u]=ordre dfs de u (dfs[u]=i <=> order[i]=u) */
  int *order;  /* order[i]=u si u est le i-ème sommet dans le parcours pré-fixe */
  int *post;   /* post[i]=u si u est le i-ème sommet dans le parcours post-fixe */
  int *weight; /* weight[u]=nombre de descendents de u, u compris */
  int *depth;  /* depth[u]=profondeur de u dans son arbre */
  int *light;  /* light[u]=plus proche ancêtre (propre) léger de u, =-1 si u racine */ 
  int *apex;   /* apex[u]=premier sommet de la branche lourde de u */
  int *heavy;  /* heavy[u]=1 ssi l'arête entre u et son parent est lourde, =0 si u racine */
} tree;


tree *new_tree(const int n){
/*
  Renvoie un objet de type "tree", un arbre enraciné à n sommets. Les
  champs sont initialisés à leurs valeurs par défaut. Si n>0, alors
  les tableaux simple de taille n sont alloués, mais pas les doubles
  tableaux comme child.
*/

  NALLOC(tree,T,1);
  T->n=imax(n,0);
  T->nroot=-1;
  T->lroot=NULL;
  T->height=NULL;
  T->root=NULL;
  T->parent=NULL;
  T->nchild=NULL;
  T->child=NULL;
  T->dfs=NULL;
  T->order=NULL;
  T->post=NULL;
  T->weight=NULL;
  T->depth=NULL;
  T->light=NULL;
  T->apex=NULL;
  T->heavy=NULL;

  if(n>0){
    ALLOC(T->lroot,n);
    ALLOC(T->height,n);
    ALLOC(T->root,n);
    ALLOC(T->parent,n);
    ALLOC(T->nchild,n);
    ALLOC(T->child,n);
    ALLOC(T->dfs,n);
    ALLOC(T->order,n);
    ALLOC(T->post,n);
    ALLOC(T->weight,n);
    ALLOC(T->depth,n);
    ALLOC(T->light,n);
    ALLOC(T->apex,n);
    ALLOC(T->heavy,n);
  }

  return T;
}


void free_tree(tree *T){
/*
  Libère un arbre T et tous ses tableaux.
*/
  if(T){
    free(T->lroot);
    free(T->height);
    free(T->root);
    free(T->parent);
    free(T->nchild);
    FREE2(T->child,T->n);
    free(T->dfs);
    free(T->order);
    free(T->post);
    free(T->weight);
    free(T->depth);
    free(T->light);
    free(T->apex);
    free(T->heavy);
    free(T);
  }
}


enum{
  TREE_PARENT_COPY = 0x0001, // duplique le tableau d'origine (sinon affecte le pointeur)

  TREE_NCHILD_FREE = 0x0002, // libère nchild[]
  TREE_CHILD_FREE  = 0x0004, // libère child[][]
  TREE_LROOT_FREE  = 0x0008, // libère lroot[]

  TREE_DFS         = 0x0010, // calcule dfs[]
  TREE_ORDER       = 0x0020, // calcule order[]
  TREE_DEPTH       = 0x0040, // calcule depth[]
  TREE_HEIGHT      = 0x0080, // calcule height[]
  TREE_POST        = 0x0100, // calcule post[]

  TREE_WEIGHT      = 0x0200, // calcule weight[]
  TREE_LIGHT       = 0x0400, // calcule light[]
  TREE_APEX        = 0x0800, // calcule apex[]
  TREE_HEAVY       = 0x1000, // calcule heavy[]
};


tree *MakeTree(int *P,const int n,const unsigned code){
/*

  NON UTILISEE, NON TESTEE

  Construit une forêt enracinée (structure tree) à partir d'une
  relation de parentée P à n sommets (tableau d'entiers P[u]=père(u)
  ou -1 s'il n'en a pas). Certaines tables de la structure sont
  initialisées ou pas suivant la valeur binaire de code. Pour tester
  les bits de code il faut utiliser l'enum TREE_xxx ci-dessus. Les
  complexités en temps et en espace sont en O(n).

  o Plus précisément, les tables parent[], nchild[], child[][] et
    lroot[] sont toujours calculées. Elles peuvent être libérées
    suivant les bits de code, sauf parent[].

  o Les tables dfs[], order[], depth[], height[] et post[] sont
    calculées ou pas suivant les bits de code, mais si l'une d'elles
    est calculée alors toutes le sont (car c'est un même parcours).

  o Enfin, les tables weight[], light[], apex[] et heavy[] peuvent
    être calculées ou pas suivant les bits de code. Ces dernières
    tables vont calculer de manière intermédiaire dfs, order, depth,
    height et post puis les supprimer (éventuellement).

  A part les tables toujours calculées, toutes les tables
  intermédiaires calculées sont supprimées sauf si le bit de code
  correspondant indique le contraire.
*/
  if(n<=0) return NULL; /* rien à faire ou problème */

  tree *T=new_tree(n); /* T->n=n>0 */
  int u,p,i;

  /* copie le tableau P ou pas */
  
  if(code&TREE_PARENT_COPY) ALLOCZ(T->parent,n,P[_i]);
  else T->parent=P;

  /* calcule le nombre de racines et le nombre de fils pour chaque
     sommets: T->nroot, T->nchild, T->root */

  ALLOCZ(T->nchild,n,0);
  for(u=0;u<n;u++){
    p=T->parent[u]; /* p=père(u) */
    if(p<0){
      T->root[u]=u;
      T->nroot++;
    }else T->nchild[u]++;
  }

  /* alloue les tableaux T->child et remet à zéro T->nchild */

  for(u=0;u<n;u++)
    if(T->nchild[u]){
      ALLOC(T->child[u],T->nchild[u]);
      T->nchild[u]=0;
    }else T->child[u]=NULL;
  
  /* remplit les tableaux T->child et rétablit T->nchild */
  /* remplit aussi la liste des racines T->lroot */

  ALLOC(T->lroot,T->nroot);
  T->nroot=0; /* on va le recalculer */
  for(u=0;u<n;u++){
    p=T->parent[u]; /* p=père(u) */
    if(p<0) T->lroot[T->nroot++]=u; /* ajoute u aux racines */
    else T->child[p][T->nchild[p]++]=p; /* ajoute u aux fils de p */
  }

  /* parcoure la forêt si nécessaire */
  /* calcule T->dfs, T->order, T->depth, T->height, T->post */

  if((code>>4)==0) goto maketree_fin;
  // vrai ssi l'un des bits de code après le 4e est mis
  // NB: <=> code&(TREE_DFS|TREE_ORDER|TREE_DEPTH|...|TREE_HEAVY)

  NALLOC(int,pile,n); /* pile */
  NALLOC(int,next,n); /* next[u]=compteur courant du prochain voisin de u à visiter */
  int sp=-1; /* sp=sommet de la pile, pile[sp]=dernier élément empilé */
  int v,dfs=0; /* dfs=date de première visite */
  p=0; /* ordre post-fixe */

  for(i=0;i<T->nroot;i++){ /* pour chaque racine faire ... */

    u=T->root[i];
    pile[++sp]=u; /* empile la racine i */
    next[u]=0; /* 1er voisin de u à visiter */
    T->dfs[u]=dfs; /* un sommet visité */
    T->order[dfs++]=u; /* ordre pré-fixe des sommets */
    T->depth[u]=0; /* profondeur du sommet u */
    T->height[i]=0; /* hauteur de la i-ème racine */

    while(sp>=0){ /* tant que la pile n'est pas vide */
      u=pile[sp]; /* u=sommet courant sur la pile */
      if(next[u]<T->nchild[u]){ /* on visite le voisin de u d'indice next[u] */
	v=T->child[u][next[u]++]; /* v=voisin de u à empiler */
	pile[++sp]=v; /* on empile le voisin */
	T->dfs[v]=dfs; /* date de première visite de v */
	T->order[dfs++]=v; /* ordre du sommet */
	T->depth[v]=T->depth[u]+1; /* hauteur de v */
	T->height[i]=imax(T->height[i],T->depth[u]);
      }else{ /* on a visité tous les voisins de u */
	T->post[p++]=pile[sp--]; /* on dépile u, ordre post-fixe des sommets */
      }
    }
    
  }
  
  free(pile);
  free(next);

  /* calcule le poids des sommets */

  if(code&(TREE_WEIGHT|TREE_LIGHT)){
    ALLOCZ(T->weight,n,1); /* tout le monde a poids 1 au départ */
    for(i=0;i<n;i++){
      u=T->post[i]; /* parcours post-fixe */
      p=T->parent[u]; /* p=père de u */
      if(p>=0) T->weight[p] += T->weight[u]; /* le père reçoit le poids de son fils */
    }
  }

  /* calcule l'ancêtre léger de chaque sommets (il faut les poids) */

  if(code&TREE_LIGHT){ // à finir
    ALLOC(T->light,n);
    for(i=0;i<n;i++){
      u=T->order[i]; /* parcours préfixe */
      p=T->parent[u]; /* p=père de u */
    }
  }

  /* calcule apex[] */

  if(code&TREE_APEX){ // à finir
    ALLOC(T->apex,n);
    ;
  }

  /* calcule heavy[] */

  if(code&TREE_HEAVY){ // à finir
    ALLOC(T->heavy,n);
    ;
  }

  // libérations si nécesaire

  if(!(code&TREE_DFS))   { free(T->dfs);      T->dfs   =NULL; }
  if(!(code&TREE_ORDER)) { free(T->order);    T->order =NULL; }
  if(!(code&TREE_DEPTH)) { free(T->depth);    T->depth =NULL; }
  if(!(code&TREE_HEIGHT)){ free(T->height);   T->height=NULL; }
  if(!(code&TREE_POST))  { free(T->post);     T->post  =NULL; }
  if(!(code&TREE_WEIGHT)){ free(T->weight);   T->weight=NULL; }
  if(!(code&TREE_LIGHT)) { free(T->light);    T->light =NULL; }
  if(!(code&TREE_APEX))  { free(T->apex);     T->apex  =NULL; }
  if(!(code&TREE_HEAVY)) { free(T->heavy);    T->heavy =NULL; }

 maketree_fin:

  // libérations si nécesaire
  
  if(code&TREE_NCHILD_FREE){ free(T->nchild);   T->nchild=NULL; }
  if(code&TREE_CHILD_FREE) { FREE2(T->child,n); T->child =NULL; }
  if(code&TREE_LROOT_FREE) { free(T->lroot);    T->lroot =NULL; }

  return T;
}


/***

  Fonctions de comparaisons pour les tris. Attention ! return x-y; ne
  marche pas toujours, même pour comparer des entiers. Par exemple,
  x=INT_MAX et y=INT_MIN provoque un dépassement avec x-y. Il est
  préférable d'utiliser return (x>y) - (x<y); qui marche toujours.

***/


int fcmp_int(const void *P,const void *Q){
/* Compare deux entiers, pour qsort(). */
  const int p=*(int*)P;
  const int q=*(int*)Q;
  return (p>q) - (p<q);
}


int fcmp_int_inv(const void *P,const void *Q){
/* Comme fcmp_int(), mais dans l'ordre inverse. */
  return fcmp_int(Q,P);
}


int fcmp_double(const void *P,const void *Q){
/* Compare deux doubles, pour qsort(). */
  const double p=*(double*)P;
  const double q=*(double*)Q;
  return (p>q) - (p<q);
}


int fcmp_point(const void *P,const void *Q){
/* Compare deux points, pour qsort(). */
  if(((point*)P)->x < ((point*)Q)->x) return -1;
  if(((point*)P)->x > ((point*)Q)->x) return 1;
  if(((point*)P)->y < ((point*)Q)->y) return -1;
  if(((point*)P)->y > ((point*)Q)->y) return 1;
  return 0;
}


int fcmp_profile(const void *P,const void *Q){
/*
  Compare deux profiles, pour qsort(). Les profiles de plus grande
  longueur sont classés avant les plus courts, ceux-ci étant plus
  discriminant.
*/
  int* A=*(int**)P;
  int* B=*(int**)Q;

  if(*A>*B) return -1; /* si long(A)>long(B), alors A<B */
  if(*A<*B) return 1; /* si long(A)<long(B), alors A>B */
  /* ici, profiles de même longueur n=A[0] */

  int u=2,n=*A; /* surtout ne pas utiliser A[1] */

  for(;u<n;u++){
    if(A[u]<B[u]) return -1;
    if(A[u]>B[u]) return 1;
  }

  return 0;
}


int fcmp_graphid(const void *P,const void *Q){
/* Compare les identifiants de deux graphes. Sert pour qsort() et
   bsearch(). Ici, P et Q sont des pointeurs de (graph*). */
  const int p=(*(graph**)P)->id;
  const int q=(*(graph**)Q)->id;
  return (p>q) - (p<q);
}


int fcmp_stretch(const void *P,const void *Q){
/*
  Compare les ratios P.x/P.y et Q.x/Q.y, dans le cas de ratios
  irréductibles. Sert pour routing_test().
*/
  const int p=((triplet*)P)->x*((triplet*)Q)->y;
  const int q=((triplet*)P)->y*((triplet*)Q)->x;
  return (p>q) - (p<q);
}


int fcmp_tabint(const void *P,const void *Q){
/*
  Compare T[*P] avec T[*Q] où T est un tableau d'entiers qui doit être
  initialisé lors d'un premier appel avec fcmp_tabint(NULL,T).
  Attention ! les valeurs du tableau initial à trier (*P et *Q) doivent
  être des indices dans T, donc dans [0,|T|[.
*/
  static int *T; // tableau global mais local à la fcmp_tabint()
  if(P==NULL){ T=(int*)Q; return 0; } // fixe le tableau T
  
  const int p=T[*(int*)P];
  const int q=T[*(int*)Q];
  return (p>q) - (p<q);
}


int fcmp_tabinteq(const void *P,const void *Q){
/*
  Comme fcmp_tabint(P,Q) qu'en cas d'égalité (soit si *P==*Q) on
  renvoie la comparaison entre *P et *Q et non 0. Donc le tri
  s'effectue d'abord selon T[], puis selon les indices de T en cas
  d'égalité.
*/
  static int *T; // tableau global mais local à la fcmp_tabint()
  if(P==NULL){ T=(int*)Q; return 0; } // fixe le tableau T
  
  const int p=*(int*)P;
  const int q=*(int*)Q;
  if(T[p]==T[q]) return p-q; 
  return T[p]-T[q];
}


/* code pour ReadRange() et InRange() */
enum{
  R_EQ,   // code =x
  R_INF,  // code <x
  R_SUP,  // code >x
  R_INT,  // code x-y
  R_TRUE  // code t
};


int ReadRange(char *s,int *R){
/*
  Lit une chaîne de caractères décrivant un intervalle de valeurs
  entières, et renvoie dans le tableau d'entiers R les valeurs et les
  codes correspondant pour que la fonction InRange(x,R) puisse
  fonctionner. En quelque sorte cette fonction prépare la fonction
  InRange(). On ignore les caractères non reconnus (pas d'erreur). On
  renvoie le nombre d'opérations décodées, c'est-à-dire le nombre de
  valeurs écrites dans le tableau R, nombre qui est aussi écrit dans
  la première entrée de R.

  Ex: s="1,-3,5-7,>30,<50" (on interprète les "," comme des "ou")
  => R={12,R_EQ,1,R_EQ,-3,R_INT,5,7,R_SUP,30,R_INF,50} (12=taille(R))

  La valeur des codes d'opérations (R_EQ, R_INF, ...) sont données par
  l'enum ci-dessus.
*/
  if(R==NULL) return 0;
  if(s==NULL){ R[0]=R_EQ; return 0; }

  int i,r,p,x,start,c;
  i=x=c=0;
  r=start=p=1;

  /* r=indice de R[] */
  /* i=indice de s[] */
  /* x=valeur entière lue */
  /* c=1 ssi le code d'opération a été détecté */
  /* start=1 ssi on va commencer à lire un entier */
  /* p=1 ou -1, signe de x */

  while(s[i]!='\0'){
    if(s[i]=='='){ R[r++]=R_EQ; c=start=p=1; }
    if(s[i]=='<'){ R[r++]=R_INF; c=start=p=1; }
    if(s[i]=='>'){ R[r++]=R_SUP; c=start=p=1; }
    if(s[i]=='-'){
      if(start) p=-p;
      else{ R[r++]=R_INT; R[r++]=x; c=start=p=1; }
    }
    if(s[i]=='t'){ x=R_TRUE; c=r=1; break; } /* t=true, pour avoir false faire "not" et "t" */
    if(s[i]=='p') PVALUE=1; /* pas de code pour "p" */
    if(s[i]==','){
      if(c==0) R[r++]=R_EQ; /* code '=' par défaut */
      R[r++]=x; c=0; start=p=1;
    }
    if(('0'<=s[i])&&(s[i]<='9')){
      if(start) start=0;
      x=x*10+p*(s[i]-'0'); /* on lit x en base 10 en tenant compte du signe p */
    }
    if(start) x=0;
    i++;
  }

  if(PVALUE==i){ x=R_TRUE;c=1; } /* si s="p", alors comme "t" */
  if(c==0) R[r++]=R_EQ;
  R[r++]=x; /* on ajoute le dernier opérande */
  R[0]=r;
  return r;
}


int InRange(const int x,const int* R){
/*
  Détermine si x appartient aux valeurs décrites par le "range" R.
  R[0] est la taille de R, R[0] compris.
*/
  int i,n,t;
  n=R[t=0]; /* n=taille(R) */
  i=1; /* commence à lire R[1] */
  CVALUE=x;

  while(i<n){
    switch(R[i++]){ /* lit le code d'opération */
    case R_EQ  : t=(x==R[i++]); break;
    case R_INF : t=(x<R[i++]); break;
    case R_SUP : t=(x>R[i++]); break;
    case R_INT : t=((R[i]<=x)&&(x<=R[i+1])); i+=2; break;
    case R_TRUE: return 1;
    default: Erreur(16); /* ne devrait jamais se produire */
    }
    if(t) break;
  }
  return t;
}


/***********************************

       ROUTINES SUR
       LES GRAPHES

***********************************/


static inline void degres_zero(graph *G){
/*
  Met à zéro tous les degrés d'un graphe. C'est très utile pour se
  servir par exemple des macros ADD_EDGE() et ADD_ARC().
*/
  memset(G->d,0,G->n*sizeof(int));
}


int NbEdges(graph *G){
/*
  Retourne le nombre d'arêtes d'un graphe symétrique G ou bien le
  champs G->m s'il est positif. Si G->m<0, alors G->m est mis à jour à
  partir de la somme des G->d[i].
*/
  int m=G->m;
  if(m<0){
    int i;
    const int n=G->n;
    for(i=m=0;i<n;i++)
      m += G->d[i];
    G->m=(m>>=1);
  }
  return m;
}


int Degree(const graph *G,const int maxi){
/*
  Renvoie le degré maximum (si maxi=1) ou minimum (si maxi=0) d'un
  graphe G. On renvoie -1 si G est nul, n'a pas de sommet ou est une
  famille de graphes.
*/
  if((G==NULL)||(G->n<=0)||(G->f>0)) return -1;
  const int n=G->n;
  int i=1,d=G->d[0];
  if(maxi) for(;i<n;i++) d=imax(d,G->d[i]);
  else for(;i<n;i++) d=imin(d,G->d[i]);
  return d;
}


void PrintGraphList(const graph *G){
/*
  Affiche le graphe G sous la forme d'une liste d'adjacence. Tient
  compte de SHIFT et de VERTEX0.
*/
  if(G==NULL){ printf("NULL\n"); return; }
  const int n=(VERTEX0<0)?G->n:VERTEX0+1;
  int u,d,i;

  for(u=(VERTEX0<0)?0:VERTEX0;u<n;u++){
    printf("%i:",u+SHIFT);
    for(i=0,d=G->d[u];i<d;i++){
      printf(" %i",G->L[u][i]+SHIFT);
    }
    printf("\n");
  }
  return;
}


void PrintGraphMatrix(graph *G){
/*
  Affiche le graphe G sous la forme d'une matrice d'adjacence complète
  ou triangulaire supérieure (en tennant compte du FORMAT, smatrix ou
  matrix). La complexité en espace est seulement de O(n).
*/
  int u,d,i,z,t;
  const int n=G->n;

  NALLOCZ(int,M,n,0);
  t=(FORMAT==F_smatrix);

  for(u=z=0;u<n;u++){
    if(t) z=u;
    for(i=0,d=G->d[u];i<d;M[G->L[u][i++]]=1);
    for(i=0;i<n;i++)
      if(i<z) printf(" ");
      else printf("%c",'0'+M[i]);
    for(i=0;i<d;M[G->L[u][i++]]=0); /* remet rapidement M[] tout à 0 */
    printf("\n");
  }

  free(M);
  return;
}


void PrintPath(graph *G,path *P)
/*
  Affiche le chemin P d'un graphe G.
  Sert pour le débugage.
*/
{
  if((G==NULL)||(P==NULL))
    printf("NULL\n");
  else{
    int i,j,u,d;
    for(i=0;i<P->n;i++)
      if(P->V[P->P[i]]!=i) break;
    if(i<P->n) goto error;
    for(u=0;u<G->n;u++)
      if((P->V[u]>=0)&&(P->P[P->V[u]]!=u)) break;
    if(u<G->n) goto error;
    printf("P->aux:");
    if(P->aux==NULL) printf(" NULL\n");
    else{
      printf("\n");
      for(i=0;i<P->n;i++){
	u=P->P[i];
	d=P->aux[u][0];
	printf("  %i:",u);
	for(j=1;j<=d;j++){
	  printf(" %i",P->aux[u][j]);
	}
	printf("\n");
      }
    }
  }
  return;
  
 error:
  printf("Chemin incohérent.\n");
  return;
}


int *SortGraph(graph *G,const int code){
/*
  Force le tri (même si G->sort=1) des listes d'adjacence d'un graphe
  G, c'est-à-dire pour chaque sommet u, G->L[u] est une liste
  d'entiers triés par ordre croissant. Le champs G->sort est mis à
  jour.  L'algorithme effectue un simple appel à qsort(). Sa
  complexité est à peu près en O(n+m*log(m/n)).

  Si code=0, on s'arrête après l'étape du tri. Sinon, on lance une
  étape de vérification du graphe: présence de multi-arêtes, de
  boucles, etc.

  Le temps de la vérification est comparable à celui du tri. Le
  résultat de la vérification est un tableau de statistiques S de
  taille fixe (déclaré en static qui ne doit pas être libéré par
  l'appelant), ayant la signification suivante:

    S[0]=nombre de boucles
    S[1]=nombre de multi-arcs
    S[2]=nombre d'arcs (avec multi-arcs et boucles)
    S[3]=nombre d'adjacence non-symétriques
    S[4]=nombre de voisins d'ID < 0
    S[5]=nombre de voisins d'ID ≥ n
    S[6]=1 ssi G est simple et non-orienté
    S[7]=degré maximum
    S[8]=degré minimum
    S[9]=nombre de sommets isolés

  À l'issue de la vérification, G->sym est mise à jour.
*/

  if(G==NULL) return NULL;
  const int n=G->n;
  int u;

  /* trie G */
  for(u=0;u<n;u++)
    QSORT(G->L[u],G->d[u],fcmp_int);

  G->sort=1;
  if((code==0)||(n==0)) return NULL;
  
  /* statistiques sur G */
  static int S[10]; /* static est important, car on fait return S */
  int v,i,d,w;
  
  memset(S,0,sizeof(S)); /* initialise les stats à 0, NB: sizeof(S)=40 */
  S[7]=S[8]=G->d[0]; /* il faut G->d<>NULL */
  
  for(u=0;u<n;u++){ /* parcoure le graphe */
    d=G->d[u]; S[7]=imax(S[7],d); S[8]=imin(S[8],d);
    S[9] += (d==0); /* un sommet isolé */
    S[2]+=d; /* ajoute le nombre de voisins */
    w=-1; /* w=voisin précédant le voisin courant v */
    for(i=0;i<d;i++){ /* pour chaque voisin */
      v=G->L[u][i];
      if(u==v) S[0]++; /* une boucle */
      if(v==w) S[1]++; /* une multi-arête */
      w=v; /* mémorise le dernier voisin rencontré */
      if(v<0){ S[4]++; continue; } /* un voisin négatif */
      if(v>=n){ S[5]++; continue; } /* un voisin trop grand */
      if(bsearch(&u,G->L[v],G->d[v],sizeof(int),fcmp_int)==NULL) S[3]++; /* un arc asymétrique */
    }
  }

  S[6]=((S[0]+S[1]+S[3]+S[4]+S[5])==0); /* vrai ssi G simple et non-orienté */
  G->sym=(S[3]==0);
  return S;
}


void PrintGraph(graph *G){
/*
  Affiche un graphe ou une famille de graphes au format standard sous
  forme compacte. Utilise WIDTH. Effet de bord: le (ou les) graphes
  sont triés par ordre croissant, et donc G->sort=1 en sortie. Si le
  graphe est asymétrique, des sommets peuvent être affichés comme
  sommets isolés alors qu'ils ne le sont pas.
*/
  if(G==NULL){ printf("NULL\n"); return; }

  int u,v,i,k,n,ligne,nk=(G->f>0)?G->f:1;
  graph *H;
  int *P;

  for(k=0;k<nk;k++){

    if(G->f>0){
      H=G->G[k];
      printf("[%i]",H->id);
    }else H=G;

    SortGraph(H,0); // ordre croissant
    n=H->n;
    ALLOCZ(P,n,0);
    i=u=ligne=0;
    v=-1;

    while(i<n){
      /* si u==i, alors u=tête d'un nouveau chemin */
      while((v<u)&&(P[u]<H->d[u])) v=H->L[u][P[u]++];
      if(v<u){ /* on a pas trouvé v => fin d'un chemin */
	if(H->d[u]==0){ /* cas des sommets isolés */
	  printf(" %i",u);
	  if(++ligne==WIDTH){ printf("\n"); ligne=0; }
	}
	u=(i+=(u==i));
	v=-1;
      }
      else{ /* u a un voisin v>u */
	if((u==i)||(ligne==0)) printf(" %i",u); /* on affiche la tête */
	printf("-%s%i",(H->sym)?"":">",v); /* on affiche -v ou ->v */
	if(++ligne==WIDTH){ printf("\n"); ligne=0; }
	u=v; /* on continu avec v */
	v=-1;
      }
    } /* fin du while */

    if(ligne>0) printf("\n"); /* newline si fini avant la fin de ligne */
    free(P);
  }

  G->sort=1; /* effet de bord */
  return;
}


void GraphRealloc(graph *G,const int *D){
/*
  Redimensionne le graphe G à G->n sommets suivant le tableau de degré
  D. On réajuste en premier les tableaux G->d et G->L pour qu'ils
  aient une taille G->n, puis on réajuste les listes d'adjacences des
  sommets de G suivant le tableau des degrés D (qui doit être de
  taille au moins G->n). Si D[u] est plus petit que G->d[u], alors la
  liste G->L[u] est tronquée. Si D[u] est plus grand que G->d[u],
  alors G->L[u] est réajusté. Le degré G->d[u] est initialisé au
  minimum de G->d[u] et D[u]. NB: le nombre d'arêtes G->m, qui a pu
  changer, est réinitialisé à -1. G->sort n'est pas changé car l'ordre
  des listes G->L n'est pas modifié.

  Pour plonger G dans un graphe complet faire:
    NALLOCZ(int,D,G->n,G->n-1);
    GraphRealloc(G,D);
    free(D);
*/
  const int n=G->n;
  int u,d;
  for(u=0;u<n;u++){
    d=D[u];
    REALLOC(G->L[u],d);
    G->d[u]=imin(G->d[u],d);
  }

  /* Il ne faut pas réajuster G->d et G->L avant la boucle for(u=...)
     car potentiellement on libère G->d et G->L. Or il est possible
     d'avoir D=G->d. */

  REALLOC(G->d,n);
  REALLOC(G->L,n);
  G->m=-1; /* le nombre d'arêtes n'est plus à jour */
  return;
}


graph *new_subgraph(const graph *G){
/*
  Renvoie un nouveau graphe R vide de manière similaire à
  R=new_graph(G->n), mais en plus dimensionne chaque R->L[u] à
  G->d[u]. On initialise aussi R->d[u] à 0. Le graphe renvoyé est donc
  un sous-graphe couvrant de G sans aucune arête.
*/
  if((G==NULL)||(G->n<0)) return NULL;

  const int n=G->n;
  graph *R=new_graph(n);
  int u;

  for(u=0;u<n;u++) ALLOC(R->L[u],G->d[u]);
  degres_zero(R);

  return R;
}


graph *new_fullgraph(const int n){
/*
  Renvoie un graphe G comme new_graph(n), mais en plus alloue G->L[u]
  de taille max{n-1,1}, et initialise G->d[u]=0 pour tous les sommets
  u. Une fois le graphe construit, on peut rédimensionner le graphe
  grâce à GraphRealloc, comme dans l'exemple:

    graph *G=new_fullgraph(n);
      ...
      ADD_EDGE(G,u1,v1);
      ADD_EDGE(G,u2,v2);
      ...
    GraphRealloc(G,G->d);
      ...
    free_graph(G);
*/

  if(n<1) return NULL;
  graph *G=new_graph(n);
  const int n1=imax(n-1,1);
  int u;
  for(u=0;u<n;u++) ALLOC(G->L[u],n1);
  degres_zero(G);
  
  return G;
}


graph *ExtractSubgraph(const graph *G,const int *T,const int k,const int code){
/*
  Construit, à partir d'un graphe G et d'une liste T de k sommets, un
  nouveau graphe S (renvoyé par la fonction) correspondant au
  sous-graphe de G induit par les sommets de T (code=1) ou de G\T (si
  code=0). Les sommets du graphe S sont dans [0,k[ (ou [0,n-k[ si
  code=0).

  On peut ainsi faire une copie C du graphe G simplement en faisant:

    graph *C=ExtractSubgraph(G,NULL,0,0);

  Effet de bord: S->pint1 est alloué si T<>NULL. Dans ce cas on
  renvoie dans S->pint1 un tableau X de taille G->n indiquant la
  renumérotation de G: pour tout sommet u de G (u dans [0,G->n[)
  S->pint1[u]=0 si u est un sommet abscent de S et S->pint1[u]=d>0 si
  u est numéroté d-1>=0 dans S. Le nombre d'arêtes S->m du graphe S
  renvoyé est à jour. L'ordre relatif des listes de G est préservé. En
  particulier, si G->sort>0, alors le sous-graphe sera aussi
  trié. G->sym est aussi copié.
*/
  
  if(G==NULL) return NULL;
  const int n=G->n;
  int u,v,d,i,s,ns,m;

  NALLOC(int,X,n);
  for(u=1-code,i=0;i<n;i++) X[i]=u;
  if(T) for(i=0;i<k;i++) X[T[i]] ^= 1;
  for(i=d=0;i<n;i++) if(X[i]) X[i]=++d; 
  /* ici X[i]=0 si i est un sommet à supprimer */
  /* ici X[i]=d (>0) si i doit être renuméroté en d-1>=0 */

  ns=(code)?k:n-k;
  graph* S=new_fullgraph(ns);

  for(s=u=m=0;u<n;u++)
    if(X[u]){ /* si u existe, s=X[u]-1 */
      d=G->d[u];
      for(i=0;i<d;i++){
	v=G->L[u][i];
	if(X[v]){ m++; ADD_ARC(S,s,X[v]-1); } /* si v existe */
      }
      s++;
    }

  /* réduit la taille des listes */
  GraphRealloc(S,S->d);

  S->pint1=X;
  S->sort=G->sort;
  S->sym=G->sym;
  S->m=(m>>1);
  return S;
}


graph *List2Graph(list *L,const int code){
/*
  Retourne un graphe G simple à partir d'un graphe défini par une
  liste L de codes (voir File2List() pour le codage précis du type
  "list"). Certaines opérations sont effectuées sur L en fonction de
  la valeur binaire de code:

  - code&1 =1: optimisation des listes du graphe (tri par ordre croissant)
           =0: sans optimisation
  - code&2 =1: auto-détection du shift dans L (pour "load file")
           =0: pas d'auto-détection du shift
  - code&4 =1: gestion d'un sous-graphe (V,NF) => code&2=0
           =0: pas de sous-graphe

  Les codes suivants servent à List2Family():

  - code&8 =1: tri de la famille suivant les identifiants (sert pour List2Family)
           =0: pas de tri de la famille
  - code&16=1: ne libère pas la liste L (sert pour List2Family)
           =0: libère la liste L
  - code&32=1: renvoie toujours un graphe, le 1er si c'est une famille
           =0: renvoie une famille si c'est une famille

  Pour calculer le graphe (et sa liste d'adjacence) on effectue
  plusieurs passes sur L: une passe pour déterminer n; une autre pour
  calculer les degrés des sommets; et une 3e pour remplir G et pour
  éventuellement libérer la liste L.
*/
  if(L==NULL) return NULL; /* si pas de cellule, ne rien faire */

  int u,v,x,n,*D;
  graph *G;
  list *p;

  u=INT_MAX;
  if(code&4){ /* si sous-graphe définit par (V,NF) */
    p=L;
    while(p){
      p->item=V[p->item]-SHIFT;
      p=p->next;
    }
    n=NF; /* on connaît n */
  }
  else{ /* sinon, on calcule n, et on lit les valeurs min (=u) et
	   valeur max (=v) de L */
    p=L; v=0;
    while(p){
      x=p->item;
      if(x<u) u=x;
      if(x>v) v=x;
      p=p->next;
    }
    if(code&2){ /* on décale les valeurs dans [0,n[ */
      p=L;
      while(p){
	p->item -= u;
	p=p->next;
      }
      n=v+1-u;
    }else{
      if((u<0)||(v<0)) Erreur(22); /* il ne devrait pas avoir ici de valeur < 0 */
      n=v+1;
    }
  }

  ALLOCZ(D,n,0);

  /* on lit les degrés (sortant) des sommets, et les met dans le
     tableau D. NB: la variable u n'est pas initialisé, car on passe
     toujours d'abord par un item de type T_NODE */
  
  p=L; x=1; /* x=1 ssi on n'est PAS dans un groupe */
  while(p){
    v=p->item;
    if(p->type==T_NODE) x=1;
    else if(p->type==T_EDGE) { D[u]++; D[v]++; }      /* u-v */
    else if(p->type==T_OPENE){ D[u]++; D[v]++; x=0; } /* u-(v */
    else if(p->type==T_ARC)    D[u]++;                /* u->v */
    else if(p->type==T_OPENA){ D[u]++; x=0; }         /* u->(v */
    if(x) u=v;
    p=p->next;
  }

  /* initialise la liste d'adjacence G. On se sert plus tard de D[u]
     pour indiquer la prochaine position libre dans G[u][]. */

  G=new_graph(n); /* G->n=n, alloue G->d et G->L */
  for(u=0;u<n;u++){
    ALLOC(G->L[u],D[u]); /* alloue une liste pour chaque sommet */
    G->d[u]=D[u]; /* G->d[u]=deg(u) */
    D[u]=0; /* prochaine position libre dans G[u] */
  }

  /* Remplit G. On met aussi à jour G->sym (orienté ou pas). On
     pourrait tester à la volée si les listes sont triées et mettre à
     jour G->sort. */
  
  p=L; x=1; /* x=1 ssi on n'est PAS dans un groupe */
  while(p){
    v=p->item;
    if(p->type==T_NODE) x=1;
    else if(p->type==T_EDGE) { G->L[u][D[u]++]=v; G->L[v][D[v]++]=u; }      /* u-v */
    else if(p->type==T_OPENE){ G->L[u][D[u]++]=v; G->L[v][D[v]++]=u; x=0; } /* u-(v */
    else if(p->type==T_ARC)  { G->L[u][D[u]++]=v; G->sym=0; }               /* u->v */
    else if(p->type==T_OPENA){ G->L[u][D[u]++]=v; G->sym=0; x=0; }          /* u->(v */
    if(x) u=v;
    p=p->next;
  }

  /* libère L si bit-4 à 1 */
  if(!(code&16)){
    p=L;
    while(p){
      L=p;
      p=p->next;
      free(L);
    }
  }
  
  free(D); /* plus besoin de D */
  if(code&1) SortGraph(G,0);
  return G;
}


graph *List2Family(list *L,const int code){
/*
  Transforme une liste en famille de graphes.  Si L représente un
  graphe simple (pas de type T_NB ou T_ID), alors un graphe simple est
  retournée (plutôt qu'une famille à un seul élément). Donc,
  List2Family() généralise List2Graph(). On utilise List2Graph() comme
  sous-routine. Pour "code" voir List2Graph().
    
  Effet de bord:
  - la famille est triée par ID croissant si code&8=1
  - la liste L est libérée si code&16=0
  - on retourne un graphe si code&32=1 (plutôt qu'une famille)
*/
  if(L==NULL) return NULL; /* liste vide */
  if(L->type!=T_NB) return List2Graph(L,code); /* si graphe */

  /* ici on a donc une famille */
  int n=L->item; /* nb de graphes dans la liste */
  list *T;
  
  if(n<=0){ /* famille vide */
    if(code&16) /* libère éventuellement L */
      while(L){
	T=L->next; /* ici L<>NULL */
	free(L);
	L=T;
      }
    return NULL;
  }

  int i,id;
  graph *F=new_graph(0);
  list *P;

  F->f=n;
  ALLOC(F->G,n); /* F->G[0..n[: tableau de n pointeurs de graphes */
  T=L; L=L->next;
  if(!(code&16)) free(T); /* on libère l'élément (n,T_NB) */

  /* ici L=début du 1er graphe de la famille */
  for(i=0;i<n;i++){ /* pour chaque graphe */
    /* ici L pointe sur un élément (id,T_ID) */
    if((L==NULL)||(L->type!=T_ID)) Erreur(13);
    id=L->item; /* identifiant du graph i */
    T=L->next;
    if((code&16)==0) free(L); /* on libère l'élément (id,T_ID) */
    P=L=T; /* P=L=T=tête courante du graphe i */
    while((L)&&(L->type!=T_ID)){ P=L; L=L->next; } /* cherche la fin du graphe i */
    /* ici le graphe i va de T à P */
    P->next=NULL; /* on coupe la liste */
    F->G[i]=List2Graph(T,code); /* T=liste du graphe i */
    F->G[i]->id=id; /* Attention! F->G[i] n'existe qu'après l'appel à List2Graph() */
    if(code&16) P->next=L; /* recolle la liste si on ne souhaite pas la libérer */
    }

  /* éventuellement trie la famille suivant les IDs */
  if(code&8) QSORT(F->G,F->f,fcmp_graphid);

  /* extrait le premier graphe */
  if(code&32){
    graph *G=ExtractSubgraph(F->G[0],NULL,0,0); /* copie le premier graphe */
    free_graph(F); /* libère complètement la famille F */
    F=G; /* F=premier graphe */
  }

  return F;
}


list *File2List(const char *file){
/*
  Lit le fichier "file" contenant un graphe (ou une famille) au format
  standard, orientés ou non, et retourne le contenu dans une liste.
  Tient compte de -shift mais pas des noms originaux (-label 1). Dans
  le cas d'une famille de graphes, il est possible de spécifier un
  "range" pour "file" avec la forme: "file:range" où "range" est une
  liste de valeurs ayant la même signification que pour "-filter F id
  value". Par exemple, "file:5" spécifie le graphe d'identifiant 5, et
  "file:5-8" est la famille contenant les graphes d'identifiant
  5,6,7,8. Notez que "-:5" est le graphe d'identifiant 5 de la famille
  lue depuis l'entrée standard.

  Chaque élément de la liste est une paire d'entiers (item,type) où
  "type" précise le rôle joué par l'entier "item". Voir l'enum pour
  une description des types.

  Si, par exemple, le graphe est "0-1 2 1-2-3 5->4" alors la liste
  retournée sera { (0,T_NODE), (1,T_EDGE), (2,T_NODE), (1,T_NODE),
  (2,T_EDGE), (3,T_EDGE), (5,T_NODE), (4,T_ARC) }.

  Si le graphe est "0-(2 3) 4 5->(6 7)" alors la liste retournée sera
  { (0,T_NODE), (2,T_OPENE), (3,T_OPENE), (4,T_NODE), (5,T_NODE),
  (6,T_OPENA), (7,T_ARC) }. Autrement dit, T_OPENE ou T_OPENA décrive
  un groupe d'arêtes ou d'arcs.

  NB: "i-j-(k ...)" est correct mais pas "i-(j ...)-k".

  La fonction est généralisée à la lecture d'une famille de graphes.
  Si le fichier contient "[5] 0-1 [8] 0->2-1" alors la liste
  contiendra { (2,T_NB), (5,T_ID), (0,T_NODE), (1,T_EDGE), (8,T_ID),
  (0,T_NODE), (2,T_ARC), (1,T_EDGE) }, où le premier élément (n,T_NB)
  signifie qu'il s'agit d'une famille de n graphes, et où (u,T_ID)
  signifie que u est l'identifiant du graphe à venir.
*/
  FILE *f;
  list *T; /* tête de la liste */
  list *L; /* élément courant */
  list *P; /* sauvegarde le dernier élément (sentinelle qui faudra supprimer) */
  int read=1; /* pour InRange(), par défaut on lit tout */
  char *r=NULL,c[2],*s;
  int range[CMDMAX]={2,R_TRUE}; /* par défaut: range toujours vrai */
  unsigned v; /* valeur lue */
  long p; /* position dans le fichier f */
  int n=0; /* nb de graphes dans la famille */
  int t=-1; /* t<0 si on est pas dans un groupe */

  T=P=L=new_list(); /* crée la liste */

  /* TO DO: si file commence par " ", alors file représente le graphe
     lui-même (pour future option -add/-del). Par exemple, file="
     5-6,7-8-0". Dans ce cas on écrit un fichier temporaire avec "5-6
     7-8-0" et on continue normalement. On détruit ensuite ce
     fichier. NB: le contenu de file est modifié. */
  /*
    char *s;
    if(file[0]==' '){
    for(v=0; v<strlen(file); v++) if(file[v]==',') file[v]=' ';
    s=strdup("/tmp");
    f=fopen(s,"rw");
    fputs(f,file);
    rewind(f); // pas la peine de fermer le fichier
    file=s;
  }
  */
  
  /* ouverture du fichier: file ou file:range */

  f=strcmp(file,"-")? fopen(file,"r"):stdin;
  if(f==NULL){ /* on a pas réussit à ouvrir file */
    fclose(f); /* il faut fermer le fichier, même si c'est NULL ! */
    if((r=strchr(file,':'))==NULL) Erreur(7); /* est-ce file:range ? */
    *r='\0'; /* coupe file en (préfixe,r=range) */
    f=strcmp(file,"-")? fopen(file,"r"):stdin;
    if(f==NULL) Erreur(7); /* on a pas réussit à ouvrir file */
    *r++ = ':'; /* rétablit le nom de fichier original de file */
    ReadRange(r,range); /* lecture du range */
  }

  /* lecture du fichier */

  while(!feof(f)){
    p=ftell(f); /* ftell() vaut toujours -1 si f=stdin (non seekable) */
    fscanf(f,"//%*[^\n]\n"); /* essaye de lire "//" */
    if(ftell(f)>=p+2) continue; /* on a lu au moins 2 caractères -> commentaire */
    fseek(f,p,SEEK_SET);
    if(fscanf(f," [%u]",&v)>0){
      read=InRange(v,range);
      if(read){ L=Insert(P=L,v,T_ID); n++; }
      continue;
    }
    if(read){
      fseek(f,p,SEEK_SET);
      if(fscanf(f,"-%u",&v)>0){ L=Insert(P=L,v,T_EDGE); continue; }
      fseek(f,p,SEEK_SET);
      if(fscanf(f,">%u",&v)>0){ L=Insert(P=L,v,T_ARC); continue; }
      fseek(f,p,SEEK_SET);
      if(fscanf(f,"-( %u",&v)>0){ L=Insert(P=L,v,T_OPENE);
	p=ftell(f);
	if(fscanf(f," %1[)]c",c)>0) t=-1;
	else{ t=T_EDGE; fseek(f,p,SEEK_SET); }
	continue;
      }
      fseek(f,p,SEEK_SET);
      if(fscanf(f,">( %u",&v)>0){ L=Insert(P=L,v,T_OPENA);
	p=ftell(f);
	if(fscanf(f," %1[)]c",c)>0) t=-1;
	else{ t=T_ARC; fseek(f,p,SEEK_SET); }
	continue;
      }
      fseek(f,p,SEEK_SET);
      if(fscanf(f," %u",&v)>0){
	p=ftell(f);
	if(fscanf(f," %1[)]c",c)>0){ L=Insert(P=L,v,t); t=-1; }
	else{ L=Insert(P=L,v,(t<0)?T_NODE:t); fseek(f,p,SEEK_SET); }
	continue;
      }
      
      /* ici on a rien trouvé: est-ce une erreur de format ? */
      fseek(f,p,SEEK_SET);
      s=fgets(c,2,f); /* lit au plus un caractère */
      if((s==NULL)||(c[0]==' ')||(c[0]=='\n')) continue; /* ok si ' ' ou '\n' */
      Erreur(28); /* mauvais format sinon */
    }

    /* ici on a rien trouvé, mais read est faux */
    fseek(f,p,SEEK_SET); /* on a rien trouvé */
    fscanf(f," %*c"); /* lit au moins un caractère */
  }
  
  fclose(f); /* on ferme le fichier */
  free(L); /* supprime le dernier élément (la sentinelle) */
  if(L==T) return NULL; /* si on a lu aucun élément */
  P->next=NULL;

  if(n>0){ /* il s'agit d'une famille */
    /* on ajoute un nouvel élément en tête de la liste indiquant le
       nombre de graphes de la famille */
    L=new_list();
    L->item=n; /* nombre de graphes de la famille */
    L->type=T_NB;
    L->next=T;
    T=L; /* nouvelle tête de liste */
  }

  return T; /* on retourne la tête */
}


graph *File2Graph(const char *file,const int code){
/*
  Renvoie un graphe (ou une famille) à partir d'un fichier. Pour
  "code" voir List2Graph() & List2Family(). La liste intermédiaire
  calculée par File2List() est toujours libérée.
*/
  graph *G=List2Family(File2List(file),code&(63-16)); /* annule le bit-4 */
  if(G==NULL) Erreur(15);
  return G;
}


/***********************************

       ROUTINES EN VRAC

***********************************/


double PosAspect(query* const Q){
/*
  Donne le coefficient par lequel les positions XPOS/YPOS vont être
  multipliées pour le format dot pour permettre une taille de sommets
  raisonable par rapport à la longueur des arêtes. On tient compte de
  Q->n et de BOXX et BOXY.
*/
  double w=C32*sqrt(Q->n); /* la largeur est en sqrt(n) */
  if((BOXX>0)&&(BOXY>0)) w /= imin(BOXX,BOXY);
  if(LABEL>0) w *= 3; /* augmente l'aspect si besoin des LABELs (et POS) */
  return w;
}


void BoundingBox(query* const Q){
/*
  Calcule XMIN,YMIN,XMAX,YMAX des tableaux XPOS/YPOS.  Il faut que
  Q->n > 0 et XPOS,YPOS <> NULL.
*/
  int i;
  XMIN=XMAX=XPOS[0];
  YMIN=YMAX=YPOS[0];
  for(i=1;i<Q->n;i++){
    XMIN=fmin(XMIN,XPOS[i]);
    YMIN=fmin(YMIN,YPOS[i]);
    XMAX=fmax(XMAX,XPOS[i]);
    YMAX=fmax(YMAX,YPOS[i]);
  }
}


#define M_2PI 6.28318530717958647692528676655900576839 // 2π, souvent utilisé


static inline double angle(const double x,const double y){
/*
  Renvoie l'angle de [0,2π[ en radian du point de coordonnées
  cartésiennes (x,y) par rapport à l'axe des abscisses (1,0). NB:
  atan(y/x) donne un angle [-π/2,+π/2] ce qui n'est pas ce que l'on
  veut. On renvoie 0 si (x,y)=(0,0).
*/
  if(x==0){
    if(y>0) return M_PI_2;
    if(y<0) return M_PI_2+M_PI;
    return 0;
  }

  // atan(y/x) renvoie un angle entre -π/2 et +π/2
  // l'angle est correct si x>0 et y>0
  // si x,y de signe contraire alors atan(y/x)=-atan(y/x)
  
  const double a=atan(y/x);

  if(x>0){
    if(y>0) return a;
    return a+M_2PI;
  }
  
  return a+M_PI;
}


int fcmp_angle(const void *P,const void *Q)
/* Compare deux angles <XP,YP> et <XQ,YQ> pour qsort2(). */
{
  const double p=angle(*(double*)P,*(((double*)P)+1));
  const double q=angle(*(double*)Q,*(((double*)Q)+1));
  return (p>q) - (p<q);
}


static inline double det(const double a1,const double a2,
			 const double b1,const double b2){
/*
  Renvoie le déterminant de deux vecteurs colonnes A=(a1,a2) et
  B=(b1,b2). Il vaut 0 si les vecteurs A et B sont colinéaires. Le
  signe de det(A,B) est celui de sin(A,B), le sinus de l'angle entre
  (OA) et (OB), O=(0,0) étant l'origine. Dit autrement, si det(A,B)>0
  alors le point B est au-dessus de la droite (OA), si det(A,B)<0 il
  est en dessous, et si det(A,B)=0, alors il est sur la droite (OA).
*/
  return (a1*b2)-(a2*b1);
}


/*
  Rappels de géométrie 2D:
  (voir aussi distgone() et rlt())

  1. Equation cartésienne d'une droite
     
     D = { (x,y) : a*x + b*y = c }
         avec (a,b) != (0,0)
	 mais a=0 ou b=0 possible

  2. Equation cartésienne d'une droite D passant par deux points
     A=(xa,ya) et B=(xb,yb)

     D: (yb-ya)*x - (xb-xa)*y = xa*yb - ya*xb = det(A,B)

  3. Intersection de deux droites cartésiennes
     D1: a1*x - b1*y = c1
     D2: a2*x - b2*y = c2

     Soient les vecteurs colonnes: A=(a1,a2), B=(b1,b2), C=(c1,c2)
     Si det(A,B)=det(A,C)=0, alors les droites sont confondues
     Si det(A,B)=0 et det(A,C)<>0, alors il n'y a pas d'intersection
     Si det(A,B)<>0, alors il y a une seule intersection (x0,y0):
     x0 = -det(B,C)/det(A,B) et y0 = -det(A,C)/det(A,B)

     Ex1: 2x - (-4)y = 20   det(A,B)=-16+28  =  12
          7x - (-8)y = 52   det(B,C)=-208+160= -48
	                    det(A,C)=104-140 = -36
	  x0=-det(B,C)/det(A,B)=48/12=4
	  y0=-det(A,C)/det(A,B)=36/12=3

     Ex2: 7x - (-5)y = 11   det(A,B)=-21+5 = -16
          1x - (-3)y =  5   det(B,C)=-25+33=   8
	                    det(A,C)=35-11 =  24
	  x0=-det(B,C)/det(A,B)=-8/-16=1/2
	  y0=-det(A,C)/det(A,B)=-24/-16=3/2

  4. det(X,-Y) = -det(X,Y)
     det(Y,X) = -det(X,Y)
*/


/* comme QSORT(), mais pour qsort2() */
#define QSORT2(T1,T2,n,f) qsort2(T1,sizeof(*(T1)),T2,sizeof(*(T2)),n,f)


void qsort2(void *T1,const int w1,void *T2,const int w2,
	    const int n,int (*fcmp)(const void*,const void*)){
/*
  Trie simultanément deux tableaux T1 et T2 non NULL chacun de n
  éléments selon la fonction de comparaison fcmp(). Ici w1 (resp. w2)
  représentent la taille d'un élément de T1 (resp. T2). Pour cela on
  construit un nouveau tableau T où chaque élément T[i] est composée
  de la paire d'éléments <T1[i],T2[i]> où T1[i] est stocké juste avant
  T2[i]. Puis on applique qsort(T,n,w1+w2,fcmp), et enfin on remet
  dans T1 et T2 les éléments ainsi triés. La fonction de comparaison
  fcmp() doit pouvoir s'appliquer à une paire <T1[i],T2[i]>, même si
  elle peut très bien s'appliquer seulement sur T1.

  Mais le plus souvent fcmp() peut s'appliquer seulement à T1, comme
  dans l'exemple suivant:

  int A[]={ 5, 8, 2, 9, 1, 6, 3, 7, 4};
  int B[]={-5,-8,-2,-9,-1,-6,-3,-7,-4};

  QSORT2(A,B,9,fcmp_int);
  -> A[]={ 1, 2, 3, 4, 5, 6, 7, 8, 9}
  -> B[]={-1,-2,-3,-4,-5,-6,-7,-8,-9}

  QSORT2(B,A,9,fcmp_int);
  -> A[]={ 9, 8, 7, 6, 5, 4, 3, 2, 1}
  -> B[]={-9,-8,-7,-6,-5,-4,-3,-2,-1}

*/
  if(n<2) return; /* ne rien faire si pas au moins 2 éléments */

  const int w=w1+w2; /* taille en octets d'une paire <T1[i],T2[i]> */
  void *T=malloc(n*w); /* taille en octets de T */
  if(T==NULL) Erreur(3); /* si problème mémoire */
  void *t,*t1,*t2;
  int i;

  /* fusionne T1 et T2 dans T */
  t=T; t1=T1; t2=T2;
  for(i=0;i<n;i++,t+=w){
    memcpy(t,   t1,w1); t1 += w1;
    memcpy(t+w1,t2,w2); t2 += w2;
  }

  /* trie T */
  qsort(T,n,w,fcmp);
  
  /* sépare T en T1 et T2 */
  t=T; t1=T1; t2=T2;
  for(i=0;i<n;i++,t+=w){
    memcpy(t1,t,   w1); t1 += w1;
    memcpy(t2,t+w1,w2); t2 += w2;
  }

  free(T);
  return;
}


void InitXY(query* const Q){
/*
  Initialise les tableaux XPOS et YPOS suivants les options
  éventuelles de -xy. Les variables suivantes, en plus de XPOS et
  YPOS, peuvent être mise à jour: Q->n, XMIN,XMAX,YMIN,YMAX
  (BoundingBox), VSIZESTD, VSIZEXY.
*/

  int i,k;
  double sx,sy,tx,ty;

  for(;;){ /* pour pouvoir faire un break; */

    if(XYtype==XY_USER) break; /* coordonnées définies par l'utilisateur */

    if(XYtype==XY_FILE){ /* charge à partir d'un fichier et met à jour Q->n */
      Q->n=LoadXY(FILEXY);
      break;
    }

    if(XYratio<=0) Erreur(6); /* paramètre incorrect */

    if(Q->n<0) Q->n=0;
    ALLOC(XPOS,Q->n);
    ALLOC(YPOS,Q->n);

    if(XYtype==XY_UNIF){ /* uniforme dans [0,1[x[0,XYratio[ */
      for(i=0;i<Q->n;i++){
	XPOS[i]=RAND01;
	YPOS[i]=XYratio*RAND01;
      }
      break;
    }
    
    if(XYtype==XY_PLAW){ /* loi puissance autour des graines choisies dans [0,1[ */
      if(XYseedk<=0) Erreur(6); // paramètre incorrect
      ALLOC(XSEED,XYseedk);
      ALLOC(YSEED,XYseedk);
      sx=sy=0; /* calcule (sx,sy), le barycentre des XYseedk graines */
      for(i=0;i<XYseedk;i++){
	XSEED[i]=RAND01;sx+=XSEED[i];
	YSEED[i]=RAND01*XYratio;sy+=YSEED[i];
      }
      sx /= XYseedk; sy /= XYseedk;
      sx -= 0.5; sy -= XYratio/2;
      for(i=0;i<XYseedk;i++){ /* centre par rapport au barycentre */
	XSEED[i] -= sx; /* enlève le barycentre puis décale de 0.5 */
	YSEED[i] -= sy;
      }
      /* on génère les points autour des graines */
      tx=sqrt(log(XYseedk+1)/XYseedk); /* rayon r, le +1 est important car abérant pour XYseedk=1 */
      for(i=0;i<Q->n;i++){
	k=random()%XYseedk;    /* on choisit la graine numéro k au hasard */
	sx=M_2PI*RAND01;  /* angle aléatoire */
	sy=tx*pow(RAND01,XYpower); /* longueur aléatoire */
	XPOS[i]=XSEED[k]+sy*cos(sx);
	YPOS[i]=YSEED[k]+sy*sin(sx)*XYratio;
      }
      break;
    }

    if(XYtype==XY_PERM){ /* permutation de [0,Q->n[ */
      NALLOC(int,P,Q->n);
      for(i=0;i<Q->n;i++) XPOS[i]=(double)(P[i]=i); // initialise aussi P[i]
      Permute(P,Q->n); // modifie P[i]
      for(i=0;i<Q->n;i++) YPOS[i]=(double)P[i];
      free(P);
      break;
    }

    if(XYtype==XY_MESH){ /* grille de paramètre Xmesh x Ymesh */
      if((Xmesh<=0)||(Ymesh<=0)) Erreur(6);
      for(i=0;i<Q->n;i++) XPOS[i]=(double)(i%Xmesh),YPOS[i]=(double)(i/Xmesh);
      break;
    }

    if(XYtype==XY_CYCLE){ /* cycle de rayon 1 et de centre (0,0) */
      const double a=M_2PI/Q->n;
      sx=0;
      for(i=0;i<Q->n;i++){
	XPOS[i]=cos(sx);
	YPOS[i]=sin(sx)*XYratio;
	sx += a;
      }
      break;
    }

    /*
      Attention ! pour générer des points aléatoires uniformes sur un
      disque unité, il faut faire: a=M_2PI*RAND01, r=sqrt(RAND01),
      puis (x,y)=(r*cos(a),r*sin(a)). Si on utilise seulement r=RAND01
      (sans le sqrt) alors les points se retrouvent plus concentrés au
      centre du disque: http://mathworld.wolfram.com/DiskPointPicking.html
    */
    if((XYtype==XY_DISK)||(XYtype==XY_CIRCLE)||(XYtype==XY_HYPER)){
      /* star-shaped polygon (ou cercle) de centre (0,0) et rayon <= 1
	 ou disque hyperbolique */
      NALLOCZ(double,A,Q->n,M_2PI*RAND01); // Q->n angles aléatoires de [0,2π[
      QSORT(A,Q->n,fcmp_double); // trie les angles (ne tient pas compte de ROUND ...)
      for(i=0;i<Q->n;i++){ // transforme coordonnées polaires en cartésiennes
	switch(XYtype){
	case XY_CIRCLE: sx=1; break;
	case XY_DISK: sx=sqrt(RAND01); break;
	case XY_HYPER: sx=exp(-RAND01*XYpower); break;
	default: sx=RAND01;
	}
	XPOS[i]=sx*cos(A[i]);
	YPOS[i]=sx*sin(A[i])*XYratio;
      }
      free(A);
      break;
    }

    if(XYtype==XY_RPOLY){ /* polygone convexe régulier */
      if(XYpoly<3) Erreur(42);
      /*
	L'Algorithme est le suivant, algorithme qui s'applique
	indépendemment à chacun des Q->n points. On considère le triangle
	défini par un seul coté du polygone, d'angle ϴ=2π/p où
	p=XYpoly. On l'oriente pour que l'axe des abscisses
	corresponde à la médiane de l'angle ϴ. Le coté du polygone est
	ainsi un segment vertical d'abscisse cos(ϴ/2) et de hauteur
	2|sin(ϴ/2)|. On peut alors tirer un point M aléatoirement
	uniforme dans ce triangle, puis on tourne le point M d'un
	angle i*ϴ avec i aléatoire dans [0,p[.
	
	Pour tirer aléatoirement un point M dans un triangle (O,V1,V2)
	dont un coin (ici O) est l'origine on peut faire comme suggéré
	dans http://mathworld.wolfram.com/TrianglePointPicking.html.
	On choisit r1,r2 uniformes dans [0,1], puis on construit le
	point M=r1*V1+r2*V2. Ce point est aléatoire dans le
	parallélogramme définit par les 4 points (O,V1,V2,V1+V2).  Si
	M n'est pas dans le triangle, alors soit on recommence, soit
	on prend le symétrique tombant dans le triangle. Le nombre
	d'essais moyen est 2.

	Ici V1=(cos(ϴ/2),sin(ϴ/2)) et V2=(cos(ϴ/2),-sin(ϴ/2)). M est
	dans le triangle si x(M)<=X(V1)=X(V2)=cos(ϴ/2). Sinon, on
	change M en son symétrique par rapport au point (cos(ϴ/2),0).
      */
      const double t=M_2PI/XYpoly; // t=ϴ=angle défini par un coté du polygone
      const double c=cos(t/2); // c=abscisse du coté vertical
      const double s=sin(t/2); // s=demi-hauteur du coté vertical
      double a,r;

      for(i=0;i<Q->n;i++){ // indépendemment pour chaque point
	a=RAND01,r=RAND01; // tire un point M uniformément dans le triangle
	sx=c*(a+r),sy=s*(a-r); // M=(sx,sy) dans le parallélogramme
	if(sx>c) sx=2*c-sx,sy=-sy; // prend le symétrique de M par rapport à (c,0)
	r=hypot(sx,sy); // (a,r)=coordonnées polaires de la rotation de M
	a=angle(sx,sy)+t*(random()%XYpoly); // a=angle avec rotation aléatoire
	XPOS[i]=r*cos(a);
	YPOS[i]=r*sin(a)*XYratio;
      }
      break;
    }

    if(XYtype==XY_CONVEX){ /* points en position convexe, algo en Q->n^2 */
      double mx,my,d,x,y,xa,ya,xb,yb,a1,b1,a2,b2,c2,h;
      int t,j,b=random()&1;

      // on suppose que les points P(0) ... P(i-1) sont déjà en
      // position convexes, l'intérieur contenant l'origine (0,0), on
      // souhaite placer un nouveau point Q=P(i)

      for(i=0;i<Q->n;i++){
	
	// calcule tx=angle aléatoire du nouveau point Q
	tx=M_2PI*RAND01;

	// pour i=0,1,2, tx est réduit à l'un des trois secteurs non
	// adjacents d'angle π/3 (les 3 cônes positifs ou négatifs des
	// 6 secteurs d'angle π/3) de sorte que l'origine sera
	// forcément à l'intérieur du triangle (b=bit aléatoire). On
	// pourrait penser à recentrer les points selon leur
	// barycentre, mais cela produit des erreurs car les points
	// peuvent maintenant être en dehors du cercle de rayon 1.
	if(i<3) tx=(2*i+b)*M_2PI/6 + fmod(tx,M_2PI/6);
	
	DEBUG(printf("tx=%.02lf %03.0lf\n",tx,360*tx/M_2PI););
	
	// cherche k dans [0,i] tq P(k-1) < Q < P(k). Si Q est le plus
	// grand des points, alors k=i. Si c'est le plus petit, k=0.
	// On fait une recherche linéaire, bien qu'on pourrait la
	// faire en log(i)
	
	for(k=0;k<i;k++){ // P(k)=(XPOS[k],YPOS[k])
	  ty=angle(XPOS[k],YPOS[k]);
	  if(ty==tx) tx=nextafter(tx,7); // évite les angles égaux. NB: 7 > 2π > tx
	  if(tx<ty) break;
	}
	DEBUG(PRINT(k););
	
	// S=[sx,sy[ segment d'angle tx où l'on cherche Q
	// D1=droite contenant Q passant par (0,0) contenant S
	// D1: a1*x - b1*y = 0

	a1=sin(tx);
	b1=cos(tx);
	
	// M=(mx,my)=point maximum pour Q selon la direction D1. De
	// manière générale, M est sur l'ellipse de demi-hauteur
	// XYratio, la demi-longueur valant 1. L'équation paramétrique
	// donne le rayon r(t) en fonction de l'angle t qui est: r(t)
	// = b/sqrt(1-e*cos(t)^2) où a=demi-longueur=1, b=demi-hauteur
	// et e=1-(b/a)^2. Si XYratio=1, alors r(t)=1.

	d=XYratio/sqrt(1-(1-XYratio*XYratio)*b1*b1); // d=r(tx)
	DEBUG(PRINTD(d););
	mx=d*b1;
	my=d*a1;

	// S = [sx,sy[
	sx=nextafter(0,1); // borne inf de S
	sy=hypot(mx,my);   // borne sup de S

	// on réfuit S en fonction des trois droites définies par les
	// 4 points successifs P(k-2), P(k-1), P(k), P(k+1)
	if(i>2){ // rien à faire si i=0,1 ou 2

	  /* Il est possible d'avoir P(k-2)=P(k+1) lorsque i=3. Pour
	     les trois droites successives [P(k-2),P(k-1)],
	     [P(k-1),P(k)] et enfin [P(k),P(k+1)], on va calculer leur
	     droite D2, puis l'intersection entre D1 et D2. On calcule
	     ensuite une nouvelle borne pour sx,sy. */

	  j=(k-2+i)%i; // au départ A=P(k-2) et B=P(k-1).
	  for(t=-1;t<=1;t++){ // on répète trois fois: t=-1,0,+1
	    xa=XPOS[j],ya=YPOS[j]; // point A=(xa,ya)
	    j=(j+1)%i; // point suivant
	    xb=XPOS[j],yb=YPOS[j]; // point B=(xb,yb)
	    // droite D2 passant par A et B
	    //  D2: a2*x - b2*y = c2
	    a2=yb-ya;
	    b2=xb-xa;
	    c2=det(xa,ya,xb,yb);
	    // (x,y)=intersection de D1 et D2
	    //  D1: a1*x - b1*y = 0
	    //  D2: a2*x - b2*y = c2
	    d=det(a1,a2,b1,b2);
	    if(d==0) continue; // droite suivante si pas d'intersection
	    x=-det(b1,b2,0,c2)/d;
	    y=-det(a1,a2,0,c2)/d;
	    // D1 est en fait une demi-droite. Est-ce que (x,y)
	    // appartient à cette demi-droite ? c'est-à-dire est-ce
	    // que les vecteurs (b1,a1) et (x,y) sont dans le même
	    // sens ou opposé ?  même sens <=> (b1*x>0)&&(a1*y>0)
	    if((b1*x>0)&&(a1*y>0)){ // si même sens, on modifie sx ou sy
	      DEBUG(printf("sx=%lf sy=%lf -> ",sx,sy););
	      h=hypot(x,y);
	      if(t) sy=fmin(sy,h); // t=-1 ou +1
	      else sx=fmax(sx,h);  // t=0
	      DEBUG(printf("sx=%lf sy=%lf h=%lf t=%i\n",sx,sy,h,t););
	    }
	  }
	  // ici sx<sy sinon les points ne sont pas en position convexe
	  DEBUG(if(sx>sy) printf("problème: sx>sy\n"););
	}
	
	// insère le point Q entre P(k-1) et P(k) en décalant
	// [P(k)...P(i-1)] vers [P(k+1)...P(i)]
	for(t=i;t>k;t--) XPOS[t]=XPOS[t-1],YPOS[t]=YPOS[t-1];

	// calcule, en fonction de tx et S=[sx,sy[, le nouveau point Q
	// qui devient P(k). On met sqrt() pour être plus proche d'une
	// distribution uniforme sur un disque (donc plus souvent
	// proche du bord du cercle que du centre)
	sx += (sy-sx)*sqrt(RAND01); // sx=point aléatoire dans S=[sx,sy[
	if(sx==0) sx=nextafter(0,1); // on force sx!=0 pour avoir P(k)≠(0,0)
	XPOS[k]=sx*cos(tx);
	YPOS[k]=sx*sin(tx);

	DEBUG(
	      PRINT(i);
	      for(t=0;t<=i;t++)
		printf("point %i: %+.02lf %+.02lf \t%+.02lf %03.lf\n",
		       t,XPOS[t],YPOS[t],angle(XPOS[t],YPOS[t]),
		       360*angle(XPOS[t],YPOS[t])/M_2PI);printf("\n");
	      );
      }
      break;
    }

    if(XYtype==XY_CONVEX2){ /* points en position convexe v2, algo en N*log(N) */
      
      // on part de points aléatoires, puis on calcule la différence
      // entre deux points consécutifs. La différence est nulle. On
      // trie ces points selon l'angle, puis on dessine de proche en
      // proche les points de l'enveloppe convexe (avec chaque fois un
      // angle croissant donc).

      for(i=0;i<Q->n;i++) XPOS[i]=RAND01,YPOS[i]=RAND01; /* points aléatoires */
      if(Q->n>0){ /* NB: il faut Q->n>0 */
	const double x0=XPOS[0],y0=YPOS[0]; // sauvegarde le 1er point
	for(i=0;i<Q->n-1;i++) XPOS[i]-=XPOS[i+1],YPOS[i]-=YPOS[i+1]; // différences
	XPOS[i] -= x0,YPOS[i] -= y0;
	QSORT2(XPOS,YPOS,Q->n,fcmp_angle); // trie les angles
	for(i=1;i<Q->n;i++) XPOS[i] += XPOS[i-1],YPOS[i] += YPOS[i-1]; // dessin
	}
      break;
    }

  }

  /* ici on devrait normalement avoir XPOS,YPOS <> NULL */
  if((XPOS==NULL)||(YPOS==NULL)) Erreur(8);

  if(XYnoiser>0) /* "noise" doit être avant "box" */
    for(i=0;i<Q->n;i++){
      sx=M_2PI*RAND01; /* angle aléatoire */
      sy=XYnoiser*pow(RAND01,XYnoisep); /* longueur aléatoire */
      XPOS[i] += sy*cos(sx); /* décale XPOS */
      YPOS[i] += sy*sin(sx); /* décale YPOS */
    }
  
  if((BOXX>0)&&(BOXY>0)){ /* "box" doit être après "noise" */
    BoundingBox(Q); /* calcule les BB */
    if(Q->n<2) sx=sy=0;
    else{
      sx=BOXX/(XMAX-XMIN);
      sy=BOXY/(YMAX-YMIN);
    }
    for(i=0;i<Q->n;i++){
      XPOS[i]=sx*(XPOS[i]-XMIN);
      YPOS[i]=sy*(YPOS[i]-YMIN);
    } // ici les BB sont obsolètes
  }
  
  if(ROUND<DBL_DIG){ /* arrondit éventuellement les coordonnées */
    sx=pow(10,ROUND);
    for(i=0;i<Q->n;i++){
      XPOS[i]=rint(XPOS[i]*sx)/sx;
      YPOS[i]=rint(YPOS[i]*sx)/sx;
    }
  }

  if(XYunique){ /* élimine les doubles, en triant les points */
    NALLOC(point,P,Q->n);
    for(i=0;i<Q->n;i++) P[i].x=XPOS[i],P[i].y=YPOS[i];
    QSORT(P,Q->n,fcmp_point); /* tri les points */
    point p=P[0]; p.x -= 1.0; /* ici le point p <> du 1er élément */
    for(i=k=0;i<Q->n;i++) /* k=nombre de points uniques */
      if(fcmp_point(P+i,&p)){ /* copie que si différent de l'élément p */
	p=P[i];
	XPOS[k]=p.x;
	YPOS[k++]=p.y;
      }
    free(P);
    if(k<Q->n){
      Q->n=k;
      REALLOC(XPOS,Q->n);
      REALLOC(YPOS,Q->n);
    }
  }

  /* on (re)calcule les BB */
  BoundingBox(Q);

  /* mise à jour de la taille des sommets */
  VSIZESTD *= XYvsize;
  VSIZEXY  *= XYvsize;

  return;
}


color *GradColor(const color *T,const int n,const int m){
/*
  Retourne un tableau de m couleurs formant un dégradé obtenu à partir
  d'un tableau T de n couleurs. Pour avoir un dégradé simple d'une
  couleur T[0] à T[1] il faut initialiser T[0],T[1] et poser n=2. Pour
  avoir un dégradé cyclique, il suffit de répéter la couleur T[0] en
  dernière position de T (et ajouter 1 à n, donc d'avoir
  T[n-1]=T[0]). Il faut dans tous les cas n>1 et m>0. On peut avoir
  m<n. Dans ce cas on prend la première couleur de T, puis la i-ème
  couleur est (i*(n-1))/(m-1).
*/
  color c1,c2;
  int i,j,k,r,q,n1,m1,dr,dg,db;

  if(T==NULL) return NULL; /* normalement ne sert à rien */

  NALLOC(color,P,m);
  c2=P[0]=T[0];
  if(m==1) return P;
  /* maintenant m >= 2 */

  m1=m-1; n1=n-1; /* valeurs utilisées souvent */

  if(m<=n){ /* cas où il y a moins de couleurs demandées que dans T */
    for(i=1;i<m;i++) /* m-1 fois */
      P[i]=T[(i*n1+m1-1)/m1]; /* le "+m-2" est pour arrondir à l'entier sup. */
    return P;
  }

  /*
    Cas m>n.  Soient B_1,B_2,...B_(n-1) les n-1>0 blocs de couleurs,
    B_i commençant juste après la couleurs T[i-1] et se terminant avec
    la couleur T[i]. On doit répartir m-1 couleurs dans ces n-1 blocs,
    la couleurs T[0] étant déjà dans P. On met alors
    floor((m-1)/(n-1)) couleurs par blocs, et une de plus pour B_i si
    i<=(m-1)%(n-1).
   */
  r=m1%n1; /* il reste r couleurs */
  q=(m1/n1)+(r>0); /* nombre de couleurs par blocs (+1 au départ si r>0) */
  for(i=j=k=1;i<n;i++){ /* on traite le bloc B_i, P[k]=prochaine couleur libre */
    c1=c2;   /* c1=T[i-1] */
    c2=T[i]; /* c2=T[i] */
    dr=c2.r-c1.r;
    dg=c2.g-c1.g;
    db=c2.b-c1.b;
    for(j=1;j<=q;j++){ /* extrapolation linéaire de q points dans ]c1,c2] */
      P[k].r=c1.r+(j*dr)/q;
      P[k].g=c1.g+(j*dg)/q;
      P[k].b=c1.b+(j*db)/q;
      k++;
    }
    if(i==r) q--; /* une couleur de moins pour B_{r+1}...B_{n-1} */
  }
  return P;
}


int graphical(const int *S,int k){
/*
  Vérifie si la suite S=(n_1,d_1,n_2,d_2,...,n_k,d_k) constituée de k
  couples est graphique ou pas, c'est-à-dire s'il existe au moins un
  graphe simple ayant exactement n_i sommets de degré d_i. On renvoie
  une valeur <0 si la séquence n'est pas graphique, et sinon on
  renvoie n=∑_i n_i, c'est-à-dire le nombre de sommets du graphe.

  L'algorithme est basé sur le test d'Erdős and Gallai (1960) qui ont
  prouvé qu'une suite de degrés (t_1,...,t_n) triée dans l'ordre
  décroissant est graphique ssi la somme des degrés des sommets est
  paire et la suite vérifie la propriété suivante (ici t_i est le
  degré du sommet i=1..n):

  (1) ∑_{i=1}^r t_i <= r*(r-1) + ∑_{i=r+1}^n min{r,t_i}

  pour chaque entier r=1..n-1 (Skiena 1990, p. 157). Cette propriété
  se généralise aux graphes orientés.  Tripathi et Vijay (2003) on
  montré qu'il suffit de vérifier (1) pour les valeurs de r où les t_i
  changent (voir aussi [DF05]).

  Si on note r_j=n[0]+...+n[j], il faut donc vérifier (1) seulement
  pour les k valeurs de r suivantes: r_0, r_1,..., r_j, ...,
  r_{k-1}. A l'étape j=0...k-1, l'équation (1) se réécrit donc en:
  
  (2) ∑_{i=0}^j n[j]*d[j] <= r_j*(r_j-1) + ∑_{i=j+1}^{k-1} min{r_j,n[i]*d[i]}

  On remarque aussi que l'équation (1) ou (2) est trivialement vraie
  dès que r*(r-1) ou r_j*(r_j-1) >= m car le terme de gauche est
  toujours <= m où m=∑_i n_i*d_i. C'est important pour éviter de faire
  un produit r*(r-1) qui peut être très grand par rapport à m et
  dépasser la capacité des entiers (et produire un faux négatif). Par
  exemple pour k=1, n_0=10^6 et d_0=3, on a m=3,000,000 et r*(r-1) qui
  peut atteindre ≃ 10^12 ... En pratique dès que r>sqrt(m), alors
  l'équation (1) sera vraie.

*/
  if((S==NULL)||(k<=0)) return -1;

  int i,j,m,v;
  for(i=0;i<2*k;i++)
    if(S[i]<0) return -1; /* une valeur <0 => pas graphique */

  /* range les valeurs de S dans n[i] et d[i] */
  NALLOCZ(int,n,k,S[2*_i]);   /* n[0] ... n[k-1]: les n_i */
  NALLOCZ(int,d,k,S[2*_i+1]); /* d[0] ... d[k-1]: les d_i */
  QSORT2(d,n,k,fcmp_int_inv); /* trie les valeurs par ordre décroissant */

  /* fusionne les valeurs de d_i identiques */
  for(j=0,i=1;i<k;i++)
    if(d[i]==d[j]) n[j]+=n[i];
    else j++;
  k=i+1; /* réajuste k */

  m=0; /* m=somme des degrés des sommets=∑_i n_i*d_i */
  v=0; /* v=nombre de sommets=∑_i n_i */

  for(i=0;i<k;i++){
    m += n[i]*d[i];
    v += n[i];
  }

  if(m&1) v=-1; /* somme de degré impaire => pas graphique */
  else{ /* somme des degrés paire */
    /* on doit vérifier, pour chaque j, que s1 <= r_j*(r_j-1) + s2 */
    int rj=0; /* rj=∑_{i=0}^j n_i */
    int s1=0; /* s1=∑_{i=0}^{j-1} n_i*d_i */
    int s2;   /* s2=∑_{i=j+1}^{k-1} min{r_j,n_i*d_i} */
    const int c=1+(int)sqrt(m); /* NB: r>=c => r(r-1)>=m */
    
    for(j=0;j<k;j++){ /* étape j=0..k-1 */
      rj += n[j]; if(rj>=c) break; /* => rj(rj-1)>=m => équation vraie */
      s1 += n[j]*d[j];
      for(i=j+1,s2=0;i<k;i++) s2 += imin(rj,n[i]*d[i]);
      if(s1>rj*(rj-1)+s2){ v=-1; break; } /* équation fausse */
    }
  }
  
  free(n);
  free(d);
  return v;
}


/***********************************

         ROUTINES POUR LES
          FONCTIONS adj()

***********************************/


static inline double Norme_dxy(const double dx,const double dy){
/*
  Calcule la distance du vecteur de coordonnées (dx,dy) selon la norme
  définie par NORM. Renvoie NORM_FAIL si aucune norme n'a été trouvée
  (ce qui en principe ne devrait jamais arriver).

  Principe pour le calcul de NORM_POLY: On considère un polygone
  convexe régulier centré sur l'origine de cercle inscrit de rayon
  unité et orienté de sorte que son coté le plus à droit soit
  vertical. Ce coté est numéroté 0, les autres successivement en
  tournant vers la gauche 1,2,...,(NORM_poly)-1. Notez que
  NORM_poly>=3. Tout d'abord on cherche le numéro i du coté du
  polygone qui est coupé par la demi-droite D partant de l'origine et
  passant par le point (x,y). Soit a l'angle de la droite D. L'angle
  inférieur du coté i juste en dessous de D est i*p0-p0/2, où p0 est
  l'angle entre deux coins consécutifs du polygone, c'est-à-dire
  i*p0-p0/2 <= a < i*p0+p0/2. On tourne (virtuellement) la figure à
  droite de i cotés, soit d'un angle de i*p0, de sorte que coté i soit
  vertical.  Il intersecte alors la droite des abscisses en son milieu
  à l'abscisse 1 précisément. Le point d'intersection entre D et ce
  coté est à distance h=1/cos(b) où b=a-i*p0 est l'angle entre D et le
  milieu du coté i puisque h*cos(b)=1. La distance recherchée est
  simplement hypot(x,y)/h = hypot(x,y)*cos(b). Notons qu'il n'est pas
  possible d'avoir cos(b)=0 puisque D ne peut être parallèle au i-ème
  coté.

  La norme hyperbolique correspond à la distance hyperbolique par
  rapport à l'origine (0,0). Pour le calcul hyperbolique entre deux
  point quelconques, il faut passer par la donction dist_ij() qui
  dépend alors des points et pas seulement de la différence des
  points.
*/
  switch(NORM){
  case NORM_L1:   return (dx+dy);       /* norme L1 */
  case NORM_L2:   return hypot(dx,dy);  /* norme L2 */
  case NORM_LMAX: return fmax(dx,dy);   /* norme Lmax */
  case NORM_LMIN: return fmin(dx,dy);   /* norme Lmin */
  case NORM_POLY:;                      /* norme polygonale */
    const double p0=M_2PI/NORM_poly;    /* NB: ici NORM_poly>=3 */
    const double a=angle(dx,dy); 
    return hypot(dx,dy)*cos(a-p0*floor(a/p0+0.5));
  case NORM_HYPER: return 2*atanh(hypot(dx,dy)); /* norme hyperbolique */
  default: return NORM_FAIL;            /* norme indéterminée */
  }
}


static inline double dist_ij(const int i,const int j){
/*
  Calcule la distance entre les points (XPOS[i],YPOS[i]) et
  (XPOS[j],YPOS[j]) selon la norme définie par NORM (voir Norme_dxy).
  Il faut faire attention que la norme n'est pas toujours symétrique
  en i et j, comme pour NORM_POLY avec p impair.

  Pour la distance hyperbolique, on ne fait pas appel à la norme qui
  n'est définie que pour un vecteur depuis l'origine, pas pour une
  différence de points. Aussi le points doivent être situés à
  l'intérieur du disque unité centré sur l'origine, sinon soit il y
  une division par 0 (si un des points est situé sur le cercle) ou la
  valeur renvoyée peut-être négative.
*/
  if(NORM==NORM_HYPER){
    const double u=XPOS[i]*XPOS[i]+YPOS[i]*YPOS[i];
    const double v=YPOS[j]*XPOS[j]+YPOS[j]*YPOS[j];
    const double dx=XPOS[i]-XPOS[j];
    const double dy=YPOS[i]-YPOS[j];
    return acosh(1+2*(dx*dx+dy*dy)/(1-u*u)/(1-v*v));
  }
  return Norme_dxy(fabs(XPOS[i]-XPOS[j]),fabs(YPOS[i]-YPOS[j]));
}


double distgone(const int u,const int v,
		const int i,const int p,const int k,const double w){
/*
  Calcule la distance P_i(u,v) qui n'est pas symétrique en u et v. Il
  s'agit de la "distance p-gone" (un polygone régulier à p cotés)
  relative à la direction i (axe d'angle i*2π/k) entre les points
  d'indice u et v, restreint au cône de visibilité d'angle w*(p-2)*π/p
  (d'angle w*π si p est infini, c'est-à-dire si p<3), avec w=0...1
  (voir aussi la définition du thetagone dans l'aide en ligne). Ici,
  k>0 est le nombre de directions. La fonction renvoie DBL_MAX si la
  distance est infinie ce qui est possible avec l'étape 1 de
  l'algorithme. L'algorithme est en O(1), notamment indépendant de p
  et k.

  Soient a_j (j=0...p-1) les sommets du p-gone P_i de rayon unité avec
  a_0=u et numérotés consécutivement en tournant dans le sense
  positif. Soit c le centre de P_i. Donc dist(u,c)=1, les a_j étant
  sur un cercle de rayon unité. On remarque que l'angle (u,a_j) est
  indépendant du rayon du p-gone, il ne dépant que de j. En fait,
  l'angle (a_j,u,a_{j+1}) vaut la moitié de l'angle (a_j,c,a_{j+1}),
  soit π/p. L'angle entre deux cotés consécutifs d'un p-gone vaut
  (p-2)*π/p.

  L'algorithme de calcul pour distgone(u,v) est le suivant:

  1. Trouver la direction j tq (u,v) soit dans le cône de visibilité
     et dans la région [(u,a_j),(u,a_{j+1})[.  Si j n'existe pas,
     alors on renvoit une distance infinie. Si p est infini, a_j est
     simplement sur la droite (u,v).

  2. On calcule l'intersection v' entre les droites (u,v) et
     (a_j,a_{j+1}). Si p est infini, v'=a_j. L'intersection existe
     forcément. Eventuellement v est sur la droite (u,a_j).

  3. distgone(u,v)=dist(u,v)/dist(u,v').

*/
  int j;
  double xu,xv,dxc,dxa,dxb,dxv;
  double yu,yv,dyc,dya,dyb,dyv;
  double hv,A,Ac,Aw;

  xu=XPOS[u],yu=YPOS[u]; /* coordonnées de u */
  xv=XPOS[v],yv=YPOS[v]; /* coordonnées de v */
  Ac=(double)i*M_2PI/k; /* angle (u,c), c=centre de P_i */
  dxc=cos(Ac),dyc=sin(Ac); /* coordonnées du centre c dans le repère u */

  dxv=xv-xu;dyv=yv-yu; /* coordonnées de v dans repère u */
  hv=hypot(dxv,dyv); /* |v-u|=dist(u,v) */
  if(hv==0) return 0; /* si u,v ont les mêmes coordonnées */

  /*
    Rappel: Si a et b sont deux vecteurs, alors le produit scalaire
    (dot product) est le réel a.b = xa*xb + ya*yb = |a|*|b|*cos(a,b),
    où |a|=hypot(xa,ya)=sqrt(xa^2 + ya^2). Notons aussi que
    |a|*|b|*sin(a,b) = det(xa,ya,xb,yb). Donc le signe de cos(a,b) est
    celui de xa*xb + ya*yb. Pour calculer sin(a,b) il faut faire une
    rotation de +π/2 au vecteur a et calculer cos(a',b) où
    a'=(-ya,xa). Donc sin(a,b) = cos(a',b) = (-ya*xb + xa*yb) /
    (|a|*|b|). Et le signe de sin(a,b) est celui de xa*yb - ya*xb =
    det(xa,ya,xb,yb).
  */

  /* Aw=demi-angle du cône de visibilité */
  Aw=w*M_PI_2; /* si p infini */
  if(p>2) Aw *= (double)(p-2)/p; /* si p est fini */

  /*
    Il faut bien sûr que (u,v) soit dans le cône de visibilité. La
    bissectrice de ce cône est l'axe (u,c) et son angle est
    w*(p-2)π/p (w*π si p infini). On note (w1,u,w2) le cône en
    question. Il faut que (u,v) soit entre (u,w1) (compris) et (u,w2)
    (non compris). Donc si sin(w1,u,v) < 0 ou si sin(w2,u,v) > 0 alors
    il n'existe pas de j (et donc on retourne une distance infinie).
  */

  A=Ac-Aw; /* A=angle (c,w1) */
  dxa=cos(A);dya=sin(A); /* coordonnées de w1 relatif à u */
  if(det(dxa,dya,dxv,dyv)<0) return DBL_MAX; /* v avant w1 */

  A=Ac+Aw; /* A=angle (c,w2) */
  dxa=cos(A);dya=sin(A); /* coordonnées de w2 relatif à u */
  if(det(dxa,dya,dxv,dyv)>=0) return DBL_MAX; /* v après ou sur w2 */

  /*
    Ici v est dans le cône de visibilité et donc la droite (uv)
    intersecte P_i en un point v'.
  */

  Ac -= M_PI; /* Ac=Ac-π */

  /* Cas p infini */
  if(p<3){
    /*
      On raisone dans le repère u.  On pose c'=(dyc,-dxc),
      c'est-à-dire une rotation de -π/2 de (u,v). On a |uc'|=1. On
      calcule l'angle A=(uv,uc'), en fait cos(A). On obtient v' en
      tournant autour de c d'un angle Ac-π+2A ...
     */
    A=acos(det(dxv,dyv,dxc,dyc)/hv); // A=acos((dxv*dyc-dyv*dxc)/hv);
    A=Ac+2*A;
    dxa=dxc+cos(A);dya=dyc+sin(A);
    return hv/hypot(dxa,dya);
  }

  /*
    Cas p fini.  On cherche j de sorte qu'en tournant dans le sens
    positif, le vecteur (u,v) soit compris entre (u,a_j) (compris) et
    (u,a_{j+1}) (non compris). La droite (u,v) intersecte le segment
    [a_j, a_{j+1}[. L'indice j recherché est l'entier tq: (j-1)*π/p <=
    angle(a_1,u,a_j) < j*π/p. Et donc, j=1+floor{acos(a_1,u,v)/(π/p)}.
  */

  const double Ap=M_2PI/p; /* valeur souvent utilisée */
  A=Ac+Ap; /* angle (c,a_1) */
  dxa=dxc+cos(A);dya=dyc+sin(A); /* coordonnées de a_1 relatif à u */

  /* Aw=cos(a_1,u,v) = (a_1-u).(v-u) / dist(a_1,u)*dist(u,v) */
  Aw=det(dxa,-dya,dyv,dxv)/(hypot(dxa,dya)*hv);
  j=(int)((acos(Aw)*(double)p)/M_PI); /* en fait, la variable j vaut "j-1" */
  A += Ap*j; /* angle (c,a_j): on part de a_1, donc on décale de j-1 cônes */
  dxa=dxc+cos(A);dya=dyc+sin(A); /* coordonnées de a_j relatif à u */
  A += Ap; /* angle (c,a_{j+1}) */
  dxb=dxc+cos(A)-dxa;dyb=dyc+sin(A)-dya; /* vecteur (a_j,a_{j+1}) */

  /*
    Calcule l'unique intersection entre la droite Dv définie par le
    vecteur (u,v) et la droite Dj définie par le vecteur
    (a_j,a_{j+1}). Dans le repère u, les équations sont:

    Dv: dyv*X - dxv*Y = 0
    Dj: dyb*X - dxb*Y = B avec B=dxa*dyb-dxb*dya

    car a_j=(dxa,dya) appartient à Dj.
    L'intersection (x0,y0) est (dans le repère u):
    en faisant dxv*Dj-dxb*Dv, on a: x0=dxv*B/(dxv*dyb-dxb*dyv)=dxv*A
    en faisant dyb*Dv-dyv*Dj, on a: y0=dyv*B/(dxv*dyb-dxb*dyv)=dyv*A
  */
  A=det(dxa,dya,dxb,dyb)/det(dxv,dyv,dxb,dyb);
  return hv/hypot(dxv*A,dyv*A);
}


double func1(const double k,const void *n){
/*
  Fonction définie par:

                   f(k,n) := 2*n/k + (k-1)*(H(k)+1)

  avec H(k) := 1+1/2+1/3+...1/k ~ ln(k)+0.577... + 0.5/k + o(1/k),
  c'est-à-dire le k-ième nombre harmonic.  Le minimum de cette
  fonction est 2*sqrt(n*ln(n*ln(n))) et est atteint pour un k ~
  0.5*sqrt(2n/ln(n/ln(n)) ce qui est toujours dans l'intervalle [1,n].
*/
#define EULER_MASCHERONI 0.5772156649015328606
  return 2.0*(*((int*)n))/k + (k-1.0) * (log(k) + EULER_MASCHERONI + 0.5/k + 1.0);
}


double Minimize(double (*f)(double,const void*),void *info,double ax,double bx,double tol){
/*
  Calcule et renvoie l'abscisse x0 de l'intervalle [ax,bx] du minimum
  de la fonction f(x,info). La valeur de tolérance "tol" indique la
  précision souhaitée, tol=0 signifiant la précision maximale. Pour la
  cherche d'un maximum, il suffit de prendre -f(x,info).
  
  L'algorithme est un mélange de la recherche selon le nombre d'or
  (afin de minimisé le nombre d'appels à f(x,info) et de
  l'interpolation quadratique (pour mieux approché la
  solution).

  Pour la recherche selon le nombre d'or, l'idée est qu'on choisit
  deux points v<w de [a,b] (au départ a=ax et b=bx), et de réduire la
  recherche à [a,w] ou à [v,b], intervalles se chevauchant. Si
  f(v)<f(w), alors le minimum ne peut être que dans [a..v..w] à cause
  du point v. On recommence donc avec [a,w] sinon avec [v,b]. Un choix
  judicieux de v et w (basé sur le nombre d'or) permet de ne calculer
  f() sur qu'un seul nouveau point.

  Pour l'interpolation quadratique, l'idée est qu'avec un point x de
  [a,b] calcule la parabole passant par f(a),f(b) et f(x), et on prend
  comme nouveau milieu le point le plus bas de la parabole dans [a,b].
  L'algorithme mélange les deux techniques.

  Voir le "free software optimize.c" du R's project dans:
  https://svn.r-project.org/R/trunk/src/library/stats/src/

  Notes sur math.h:
  - nextafter(x,y)=plus petit double après x en direction de y 
  - fma(x,y,z)=x*y+z
  - fdim(x,y)=x-y si x>y, et 0 sinon
*/

#define K_MINIMUM        0.3819660112501051518 // = (3-sqrt(5))/2

  double a,b,d,e,p,q,r,u,v,w,x;
  double fu,fv,fw,fx,xm,tol1,tol2,tol3;

  static double eps=-1; /* on calcule eps qu'une seul fois */
  if(eps<0) eps=sqrt(nextafter(0,1)); /* racine carrée du plus petit double > 0 */
  if(tol<=0) tol=1E-10; /* en dessous, cela ne marche pas toujours !?! */

  a=ax,b=bx;
  w=v=x=fma(K_MINIMUM,b-a,a); // K*(b-a)+a
  fw=fv=fx=f(x,info);

  tol1=nextafter(1,2); /* le plus petit double > 1 */
  tol3=tol/3;
  d=e=0;

  for(;;){
    xm=(a+b)/2; // xm=milieu de [a,b]
    tol1=eps*fabs(x)+tol3;
    tol2=2*tol1;

    /* critère d'arrêt */
    if(fabs(x-xm)<=tol2-(b-a)/2) break;

    p=q=r=0;
    if(fabs(e)>tol1){ /* fit parabola */
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=(q-r)*2;
      if(q>0) p=-p; else q=-q;
      r=e;
      e=d;
    }

    if((fabs(p)>=fabs(q*r/2))||(p<=q*(a-x))||(p>=q*(b-x))){
      /* étape: recherche nombre d'or */
      e=(x<xm)? b-x : a-x;
      d=K_MINIMUM*e;
    }
    else{ /* étape: interpolation quadratique */
      d=p/q;
      u=x+d; 
      /* u ne doit pas trop près de a ou b */
      if((u-a<tol2)||(b-u<tol2)){ d=tol1; if(x>=xm) d=-d; }
    }
    
    /* on va évaluer f en u, qui ne doit pas être trop près de x */
    if(fabs(d)<tol1) u=(d>0)? x+tol1 : x-tol1;
    else u=x+d;

    fu=f(u,info); // calcul de f(u)

    /* met à jour a,b,v,w,x puis recommence */
    if(fu<=fx){
      if(u<x) b=x; else a=x;
      v=w; fv=fw;
      w=x; fw=fx;
      x=u; fx=fu;
    }else{
      if(u<x) a=u; else b=u;
      if((fu<=fw)||(w==x)){
	v=w; fv=fw;
	w=u; fw=fu;
      }else
	if((fu<=fv)||(v==x)||(v==w)){ v=u; fv=fu; }
    }
    
  }
  
  return x;
}


enum{
  DYCK_WALK,
  DYCK_WORD,
  DYCK_TREE,
  DYCK_KTREE,
};


int* Dyck(int *R,const int n,const int k,const int code){
/*
  Construit de manière aléatoirement uniforme un mot de k-Dyck composé
  de n segments montant (codé par la lettre k), chacun de longueur k,
  et de kn pas descendant (codé par 0). On doit avoir n>0 et k>0. Le
  résultat est donc en général un mot de n+kn lettres composés de n
  valeurs k et de kn valeurs 0. La propriété de ce mot est que tout
  préfixe possède au moins autant de pas montant que descendant. Le
  mot classique de Dyck et les arbres binaires sont obtenus pour k=1,
  une suite de 0 et de 1.

  Le résultat est renvoyé et mis dans R, qui est alloué si R=NULL.
  Suivant la valeur de code on renvoit des variantes du mot de
  Dyck. Attention, la taille de R, fonction de n et k, varie avec code
  (voir ci-dessous).

    code    taille  résultat
  --------------------------------------------------------------------
  DYCK_WALK  n+nk   mot ayant n valeurs k et kn valeurs 0
  DYCK_WORD  n+nk   même mot mais décalé sur un minimum
  DYCK_TREE  nk+1   arbre DFS à nk+1 noeuds
  DYCK_KTREE n+nk+1 arbre (k+1)-ary à n+nk+1 noeuds dont n de deg. k+1
  --------------------------------------------------------------------

  On peut obtenir un mot de longueur p ayant q valeurs non nulles
  aléatoire uniforme avec DYCK_WALK en choisissant n=q et k=p/q-1.

  Exemple pour n=5 et k=1:
  code=DYCK_WALK -> 1101001010
                                            3   4           3   4
                       /\/\                / \ / \           \ /
                    /\/    \/\        1   2   2   2   5     1 2 5
  code=DYCK_WORD -> 1011010010       / \ /         \ / \     \|/
  code=DYCK_TREE -> [-1,0,0,2,2,0]  0   0           0   0     0

  L'algorithme principal, en O(kn), consiste à tirer aléatoirement n
  valeurs k et k*n valeurs 0. On tire k avec une probabilité
  proportionelle au nombre de valeurs k restant à tirer sur le nombre
  total de valeurs restant à tirer. C'est un tirage aléatoire
  uniforme.

  L'arbre DFS (DYCK_TREE) est construit de sorte que le mot de Dyck
  forme le parcours DFS de l'arbre: 1 on descend le long d'une arête,
  et 0 on en revient. NB: Si k>1, on fait comme si on avait k pas de
  1. Les propriétés de cet arbre sont nombreuses. Par exemple:
  - R[0] = -1, car 0 est la racine.
  - nk = dernière feuille de l'arbre, R[nk] est donc sont père
  - si u-R[u]>1 alors u démarre une nouvelle branche

  L'arbre DYCK_KTREE est un DFS modifé: on posant d'abord tous les
  fils avant la récursion. Voici quelques exemples.

  Exemple pour n=2 et k=2: arbre ternaire à 2 sommets internes

                         /\           4 5 6
                      /\/  \           \|/
                     /      \         1 2 3
  code=DYCK_WORD  -> 2 02 000          \|/
  code=DYCK_KTREE -> [-1,0,0,0,2,2,2]   0

  Exemple pour n=3 et k=1: arbre binaire à 3 sommets internes
  code=DYCK_WORD  -> 110100
  code=DYCK_KTREE -> [-1,0,0,1,1,4,4]

        5   6
         \ /
      3   4
       \ /
        1   2
         \ /
          0

  Exemple pour n=3 et k=2: arbre ternaire à 3 sommets internes
  code=DYCK_WORD  -> 202200000
  code=DYCK_KTREE -> [-1,0,0,0,2,2,2,4,4,4]

      7 8 9
       \|/
        4 5 6
         \|/
        1 2 3
         \|/
          0

  Pour construire l'arbre (k+1)-ary (DYCK_KTREE) on procède selon un
  DFS modifié (on pose les fils avant la récursion), en lisant
  succéssivement les n+kn valeurs du mot de k-Dyck. Plus précisément,
  depuis le sommet courant u: si on lit k, on ajoute k+1 fils à u, le
  nouveau sommet courant devenant le 1er fils de u. Si on lit 0, on
  remonte les ancêtres de u jusqu'à trouver un ancêtre v dont le fils
  f menant à u ne soit pas le dernier fils (le (k+1)-ième). Le nouveau
  sommet courant est alors f+1. Pour déterminer si f est le dernier
  fils (ou pas) il suffit de tester si f%(k+1)=0 ou pas.
*/

  const int m=n+k*n; /* longueur du mot de Dyck */
  int *B; /* mot de Dyck */
  int i;

  if(R==NULL){
    if(code==DYCK_WALK)  ALLOC(R,m);
    if(code==DYCK_WORD)  ALLOC(R,m);
    if(code==DYCK_TREE)  ALLOC(R,m-n+1);
    if(code==DYCK_KTREE) ALLOC(R,m+1);
  }
  if(code==DYCK_WALK) B=R; else ALLOC(B,m);

  /* DYCK_WALK: construit B */
  /* commun à DYCK_WORD et DYCK_TREE */
  
  int t=n; /* t=nombre de k à tirer */
  int s=m; /* s=nombre total de valeur à tirer */

  for(i=0;s>0;s--,i++){
    B[i]=((random()%s)<t); /* k avec proba proportionnel à t */
    if(B[i]) B[i]=k,t--; /* enlève un k */
  }
  if(code==DYCK_WALK) return R; /* NB: ici R=B */
  
  /* cherche la position r dans B de la hauteur minimum */
  /* commun à DYCK_WORD et DYCK_TREE */

  int r=-1;
  int h=0;
  
  for(i=s=0;i<m;i++){
    s += B[i]? k : -1; /* s=hauteur courante = +k ou -1 */
    if(s<h) h=s,r=i;   /* h=hauteur minimum */
  }
  r=(r+1)%m; /* la position r recherchée */
  
  /* DYCK_WORD: décale le mot de B vers R */
  
  if(code==DYCK_WORD){
    for(i=0;i<m;i++) R[i]=B[(i+r)%m];
    goto fin_dyck;
  }

  /* DYCK_TREE: décale et construit un arbre DFS, R[v]=parent du sommet v */
  /* Si k>1, on fait comme si on avait k pas de 1 */

  if(code==DYCK_TREE){
    int u=0; /* u=dernier sommet visité (=racine) */
    int v=1; /* v=prochain sommet à visiter (=sommet 1) */
    R[0]=-1; /* père de la racine */
  
    for(i=0;i<m;i++) /* parcoure toutes les valeurs de B */
      if(B[t=(i+r)%m])
	while(B[t]){ /* tantque c'est > 0, on monte */
	  R[v]=u; /* père de v=dernier sommet visité */
	  u=v++; /* met à jour le dernier sommet visité */
	  B[t]--;
	}
      else u=R[u]; /* si c'est un 0, on descend */
    goto fin_dyck;
  }

  /* DYCK_KTREE: on a n noeuds internes de degré k+1 */

  if(code==DYCK_KTREE){
    int u=0,s=1; /* u=sommet courant, indexe dans R */
    R[0]=-1;     /* père de la racine */
    for(i=0;i<m;i++) /* parcoure toutes les valeurs de B */
      if(B[(i+r)%m]){
	for(t=0;t<=k;t++) R[s++]=u; /* si on lit k, on pose k+1 fils */
	u=s-t; /* le nouveau sommet courant est le 1er fils = s-k-1=s-t */
      }
      else{
	while(u%(k+1)==0) u=R[u]; /* si on lit 0, on remonte tous les derniers fils */
	u++; /* le nouveau sommet courant est u+1 */
      }
    goto fin_dyck;
  }

 fin_dyck:
  free(B);
  return R;
}


int NextPermutation(int *P,const int n,const int *C){
/*
  Génère, à partir d'une permutation P, la prochaine dans l'ordre
  lexicographique suivant les contraintes définies par le tableau C
  (voir ci-après). Mettre C=NULL s'il n'y a pas de contrainte
  particulière. On renvoie 1 si la prochaine permutation à pu être
  déterminée et 0 si la dernière permutation a été atteinte. Dans ce
  cas la permutation la plus petite selon l'ordre lexicographique est
  renvoyée. On permute les éléments de P que si leurs positions sont
  entre C[j] et C[j+1] (exclu) pour un certain indice j. On peut
  initialiser P avec ALLOCZ(P,k,_i) ou si le tableau P existe déjà
  avec NextSet(P,-1,k). On l'utilise comme ceci (ici sans contrainte
  C):

  NextSet(P,-1,n);
  do{
    PRINTT(P,n);
    ...;
  }while(NextPermutation(P,n,NULL));

  Ex: si C={2,3,5}. Les permutations sous la contrainte C sont:
  (on peut permuter les indices {0,1}{2}{3,4})

                 0 1 2 3 4 (positions dans P)
	      P={a,b,c,d,e}
	        {b,a,c,d,e}
		{a,b,c,e,d}
		{b,a,c,e,d}
  
  Evidemment, il y a beaucoup moins de permutations dès que le nombre
  de contraintes augmente. Par exemple, si C contient k intervalles de
  même longueur, alors le nombre de permutations sera de (n/k)!^k au
  lieu de n!. Le rapport des deux nombre est d'environ k^n.

  Concrêtement, pour:
  - n=9 et k=3, on a 216 permutations au lieu de 362.880 (k^n=19.683)
  - n=12 et k=3, on a 13.824 permutations au lieu de 479.001.600 (k^n=531.441)

  Le dernier élément de C doit être égale à n-1 (sentinelle), le
  premier étant omis car toujours = 0. Donc C est un tableau à au plus
  n éléments. Si C=NULL, alors il n'y a pas de contrainte
  particulière, ce qui est identique à poser C[0]=n.

  On se base sur l'algorithme classique (lorsque C=NULL, sinon on
  l'applique sur l'intervalle de positions [C[j],C[j+1][):

  1. Trouver le plus grand index i tel que P[i] < P[i+1].
     S'il n'existe pas, la dernière permutation est atteinte.
  2. Trouver le plus grand indice j tel que P[i] < P[j].
  3. Echanger P[i] avec P[j].
  4. Renverser la suite de P[i+1] jusqu'au dernier élément.

*/
  int i,j,a,b,c,T[1];

  if(C==NULL){
    T[0]=n;
    C=T;
  }

  b=C[i=j=0]; /* j=indice de la prochaine valeur à lire dans C */
  c=-1;

  /* étape 1: on cherche l'intervalle [a,b[ avec i tq P[i]<P[i+1] */
 etape1:
  for(a=i;i<b-1;i++) if(P[i]<P[i+1]) c=i; /* on a trouvé un i tq P[i]<P[i+1] */
  if(c<0){ /* le plus grand i tq P[i]<[i+1] n'existe pas */
    for(i=a;i<b;i++) P[i]=i; /* on réinitialise P[a]...P[b-1] */
    if(b==n) return 0; /* alors on a fini d'examiner C */
    b=C[++j]; /* b=nouvelle borne */
    goto etape1;
  }
  i=c; /* i=le plus grand tq P[i]<P[i+1] avec a<=i,i+1<b */

  /* étape 2: cherche j=le plus grand tq P[i]<P[j] */
  for(j=i+1;j<b;j++) if(P[i]<P[j]) c=j;
  j=c;

  /* étape 3: échange P[i] et P[j] */
  SWAP(P[i],P[j],c);

  /* étape 4: renverse P[i+1]...P[b-1] */
  for(++i,--b;i<b;i++,b--) SWAP(P[i],P[b],c);

  return 1;
}


int NextSet(int *S,const int n,const int k){
/*
  Calcule les sous-ensembles de k entiers de [0,n[. Si n<0, alors S
  est initialisé au plus petit ensemble possible, soit S={0,1,2,
  ...,k-1}. L'idée est de maintenir une sorte de compteur S qui
  représente le sous-ensemble courant et qui va passer par tous les
  sous-ensembles possibles. Les éléments sont toujours rangés dans
  l'ordre croissant. On renvoie 1 si on a pu construire le prochain
  sous-ensemble, et 0 si S était le dernier sous-ensemble. Dans ce cas
  l'ensemble le plus petit est écrit dans S. On l'utilise comme ceci:

  NextSet(S,-1,k);
  do{
    ...;
    PRINTT(S,k); // traitement de S
    ...;
  }while(NextSet(S,n,k));

  On peut facilement en déduire un algorithme pour générer des
  multi-ensembles, comme par exemple tous les multi-ensembles du type
  [2,0,0,2,1,0] comprennant 3 fois 0, 1 fois 1 et 2 fois 2:
  NextMultiSet(S,C,k) avec C=[3,1,2] (voir la fonction ggosset()).

  La stratégie pour "incrémenter" S est la suivante : on essaye
  d'incrémenter S[i] (au départ i=0) tout en restant strictement
  inférieur à l'élément suivant S[i+1]. Si cela marche on a trouvé le
  prochain sous-ensemble, sinon on pose S[i]=i et on recommence avec
  indice i suivant. Si S[k-1] atteint n c'est que S était le dernier
  sous-ensemble.

  L'algorithme est en O(k) dans le pire des cas, mais de manière
  amortie c'est beaucoup moins car on incrémente moins souvent S[j]
  que S[i] si j>i.
*/
  int i=0,j,s;

  if(n<0){
    for(;i<k;i++) S[i]=i;
    return 1;
  }
  
  while(i<k){
    s=++S[j=i++];
    if(i==k){ if(s<n) return 1; }
    else{ if(s<S[i]) return 1; }
    S[j]=j;
  }
  
  return 0;
}


int NextArrangement(int *S,int *P,const int n,const int k){
/*
  Permet de générer tous les arrangements de k entiers de
  [0,n[. L'arrangement est représenté par les tableaux S et P de k
  entiers. S représente un ensemble de k entiers et P une permutation
  de [0,k[. Ainsi, l'arrangement A=(4,2,7,3) est représenté par
  S=(2,3,4,7) et P=(2,0,3,1). Autrement dit A[i]=S[P[i]] pour tout
  i=0...k-1.

  L'idée est d'essayer d'incrémenter le double compteur S,P. On
  incrémente P en premier avec NextPermutation(). Si on est arrivé à
  la fin de P, on incrémente S avec NextSet(). Si n<0, alors S et P
  sont initialisés au plus petit arrangement possible, soit S = P =
  (0,1,2, ...,k-1). On renvoie 1 si on a pu trouver le prochain
  arrangement, 0 si c'était le dernier arrangement possible. Dans ce
  cas l'arrangement le plus petit est écrit dans S,P.
*/
  if(n<0){
    int i;
    for(i=0;i<k;i++) S[i]=P[i]=i;
    return 1;
  }
  if(!NextPermutation(P,k,NULL)) return NextSet(S,n,k);
  return 1;
}


int *NextPart(int *S,const int n,int s,int *C){
/*
  Permet de générer toutes les suites S de n>0 entiers >=0 dont la
  somme fait s et dont la i-ème part S[i] ne dépasse pas C[i]. Il faut
  que s <= ∑_{i=0}^(n-1) C[i], n=1 est possible, de même que C[i]=s.

  Initialement S est la suite courante de somme s et on renvoie dans S
  la prochaine suite (la fonction renvoie aussi S). On renvoie NULL si
  on a atteint la dernière suite, et on remplit S avec le première
  suite. Si S=NULL, alors S est allouée et initialisée à la première
  suite. La première suite de somme s est obtenue en remplissant
  autant que possible les parts S[n-1],S[n-2],...

  L'algorithme est le suivant:
   1. on détermine le plus grand indice j tq S[j]>0
   2. si j=0, alors on a finit: on va à l'étape 6 avec x=s+1 et i=-1
   3. on détermine le plus grand indice i<j tq S[i] peut être augmenté
   4. on calcule x = ∑_{j=i+1}^(n-1) S[i]
   5. on incrémente S[i]
   6. on remplit S[i+1]...S[n-1] avec la première suite de somme x-1

  On l'utilise comme ceci:

  int s=n/2;
  NALLOCZ(int,C,n,1); // compteur binaire avec au plus n/2 valeurs 1
  int *S=NextPart(NULL,n,s,C); // initialisation de S à 0
  do{ PRINTT(S,n); // traitement de la partition S
      ...;
  }while(NextPart(S,n,s,C));

  Exemple: s=n=5

  C = 1 2 2 1 1
  S = 0 1 2 1 1
      0 2 1 1 1
      0 2 2 0 1
      0 2 2 1 0
      1 0 2 1 1
      1 1 1 1 1
      1 1 2 0 1
      1 1 2 1 0
      1 2 0 1 1
      1 2 1 0 1
      1 2 1 1 0
      1 2 2 0 0
*/
  int x,i,j,r;

  i=0;
  r=(S==NULL);
  if(r) ALLOC(S,n);
  else i=n-1;
  
  /* calcule le plus grand indice i tq S[i]>0 */
  while((i>0)&&(S[i]==0)) i--;
  x=S[i--]; /* rem: si i=0, alors i -> -1 */

  /* calcule le plus grand indice j<i tq S[j] peut être augmenté */ 
  while((i>=0)&&(S[i]==C[i])) x += S[i--];

  if(i>=0){ S[i]++; s=x-1; } /* si cet indice n'existe pas i<0 => FIN */
  
  /* écrit la première suite de somme s dans S[i+1]...S[n-1] */
  for(j=n-1;j>i;j--){
    x=imax(s,0);
    x=imin(C[j],x);
    S[j]=x;
    s -= x;
  }

  /* on retourne S sauf si i<0 et r=0 (<=> FIN ) */
  return ((i<0)&&(!r))? NULL : S;
}


int *NextIntPartition(int *S,const int n, int *t){
/*
  Permet de générer toutes les partitions d'un entier n>0 dans l'ordre
  lexicographique croissant. L'entier retourné dans *t est le nombre
  de parts de la partition renvoyée. Il est possible de mettre t=NULL
  et dans ce cas le nombre de parts n'est pas renvoyé.  Si S=NULL,
  alors S est allouée (de taille n) et initialisée à la première
  partition de n, soit S={1,...,1}. Dans tous les cas, S est renvoyée.
  Si S valait la dernière partition de n, soit S={n}, alors NULL est
  renvoyé et *t=1.

  On l'utilise comme ceci:

  int t; // t=pour connaître le nombre de parts
  int *S=NextIntPartition(NULL,n,&t); // initialisation de S
  do{ PRINTT(S,t); // traitement de la partition S en t parts
      ...;
  }while(NextIntPartition(S,n,&t));

  S = 1 1 1 1 1 1
  S = 1 1 1 1 2
  S = 1 1 1 3
  S = 1 1 2 2
  S = 1 1 4
  S = 1 2 3
  S = 1 5
  S = 2 2 2
  S = 2 4
  S = 3 3
  S = 6

  La complexité en temps est constante en moyenne. L'algorithme est
  basé sur http://jeromekelleher.net/category/combinatorics.html.
*/
  static int k=-1; // k en statique pour aller plus vite
  if((S==NULL)||(k<0)){ // première fois ?
    if(S==NULL) ALLOC(S,n); // alloue S
    S[k=1]=n; S[0]=0; // pour initialiser S
  }
  
  if(k==0) return NULL; // dernière part atteinte
  int y=S[k--]-1;
  int x=S[k]+1;
  while(x<=y) S[k++]=x, y-=x;
  S[k]=x+y; 
  if(t) *t=k+1; // nombre de parts
  return S;
}


int *RandomIntPartition(int *S,const int n){
/*
  Produit une partition aléatoire uniforme de l'entier n>0. Le tableau
  S doit être de taille n au moins. Si S=NULL, il est alloué et
  renvoyé. En retour S[i] représente le nombre de fois où la valeur
  i+1 est dans la partition. La complexité moyenne de l'algorithme est
  en O(n^1.25).

  L'algorithme à rejet ci-dessous calcule une partition S =
  (Z_1,...,Z_n) de n, Z_i étant le nombre de fois où l'entier i>0
  apparaît dans la partition:
  
  1. Soit X=-π/√(6n).

  2. Générer n-1 variables indépendantes Z_2,...,Z_n, où Z_i suit une
     loi géométrique de paramètre p_i=1-exp(X)^i, donnant le nombre
     d'échecs avant le 1er succès de probabilité p_i. Pour générer une
     telle variable il suffit de faire Z_i = floor{ ln(U) / ln(1-p_i)
     } = floor{ ln(U) / (i*X) } où U=RAND01. Poser Z_1 = n - ∑_{i=2}^n
     i*Z_i.

  3. Si Z_1<0 ou ln(RAND01)>Z_1*X, recommencer en 2.

  D'après [AD15] le nombre de répétitions moyen est au plus
  2*π*(6n)^(1/4) < 10*n^0.25. On peut faire O(1) répétitions, mais
  c'est bien plus compliqué. En fait, l'algorithme pourrait être
  accéléré car, avec grande probabilité (cf. Erdos-Lehner 1941), la
  plus grande part est proche de 2c*√n*log(n) où c=π/√6 et donc Z_i=0
  lorsque i>>√n*log(n). Donc on pourrait utiliser un seul random pour
  déterminer à partir de quelle valeur de i on met tous les Z_i=0. Il
  est montré dans [AD15, pp. 22] qu'il est possible de faire O(√n)
  appels à random pour générer tous les Z_i.

  Une méthode qui génère aussi une partition aléatoire, mais pas
  uniformément, consiste à tirer un tableau S aléatoire de taille n
  dont les valeurs sont entières et >=0 et dont la somme fait n ce qui
  peut se faire avec la ligne de code suivante: ALLOCZ(S,n,0);
  for(i=0;i<n;i++) S[random()%n]++;
*/
  if(S==NULL) ALLOC(S,n);
  const double x=-M_PI/sqrt(6*n);
  int i,s; // s=S[0]=Z_1
  double p;

  do{
    s=n; i=1; p=x;
    while((s>=0)&&(i<n)){
      p += x; // au début p=2*x
      S[i]=(int)(log(RAND01)/p); // S[i] = Z_{i+1} = ln(U)/((i+1)*x)
      s -= (i+1)*S[i];
      i++;
    }
  }while((s<0)||(log(RAND01)>x*s));

  S[0]=s;
  return S;
}


int SetCmp(int *T1,int *T2,int n1,int n2)
/*
  Compare deux tableaux d'entiers T1 et T2 de taille n1 et n2 triés
  par ordre croissant. Les tableaux peuvent être de taille nulle. La
  valeur renvoyée est un entier interprété en binaire comme suit:

  bit-0: 1 ssi T1 intersecte T2 (possède au moins 1 élément commun)
  bit-1: 1 ssi T1 égale T2
  bit-2: 1 ssi T1 est strictement inclu dans T2
  bit-3: 1 ssi T2 est strictement inclu dans T1

  Les valeurs possibles sont donc: 0,1,2,3,4,5,8,9 (=[0,9]\{6,7})
  La complexité est O(n1+n2).
*/
{
  if(n1==0) return (n2==0)? 2:4;
  if(n2==0) return 8;
  /* ici T1 et T2 contiennent au moins 1 élément */

  if((T1[n1-1]<T2[0])||(T2[n2-1]<T1[0])) return 0; /* cas trivial de disjonction */

  int i1,i2,r;
  i1=i2=0;
  r=14; /* tous les bit à 1 sauf b0 */

  while((i1<n1)&&(i2<n2)){
    if(T1[i1]==T2[i2]){
      i1++; i2++; /* T1 et T2 ont un élément en commun */
      r |= 1; continue; /* met b0 */
    }
    r &= 13; /* annule b1 (15-2) car T1<>T2 */
    if(T1[i1]<T2[i2]){
      i1++; /* T1 contient des éléments qui ne sont pas dans T2 */
      r &= 11; /* annule b2 (15-4) car T1 ne peux pas contenir T1 */
    }else{
      i2++; /* T2 contient des éléments qui ne sont pas dans T1 */
      r &= 7; /* annule b3 (15-8) car T2 ne peux pas contenir T1 */
    }
  }

  if(i1<n1) r &= 9; /* annule b2 et b1 (15-4-2) */
  if(i2<n2) r &= 5; /* annule b3 et b1 (15-8-2) */
  if(r&2)   r &= 3; /* annule b3 et b2 (15-8-4) */

  return r;
}


int SetSearch(const int u,const int *T,const int n,const int sort){
/*
  Cherche l'entier u dans le tableau T de taille n. Si u est dans T,
  on renvoie son indice sinon on renvoie -1. Si sort=1, alors on
  suppose T trié par ordre croissant et la recherche est dichotomique,
  sinon la recherche est linéaire.
*/
  if(sort){ /* recherche dichotomique */
    int *t=bsearch(&u,T,n,sizeof(int),fcmp_int);
    return t? (int)(t-T) : -1;
  }
  /* recherche linéaire */
  int i;
  for(i=0;i<n;i++) if(u==T[i]) return i;
  return -1;
}


int Binom(const int n,const int k){
/*
  Calcule l'entier B={n choose k}.  L'algorithme utilisé ici est en
  O(k). Il n'utilise que des multiplications et divisions entières sur
  des entiers en O(B), sans aucun tableau.

  L'algorithme classique issu du Triangle de Pascal, qui lui est en
  O(n*k), utilise un espace en O(k) (tableaux d'entiers en O(B)). Par
  contre il n'utilise que des additions sur des entiers en O(B).

  Principe: B = n x (n-1) x ... x (n-k+1) / k x (k-1) x ... x 1

  On réecrit B en (((((n/1) x (n-1)) / 2) x (n-2)) / 3) ...
  c'est-à-dire en multipliant le plus grand numérateur et en divisant
  par le plus petit numérateur. Chaque terme est ainsi un certain
  binomial, et donc toujours un entier.

  Catalan(n) = Binom(2n,n)/(n+1). On peut aussi calculer ce nombre
  avec Catalan(0)=1 et Catalan(n) = (2*(2n-1)*Catalan(n-1)/(n+1).
*/
  int B,s,i;
  for(B=s=n,i=2;i<=k;i++) B=(B*(--s))/i;
  return B;
}


/* code pour la fonction SortInt() */
enum{
  SORT_INC,     /* tri croissant selon les valeurs */
  SORT_DEC,     /* tri décroissant selon les valeurs */
  SORT_INC_RANK,/* donne le rang des valeurs de T dans l'ordre croissant */
  SORT_DEC_RANK,/* donne le rang des valeurs de T dans l'ordre décroissant */
  SORT_FREQv,   /* donne la fréquence des valeurs de T */
  SORT_FREQe,   /* donne la fréquence des éléments de T */
  SORT_INDEXi,  /* donne l'indice des éléments de T dans l'ordre croissant */
  SORT_INDEXe   /* donne les éléments de T dans l'ordre croissant */
};


int *SortInt(int *T,int *R,const int n,const int a,int *m,const int code){
/*
  Trie par ordre croissant un tableau T non NULL de taille n>0 dont
  les valeurs sont des entiers de [a,a+*m[ si m<>NULL. Sinon, les
  valeurs de T sont supposées être dans [a,a+n[. Pour simplifier le
  texte qui suit, je note m pour dire *m.

  La complexité en temps est O(n+m), ce qui est mieux que qsort(). Le
  tableau T n'est pas modifié. Le résultat est rangé dans le tableau R
  qui doit être de taille au moins n, et est aussi renvoyé par la
  fonction. Il R=NULL, alors R est alloué et retourné par la
  fonction. Pour être compétitif avec qsort() pour m=n il faut,
  théoriquement, n>32.

  L'opération de tri dépend de "code" qui est interprétée comme suit
  (la complexité est le nombre d'étapes dans les boucles "for"):

    v = valeur d'un élément de T, v dans [a,a+m[
    d = nombre de valeurs v distinctes de T, d dans [1,min(n,m)]
    e = indice du tableau T, e dans [0..n[
    i = indice du tableau R, i dans [0..n[
    r = rang dans un tableau, r dans [0..d[

  - code=SORT_FREQv: renvoie un tableau F[0..m[ de fréquence des
    valeurs de T, c'est-à-dire que F[i] est le nombre de fois où la
    valeur i+a apparaît dans T. Le tableau R et la variable m ne sont
    pas modifiés. Complexité: m+n.

    Ex: T = [12 11 12 13 12 15 16 13]  avec n=8, a=10, m=7
        F = [0 1 3 2 0 1 1] de taille m

  - code=SORT_FREQe: renvoie dans R[0..n[ un tableau de fréquence des
    éléments de T, c'est-à-dire où R[e] est le nombre de fois où T[e]
    apparaît dans T. La variable m n'est pas modifiée. Complexité:
    m+2n.

    Ex: T = [12 11 12 13 12 15 16 13]  avec n=8, a=10, m=7
        R = [ 3  1  3  2  3  1  1  2] de taille n

  - code=SORT_INC ou SORT_DEC: renvoie dans R[0..n[ le tableau T trié
    par ordre croissant (ou décroissant). La variable m n'est pas
    modifiée. Complexité: 2m+2n.

    Ex: T = [12 11 12 13 12 15 16 13]  avec n=8, a=10, m=7
        R = [11 12 12 12 13 13 15 16] de taille n

  - code=SORT_INC_RANK ou SORT_DEC_RANK: renvoie dans R[0..n[ un
    tableau de rangs où R[e] est le rang r dans [0..d[ de l'élément
    T[e] dans la version triée dans l'ordre croissant (ou décroissant)
    de T. Le tableau R est modifié et on renvoie d dans m. Complexité:
    3m+2n.

    Ex: T = [12 11 12 13 12 15 16 13]  avec n=8, a=10, m=7
        R = [ 1  0  1  2  1  3  5  2] et m=5

  - code=SORT_INDEXi: renvoie dans R[0..n[ un tableau d'indices où
    R[i]=e est l'élément de T en position i dans la version triée par
    ordre croissant de T. Pour obtenir un tri de T il suffit de lister
    T[R[i]] pour i=0..n-1. La variable m n'est pas modifiée.
    Complexité: 2m+2n.

    Ex: T = [12 11 12 13 12 18 15 11]  avec n=8, a=10, m=9
        R = [ 1  7  0  2  4  3  6  5]

  - code=SORT_INDEXe: renvoie dans R[0..n[ un tableau d'indices où
    R[e] est la position de T[e] dans la version triée par ordre
    croissant de T. La variable m n'est pas modifiée. Complexité:
    2m+2n.

    Ex: T = [12 11 12 13 12 18 15 11]  avec n=8, a=10, m=9
        R = [ 2  0  3  5  4  7  6  1]
*/
  int i,r,t;
  int k=(m==NULL)? n:*m;

  /* initialise F[0..m[ */
  NALLOCZ(int,F,k,0); /* coût: m */
  
  /* calcule F[i]=fréquence de la valeur v=i+a dans T */
  for(i=0;i<n;i++) F[T[i]-a]++; /* coût: n */

  if(code==SORT_FREQv) return F;

  /* alloue R, si nécessaire */
  if(R==NULL) ALLOC(R,n);

  if(code==SORT_FREQe){
    for(i=0;i<n;i++) R[i]=F[T[i]-a]; /* coût: n */
    free(F); return R;
  }

  if(code==SORT_INC){ /* R=tableau T trié, ordre croissant */
    for(i=r=0;i<k;i++) for(t=F[i];t>0;t--) R[r++]=i+a; /* coût: m+n */
    free(F); return R;
  }

  if(code==SORT_DEC){ /* R=tableau T trié, ordre décroissant */
    for(i=0,r=n;i<k;i++) for(t=F[i];t>0;t--) R[--r]=i+a; /* coût: m+n */
    free(F); return R;
  }

  /* calcule F[i]=nombre de valeurs de T qui sont < i+a.     
     Ex:  T = [2 1 2 3 2 8 5]  avec n=7, m=9, a=0
	  F = [0 0 1 4 5 5 6 6 6]
  */
  for(i=r=0;i<k;i++){ t=F[i]; F[i]=r; r += t; } /* coût: m */

  if(code==SORT_INDEXi){
  for(i=0;i<n;i++) R[F[T[i]-a]++]=i; /* coût: n */
    free(F); return R;
  }

  if(code==SORT_INDEXe){
    for(i=0;i<n;i++) R[i]=F[T[i]-a]++; /* coût: n */
    free(F); return R;
  }

  /* pour SORT_INC_RANK ou SORT_DEC_RANK */
  /* calcule dans F[i] le rang dans [0,d[ du nb de valeurs < i+a */
  
  for(t=r=-1,i=0;i<k;i++){ /* coût: m */
    if(F[i]!=t) { t=F[i]; r++; }
    F[i]=r;
  }
  if(m!=NULL) *m=r+1; /* m=nb de valeurs différentes de T */

  if(code==SORT_INC_RANK){
    for(i=0;i<n;i++) R[i]=F[T[i]-a]; /* coût: n */
    free(F); return R;
  }
  
  if(code==SORT_DEC_RANK){
    for(i=0;i<n;i++) R[i]=r-F[T[i]-a]; /* coût: n */
    free(F); return R;
  }

  free(F);
  free(R);
  Erreur(21);
  return NULL;
}


/***********************************

       BFS, DFS, ...
       (DIJKSTRA)

***********************************/


/*
  Structure de données pour le BFS
*/
typedef struct{
  int root;  /* racine ou source du BFS. */
  int radius;/* eccentricité du sommet source. */
  int *D;    /* D[u]=distance entre u et root. Les sommets u avec
	        D[u]=-1 sont à une distance infinie de root (situation
	        initiale par défaut). En entrée, les sommets avec
	        D[u]=-2 sont considérés comme inexsitant dans le
	        graphe. Si D=NULL, D est alloué puis initialisé à
	        -1. */
  int *P;    /* P[u]=père de u dans un arbre BFS de racine root, avec
	        P[root]=-1. Seuls les sommets u<>root avec D[u]>=0 ont
	        un père défini. Si P=NULL, alors P est alloué. Il
	        n'est pas nécessaire de l'initialiser. */
  int n;     /* nombre de sommets parcourus = nombre d'éléments dans
		la file */
  int *file; /* contenu de la file = liste des n sommets parcourus. La
		taille de ce tableau est toujours le nombre de sommets
		du graphe. */
  int cycle; /* longueur du plus petit cycle passant par la source
		rencontré lors du parcours. Cette valeur (>2) est
		indépendente du parcours spécifique des voisins (ie de
		tout BFS). Elle ne dépend que de la structure
		non-étiquetée du graphe et de la source. Si cycle<0,
		alors cette longueur est infinie.  On peut déterminer
		la maille du graphe en prennant le minimum obtenu pour
		chacun des sommets du graphe comme source. Cette
		valeur n'est pas correcte si le graphe est orienté. */
  int clean; /* permet d'initialiser à -1 les sommets parcourus dans
		le précédent bfs(). Par défaut, clean=0, et le tableau
		D n'est pas initialisé (sauf si D=NULL, dans ce cas il
		est initialisé complètement à -1). Si clean=2, alors
		on remet D[u] à -1 seulement pour les sommets u de
		file. C'est plus rapide si bfs() est utilisé plusieurs
		fois avec le même paramètre et tous les sommets ne
		sont pas visités (vmax ou hmax >0). Si clean=1, alors
		on force l'initialisation complète de D à -1, puis on
		passe clean=2 (pour les appels suivants). */
  int tf;    /* tête de file, ie nombre d'éléments parcourus dont tous
		les voisins ont été enfilés (tf<=n). */
  int cont;  /* cont=1 si on poursuit le bfs() précédant là où on
		s'était arrêté dans le cas d'un arrêt par hmax. Par
		défaut cont=0. Après le 1er bfs() où cont=1, cont
		passe à 2. Cela sert à augmenter progressivement la
		hauteur jusqu'à obtenir une certaine condition. */
  int vmax;  /* arrête le parcours lorsque vmax sommets on été parcourus */
  int hmax;  /* arrête le parcours lorsqu'un sommet de hauteur > hmax est atteint */
} param_bfs;


param_bfs *new_param_bfs(void){
/*
  Crée et initialise une structure pour la fontion bfs(). C'est la
  fonction bfs() qui doit se charger, éventuellement, d'allouer les
  tableaux P et D (qui sont ici initialisé à NULL). Le tableau D peut
  être utilisé pour effacer des sommets.
*/
  NALLOC(param_bfs,X,1);
  X->D=X->P=X->file=NULL;
  X->radius=X->n=X->clean=X->cont=X->tf=0;
  X->cycle=X->root=X->vmax=X->hmax=-1;
  return X;
}


void free_param_bfs(param_bfs *X){
  if(X==NULL) return;
  free(X->D);
  free(X->P);
  free(X->file);
  free(X);
}


param_bfs *bfs(const graph *G,const int source,param_bfs *X){
/*
  Effectue un parcours en largeur (BFS) d'un graphe G (orienté ou non)
  depuis le sommet source. Les résultats du parcours (comme les
  distances à la source, le père, le nombre de sommets parcourus,
  etc.) sont stockées dans la variable X qui est renvoyée. Si X=NULL,
  X est d'abord allouée.

  On peut également réaliser un BFS seulement sur un sous-ensemble de
  sommets (sous-graphe), limiter le parcours à une profondeur donnée
  (en fixant X->hmax), ou à un nombre de sommets parcourus (en fixant
  X->vmax). Il est important que la table X->D soient correctement
  initialisée pour les sommets à parcourir. Les tables X->P et X->file
  ne sont jamais initialisés. Si elles sont NULL, elles sont
  simplement allouées. Seules les valeurs des X->n sommets parcourus
  (ceux dans X->file) sont garanties d'être correctes. A noter que le
  parcours est lancé depuis la source, même si source est dans T.

  On l'utilise comme ceci:

    param_bfs *X=bfs(G,s,NULL); // BFS depuis le sommet s dans G
    ... // X->D[u]=distance entre s et u
    ... // X->P[u]=père de u, ou -1 s'il n'existe pas
    ... // X->cycle=longueur du plus petit cycle passant par s
    ... // X->n=nombre de sommets parcourus
    ... //
    free_param_bfs(X); // libère la variable crée X

  Pour réaliser un BFS d'un sous-graphe G\T de G (évitant les k
  sommets du tableau T):

    param_bfs *X=new_param_bfs();   // par défaut X->clean=0
    ALLOCZ(X->D,G->n,-1);           // alloue X->D et met tout à -1
    for(i=0;i<k;i++) X->D[T[i]]=-2; // les sommets à enlever doivent être à -2
    bfs(G,s,X);
    ... // X->D[u]=distance entre s et u dans G\T
    ... // X->P[u]=père de u, ou -1 s'il n'existe pas
    ... // X->cycle=longueur du plus petit cycle dans G\T passant par s
    ... // X->n=nombre de sommets parcourus dans G
    ... //
    free_param_bfs(p);

  Pour faire des appels multiples à bfs() et optimiser le temps total:

    param_bfs *X=new_param_bfs();
    X->clean=1;  // initialisation complète puis partielle de X->D
    ...
    X->vmax=100; // pour parcourir au plus 100 sommets
    bfs(G,u,X);  // initialise complètement X->D à -1 avant le bfs
    ...
    X->hmax=3;   // pour faire un bfs à distance au plus 3 de v
    bfs(G,v,X);  // X->D sera initialisé en initialisant seulement les sommets
    ...          // parcourus au bfs() précédant (donc au plus 100 sommets)
    ...
    bfs(G,w,X);  // initialisation partielle de X->D
    ...
    free_param_bfs(X);

  La complexité est proportionnel au nombre d'arêtes dans la boule des
  sommets parcourus à condition que X->D ne soit pas NULL, puisque
  sinon une initialisation de X->D en O(n) sera effectuée. C'est un
  point important si on lance plusieurs BFS partiels à partir de
  sommets d'un même graphe. Pour être efficace il faut, à partir du 2e
  bfs, rétablir X->D[u]=-1 pour tous les sommets u de la X->file. On
  peut le réaliser en mettant X->clean=1. La complexité pour chaque
  appel (à part le premier qui initialise complètement X->D à -1)
  reste en le nombre d'arêtes de la boule des sommets parcourus.

  ALGORITHME:
  - prendre le sommet u en tête de file
  - pour tout ses voisins v non marqués:
    - enfiler v
    - marquer v
  
  Si D<>NULL, alors D n'est pas initialisé (sauf si clean>0). Dans ce
  cas, il faut que D[u]=-1 pour tout sommet u, sauf pour les sommets
  que l'on souhaite ne pas visiter où il faut D[u]<>-1 (typiquement,
  D[u]=-2). On peut initialiser partiellement ou complètement D avec
  clean=1 ou clean=2.

  Pour déterminer X->cycle, on calcule la longueur L du tout premier
  cycle crée. On remarque que X->cycle peut être L ou L-1.  Il est L-1
  si plus tard on rencontre deux sommets voisins sur ce même niveau.
*/
  const int n=G->n;
  if((source<0)||(source>=n)) Erreur(48); // au cas où
  int i,u,v,d,ff,tf,h;

  if(X==NULL) X=new_param_bfs(); /* NB: X->n=X->clean=0, X->P=X->file=NULL */
  if(X->P==NULL) ALLOC(X->P,n); /* alloue tableau si P==NULL */
  if(X->file==NULL) ALLOC(X->file,n); /* alloue la file */
  
  if(X->D==NULL) ALLOCZ(X->D,n,-1); /* alloue et initialise D en O(n) */
  else{ /* initialisation partielle ou pas de D */
    if(X->cont<2){
      if(X->clean==1) for(u=0;u<n;u++) X->D[u]=-1; /* initialisation complète, et les -2 ??? */
      if(X->clean==2) for(i=0;i<X->n;i++) X->D[X->file[i]]=-1; /* initialisation partielle */
    }
  }
  if(X->clean==1) X->clean=2; /* la prochaine fois, initialisation partielle de X->D */

  if(X->cont==2){
    tf=X->tf;
    ff=X->n;
  }else{
    tf=0; /* tf=tête de file, pointe sur la tête */
    ff=0; /* ff=fin de file, pointe sur prochain élément libre */
    if(X->cont==1) X->cont=2; /* la prochaine fois on continue */
    X->root=source; /* la racine est la source */
    X->P[source]=-1; /* pas de père pour la source */
    X->D[source]=0; /* distance=0 pour la source */
    X->file[ff++]=source; /* enfile le 1er sommet (=source), même s'il est supprimé */
    X->cycle=(G->sym)? n+1 : 0; /* X->cycle non défini (=0) si G orienté */
  }
  
  const int hmax=(X->hmax==-1)? n : X->hmax; /* si X->hmax non défini */
  const int vmax=(X->vmax==-1)? n : X->vmax; /* si X->vmax non défini */
  h=(G->sym)? 1+(X->cycle/2) : 0;  /* h=hauteur à partir de laquelle le plus court
				        cycle ne peut plus apparaître */

  while(tf<ff){
    u=X->file[tf]; /* défile la tête de file */
    if(X->D[u]>=hmax) break; /* fin si on a atteint une hauteur >=
				hmax. Dans ce cas, tous ceux de
				hauteur <= hmax on été enfilés. */
    tf++;
    for(i=0,d=G->d[u];i<d;i++){ /* pour tout voisin v de u */
      v=G->L[u][i];
      if(X->D[v]==-1){ /* si v voisin non marqué, si =-2 on saute le sommet */
	X->P[v]=u; /* le père de v est u */
	X->D[v]=X->D[u]+1; /* hauteur(v)=1+hauteur(père(v)) */
	X->file[ff++]=v; /* enfile v */
	if(ff>vmax){ tf=ff; i=d; } /* fin si parcouru vmax sommets */
      }else /* si v a déjà été visité ou s'il doit ne pas être visité (=-2) */
	if((X->D[u]<h)&&(v!=X->P[u])&&(X->D[v]!=-2)){ /* sinon X->cycle ne peut plus être améliorée */
	  h=X->D[u]+1; /* pas au delà de X->D[u] */
	  X->cycle=imin(X->cycle,h+X->D[v]);
	}
    }
  }

  if(X->cycle>n) X->cycle=-1; /* si > n, alors pas trouvé de cycle -> -1 */
  X->n=ff;  /* nb de sommets parcourus */
  X->tf=tf; /* nb de sommets parcourus dont les voisins ont été enfilés */
  X->radius=X->D[u]; /* hauteur du dernier sommet défilé (le plus éloigné) */
  
  /* c'est une mauvaise idée de faire ici un REALLOC(X->file,ff)
     car lors d'appels suivants avec la même structure, X->file n'aura
     plus la bonne taille ! */

  return X;
}


/*
  Structure de données pour le DFS
*/
typedef struct{
  int nc; // nc=nombre de composantes connexes du graphe
  int na; // na=nombre de sommets d'articulation (ou cut-vertex)
  int *C; // C[u]=couleur de la composante de u, entier de [0,nc[ ou sommet supprimé
  int *P; // P[u]=parent de u dans une forêt couvrante, P[racine]=-1
  int *R; // R[i]=i-ème sommet racine dans la forêt couvrante, i dans [0,nc[
  int *d; // d[u]=date de début de visite du sommet u, entier de [0,n[
  int *A; // A[u]=vrai ssi u est un sommet d'articulation, u dans [0,n[
  int *H; // H[u]=profondeur de u dans l'arborescence, H[racine]=0
} param_dfs;


param_dfs *new_param_dfs(const int n){
/*
  Attention ! le tableau X->C n'est pas alloué par cette fonction,
  même si n>0. Cela doit être fait par l'utilisateur.
*/

  NALLOC(param_dfs,X,1);

  X->C=X->P=X->R=X->d=X->A=X->H=NULL;
  X->nc=X->na=0;

  if(n>0){
    // X->C est alloué par l'utilisateur ou par dfs()
    ALLOC(X->P,n);
    ALLOC(X->R,n);
    ALLOC(X->A,n);
    ALLOC(X->d,n);
    ALLOC(X->H,n);
  }
  return X;
}


void free_param_dfs(param_dfs *X){
  if(X==NULL) return;
  free(X->C);
  free(X->P);
  free(X->R);
  free(X->A);
  free(X->d);
  free(X->H);
  free(X);
}


param_dfs *dfs(const graph *G,int source,param_dfs *X){
/*
  Effectue un parcours en profondeur de toutes les composantes
  connexes du graphe G depuis le sommet source. Version non
  récursive. On détermine également tous les sommets d'articulations
  (voir la définition de param_dfs pour lire le résultat), ainsi que
  les composantes connexes. On l'utilise comme suit:

  param_dfs *p=dfs(G,s,NULL); // DFS dans G depuis s
  ...
  free_param_dfs(p);

  ou alors, pour un DFS dans G évitant les sommets de T:

  param_dfs *p=new_param_dfs(G->n);
  ALLOCZ(p->C,G->n,-1);
  for(i=0;i<G->n;i++) p->C[T[i]]=-2;
  dfs(G,s,p);
  ...
  free_param_dfs(p);

  Si p->C[u]=-2, alors le sommet u n'est pas parcouru (il est
  virtuellement supprimé de G). Les autres sommets v (non supprimés)
  doivent avoir p->C[v]=-1. Si p->C==NULL (par défaut), alors ce
  tableau est alloué et initialisé à -1. Il sera libéré par
  free_param_dfs(p).
*/
  
  if(G==NULL) return NULL;
  const int n=G->n;
  if((source<0)||(source>=n)) Erreur(48); /* au cas où */

  int u,i,d,v,k,t,r=0;
  int tete,nc,na,b;

  if(X==NULL){ r=1; X=new_param_dfs(n); }
  if(X->C==NULL) ALLOCZ(X->C,n,-1);
  for(i=0;i<n;i++) X->A[i]=0;

  nc=na=0;
  NALLOC(int,pile,n);  /* pile */
  NALLOC(int,next,n);  /* next[u]=prochain voisin de u à visiter */
  NALLOC(int,level,n); /* level[u]=... */
  t=tete=-1;

  for(k=0;k<n;k++,source=(source+1)%n)
    /* on essaye tous les sommets à partir de source */
    if(X->C[source]==-1){ /* si ==-2 ou >=0 alors on saute le sommet */
      pile[++tete]=source;
      next[source]=0; /* premier voisin à visiter */
      X->P[source]=-1;
      X->H[source]=0;
      X->R[nc]=source;

      while(tete>=0){ /* tant que la pile n'est pas vide */
	u=pile[tete]; /* u=sommet courant */
	i=next[u]; /* indice du prochain voisin de u à visiter */
	if(i==0){
	  X->C[u]=nc; /* couleur de la composante courante */
	  level[u]=X->d[u]=++t; /* date de début de visite */
	}
	d=G->d[u]; /* degré de u */
	b=1; /* sentiennelle pour savoir si on a trouvé un v */
	while(i<d){ /* on cherche le prochain voisin v de u non visité */
	  v=G->L[u][i++]; /* pour tous les voisins v de u */
	  if(X->C[v]==-1){ /* si v n'a jamais été visité */
	    if((u==source)&&(t>X->d[u])&&(!X->A[u])) na++,X->A[u]=1; /* u=cut-vertex */
	    X->P[v]=u; /* père(v)=u */
	    X->H[v]=X->H[u]+1; /* mise à jour de la profondeur */
	    pile[++tete]=v; /* on empile v */
	    next[v]=b=0; /* le prochain voisin de v à visiter */
	    next[u]=i; /* le prochain voisin de u à visiter */
	    break;
	  } else /* v existe et a déjà été visité */
	    if((X->C[v]>=0)&&(v!=X->P[u])) /* si (u,v) est un arc de retour */
	      level[u]=imin(level[u],X->d[v]);
	}
	if(b){ --tete; /* il n'y a plus de voisin v: on dépile u pour toujours */
	  if((v=(X->P[u]))>=0){ /* si u n'est pas une racine, il a un père v */
	    level[v]=imin(level[v],level[u]); /* met à jour level du père de u */
	    if((v!=source)&&(level[u]>=X->d[v])&&(!X->A[v])) na++,X->A[v]=1; /* v=cut-vertex */
	  }
	}
      } /* fin du while(i<d) */

      nc++; /* ici nc=nb de composantes visitées */
    }

  X->nc=nc;
  X->na=na;
  free(pile);
  free(next);
  free(level);

  /* on réduit le tableau ->R au minimum que si dfs() l'a
     alloué. Sinon, il ne faut pas toucher aux pointeurs de X */
  if(r) REALLOC(X->R,nc);

  return X;
}


typedef struct{
  int n;        // nombre de sommets du graphe
  int source;   // source
  int *parent;  // tableau des parents: parent[u], -1 si pas de père
  double *dist; // tableau de distance: dist[u]

  // variables internes pour Bellman_Ford()
  int *stack1;
  int *stack2;
  int *vec1;
  int *vec2;
} param_bellman;


param_bellman *new_param_bellman(const int n){
  NALLOC(param_bellman,p,1);
  p->n=n;
  p->source=-1;
  p->parent=p->stack1=p->stack2=p->vec1=p->vec2=NULL;
  p->dist=NULL;
  if(n>0){
    ALLOC(p->parent,n);
    ALLOC(p->dist,n);
    ALLOC(p->stack1,n);
    ALLOC(p->stack2,n);
    ALLOC(p->vec1,n);
    ALLOC(p->vec2,n);
  }
  return p;
}


void free_param_bellman(param_bellman *p){
  if(p){
    free(p->parent);
    free(p->dist);
    free(p->stack1);
    free(p->stack2);
    free(p->vec1);
    free(p->vec2);
    free(p);
  }
  return;
}


param_bellman *Bellman_Ford(graph *G,const int source,param_bellman *p){
/*
  Calcule les distances du sommet source à tous les autres de G selon
  l'algorithme de Bellman-Ford. On retourne un tableau ->parent et
  ->dist. Il permet d'avoir des poids négatifs sur les arcs, cependant
  il ne doit pas avoir de cycle absorbant. Renvoie p, et si p=NULL on
  l'alloue. On utilise en interne 4 tableaux de taille n, inclus dans
  la variable p.
  
  A faire: gérer les sommets éteints.
*/  
  if(G==NULL) return NULL;
  if(G->W==NULL) return NULL;
  const int n=G->n;
  if((source<0)||(source>=n)) return NULL;
  if(p==NULL) p=new_param_bellman(n);

  int u,v;
  double w;
  int i,i1,i2,d;
  int *tmp;

  /* initialisations */
  for(u=0;u<n;u++){
    p->dist[u]=DBL_MAX; /* +∞ */
    p->parent[u]=-1; /* pères non définis, y compris pour la source */
    p->vec1[u]=p->vec2[u]=1; /* vec1[u]=0 ssi u est dans p->stack1, idem pour vec2 */
  }
  p->dist[source]=0; /* distance nulle pour la source */
  i1=i2=0;
  p->stack1[i1++]=p->source=source; /* empile la source dans stack1 */
  p->vec1[source]=0; /* source est dans la stack1 */

  while(i1>0){
    while(i1>0){
      u=p->stack1[--i1]; /* dépile u */
      p->vec1[u]=1; /* u n'est plus dans la stack1 */
      d=G->d[u]; /* d=degré de u */
      for(i=0;i<d;i++){ /* pour tous les voisins v de u */
	v=G->L[u][i];// if((p->inc)&&(v<u)) continue;
	w=p->dist[u]+G->W[u][i];
	if(w<p->dist[v]){
	  p->parent[v]=u;
	  p->dist[v]=w;
	  if(p->vec2[v]){ /* empiler v dans stack2 s'il n'y est pas */
	    p->stack2[i2++]=v;
	    p->vec2[v]=0;
	  }
	}
      }
    }
    /* échange stack1 et stack2 */
    SWAP(i1,i2,i);
    SWAP(p->stack1,p->stack2,tmp);
    SWAP(p->vec1,p->vec2,tmp);
  }

  return p;
}


/*
  Dijkstra(G,W,s):

  1. Init(G,s)
  2. S={}
  3. Q=V(G)
  4. while Q<>{}
  5.  do u=Extract-Min(Q)
  6.     S=S u {u}
  7.     for each v in Adj[u]: Relax(u,v,W)

  Init(G,s):
  1. for each vertex v in V(G):
  2.   do d[v]=+∞
  3.      pi[v]=NULL
  4. d[s]=0

  Relax(u,v,W):
  1. if d[v]>d[u]+W[u,v]
  2.   then d[v]=d[u]+W[u,v]
  3.        pi[v]=u

  Extract-Min(Q):
  extract vertex u with smallest d[v]
*/


enum{
  IW_UNW,  // graphe valué uniformément
  IW_GEO,  // graphe euclidien
  IW_POS,  // IW_UNW ou IW_GEO suivant la valeur de POS
  IW_R01,  // valuation aléatoires dans [0,1]
};


int InitWeights(graph *G,int code){
/*
  Initialise le poids des arcs du graphe G, c'est-à-dire le tableau
  G->W, qui est éventuellement alloué. L'initialisation dépend de code
  (cf. l'enum IW_xxx). Retourne 1 ssi tout c'est bien passé, 0 sinon.
*/
  if(G==NULL) return 0;
  if(code==IW_POS) code=(POS)? IW_GEO : IW_UNW;
  if((code==IW_GEO)&&(G->xpos==NULL)) return 0;

  const int n=G->n;
  int u,v,i,d;
  if(G->W==NULL) ALLOCZ(G->W,n,NULL);

  for(u=0;u<n;u++){
    d=G->d[u];
    free(G->W[u]);
    ALLOC(G->W[u],d);
    for(i=0;i<d;i++){
      v=G->L[u][i];
      switch(code){
      case IW_GEO: G->W[u][i]=hypot(G->xpos[u]-G->xpos[v],G->ypos[u]-G->ypos[v]); break;
      case IW_UNW: G->W[u][i]=1; break;
      case IW_R01: G->W[u][i]=RAND01; break;
      default: return 0; // problème
      }
    }
  }
  
  return 1;
}

/***********************************

       ISOMORPHISM, SUBGRAPH,
       MINOR, INDUCEDSUBGRAPH,
       PATHS, PS1, TREEWIDTH ...

***********************************/


int *Isomorphism(graph *G,graph *H){
/*
  Renvoie un tableau P <> NULL ssi G est isomorphe à H. Si tel est le
  cas, le tableau P indique le morphisme de H vers G. Après l'appel,
  le graphe H est modifié: ses listes sont triées si H->sort=0 est
  faux (sauf si G=H - même pointeur). Le graphe G n'est par contre pas
  modifié. Dans H->int1 est retourné le nombre de tests effectués.
  Moins le graphe possède de symétrie, plus faible est le nombre de
  tests (et rapide est la fonction).

  On applique l'algorithme suivant. Pour chacun des deux graphes et
  chaque sommet u, on calcule son "profile" un tableau noté
  profile[u]: profile[u][i+2] = le nombre de sommets à distance i de
  u, en commençant à partir de i=1. Donc, profile[u][3] est le degré
  de u. Ceci est calculé par un simple BFS, les indices 0,1,2 étant
  réservés. On utilise profile[u][0] pour coder la taille du tableau
  profile[u], profile[u][1] pour stocker le nom de u (plutôt que le
  nombre de sommet à distance 0 qui est toujours 1) et profile[u][2]
  pour stocker la longueur du plus petit cycle passant par u. Cette
  dernière information étant calculée "gratuitement" par le BFS.

  On ordonne ensuite les sommets des graphes suivant les profiles des
  sommets avec qsort(). On ne renumérote pas les sommets dans le
  graphe, mais plutôt on donne un ordre: c'est possible avec qsort()
  car le nom original u est dans profile[u][1] le sommet u. Si bien
  que profile[i][1] sera le sommet u de rang i. On priviligie les
  profiles le grande taille (que l'on classe en premier) ce qui est
  plus discriminant. Les isomorphismes (=permutations) à tester ne
  concernent que les sommets de même profile. On construit les
  contraintes dans un tableau C, C[j] indiquant que les sommets de
  rang C[j-1] à C[j] (exclu) ont même profile, et sur lequel se base
  NextPermutation().

  Sur le même principe, on pourrait imaginer un profile plus complexe,
  comme suit: à chaque distance i et sommet u, on calcule le graphe
  G[u][i] induit par les sommets à distance i de u. On peut alors
  calculer le profile de chaque sommet de G[u][i] et ordonner les
  sommets selon celui-ci.
*/
  
  H->int1=0; /* par défaut, 0 tests */
  if((G==NULL)||(H==NULL)) return NULL;
  if(G->n!=H->n) return NULL;

  int *P; /* isomorphisme final */
  const int n=G->n;

  if(G==H){ /* isomorphisme trivial si même emplacement mémoire */
    ALLOCZ(P,n,_i);
    return P;
  }

  param_bfs *param=new_param_bfs(); /* pour le BFS */
  int **profile,**profileG=NULL,**profileH;
  int *R,*C; /* permutation et contraintes (sur les rangs) */
  int u,v,t,r,i;
  graph *M;

  for(M=G;;){ /* on fait la même chose pour M=G puis M=H */
    ALLOC(profile,n); /* profile[u] */
    for(u=0;u<n;u++){ /* faire un BFS pour tous les sommets u */
      bfs(M,u,param); /* le premier BFS va allouer param->D et param->P */
      t=3+param->radius; /* taille du tableau profile[u] */
      ALLOCZ(profile[u],t,0); /* initialise le profile */
      for(v=0;v<n;v++){
	i=param->D[v];
	if(i>0) profile[u][i+2]++; /* compte les sommets à distance i>0 de v */
	param->D[v]=-1; /* réinitialise les distances pour les BFS suivants */
      }
      profile[u][0]=t; /* taille du profile */
      profile[u][1]=u; /* nom du sommet, pour qsort() */
      profile[u][2]=param->cycle; /* maille */
    }
    QSORT(profile,n,fcmp_profile); /* trie les profiles */

    if(M==H){ profileH=profile; break; } /* on s'arête si M=H */
    profileG=profile; /* on refait la boucle pour H */
    M=H;
  }
  free_param_bfs(param);

  /* on verifie que profileG "=" profileH */
  for(u=0;u<n;u++)
    if(fcmp_profile(profileG+u,profileH+u)){
      P=NULL;
      goto fin_noniso;
    }

  /* calcule les contraintes */
  /* ici les profiles de G et H sont identiques */
  ALLOC(C,n);
  R=profile[0]; /* R=profile du premier sommet. NB: profile=profileH. */
  for(u=t=0;u<n;u++){
    if(fcmp_profile(&R,profile+u)){ /* si profiles différent */
      R=profile[u];
      C[t++]=u;
    }
  }
  C[t]=n;

  ALLOC(P,n);
  ALLOCZ(R,n,_i); /* initialise l'isomorphisme sur les rangs des profiles */
  if(!H->sort) SortGraph(H,0); /* on trie H pour la recherche dichotomique */
  H->int1=0; /* compte le nb de tests */

  /* vérifie, pour chaque permutation P, que P(G)=H */

  do{
    H->int1++;

    /* calcule P en fonction de R */
    for(r=0;r<n;r++) P[profileG[r][1]]=profileH[R[r]][1];
    for(r=0;r<n;r++){ /* on commence par les profiles longs: r=0,1,2, ... */
      u=profileG[r][1]; /* v=P[u] est le sommet correspondant à u dans H */
      for(i=0,t=G->d[u];i<t;i++) /* on regarde si chaque voisin de u est dans H->L[v] */
	if(bsearch(P+(G->L[u][i]),H->L[P[u]],t,sizeof(int),fcmp_int)==NULL){
	  /* alors élément non trouvé */
	  i=t;r=n; /* prochaine permutation à tester */
	}
    } /* ici r=n (trouvé) ou n+1 (pas encore trouvé) */
    if(r==n) goto fin_iso; /* on a trouvé un isomorphisme P */
  }
  while(NextPermutation(R,n,C));

  /* si on arrive ici, c'est qu'on a pas trouvé l'isomorphisme */
  free(P);
  P=NULL;

 fin_iso:
  free(R);
  free(C);

 fin_noniso:
  FREE2(profileG,n);
  FREE2(profileH,n);

  return P;
}


edge *ListEdges(graph *G){
/*
  Construit la liste des arêtes de G, chaque arête uv ne figure qu'une
  seule fois. On a, pour tout i, E[i].u < E[i].v. Le champs G->m est
  mis à jour.
*/
  const int m=NbEdges(G);
  const int n=G->n;
  int u,v,i,j,d;

  NALLOC(edge,E,m);
  for(u=j=0;u<n;u++){
    for(i=0,d=G->d[u];i<d;i++){
      v=G->L[u][i];
      if(u<v) E[j].u=u,E[j++].v=v;
    }
  }
  return E;
}


graph *Subgraph(graph *G,graph *H){
  /*
    Détermine si G est un sous-graphe de H s'ils ont même nombre de
    sommets. On renvoie un sous-graphe S de H isomorphe à G, et on
    renvoie NULL si H ne possède pas G comme sous-graphe ou si G et H
    n'ont pas le même nombre de sommets.

    Effets de bord:
    - les listes de G sont triées et donc G->sort=1,
    - H->int1 contient le nombre total de tests effectués,
    - S->pint1 contient l'isomorphisme de S vers G, et donc de H vers G.

    L'algorithme est le suivant: on teste d'abord si la séquence des
    degrés de H est bien supérieure à celle de G (ceci prend un temps
    O(n)). Si c'est le cas, on effectue, pour tous les sous-graphes S
    de H qui ont autant d'arêtes que G, un test d'isomorphisme entre S
    et G grâce à isomorphisme(S,G).
  */
  const int n=H->n;
  H->int1=0; /* nb de tests = 0 par défaut */
  if(n!=G->n) return NULL; /* pas le même nombre de sommets */

  /* on trie en O(n) les deux listes */
  int* const Eh=SortInt(H->d,NULL,n,0,NULL,SORT_INC);
  int* const Eg=SortInt(G->d,NULL,n,0,NULL,SORT_INC);
  int i;

  /* on s'arrête si, pour un rang i donné, degH(i)<degG(i) */
  for(i=0;i<n;i++) if(Eh[i]<Eg[i]) break;
  free(Eh);
  free(Eg);
  if(i<n) return NULL; /* G ne peut pas être un sous-graphe de H */

  int mg=NbEdges(G);
  int mh=NbEdges(H);
  graph *S=new_subgraph(H); /* S=sous-graphe de H, alloue S->L et dimensionne S->d */
  edge* const E=ListEdges(H); /* liste des arêtes de H: e_j=(u,v) -> E[j].u=u et E[j].v=v */
  int *B; /* les arêtes de S, sous-ensemble d'indices d'arêtes de H */
  int *P; /* isomorphisme S->G */
  int u,v,j,d;
  
  /* initialise le sous-ensemble d'arêtes B de H avec mg arêtes */
  ALLOC(B,mg);
  NextSet(B,-1,mg);
  d=0; /* d=compteur pour le nombre total de tests */
  
  do{
    
    /* remplit S avec les arêtes de B */
    degres_zero(S); /* position libre pour le sommet u de S */
    for(i=0;i<mg;i++){ j=B[i]; /* j=numéro de la i-ème arête de B */
      u=E[j].u; v=E[j].v; /* l'arête j de E est (u,v) */
      ADD_EDGE(S,u,v); /* ajoute l'arête u,v */
    }
    
    /* il vaut mieux que G soit le 2e paramètre, car il va être trié
       la première fois par Isomorphism(), puis plus jamais grâce au
       test de G->sort, alors que S serait trié à chaque appel */
    P=Isomorphism(S,G);
    d += 1+G->int1; /* on ajoute 1 pour donner le nombre d'ensembles testés */
  }while((P==NULL)&&(NextSet(B,mh,mg)));
  
  H->int1=d; /* nombre total de tests */

  if(P==NULL){ /* on a pas trouvé de sous-graphe de H isomorphe à G */
    free_graph(S);
    S=NULL;
  }
  else S->pint1=P; /* S isomorphe à G, sous-graphe de H */

  free(E);
  free(B);
  return S;
}


graph *MatrixToGraph(int **M,const int n){
/*
  Renvoie le graphe correspondant à la matrice d'adjacence n x n où
  seule la partie inférieure est utilisée. Les listes du graphe de
  retour sont triées et les champs ->m et ->sort sont mise à jour.
*/
  if((M==NULL)||(n<=0)) return NULL;
  int u,v,m;
  graph *G=new_fullgraph(n);

  for(u=m=0;u<n;u++)
    for(v=0;v<u;v++)
      if(M[u][v]){
	ADD_EDGE(G,u,v); /* ajoute u-v et v-u */
	m++; /* une arête de plus */
      }
  
  /* réduit les listes */
  GraphRealloc(G,G->d);

  G->m=m;
  G->sort=1;
  return G;
}


graph *GraphOfColor(graph *G,int *col,const int k){
/*
  Renvoie un graphe C dont les sommets sont les entiers de [0,k[ (les
  valeurs du tableau col[] qui doit avoir une taille G->n) et les
  arêtes les paires uv telle qu'il existe une arête xy de G avec
  col[x]=u et col[y]=v. La valeur C->m est mise à jour, et les listes
  de C sont triées (C->sort=1).
*/
  if((k<0)||(col==NULL)||(G==NULL)) return NULL;

  const int n=G->n;
  int u,v,cu,cv,i,d;
  graph *C; /* le graphe des couleurs renvoyé */

  /* matrice d'adjacence inférieure à 0 */
  NALLOCMAT(int,M,k,k-1); /* matrice d'adjacence du graphe des couleurs */
  for(u=0;u<k;u++)
    for(v=0;v<u;M[u][v++]=0);

  for(u=0;u<n;u++) /* parcourt G et remplit la partie inférieure de M */
    for(i=0,d=G->d[u];i<d;i++){
      v=G->L[u][i];
      if(u<v){ /* si cu=cv, on ne fait rien */
	cu=col[u];
	cv=col[v];
	if(cu>cv) M[cu][cv]=1;
	if(cv>cu) M[cv][cu]=1;
      }
    }
  
  C=MatrixToGraph(M,k);
  FREE2(M,k);
  return C;
}


// Routine de fusion pour UNION-FIND
#define UNION(x,y,c,r)				\
  do{			       			\
    if(r[x]>r[y]) c[y]=x;	       		\
    else{ c[x]=y; if(r[x]==r[y]) r[y]++; }     	\
  }while(0)


int FindSet(const int x,int *p){
/* Routine FIND pour le problème d'UNION-FIND. Donne le représentant
   de x dans le tableau p ou encore la racine de x dans la forêt
   couvrante p. Le temps total de m requêtes à FindSet(), pour un
   tableau p de taille n est O(m*𝛼(n)) où 𝛼(n) est la fonction inverse
   d'Ackerman, 𝛼(n)<=4 en pratique lorsque que pour la fusion on
   utilise les heuristiques du rang (réalisée par UNION) et de la
   compression de chemin (réalisée par FIND).

   Ex: on utilise typiquement FindSet() pour savoir si, lorsqu'on
   ajoute une arête à un graphe, on crée un cycle ou pas. Si on ne
   crée par de cycle, il faut appeler Union(x,y,rang,couleur)

   NALLOCZ(int,couleur,n,_i);
   NALLOCZ(int,rang,n,0);

   for(u=...)
     for(v=...){ // pour chaque arête {u,v}
       couleur[u]=x=FindSet(u,couleur); // x=représentant de u
       couleur[v]=y=FindSet(v,couleur); // y=représentant de v
       if(x==y){ // il y a un cycle
       ...;
       }else{ // il n'y a pas de cycle
       ...;
       // fusionne la composante de u et de v
       UNION(x,y,couleur,rang);
       }
    }
*/
  if(x!=p[x]) p[x]=FindSet(p[x],p);
  return p[x];
}


int *Minor(graph *H,graph *G)
/*
  Détermine si H est un mineur de G. Si c'est le cas, un tableau T est
  renvoyé, sinon NULL est renvoyé. Le tableau T code un modèle du
  mineur de H dans G. Plus précisément, T[u]=c, où u est un sommet de
  G, est le nom du sommet c de H si bien que l'ensemble des sommets u
  de G tel que T[u]=c forme le super-noeud c.

  L'algorithme est le suivant: on effectue, pour tous les ensembles de
  contractions d'arêtes de G produisant un mineur avec autant de
  sommets que H, un test de sous-graphe grâce à Subgraph().
*/
{
  graph *C; /* graphe des couleurs = graphe contracté */
  graph *S=NULL; /* sous-graphe éventuellement isomorphe à C */
  int *B; /* sous-ensemble (d'indices) d'arêtes de G à contracter */
  int *couleur; /* les sommets de même couleurs seront contractés */
  int *rang; /* pour union-find rapide */
  edge e;
  
  int nh=H->n;
  int ng=G->n;
  int c=ng-nh; /* c=nb de contractions à effectuer */
  int t;
  H->int1=t=0; /* initialise le nb de tests */
  if(c<0) return NULL; /* pas de mineur */

  edge *E=ListEdges(G); /* E=liste des arêtes de G, met à jour G->m */
  int mg=G->m;
  if(c>mg){ /* fini s'il faut contracter plus d'arêtes qu'il y a dans G */
    free(E);
    return NULL;
  }

  int i,u,v,x,y;
  int test=((c<<1)<ng); /* vrai ssi on liste suivant les arêtes ou suivant les sommets */
  ALLOC(B,c); NextSet(B,-1,c); /* B=premier sous-ensemble de c arêtes de E */
  ALLOC(couleur,ng); /* couleur des sommets */
  ALLOC(rang,ng); /* rang des sommets, sert pour UNION-FIND */

  /*
    On pourrait générer des sous-ensembles acycliques d'arêtes
    directement, en combinant NextSet() et le test d'acyclicité.
   */

  do{

    t++; /* on compte 1 test par ensemble B testés */

    /* initialise la couleur et le rang des sommets */
    for(u=0;u<ng;u++){ couleur[u]=u; rang[u]=0; }

    /* on teste rapidement (avec UNION-FIND) si B contient un cycle */
    for(i=0;i<c;i++){ e=E[B[i]]; /* e=i-ème arête de B */
      u=e.u; couleur[u]=x=FindSet(u,couleur);
      v=e.v; couleur[v]=y=FindSet(v,couleur);
      if(x==y) break; /* y'a un cycle, sinon on fait UNION */
      UNION(x,y,couleur,rang);
    }

    if(i==c){ /* si B est acyclique, on fait un test de sous-graphe */
      if(test)
	/* on met à jour la couleur de chaque sommet. Suivant les
	   valeurs respectives de c et ng (test) on met à jour soit
	   suivant les arêtes ou suivant les sommets. */
	for(i=0;i<c;i++){
	  e=E[B[i]]; /* e=i-ème arête de B */
	  u=e.u; couleur[u]=FindSet(u,couleur);
	  v=e.v; couleur[v]=FindSet(v,couleur);
	}
      else
	for(u=0;u<ng;u++) couleur[u]=FindSet(u,couleur);

      /* on recadre les couleurs dans [0,c[. Complexité: 4ng */
      for(i=0;i<ng;rang[i++]=0);
      for(i=0;i<ng;i++) rang[couleur[i]]++; /* rang=fréquence */
      for(i=u=0;i<ng;i++) /* repère les couleurs manquantes */ 
	if(rang[i]==0) u++; else rang[i]=u;
      /* ici rang[i]=nb de zéros (=valeurs manquantes) dans rang[0..i[ */
      for(i=0;i<ng;i++) couleur[i] -= rang[couleur[i]];

      C=GraphOfColor(G,couleur,nh);
      S=Subgraph(H,C); /* avant l'appel, S=NULL nécessairement */
      t += C->int1;
      free_graph(C); /* on a plus besoin du graphe des couleurs */
    }

  }while((S==NULL)&&(NextSet(B,mg,c)));

  H->int1=t;
  free(B);
  free(E);

  /* on a rien trouvé */
  if(S==NULL){
    free(couleur);
    free(rang);
    return NULL;
  }
  
  /* on a trouvé un mineur, on construit le modèle dans rang[] */
  for(u=0;u<ng;u++) rang[u]=S->pint1[couleur[u]];
  free_graph(S);
  free(couleur);
  return rang;
}


int *InducedSubgraph(graph *H,graph *G)
/*
  Indique si H est un sous-graphe induit de G.  La fonction renvoie un
  ensemble X de sommets tel que G[X] est ismomorphe à H. Evidemment
  |X|=H->n. On renvoie dans G->int1 le nombre de tests effectués, et
  dans G->pint1 l'isomorphisme entre H et G[X].

  L'algorithme consiste à générer tous les ensembles X possibles de
  |V(H)| sommets et à tester l'isomorphisme entre G[X] et H.
 */
{
  if((G==NULL)||(H==NULL)) return NULL;
  int ng=G->n,nh=H->n;
  if(nh>ng) return NULL;

  graph *S;
  int *P,t=0;
  NALLOC(int,X,nh);
  NextSet(X,-1,nh); /* premier sous-ensemble */

  do{
    t++;
    S=ExtractSubgraph(G,X,nh,1);
    P=Isomorphism(S,H);
    t += H->int1;
    free_graph(S);
  }while((P==NULL)&&(NextSet(X,ng,nh)));

  G->int1=t;
  free(G->pint1); /* pour éviter les fuites mémoires */
  G->pint1=P;
  if(P==NULL){ free(X); X=NULL; }
  return X;
}


int NextPath(graph *G,path *P,int j){
/*
  Cette fonction (récursive) permet de trouver tous les chemins
  simples entre deux sommets. Plus précisément, on met à jour le
  chemin P de sorte qu'on trouve un nouveau chemin allant du j-ème
  sommet de P au dernier tout en évitant la première partie du chemin
  allant du premier sommet de P au j-ème. Si un tel chemin a pu être
  trouvé on renvoie 1, sinon on renvoie 0. On met j=-1 s'il s'agit de
  la création d'un chemin de P->P[0] à P->P[1].

  On l'utilise comme ceci:

    path *P=new_path(G,NULL,G->n); // crée un chemin vide P mais avec G->n sommets possibles
    P->P[0]=x; // origine du chemin
    P->P[1]=y; // destination du chemin, y=x possible
    if(NextPath(G,P,-1)==0){ ... } // crée le premier chemin, ==0 s'il n'y en a pas
    do{ // traitement du chemin P
      ...
    }while(NextPath(G,P,0); // tantqu'il existe un nouveau chemin P
    free_path(P); // libère le chemin P

  Plus précisément, étant donnés un chemin P=x-...-v_j-...-y du graphe
  G et un sommet v_j du chemin (v_j=j-ème sommet du chemin), la
  fonction tente de compléter P par un autre chemin de v_j à y évitant
  x-...-v_(j-1). Si ce chemin à été trouvé, alors la partie de v_j à y
  de P est mise à jour et on renvoie 1. Sinon, le chemin est coupé
  après v_j et on renvoie 0. Dans tous les cas P est un chemin à jour
  de G. Si j<0, alors on initialise P par un chemin allant de
  x=P->P[0] à y=P->P[1].

  Algorithme: on essaye d'améliorer en premier le sous-chemin de
  v_{j+1} à y (récursivement). Si cela n'est pas possible, on calcule
  un nouveau chemin de v_j à y passant par un autre voisin v de v_j
  (autre que v_{j+1}) et évitant le chemin x-...-v_j. On passe en
  revue ainsi tous les voisins v de v_j qui ne sont pas dans
  P[0]...P[j]. Si aucun des voisins ne possède un tel chemin, c'est
  qu'il y en a plus et on retourne 0.

  Comme il faut tester les voisins v de v_j qu'une seule fois (pour un
  chemin P[0]...P[j] donné), on utilise le champs aux[v_j][i] de P qui
  donne le i-ème voisin libre de v_j avec la convention que
  aux[v_j][0] est le nombre de voisins encore possibles.

  Effet de bord: P est mis à jour.
*/
  if((P==NULL)||(G==NULL)) Erreur(-1); /* ne devrait jamais arriver */
  param_bfs *p;
  int i,x,y,u,v,n;

  if(j<0){ /* initialisation du premier chemin */
    n=G->n;
    if(P->aux==NULL) ALLOCMAT(P->aux,n,n);
    for(u=0;u<n;P->V[u++]=-1); /* vide le chemin */
    x=P->P[0];
    y=P->P[1];
    if((x<0)||(y<0)||(x>=n)||(y>=n)){ P->n=0; p=NULL; goto fin_0; } /* sommets inexistant */
    p=bfs(G,x,NULL); /* calcule le chemin */
    i=p->D[y];
    if(i<0){
      P->V[x]=0; /* x est en première position dans P */
      P->n=1;
      goto fin_0;
    }
    /* on initialise aux[x] et aux[y] qu'une seule fois */
    P->aux[y][0]=0;
    j=-1; /* pour que la longueur soit correcte */
    goto fin_1;
  }

  n=P->n;

  if(j+1==n) return 0; /* si x=y, alors pas de prochain chemin */
  if(NextPath(G,P,j+1)) return 1; /* c'est le NextPath à partir du voisin de v_j */

  /* Ici on ne peut pas augmenter v_(j+1)...y. Donc, il faut calculer
     un premier chemin depuis le prochain voisin v disponible de u=v_j
     et y qui évite x=P[0]-...-P[j]. Si cela ne marche pas avec le
     voisin v, il faut essayer le suivant, etc. */

  /* efface depuis P la fin du chemin P[j+1]...P[n-1] */
  /* pour ne garder que P[0]...P[j] */
  for(i=j+1;i<n;i++) P->V[P->P[i]]=-1;

  /* rem: ici t<d */
  p=new_param_bfs();
  ALLOC(p->D,G->n);
  y=P->P[n-1];
  u=P->P[j];
  i=-1;

  while((P->aux[u][0])&&(i<0)){ /* tant qu'il y a encore des voisins de u non testés */
    v=P->aux[u][(P->aux[u][0])--]; /* lit et enlève le dernier voisin dispo */
    if(P->V[v]>=0) continue; /* si le voisin est dans P[0] ... P[j] on le saute */
  
    /* initialise p->D: on doit le faire à chaque appel */
    for(i=0;i<G->n;i++) p->D[i]=-1;
    /* on enlève P[0]...P[j] du graphe pour le BFS */
    for(i=0;i<=j;i++) p->D[P->P[i]]=-2;

    /* calcule un chemin de v à y dans G\(P[0]-...-P[j]) */
    bfs(G,v,p);

    /* a-t-on trouvé un chemin ? */
    i=p->D[y]; /* si i>=0, c'est oui */
  }

  /* ici i=distance de u=v_j à y */
  if(i<0){ /* on a pas trouvé de chemin */
  fin_0:
    free_param_bfs(p);
    return 0;
  }

  /* on a trouvé un chemin, on met à jour P */
 fin_1:
  P->n=i+j+2; /* nouvelle longueur */
  /* ajoute à P[0]...P[j] le chemin de P[j+1] à y en partant de y */
  while(i>=0){
    P->P[i+j+1]=y;
    P->V[y]=i+j+1;
    y=p->P[y];
    i--;
  }

  /* initialise aux[u] pour tous les sommets u de P[j+1] à y non
     compris. Pour P[j] on le fait au moment de la lecture de v dans
     aux[u], et pour y on le fait une seule fois à la création
     (avec aux[y][0]=0) */

  for(i=j+1,n=P->n-1;i<n;i++){ /* si j=-1, alors on fait pour x aussi */
    u=P->P[i];
    P->aux[u][0]=0;
    for(v=0;v<G->d[u];v++){
      j=G->L[u][v]; /* j=voisin de u */
      /* si pas dans P on l'ajoute comme voisin possible */
      if(P->V[j]<0) P->aux[u][++P->aux[u][0]]=j;
    }
  }

  free_param_bfs(p);
  return 1;
}


#define LCONF 9  /* nb de bits pour coder un sommet dans le graphe des conflits */
#define CONFMAX (1<<LCONF) /* = 2^LCONF = nb max de sommets du graphe des conflits */
#define CONFMASK (CONFMAX-1) /* = CONFMAX-1 = mask pour un sommet */
#define CONFC 2 /* constante pour modifier le code des noeuds à 0. Il faut CONFC>1 */
#define NCMAX 512 /* taille maximum de la liste de composantes maximales. NCMAX<CONFMAX */

/* ensemble de variables pour gérer le graphe des conflits */
typedef struct{
  int n;   /* nb de sommets du graphe d'origine, sert pour la règle 5 */
  int nbi; /* nb de noeuds avec code=-1 (indéterminée) */
  int nbzi;/* nb de noeuds de code 0 indépendants */
  int outmem; /* 1 ssi il y a eut un dépassement de CONFMAX ou NCMAX */
  int paire[CONFMAX]; /* 1er noeud de la paire xy */
  int path[CONFMAX];  /* 1er noeud du chemin P, ne sert que pour PrintConflit() */
  int code[CONFMAX];  /* code du noeud i. Valeurs possibles: -1,0,1 */
  int nbc[CONFMAX];   /* nb de noeuds de la paire i dont le code est <1 */
  int *comp[CONFMAX];  /* liste des noeuds de la composante de i */
  int *cmax[NCMAX+1]; /* liste des composantes maximales */
  int tmax[NCMAX+1]; /* taille des composantes maximales */
  int ncmax; /* taille de la liste cmax */
  int nodemax; /* noeud correspondant à la composante maximale de la liste cmax */
  int valmax; /* valeur (VRAI/FAUX) des composantes maximales */
  int tcomp[CONFMAX]; /* taille de comp[i] */
  int x[CONFMAX],y[CONFMAX]; /* paire du noeud i, ne sert que pour PrintConflit() */
  graph *G; /* graphe des conflits */
  /* Attention! les noeuds des listes d'adjacence de G peuvent être >
     CONFMAX et donc aussi > G->n. On s'en sert pour coder un type
     d'arête. Donc si u est un noeud de G et v=G->L[u][i], alors le
     i-ème voisin de u est (v&CONFMASK) et le type de l'arête uv est
     (v>>LCONF). */
} conflit;


void PrintConflit(conflit *c)
/*
  Affiche le graphe des conflits, et rien si c->G->n=0.
*/
{
  if(c->G->n==0) return;

  int u,v,i,t,d;
  char T1[]="|=><";
  char T2[]="-01*";
  char T3[]="_/.";

  // UTF-8:
  // char T3[]="━┛o";
  // \xe2\x94\x80: ─
  // \xe2\x94\x81: ━
  // \xe2\x94\x98: ┘
  // \xe2\x94\x99: ┙
  // \xe2\x94\x9a: ┚
  // \xe2\x94\x9b: ┛

  printf("code (x,y) [comp] node-path-paire: node[type]\n");
  for(u=0;u<c->G->n;u++){
    t=c->code[u];
    if(c->nbzi>0){ /* si on a des zéros indépendants */
      if(t==0) t=2; /* c'est un zéro indépendant */
      else if(t==CONFC) t=0; /* c'est un zéro marqué -> remet sa valeur d'origine */
    }
    printf(" %c ",T2[t+1]);
    if(c->paire[u]==u)
      printf("(%i,%i) [",c->x[u],c->y[u]);
    else printf("\t [");
    for(v=0;v<c->tcomp[u];v++){
      printf("%i",c->comp[u][v]);
      if((c->n>9)&&(v<c->tcomp[u]-1)) printf(",");
    }
    printf("]\t%s%i",(v<5)?"\t":"",u);
    if(u>c->path[u]) printf("%c \t:",T3[1]); else{
      printf("%c",T3[0]);
      if(u>c->paire[u]) printf("%c\t:",T3[1]); else{
	printf("%c%c\t:",T3[0],T3[2]);
      }
    }
    for(i=0,d=imin(c->G->d[u],WIDTH);i<d;i++){
      v=c->G->L[u][i];
      t=(v>>LCONF);
      printf(" %i%c",v&CONFMASK,T1[t]);
    }
    if(c->G->d[u]>WIDTH) printf("...");
    printf("\n");
  }
  printf("#nodes in conflict graph: %i\n",c->G->n);
  printf("#heavy indep. components: %i\n",c->nbzi);
  printf("#unspecified values: %i\n",c->nbi);
  if(c->outmem){
    printf("!!! Out of Memory !!!\n");
    if(c->ncmax>=NCMAX)
      printf("#nodes in conflit graph exceeded (%i)\n",CONFMAX);
    else
      printf("#maximal components exceeded (%i)\n",NCMAX);
  }
  return;
}


void ps1_delxy(conflit *c,const int w)
/*
  Supprime du graphe des conflits c la dernière paire créee, et met
  ainsi à jour c. Ici w est l'indice du premier noeud de la dernière
  paire dans le graphe des conflits. Il faut supprimer tous les noeuds
  u >= w et tous les arcs entre les noeuds u et v<w.

  La liste des composantes maximales de la dernière paire sera effacée
  à la fin de la paire courante (voir le code après nextxy:).
*/
{
  int u,i,v,j;

  for(u=w;u<c->G->n;u++){
    for(i=0;i<c->G->d[u];i++){ /* pour tous les arcs v->u */
      v=c->G->L[u][i]&CONFMASK; /* v est le i-ème voisin de u */
      if(v<w) /* seulement pour les noeuds v avant w */
	for(j=0;j<c->G->d[v];j++)
	  /* on cherche et supprime u de la liste de v */
	  if((c->G->L[v][j]&CONFMASK)==u)
	    c->G->L[v][j]=c->G->L[v][--(c->G->d[v])]; /* supprime u */
    }
    free(c->G->L[u]); /* supprime la liste de u */
    free(c->comp[u]); /* supprime la composante de u */
  }

  c->G->n=w; /* met à jour le nombre de noeuds du graphe des conflits */
  return;
}


int ps1_addmax(int *C,const int t,conflit *c)
/*
  Ajoute une composante C (de taille t) à la liste des composantes
  maximales déjà rencontrées, et maintient cette propriété. Il faut
  que la liste (c->cmax) soit de taille suffisante pour accueillir une
  nouvelle composante. La fonction renvoie VRAI ssi la composante C a
  été ajoutée à la liste.

  Algorithme:
    1. Pour toute composante X de L (=la liste):
       1.1. Si C est inclue ou égale à X alors FIN
       1.2. Si C contient X alors supprimer X de L
    2. Ajouter C à L

  On remarque aussi qu'il n'est pas possible de rencontrer plusieurs
  fois la même composante C avec des valeurs différentes. Si cela
  arrive, c'est pour deux chemins différents, disons Q1 et Q2. Soit X
  les voisins de C qui sont à la fois dans Q1 et Q2.  X se décompose
  en segments, chacun étant un sous-chemin de Q1 inter Q2. Tout sommet
  voisin de C doit être dans X sinon C aurait un sommet de trop. Donc
  tout chemin de G contenant X tel que C est une composante de G\P,
  doit être parallèle à Q1 et Q2. En particulier Q1 et Q2 sont
  parallèles ... et alors [A FINIR] ? Bon, bah en fait c'est
  possible. Q1 peut se réduire à une arête xy et Q2 à xzy (un sommet
  de plus). Et avec Q1 c'est vrai, mais avec Q2 cela devient faux. On
  a des contre-exemple avec des sommets (z) de degré deux.
*/
{
  int **L=c->cmax,n=c->ncmax; /* L=liste et n=|L|*/
  int *T=c->tmax; /* T=taille de chaque composante */
  int i,r;

  /* on passe en revue chaque composante la liste L */
  for(i=0;i<n;i++){
    r=SetCmp(C,L[i],t,T[i]); /* compare C et L[i] */
    if(r&6) return 0; /* C est strictement contenu (&4) ou égale (&2) à L[i] */
    if((r&8)&&(n>0)){ /* L[i] est contenu (strictement) dans C et L non vide */
      free(L[i]); /* libère ce tableau de sommets précédemment alloué par un ps1_addmax() */
      /* si L[i] dernier élément de L (i=n-1), alors il n'y a plus rien à faire */
      n--;
      if(i<n){ /* si c->cmax[i] pas dernier élément, il faut déplacer le dernier */
	L[i]=L[n]; L[n]=NULL; /* pour éviter plus tard un double free. NB: i<>n */ 
	T[i]=T[n]; /* nombre de sommets dans L[i] */
	i--; /* pour recommencer avec ce nouveau L[i] */
      }
    }
  }

  /*ajoute C à la fin de L */
  ALLOCZ(L[n],t,C[_i]); /* alloue et copie C dans L[i] */
  T[n]=t; /* taille de C */
  c->ncmax=n+1; /* nouvelle taille de L */
  return 1; /* on a ajouté C à L */
}


int ps1_push(int x,const int v,conflit *c)
/*
  Affecte la valeur v (=0 ou 1) dans le noeud x et propage, en
  appliquant les règles décrites ci-après, à tous ses voisins et
  récursivement à tout le graphe. Renvoie 1 si une contradiction
  apparaît, ou bien si tous les noeuds de la paire de x ont pour code
  1. Sinon, on renvoie 0 (terminaison standard).

  Attention! Il faut appliquer cette fonction que si la paire de x est
  complète. Donc il ne faut pas l'appliquer si on vient de créer le
  noeud x, car on n'est pas sûr de ne pas faire un "Out of Memory" un
  peu plus tard sur la même paire.

  Effet de bord: met à jour plusieurs champs de la variable c, mais ne
  modifie pas le graphe c->G.

  Règles pour l'arc x->y:

  x=noeud x
  y=noeud y
  v=valeur 0 ou 1 qui vient d'être écrite dans x
  t=type de l'arête: 0=(Tx|Ty), 1=(|Tx|=|Ty|), 2=(Tx<Ty), 3=(Tx>Ty)
  c=graphe des conflits

  (Attention! "Tx|Ty" signifie que les composantes Tx et Ty sont
  disjointes.)

  règle 1: si Tx|Ty (disjoint) et v=0, alors écrire 1 dans y
  règle 2: si Tx=Ty, alors écrire v dans y
  règle 3: si Tx<Ty et v=0, alors écrire 0 dans y
  règle 4: si Tx>Ty et v=1, alors écrire 1 dans y
  règle 5: si Tx|Ty, v=1 et |Tx|+|Ty|=n, alors écrire 0 dans y

  On applique également la "règle du dernier -1", à savoir si la paire
  de x, après avoir écrit v, possède exactement 1 seul noeud de valeur
  -1, alors on écrit 0 dans ce noeud.

  La "règle du max" et la "règle de l'influence des voisins" sont
  appliquées plus directement par la fonction ps1() lors de la
  création d'une paire. Elles n'ont pas lieu d'être lors de la
  propagation.
*/
{
  int i,d,y,t;

  if(c->code[x]==v) return 0; /* rien à faire */
  if(c->code[x]==1-v) return 1; /* contradiction */
  /* ici on a code[x]==-1 */
  c->code[x]=v; /* écrit la valeur v */
  c->nbi--; /* et une valeur indéterminée en moins ! */
  t=(c->nbc[c->paire[x]] -= v); /* diminue seulement si on a écrit 1 */
  if(t==0) return 1; /* la paire de x est bonne, elle contient que des 1 ! */

  /* applique les règles 1-5 à tous les arcs sortant de x */
  for(i=0,d=c->G->d[x];i<d;i++){ /* pour tous les voisins de x */
    y=c->G->L[x][i]; /* y=i-ème voisin de x */
    t=(y>>LCONF); /* t=type de l'arc x->y: 0,1,2,3 */ 
    y &= CONFMASK; /* y=numéro du noeud voisin de x */

    /* applique les règles */
    switch(t){
    case 0:
      if((v==0)&&(ps1_push(y,1,c))) return 1; /* règle 1 */
      if((v==1)&&(c->tcomp[x]+c->tcomp[y]==c->n)&&(ps1_push(y,0,c))) return 1; /* règle 5 */
      break;
    case 1: if(ps1_push(y,v,c)) return 1; break; /* règle 2 */
    case 2: if((v==0)&&(ps1_push(y,0,c))) return 1; break; /* règle 3 */
    case 3: if((v==1)&&(ps1_push(y,1,c))) return 1; break; /* règle 4 */
    }
  }

  /* règle du dernier -1 ? */
  if(c->nbc[c->paire[x]]==1){
    /* on cherche alors l'unique noeud x de la paire courante qui est
       de code < 1.  NB: ce noeud existe forcément */
    x=c->paire[x]; /* x=premier noeud de la paire de x */
    while(c->code[x]==1) x++; /* passe au suivant si le code est 1 */
    return ps1_push(x,0,c); /* écrit 0 dans le noeud x */
  }

  return 0;
}


/* pour le débugage de PS1() */
/* LEVEL=niveau de récursion, POS=numéro de ligne */
int LEVEL=0;
#define PRINTS do{			       	\
    int _i;			       		\
    printf("%03i:%02i  ",++POS,LEVEL);		\
    for(_i=0;_i<3*LEVEL;_i++) printf(" ");	\
  }while(0)


int PS1(graph *G,path *P,const int version){
/*
  P est un chemin d'un graphe G, et G\P est formé d'une seule
  composante connexe. Il faut G et P <> NULL (initialisés avec
  new_graph() et new_path()). Renvoie 1 si P peut "séparer" G (voir
  explications ci-dessous). Il y a une optimisation avec le graphe des
  conflits si version>0. On utilise cette fonction comme ceci:

  path *P=new_path(G,NULL,G->n); // P=chemin vide, sans sommet
  int r=PS1(G,P); // r=0 ou 1
  free_path(P);
  
  Effet de bord: G->int1 retourne le nombre de tests (nombre de
  paires, nombre de chemins testés, et nombre de passes dans le graphe
  des conflits). Les autres champs de G et de P ne sont pas modifiés.

  Le paramètre "version" indique la variante du test:
  - version=0: sans le graphe des conflits
  - version=1: avec le graphe des conflits
  - version=2: comme version=1 mais sans le graphe des conflits lors de la récursivité
  - version=3: comme version=1 mais avec l'écriture de valeurs dans le graphe des conflits

  Améliorations possibles:

  - Si G n'est pas biconnexe, alors on pourrait tester si toutes ses
    composantes biconnexes sont bien évaluée à vraie. Si P
    n'intersecte pas une composante biconnexe B de G, alors il faut
    évaluer PS1(B,{}).

  - G\P étant connexe, on pourrait déjà supprimer les sommets de P qui
    forment un segment de P contenant une extrémité de P et qui n'ont
    pas de voisin dans G\P. En particulier, si une des extrémités de P
    est de degré 1, on peut la supprimer.

  - Privilégier les paires de sommets qui ne sont pas adjacents (pour
    diminuer la taille les composantes et avoir des chemins plus
    longs). Plus généralement, on pourrait privilégier les paires de
    sommets les plus distants. Ceci dit, on ne gagne probablement pas
    grand chose, et une renumérotation aléatoire des sommets devrait
    suffir pour ne pas traité les paires de sommets voisins en
    priorité.

  - On pourrait tester des cas simples pour G: arbre (tester si m=n-1,
    on sait que G est connexe), clique (tester si m=n(n-1)/2: si n<=4
    sommets alors vraie, sinon faux). (Ces tests sont déjà
    implémentés). Plus dur: G est outerplanar. En fait, si G est un
    arbre de cycles (chaque composante connexe est un cycle ou un K4),
    alors le test est vrai. C'est donc vrai en particulier si
    m=n. Pour tester si G est un arbre de cycle, il suffit de faire un
    DFS, puis de vérifier que si (u,x) et (u,y) sont deux arêtes qui
    ne sont pas dans l'arbre du DFS, alors x et y ne sont pas ancêtres
    l'un de l'autre (??? Pourquoi ???).

  - Pour tous les chemins possibles testés récursivement pour G, ne
    tester effectivement que ceux qui correspondent à des
    sous-ensembles de sommets différents puisque le résultat sera le
    même (mêmes composantes connexes). Pour cela, il faut gérer une
    table pour mémoriser les sous-ensembles testés. Notons que si un
    chemin est induit alors cela ne sert à rien de le mémoriser, il ne
    pourra jamais être rencontré de nouveau. On pourrait coder un
    sous-ensemble de n sommets par un entier sur n bits (n < 32 ou 64
    donc). La recherche/insertion pourrait être une recherche dans un
    arbre binaire de recherche.

 EXPLICATIONS:

  Soit P et Q deux chemins d'un graphe G. On dit que Q est parallèle à
  P, noté Q//P, s'il n'y a pas d'arête entre P\Q et G\(Q u P).

  Quelques propriétés:

  - La relation // n'est pas symétrique.
  - Tout chemin est parallèle au chemin vide.
  - Tout chemin est parallèle à lui-même.
  - La relation // n'est pas transitive.

  Pour la dernière proposition, on peut vérifier en fait que si R//Q
  et Q//P, alors on a R//P sauf s'il existe une arête uv telle que u
  in (P inter Q)\R et v in Q\(R u P).

  Soit P et Q deux chemins de G. On dit que Q dérive de P, noté Q///P,
  s'il existe une suite P_0,P_1,...,P_n de chemins de G avec P_0=Q et
  P_n=P tels que P_{i-1}//P_i pour chaque i=1..n.

  Quelques propriétés:

  - Si Q//P, alors Q///P. En particulier, Q///{} puisque Q//{}.
  - On peut avoir R//Q//P, soit R///P, sans avoir R//P (cf. ci-dessus).
  - Si R///Q et Q///P, alors R///P.

  On dit que Q///P dans un graphe valué (G,w) si tous les chemins
  P_0,...,P_n (en particulier Q et P) sont des plus courts chemins
  selon w.

  Soit P un chemin de G et w une valuation de G. On définit le
  potentiel pour P selon w comme score(P,w) := max_C { w(C)*|V(G)|+|C|
  } où le maximum est pris sur toute composante connexe C de G\P.

  Lemme 1. Supposons que G\P a une seule composante connexe, et soit Q
  un chemin de G parallèle à P différent de P. Alors, pour chaque
  valuation w de G, soit P est un demi-séparateur de G ou bien
  score(Q,w) < score(P,w).

  Preuve. Supposons que P ne soit pas un demi-séparateur de G pour la
  valuation w. Soit C la composante de G\P, et posons n = |V(G)|. Par
  définition, score(P,w) = w(C)*n + |C|. On a w(C) > w(G)/2, et donc
  w(P) < w(G)/2. Il suit que w(P) < w(C). Soit C' une composante de
  G\Q telle que w(C')*n+|C'| = score(Q,w). Comme Q est parallèle à P,
  soit C' est contenue dans C soit C' est contenue dans P. En effet,
  si C' intersecte C et P, alors C' contient une arête uv avec u in C
  et v in P. Bien sûr uv not in Q. Cela contredit le fait qu'il existe
  pas d'arête entre P\Q et G\(QuP).

  Si C' est contenue dans P, alors w(C') <= w(P) < w(C). Il suit que
  w(C') <= w(C)-1, soit w(C')*n <= w(C)*n - n. Clairement |C'| < n +
  |C|. D'où w(C')*n+|C'| < w(C)*n+|C|, soit score(Q,w) < score(P,w).

  Si C' est contenue dans C, alors w(C') <= w(C) et |C'| < |C| car
  Q<>P. Il suit que w(C')*n + |C'| < w(C)*n + |C|, soit score(Q,w) <
  score(P,w).

  Dans les deux cas nous avons prouvé que score(Q,w) < score(P,w).
  QED

  Soit G un graphe et P un chemin de G tel que G\P est composé d'une
  seule composante connexe (en particulier G est connexe). On définit
  PS1(G,P) le prédicat qui est VRAI ssi pour toute pondération w de G
  telle que P est un plus court chemin il existe dans (G,w) un chemin
  demi-séparateur qui dérive de P.

  Lemme 2. G est dans PS1 ssi PS1(G,{}) = VRAI.

  Preuve. En effet, en réécrivant la définition de PS1(G,{}) on déduit
  que PS1(G,{}) est VRAI ssi pour toute pondération w de G il existe
  dans (G,w) un chemin demi-séparateur qui dérive du chemin vide (tout
  chemin dérive du chemin vide). Notons que c'est nécessairement un
  plus court chemin de (G,w). C'est précisemment la définition de la
  classe PS1. QED

  L'objectif est d'avoir un test noté ps1(G,P) qui implémente
  PS1(G,P), disons qu'il s'en approche. La propriété souhaitée est que
  si ps1(G,P) est VRAI, alors PS1(G,P) aussi. En particulier, si
  ps1(G,{}) est VRAI, alors G est dans PS1.

  Algorithme pour le test ps1(G,P):

  On renvoie VRAI s'il existe une paire x,y de sommets où y n'est pas
  dans P telle que tout chemin Q de x à y:
  1. Q est parallèle à P, et
  2. pour toute composante C de G\(QuP), ps1(CuQ,Q)=VRAI.

  Lemme 3. Si ps1(G,P)=VRAI, alors PS(G,P)=VRAI.

  Preuve [A FINIR]. Par induction sur le nombre de sommets hors de
  P. Soit C est l'unique composante G\P. Si |C|=0, alors P est un
  demi-séparateur de G et donc PS(G,P) est VRAI. Supposons le lemme
  vrai pour tout entier < |C|.

  On suppose donc que ps1(G,P) est VRAI. Soit x,y la paire de sommets
  telle que y n'est pas dans P et où tous les chemins de x à y sont
  parallèles à P. En particulier, pour chaque valuation w, où P est un
  plus court chemin, tout plus court chemin Q selon w entre x et y est
  parallèle à P. Comme Q est différent de P (à cause du choix de y),
  on peut appliquer le lemme 1, et donc soit P est un demi-séparateur,
  soit score(Q,w) < score(P,w). On peut supposer qu'on est pas dans le
  premier cas, c'est-à-dire que P n'est pas un demi-séparateur de G,
  puisque sinon PS(G,P)=VRAI et le lemme est prouvé.

  Si Q est un demi-séparateur pour G, alors PS(G,P) est VRAI puisque Q
  est parallèle à P. Supposons que Q n'est pas un demi-séparateur pour
  G, et soit C' la composante de G\Q telle que w(C')>w(G)/2.

  On peut appliquer l'induction sur (C'uQ,Q) car comme Q<>P,
  |C'|<|C|. Posons G'=C'uQ. D'après le test ps1(G,P), ps1(G',Q) est
  VRAI. Donc par induction PS(G',Q)=VRAI et G' contient un chemin
  demi-séparateur pour la valuation w, disons P', parallèle à
  Q. Montrons d'abord que dans G, P' est parallèle à P.

  ...

  w(C')<w(G)/2 ...

  QED

  Lemme 4. Si ps1(G,P)=VRAI, alors soit il existe une paire x,y de
  sommets de G telle que tout chemin de x à y contient P, ou bien il
  n'existe pas de sommet de P ayant trois voisins inclus dans un cycle
  de G\P. [PAS SÛR, A VÉRIFIER]

  Preuve. [A FINIR] Soit x,y une paire de sommets de G avec y pas dans
  P telle que tous les chemins de x à y soient parallèles à P. Soient
  Q un tel chemin. Supposons que Q ne contient pas P. ...  QED

  Dit autrement, si on n'a pas cette propriété, alors ps1(G,P)=FAUX et
  il est inutile de tester tous les chemins Q possibles.  Remarquons
  que si G est 2-connexe, alors il ne peut avoir de paire x,y de
  sommets (avec y pas dans P) où tout chemin de x et y contient P. A
  montrer: ps1(G,P)=VRAI ssi tout les composantes 2-connexes G' de G
  on a ps1(G',P inter G)=VRAI ...

  [A VOIR]

  - u := ps1(CuQ,Q)
  - ajoute (C,u) à la liste L des composantes maximales (voir ps1_addmax)
  - si u = FAUX, ajouter un nouveau noeud au graphe des conflits (GC)
  - recommencer avec le chemin Q suivant de mêmes extrémités xy (s'il existe)
  - s'il n'y a plus de tel chemin:

    On essaye d'appliquer les trois règles (max, influence, dernier):
    - la règle du max en tenant compte de la liste L.
      si |L|<>1, alors on ne peut pas appliquer la règle du max
      sinon, L={(C,u)}
        si u = VRAI, on peut supprimer la paire xy de GC
        sinon, alors on peut appliquer la règle du max
    - la règle d'influence des voisins
    - la règle du dernier -1

    Si l'application d'une de ces règles (avec ps1_push) produit une
    contradiction dans GC, alors on a trouvé une bonne paire xy, et on
    renvoie VRAI.  Sinon, s'il existe une autre paire xy on recommence
    avec cette pnouvelle paire.

  A la fin du traitement de toutes les paires:

  - soit il reste des indéterminées (-1), et il faut les éliminer. On
    les force d'abord à 0. S'il y a une contradiction (en faisant des
    ps1_push), on en déduit que la valeur doit être 1 (et on fait
    ps1_push). Sinon, on force la valeur à 1. Si y a une
    contradiction, on déduit que la valeur doit être 0 (et on fait
    ps1_push). Sinon, s'il y a une valeur initialement indéterminée
    qui passe à la même valeur pour le forcage à 0 et à 1, on élimine
    cette indéterminée (et on fait ps1_push). On fait ainsi pour
    chaque indéterminée. Chaque fois qu'une indéterminée est éliminée,
    on recommence la recherche sur l'ensemble des indéterminées (et
    pas seulement sur celles qui restent). Si on arrive ainsi à
    éliminer toutes les indéterminées on peut passer à la
    suite. Sinon, on ne peut pas faire grand chose à part essayer tous
    les système possibles ... Voir MINISAT+ (système pseudo booléens)
    et Sugar qui transforme du système linéaire ou CSP en SAT.

  - soit il n'y a plus d'interminées. Dans ce cas on peut déduire un
    système d'équations linéaire indépendantes où les inconnues sont
    les poids des sommets. Si le système n'a pas de solution, alors
    renvoyer VRAI. Sinon, la solution trouvée peut renseigner sur une
    valuation possible pour prouver que G est éventuellement pas dans
    PS1(G,P).

*/
  G->int1=1; /* compte les tests */

  DEBUG(
	LEVEL++;
	int u;int v;
	if(G==GF) LEVEL=POS=0;
	printf("\n");
	PRINTS;printf("version=%i G=%p n=%i ",version,G,G->n);PRINTT(P->P,P->n);
	PRINTS;printf("G=");
	for(u=0;u<G->n;u++){
	  printf("%i:[",u);
	  for(v=0;v<G->d[u];v++){
	    printf("%i",G->L[u][v]);
	    if(v<(G->d[u])-1) printf(" ");
	  } printf("] ");
	} printf("\n");
	);
  
  const int n=G->n;

  /* Ici on élimine un certain nombre de cas faciles à tester.
     Attention ! vérifier que G est dans PS1 dans ces cas là ne suffit
     pas (ce qui revient à vérifier que PS1(G,{})=VRAI). Il faut être
     certain qu'on a en fait PS1(G,P)=VRAI. */

  if(n-(P->n)<3){
    DEBUG(PRINTS;printf("-> n-|P|<3, moins de 3 sommets hors P\n"););
    DEBUG(LEVEL--;);
    return 1;
  }
  /* Par hypothèse P est léger, donc la composante hors de P est
     lourde. Il suffit de prendre un chemin entre ces au plus deux
     sommets hors de P. */

  /* lors d'un appel récursif, le nombre d'arêtes G->m est déjà mis à
     jour car G provient alors du résultat de ExtractSubgraph(). Donc
     NbEdges(G) est calculé en temps constant dès le 2e appel. */

  if((!P->n)&&(n<6)&&(NbEdges(G)<10)){
    DEBUG(PRINTS;printf("-> |P|=0 et n<6 et pas K₅\n"););
    DEBUG(LEVEL--;);
    return 1;
  }

  /* ici on a au moins 3 sommets de G qui sont hors de P, et P est non
     vide. Alors, comme P contient au moins deux somemts, cela fait
     que G possède au moins 5 sommets. */

  if(NbEdges(G)==((n*(n-1))>>1)){
    DEBUG(PRINTS;printf("-> clique avec > 2 sommets hors P\n"););
    DEBUG(LEVEL--;);
    return 0;
  }
  /* Ici G est clique avec n>4 et n-|P|>2. Dans le cas où tous les
     poids sont à 1 (sommets et arêtes) alors, il n'existe aucun
     chemin Q parallèle P permettant de progresser. Notons qu'une
     clique G avec n sommets donne PS1(G,P)=VRAI s'il y a 0, 1 ou 2
     sommets dans G\P. */

  if(NbEdges(G)<=n){
    DEBUG(PRINTS;printf("-> m<=n\n"););
    DEBUG(LEVEL--;);
    return 1;
  }
  /* Dans ce cas, il s'agit d'un arbre avec un cycle. Soit C la
     composante lourde de G\P, et x une des extrémités de P. La
     composante C est connectée à P par au plus deux sommets, disons u
     et v (u=v possible) et u le plus près de x. Soit y le voisin de v
     dans C. On définit alors P' comme le plus court chemin de G
     allant de x à y. Il est facile de voir que P' longe P depuis x
     puis entre dans C soit par u soit par v. Dans les deux cas les
     sommets de C ne peuvent être connectés à aucun sommet de P\P',
     les deux seules arêtes connectant C à P étant détruite par P'. */

  /* ici G possède au moins 5 sommets et 6 arêtes. */

  int x,y,i,u,v,w,d;
  path *Q=new_path(G,NULL,n);
  path *R=new_path(G,NULL,n);
  param_dfs *p=new_param_dfs(n); /* p->C n'est pas alloué */
  graph *C;

  ALLOC(p->C,n);   /* pour le DFS avec sommets supprimés */
  NALLOC(int,T,n); /* pour la composante C de G */
  NALLOC(int,M,n); /* pour la compatiblité des chemins P et Q */

  conflit c; /* ensembles de variables pour gérér le graphe des conflits */
  int npaire,npath,goodxy;
  /* npaire = 1er noeud dans le graphe des conflits de la paire courante */
  /* npath = 1er noeud dans GC du chemin courant pour la paire courante */
  /* goodxy = 1 ssi la paire xy est bonne, ie. tous les chemins et comp. sont ok */

  if(version>0){
    c.n=n; /* nombre de sommets du graphe G */
    c.G=new_graph(CONFMAX); /* graphe des conflits, alloue G->d et G->L */
    c.G->n=c.nbi=c.nbzi=c.ncmax=c.outmem=0; /* c.G->n=nb de noeuds déjà crées */
  }

  /* pour toutes les paires de sommets x,y de G avec y pas dans P */

  for(x=0;x<n;x++)
    for(y=0;y<n;y++){

      if(P->V[y]>=0) continue; /* y ne doit pas être dans P */
      if((P->V[x]<0)&&(y<=x)) continue; /* si x et y pas dans P, ne tester que le cas x<y */

      /* ici on a soit:
	 1) x dans P et y pas dans P, ou bien
	 2) x et y pas dans P et x<y
      */

      (G->int1)++; /* +1 pour chaque paire testée */
      goodxy=1; /* par défaut on suppose la paire xy comme bonne */

      /* calcule le 1er chemin entre x et y */
      Q->P[0]=x; Q->P[1]=y; /* initialise les extrémités du chemin Q */
      if(!NextPath(G,Q,-1)) goto fin_ps1; /* fin si pas de chemin x->y, impossible si G connexe */
      if((version>0)&&(!c.outmem)){
	npaire=c.G->n; /* initialise le 1er noeud de la paire courante */
	c.nbc[npaire]=0; /* nombre de valeurs < 1 pour la paire courante */
	/* on efface la liste des composantes maximales, s'il y en avait */
	for(u=0;u<c.ncmax;u++) free(c.cmax[u]);
	c.ncmax=0;
      }

      DEBUG(PRINTS;printf("Paire (%i,%i) G=%p n=%i\n",x,y,G,n););

      do{ /* pour tous les chemins Q entre x et y */
	(G->int1)++; /* +1 pour chaque sommet testé */

	DEBUG(PRINTS;printf("Essai du chemin ");PRINTT(Q->P,Q->n););

	/* On vérifie que Q est parallèle avec P. Il faut qu'aucun
	   sommet de P\Q n'ait de voisin en dehors de P ou de Q. */

	for(i=0;i<P->n;i++){
	  if(Q->V[u=P->P[i]]>=0) continue; /* on est aussi dans Q */
	  /* ici on est dans P\Q */
	  for(v=0;v<G->d[u];v++){ /* vérifie chaque voisin v de u */
	    if(P->V[G->L[u][v]]>=0) continue; /* si v est dans P */
	    if(Q->V[G->L[u][v]]>=0) continue; /* si v est dans Q */
	    DEBUG(PRINTS;printf("-> chemin Q non parallèle à P\n"););
	    goodxy=0; /* cette paire n'est pas bonne */
	    goto nextxy; /* aller à la prochaine paire */
	  }
	}

	/* ici Q est parallèle à P */

	DEBUG(PRINTS;printf("-> chemin Q parallèle à ");PRINTT(P->P,P->n););

	/* on vérifie que pour chaque composante C de G\Q,
	   PS1(CuQ,Q)=VRAI. Notons que les sommets de P\Q sont
	   léger. Il ne faut donc pas les considérer. Il faut donc
	   prendre les composantes de G\(QuP).*/

	/* on enlève de G les sommets de QuP pour calculer les
	   composantes de G\(QuP) */
	for(u=0;u<n;u++) /* initialise le dfs */
	  if((P->V[u]<0)&&(Q->V[u]<0)) p->C[u]=-1; else p->C[u]=-2;
	dfs(G,0,p); /* calcule les composantes de G\(QuP) */
	DEBUG(PRINTS;printf("#composantes dans G\\(QuP) à tester: %i\n",p->nc););
	/* si aucune composante inutile de tester récursivement ce
	   chemin, G\(QuP) est vide. On peut passer au prochain chemin */
	if(p->nc==0) goto nextQ;
	/* ici, il y a au moins une composante dans G\(QuP) */

	d=R->n=Q->n; /* d=nombre de sommets du chemin Q */
	for(u=0;u<d;u++) T[u]=Q->P[u]; /* T=liste des sommets de Q */
	if((version>0)&&(!c.outmem)) npath=c.G->n; /* initialise le chemin Q au noeud courant */

	DEBUG(PRINTS;printf("Q = ");PRINTT(T,d););

	/* pour chaque composante de G\Q */

	for(i=0;i<p->nc;i++){  /* T[d..v[=i-ème composante de G\(QuP) u Q */
	  for(u=0,v=d;u<n;u++) if(p->C[u]==i) T[v++]=u;
	  /* T[0..v[=sommets de Q u C */
	  /* T[0..d[=sommets de Q, T[d..v[=sommets de la composante */
	  /* NB: les sommets de T[d..v[ sont dans l'ordre croissant */

	  DEBUG(PRINTS;printf("C\\Q=CC(%i) = ",i);PRINTT(T+d,v-d););

	  C=ExtractSubgraph(G,T,v,1); /* crée C=G[T] avec v=|T| sommets */
	  /* Attention! Q n'est pas un chemin de C, mais de G. On crée
	     donc un nouveau chemin R dans C qui est équivalent à Q */
	  for(u=0;u<v;u++) R->V[u]=-1; /* Attention! boucle avec v=C->n sommets */
	  for(u=0;u<d;u++) R->P[u]=C->pint1[Q->P[u]]-1;
	  for(u=0;u<d;u++) R->V[R->P[u]]=u;

	  DEBUG(PRINTS;printf("Q = ");PRINTT(R->P,R->n););

	  /* appel récursif avec une nouvelle version 0 ou 1:
	     - si version=0: sans le graphe des conflits -> 0
	     - si version=1: avec le graphe des conflits -> 1
	     - si version=2: comme version=1 mais sans le graphe
	       des conflits lors de la récursivité -> 0
	     - si version=3: comme version=1 mais avec l'écriture
	       de valeurs dans le graphe des conflits -> 1
	  */
	  u=PS1(C,R,version%2);

	  DEBUG(PRINTS;printf("PS1(CuQ,Q)=%i\n",u););
	  DEBUG(if(c.outmem){PRINTS;printf("PROBLEME MEMOIRE !\n");});

	  G->int1 += C->int1; /* met à jour le nombre de tests (ajoute au moins 1) */
	  free_graph(C); /* libère C qui ne sert plus à rien */

	  /* à faire que si graphe des conflits et pas eut de problème mémoire */
	  if((version>0)&&(!c.outmem)){
	    
	    /* on vérifie si la même composante n'existe pas déjà dans
	       la paire de c.G->n. Si c'est le cas, on passe à la
	       prochaine composante. NB: la composante de c.G->n est
	       dans T[d...v[ et sa taille vaut v-d. */

	    for(w=npaire;w<c.G->n;w++)
	      if(SetCmp(c.comp[w],T+d,c.tcomp[w],v-d)&2)
		goto nextC; /* si même composante, alors prochaine composante */
	    /* ici w=c.G->n */
	    /* si u=0, w sera le nouveau noeud du graphe des conflits */

	    /* ajoute T[d..v[ à la liste des composantes maximales (c.cmax) */
	    if(ps1_addmax(T+d,v-d,&c)){ c.nodemax=w; c.valmax=u; }
	    /* On stocke dans c->nodemax le noeud de la dernière
            composante ajoutée à la liste des composantes maximales.
            Si à la fin du traitement de la paire xy la liste a
            exactement une seule composante, alors c->nodemax sera
            forcément le noeud de la dernière composante ajoutée et
            c->valmax sa valeur (VRAIE/FAUSSE). */
	    if(c.ncmax>=NCMAX){ /* dépassement du nb de composantes maximales */
	      c.outmem=1; /* Ouf of Memory */
	      /* la paire xy n'est pas complète. On l'enlève, sinon
	         l'analyse des conflits pourrait ne va pas être correcte */
	      ps1_delxy(&c,npaire); /* supprime la dernière paire créee */
	      goodxy=0; /* cette paire n'est pas bonne */
	      goto nextxy; /* change de paire */
	    }
	  }
	  /* ici outmem=0 */

	if(u) goto nextC; /* si u=VRAI, on a finit le traitement de la composante */
	goodxy=0; /* une composante n'est pas bonne pour la paire xy */
	if(version==0) goto nextxy; /* si pas graphe de conflits, alors changer de paire */

	/* ici u=FAUX et version>0 */
	/* on va donc essayer d'ajouter le noeud w au graphe des conflits */
	/* ici pas de problème de mémoire, donc on ajoute le noeud w */

	(c.G->n)++; /* un noeud de plus */
	if(c.G->n==CONFMAX){ /* Out of Memory */
	  c.outmem=1;
	  ps1_delxy(&c,npaire); /* supprime la dernière paire créee, qui n'est pas bonne */
	  goto nextxy; /* ici goodxy=0 */
	}
	/* ici v=|T|, d=|Q|, w=nouveau noeud à créer */
	
	/* Attention! ne pas utiliser ps1_push(w,...) juste après la
	   création de w, car on est pas certain que la paire de w
	   sera complète ... On ne peut faire ps1_push(w,...)
	   seulement après le while(NextPath...). Idem pour
	   l'application de la règle du max ! */

	/* initialise le nouveau noeud w */

	DEBUG(PRINTS;printf("Ajoût du noeud %i dans le graphe des conflits\n",w););

	c.x[w]=x; c.y[w]=y; /* mémorise la paire x,y (pour l'affichage) */
	c.path[w]=npath; c.paire[w]=npaire; /* w rattaché à npath et à npaire */
	c.code[w]=-1; c.nbi++; /* un noeud de plus à valeur = -1 */
	(c.nbc[npaire])++; /* un noeud de plus à valeur < 1 pour cette paire */
	c.G->d[w]=0; /* initialise le degré de w */
	ALLOC(c.G->L[w],CONFMAX); /* crée la liste des voisins de w */
	ALLOCZ(c.comp[w],v-d,T[d+(_i)]); /* crée & copie la composante de w */
	c.tcomp[w]=v-d; /* taille de la composante de w */
	
	/* y'a-t'il des arcs vers les noeuds précédant w ? */
	/* calcule les arêtes sortantes de w et de leurs types */
	
	for(u=0;u<w;u++){ /* pour chaque noeud < w, v=type de l'arc */
	  if(u>=npath) v=0; /* test rapide: si u>=napth, alors clique */
	  else{ /* sinon calcule l'intersection des composantes */
	    v=SetCmp(c.comp[u],c.comp[w],c.tcomp[u],c.tcomp[w]);
	    /* ici les valeurs possibles pour SetCmp: 0,1,3,5,9 car T1 et T2 != {} */
	    if(v==1) v=-2; /* si v=1, T1 intersecte T2, et donc pas d'arête u-w */
	    /* avant: v: 0=disjoint, 3=égalité, 5=T1 sub T2, 9=T2 sub T1 */
	    v >>= 1; v -= (v==4); /* rem: si avant v=-2, alors après v=-1 */
	    /* après: v: 0=disjoint, 1=égalité, 2=T1 sub T2, 3=T2 sub T1 */
	  } /* si v<0, alors pas d'arête */
	  if(v>=0){ /* ajoute les arcs u->w et w->u */
	    ADD_ARC(c.G,u,(v<<LCONF)|w); /* u->w */
	    if(v>1) v ^= 1; /* si v=2, alors v->3, et si v=3, alors v->2 */
	    ADD_ARC(c.G,w,(v<<LCONF)|u); /* w->u (asymétrie pour v=2,3) */
	  }
	} /* fin du for(u=...) */
	
	nextC:; /* prochaine composante (next i) */
	} /* fin du for(i=0...) */

      nextQ:; /* prochain chemin Q */
      }while(NextPath(G,Q,0));
      
      /* si ici goodxy=1, c'est que pour la paire xy, tous les chemins
	 entre xy sont parallèles à P et qu'on a jamais vue de
	 composante à FAUX. Dans ce cas on a trouvé une bonne paire et
	 on a fini. */

      if(goodxy) goto fin_ps1; /* termine l'algo. */
      if(version==0) goto nextxy; /* si pas graphe des conflits */
      /* ici version>0 et goodxy=0 */
      
      /* Ici on a examiné tous les chemins Q pour la paire xy et on
      crée au moins un noeud dans le graphe des conflits. Ici tous les
      noeuds crées ont comme valeur -1. Il reste à essayer d'appliquer
      trois règles:

	 1) la règle du max: on écrit 0 s'il existe une composante de
	 xy contenant toutes les autres (y compris celle évaluées
	 récursivement à VRAI). Pour cela il faut avoir exactement une
	 seule composante maximale dans c.cmax qui doit être à
	 FAUX. Si elle est à VRAI il faut supprimer la paire xy (car
	 la plus lourde est VRAI, donc toutes les autres sont légères
	 !) et passer à la suivante.

	 2) la règle d'influence des voisins: vérifier si des voisins
	 v de chaque noeud u (de la paire xy) ne peuvent pas
	 influencer u. Par exemple, si la composante de v est la même
	 que celle de u il faut alors écrire cette valeur dans le code
	 de u.

	 3) la règle du dernier -1: si le nb de valeurs < 1 est
	 exactement 1, alors il faut écrire 0 dans ce noeud là. Il
	 faut, comme pour la règle d'influence des voisins, balayer
	 tous les sommets u de la paire xy.

	 NB: on peut traiter dans un même parcours les deux dernières
	 règles.
      */

      v=c.paire[c.G->n-1]; /* v=1er noeud de la dernière paire */
      /* NB: c.outmem=0 et c.G->n > 0 */

      /* règle du max */
      if(c.ncmax==1){ /* sinon pas de règle du max */
	if(c.valmax){ /* supprime la paire xy */
	  ps1_delxy(&c,v);
	  goto nextxy;
	}
	/* applique la règle du max: on écrit 0 dans la composante maximale */
	if(ps1_push(c.nodemax,0,&c)){ goodxy=1; goto fin_ps1; } /* fin si contradiction */
      }

      /* règle d'influence des voisins + règle du dernier -1 */
      for(u=v;u<c.G->n;u++){ /* pour tout noeud u de la dernière paire */
	if(c.code[u]>=0) continue; /* plus rien à faire */
	/* NB: si code[u]=0 ou 1, c'est qu'on a forcément déjà fait
	   un ps1_push() sur u. Et donc la règle d'influence des
	   voisins n'a pas a être testée sur u */
	else /* code=-1 */
	  if(c.nbc[c.paire[v]]==1){ /* règle du dernier -1 */
	    /* tous les sommets sont à 1, sauf u qui est ici a -1 */
	    if(ps1_push(u,0,&c)){ goodxy=1; goto fin_ps1; } /* fin si contradiction */
	    else break; /* on peut sortir du for(u...), on a tout testé */
	  }
	
	for(i=0,d=c.G->d[u];i<d;i++){ /* scan tous les voisins v de u */
	  v=c.G->L[u][i]&CONFMASK; /* v=i-ème voisin de u */
	  w=c.code[v]; if(w<0) continue; /* on ne déduit rien */
	  /* efface le code de v pour pouvoir le réécrire avec
	     ps1_push().  Il faut prendre des précaution pour que c
	     soit cohérent (il faut mettre à jour .nbi et .nbc). NB:
	     l'ancien code de v est w=0 ou 1, celui de u est aussi 0 ou 1 */
	  c.code[v]=-1;
	  c.nbi++; /* met à jour le nombre d'indéterminées */
	  c.nbc[c.paire[v]] += w; /* augmente s'il y avait 1 */
	  /* réécrit le code pour v */
	  if(ps1_push(v,w,&c)){ goodxy=1; goto fin_ps1; } /* fin si contradiction */
	  /* NB: la propagation de cette écriture va tester l'arc
	     v->u et donc modifier éventuellement le code de u */
	}
      }

    nextxy:;
    } /* fin du for(x,y=...) */
  
  /* Ici on a testé toutes paires xy et aucune n'est bonne */
  /* NB: goodxy=0 */

  if((version>0)&&(!c.outmem)){

    /* Traitement final du graphe des conflits. On pourrait ici
       essayer d'éliminer les derniers -1. Seulement, il semble
       qu'aucune des règles ne peut plus être appliquées. */

  loop_ps1:
    do{
      if(c.nbi==0) break; /* fini: il ne reste plus de -1 à modifier */
      x=c.nbi; /* mémorise le nb total de valeurs à -1 */
      //
      // que faire ? essayer toutes les valeurs possibles pour les -1 ?
      //
      // G->int1++;
    }
    while(c.nbi<x); /* tant qu'on a enlevé au moins un -1, on recommence */
    
    /* ps1x */
    if(version==3){ /* NB: si appel récursif, alors version!=3 */
      i=MEM(CPARAM,0,int);
      if(i>0){ /* lit les valeurs à l'envers */
	u=MEM(CPARAM,(2*i-1)*sizeof(int),int);
	v=MEM(CPARAM,(2*i)*sizeof(int),int);
	MEM(CPARAM,0,int)--;
	if(ps1_push(u,v,&c)){ goodxy=1; goto fin_ps1; } /* fin si contradiction */
      goto loop_ps1;
      }
    }
    
    /*
      on cherche les noeuds "0" indépendants, correspondant à des
      composantes lourdes (donc à des inégalités).
      
      Algorithme: on ne sélectionne que les noeuds de code=0 et qui
      n'ont pas de voisins v de code aussi 0 avec v<u (inclus). Ils
      doivent donc être de code 0 et minimal par rapport aux autres
      voisins de code 0. Si un tel noeud existe, on marque tous ces
      voisins de code 0 qui ne pourront plus être sélectionnés.
    */

    c.nbzi=0; /* =nb de zéros indépendants */
    for(u=0;u<c.G->n;u++){ /* pour chaque noeud du graphe des conflits */
      if(c.code[u]) continue; /* saute le noeud si pas 0 */
      for(i=0,d=c.G->d[u];i<d;i++){ /* scan tous les voisins u (qui est forcément à 0) */
	v=c.G->L[u][i]; w=c.code[v&CONFMASK];
	if(!((w==0)||(w==CONFC))) continue; /* on ne regarde que les voisins de code=0 */
	if((v>>LCONF)==3) break; /* ici i<d */
      }
      if(i==d){ /* on a trouvé un noeud u de code 0 et sans aucun arc u->v avec v<u */
	c.nbzi++; /* un de plus */
	for(i=0,d=c.G->d[u];i<d;i++){ /* modifie le code de tous les voisins à 0 de u */
	  v=c.G->L[u][i]&CONFMASK; /* v=voisin de u */
	  if(c.code[v]==0) c.code[v]=CONFC; /* on ne modifie pas deux fois un voisin */
	}
      }else c.code[u]=CONFC; /* u n'est pas bon, on ne laisse pas la valeur 0 */
    }

  }
  /* ici on a fini le traitement du graphe des conflits */

 fin_ps1: /* termine l'algo avec goodxy=0 ou 1 */

  if(version>0){ /* affichage et libération du graphe des conflits */
    if((!goodxy)&&(G==GF)) PrintConflit(&c); /* G=GF => affiche que le 1er niveau de récurrence */
    /* efface la liste des composantes maximales */
    for(u=0;u<c.ncmax;u++) free(c.cmax[u]);
    /* efface le graphe des conflits */
    for(u=0;u<c.G->n;u++) free(c.comp[u]);
    free_graph(c.G); /* c.G->L après c.G->n n'ont pas été alloués ou sont déjà libérés */
  }

  /* efface les tableaux alloués */
  free_path(Q);
  free_path(R);
  free_param_dfs(p);
  free(T);
  free(M);

  DEBUG(LEVEL--;);
  return goodxy;
}


/***********************************

       ROUTINES POUR LES
        ROUTING SCHEMES

***********************************/


int nca_bfs(int u,int v,const param_bfs *X){
/*
  Calcule le plus petit ancêtre commun entre u et v dans l'arbre donné
  par le bfs X. Si la forêt n'est pas connexe, alors on renvoie -1.

  Algorithme (de complexité dist(u,v)): on remonte d'abord le sommet
  le plus profond. Lorsque les deux sommets sont à la même profondeur
  on teste si c'est les mêmes et si oui on a trouver le nca, ou alors
  on les remontes tous les deux.
*/
  if(X->D[u]<X->D[v]){ int w; SWAP(u,v,w); }

  /* ici u est plus profond que v (ou à même profondeur) */
  while(X->D[u]!=X->D[v]) u=X->P[u];

  /* ici u et v sont à la même profondeur */
  while(u!=v) u=X->P[u],v=X->P[v];
	
  return u;
}


static inline int dist_nca(const int u,const int v,const int w,const int *D){
/*
  Calcule la distance de u à v dans un arbre via leur ancêtre commun
  w, en entré on doit avoir que D[u] est la distance de u à la racine
  de l'arbre. Par rapport à dist_bfs() il est nécessaire ici que w
  soit déterminé, en particulier que u et v soit dans le même arbre.
*/
  return D[u] + D[v] - (D[w]<<1);
}


static inline int dist_bfs(const int u,const int v,const param_bfs *X){
/*
  Calcule la distance de u à v dans l'arbre (en fait la forêt) donné
  par le bfs X. Si les sommets sont dans des composantes connexes
  différentes, alors on renvoie -1. C'est basé sur le calcule du plus
  petit ancêtre commun.
*/
  int w=nca_bfs(u,v,X);
  if(w<0) return -1;
  return dist_nca(u,v,w,X->D);
}


char *TopChrono(const int i){
/*
  Met à jour le chronomètre interne numéro i (i=0..CHRNONMAX-1) et
  renvoie sous forme de char* le temps écoulé depuis le dernier appel
  à la fonction pour le même chronomètre. La précision dépend du temps
  mesuré. Il varie entre la seconde et le 1/1000 de seconde. Plus
  précisément le format est le suivant:

  1d00h00'  si le temps est > 24h (précision: 1')
  1h00'00"  si le temps est > 60' (précision: 1s)
  1'00"0    si le temps est > 1'  (précision: 1/10s)
  1"00      si le temps est > 1"  (précision: 1/100s)
  0"000     si le temps est < 1"  (précision: 1/1000s)

  Pour initialiser et mettre à jour tous les chronomètres (dont le
  nombre vaut CHRONOMAX), il suffit d'appeler une fois la fonction,
  par exemple avec TopChrono(0). Si i<0, alors les pointeurs alloués
  par l'initialisation sont désalloués. La durée maximale est limitée
  à 100 jours. Si une erreur se produit (durée supérieure ou erreur
  avec gettimeofday()), alors on renvoie la chaîne "--error--".
  
  On n'utilise pas clock() qui se limite à 72 minutes. Pour des durées
  en micro secondes mais plus beaucoup plus longues, il faut utiliser
  gettimeofday() comme ceci:

  gettimeofday(&t0,NULL);
  ...
  gettimeofday(&t1,NULL);
  long long t = (t1.tv_sec-t0.tv_sec)*1000000LL + t1.tv_usec-t0.tv_usec;
*/
  if(i>=CHRONOMAX) Erreur(26);
  
  /* variables globales, locale à la fonction */
  static int first=1; /* =1 ssi c'est la 1ère fois qu'on exécute la fonction */
  static char *str[CHRONOMAX];
  static struct timeval last[CHRONOMAX],tv;
  int j;

  if(i<0){ /* libère les pointeurs */
    if(!first) /* on a déjà alloué les chronomètres */
      for(j=0;j<CHRONOMAX;j++)
	free(str[j]);
    first=1;
    return NULL;
  }

  /* tv=temps courant */
  int err=gettimeofday(&tv,NULL);

  if(first){ /* première fois, on alloue puis on renvoie TopChrono(i) */
    first=0;
    for(j=0;j<CHRONOMAX;j++){
      str[j]=malloc(10); // assez grand pour "--error--", "99d99h99'" ou "23h59'59""
      last[j]=tv;
    }
  }

  /* t=temps en 1/1000" écoulé depuis le dernier appel à TopChrono(i) */
  long t=(tv.tv_sec-last[i].tv_sec)*1000L + (tv.tv_usec-last[i].tv_usec)/1000L;
  last[i]=tv; /* met à jour le chrono interne i */
  if((t<0L)||(err)) t=LONG_MAX; /* temps erroné */
  
  /* écrit le résultat dans str[i] */
  for(;;){ /* pour faire un break */
    /* ici t est en millième de seconde */
    if(t<1000L){ /* t<1" */
      sprintf(str[i],"0\"%03li",t);
      break;
    }
    t /= 10L; /* t en centième de seconde */
    if(t<6000L){ /* t<60" */
      sprintf(str[i],"%li\"%02li",t/100L,t%100L);
      break;
    }
    t /= 10L; /* t en dixième de seconde */
    if(t<36000L){ /* t<1h */
      sprintf(str[i],"%li'%02li\"%li",t/360L,(t/10L)%60L,t%10L);
      break;
    }
    t /= 10L; /* t en seconde */
    if(t<86400L){ /* t<24h */
      sprintf(str[i],"%lih%02li'%02li\"",t/3600L,(t/60L)%60L,t%60L);
      break;
    }
    t /= 60L; /* t en minute */
    if(t<144000){ /* t<100 jours */
      sprintf(str[i],"%lid%02lih%02li'",t/1440L,(t/60L)%24L,t%60L);
      break;
    }
    /* error ... */
    sprintf(str[i],"--error--");
  }
  
  return str[i];
}


/* type pour une fonction f(u,v,T) renvoyant la longueur d'une route
   de u vers v étant données la table de routage globale T */
typedef int(*rt_length)(int,int,void*);


/* structure de données pour la table de routage d'un sommet */
typedef struct{
  int n;      /* taille des listes */
  int *node;  /* liste de sommets */
  int *dist;  /* distances */
  int radius; /* rayon d'une boule */
  int vpd;    /* voisin par défaut */
} table;


/*

pour avoir plusieurs noms différents pour une même structre

typedef union{
  struct{
    int n;
    int *dist;
  };
  struct{
    int n;
    int *color;
  };
} tableau;

  tableau *T=malloc(sizeof(*T));
  T->n=0;
  T->dist=NULL;
  printf("->y=%p\n",T->color);

*/

table *new_table(const int n){
/*
  Crée un objet de type table où les champs sont initialisés à leur
  valeurs par défaut, les pointeurs étant initialisés à NULL. Si n>0,
  alors les pointeurs (->node et ->dist) de taille n sont alloués, et
  le champs ->n est initialisé à n.
*/
  NALLOC(table,X,1);
  X->n=0;
  X->node=NULL;
  X->dist=NULL;
  X->radius=-1;
  X->vpd=-1;

  if(n>0){
    X->n=n;
    ALLOC(X->node,n);
    ALLOC(X->dist,n);
  }
  
  return X;
}


void free_table(table *X){
  if(X==NULL) return;
  free(X->node);
  free(X->dist);
  free(X);
  return;
}


/* tables de routage globale pour rs_cluster() */
typedef struct{
  int n; // pour libérer les n tables
  table **B; // les boules
  table **S; // les spanning trees depuis les landmarks
  table **R; // les sommets en charge des landmarks
  table **W; // tables des voisins dans le cluster
  int *H; // hahs des sommets
  int *C; // pour routage vers voisins de couleur i
  int center; // centre du cluster (pour VARIANT=1)
} rs_cluster_tables;


void free_rs_cluster_tables(rs_cluster_tables *X){
/*
  Libère les tables créées par rs_cluster().
*/
  if(X==NULL) return;
  int u;
  for(u=0;u<X->n;u++){
    free_table(X->B[u]);
    free_table(X->S[u]);
    free_table(X->R[u]);
    free_table(X->W[u]);
  }
  free(X->B);
  free(X->S);
  free(X->R);
  free(X->W);
  free(X->H);
  free(X->C);
  free(X);
  return;
}


/* tables de routage globale pour rs_dcr() */

typedef struct{ // structure "boule-contiguë" spécifique au schéma AGMNT
  int d;    // distance entre u et v
  int v;    // sommet v=CONT[u][i]
  int s;    // landmark, <0 pour dire via boule-contiguë
  int w;    // nca entre u et v (si via landmark)
} contigue;

typedef struct{
  table **B; // B[u]=tables de voisinage indexées par les sommets
  table **W; // W[u]=tables des représentants indexées par les couleurs 
  param_bfs **S; // S[u]=bfs(u) pour u landmark (sinon =NULL)
  int *H; // H[u]=hash du sommet u
  int *C; // C[u]=couleur du sommet u
  int **dist; // distances partielles (issues des landmarks)
  int n; // pour libérer les n tables
  // pour AGMNT (=NULL sinon)
  contigue **CONT;
  int *F; // F[C[u]]=|CONT[u]|=#sommets de hash C[u]
} rs_dcr_tables;


int fcmp_contigue(const void *P,const void *Q)
/*
  Compare les champs .v d'une structure "contigue" (pour AGMNT).
*/
{
  const int p=(*(contigue*)P).v;
  const int q=(*(contigue*)Q).v;
  return (p>q) - (p<q);
}


void free_rs_dcr_tables(rs_dcr_tables *X){
/*
  Libère les tables créées par rs_dcr().  On a pas besoin de libérer
  les pointeurs de X->dist, car ce sont ceux de X->S->D.
*/
  if(X==NULL) return;
  int u;
  for(u=0;u<X->n;u++){
    free_table(X->B[u]);
    free_table(X->W[u]);
    free_param_bfs(X->S[u]); // libère aussi X->dist[u]
  }
  free(X->B);
  free(X->W);
  free(X->S);
  free(X->H);
  free(X->C);
  free(X->dist);
  FREE2(X->CONT,X->n); // pour AGMNT
  free(X->F); // pour AGMNT
  free(X);
  return;
}


/* tables de routage pour rs_tzrplg() */
typedef struct{
  table **B; // boules
  table **L; // landmarks
  int* label; // les étiquettes des sommets
  int n; // pour libérer les n tables
} rs_tzrplg_tables;


void free_rs_tzrplg_tables(rs_tzrplg_tables *X){
/*
  Libère les tables créees par rs_cluster().
*/
  if(X==NULL) return;
  int u;
  for(u=0;u<X->n;u++){
    free_table(X->B[u]);
    free_table(X->L[u]);
  }
  free(X->B);
  free(X->L);
  free(X);
  return;
}


/* tables de routage pour rs_bc() */
typedef struct{
  int nbfs; // taille de L
  param_bfs **L; // liste des BFS;
  param_bfs **Lu; // liste des BFS indexé par les sommets
  int **dist; // distances partielles (issues des landmarks)
} rs_bc_tables;


void free_rs_bc_tables(rs_bc_tables *X){
/*
  Libère les tables créees par rs_bc().
*/
  if(X==NULL) return;
  int u;
  for(u=0;u<X->nbfs;u++)
    free_param_bfs(X->L[u]); // libère X->Lu et aussi X->dist

  free(X->L);
  free(X->Lu);
  free(X->dist);

  free(X);
  return;
}


/* tables de routage pour rs_hdlbr() */
typedef struct{
  table **B; // boules
  table **L; // landmarks
  int* Core; // liste des landmarks
  int* H;    // hash des sommets
  int n;     // pour libérer les n tables
  int core_size;
} rs_hdlbr_tables;


void free_rs_hdlbr_tables(rs_hdlbr_tables *X){
/*
  Libère les tables créees par rs_hdlbr().
*/
  if(X==NULL) return;
  int u;
  for(u=0;u<X->n;u++){
    free_table(X->B[u]);
    free_table(X->L[u]);
  }
  free(X->B);
  free(X->L);
  free(X->Core);
  free(X->H);

  free(X);
  return;
}


/* écart type calculé à partir de la somme s1 et de la somme s2 du
   carré de n valeurs. Si n<=0, alors on retourne NaN. On utilise la
   formule écart_type(X) := sqrt(E(X^2)-E(X)^2). Attention d'utiliser
   (s1/n)*(s1/n) plutôt que (s1*s1)/(n*n) à cause du dépassement
   arithmétique du type int de n*n. NB: nan(NULL) provoque un warning
   sur certains systèmes.
*/
#define ECARTYPE(s1,s2,n) ((n>0)?sqrt(s2/(double)n-(s1/(double)n)*(s1/(double)n)):nan("NaN"))


/* affiche le min/max et la moyenne et l'écart type d'un tableau
   d'entiers vérifiant une certaine condition. Si la condition n'est
   jamais vérifiée (moyenne sur aucun terme), alors les min/max et
   moyenne sont évalués à 0 et on ajoute le texte "(undefined)".
*/

// Ex:
//     MINMAXMOY(T[_i],n,T[_i]>0,"size of T");
//       Calcule la moyenne des termes T[_i] tq T[_i]>0
//       pour i=0..n-1. Le texte sert pour l'affichage.
//
#define MINMAXMOY(term,n,condition,string)				\
  do{									\
    int _i,m0,m1,cpt;long t,sum=0L,sum2=0L;				\
    for(_i=cpt=0;_i<(n);_i++)						\
      if(condition){							\
        t=(long)term,sum+=t,sum2+=t*t;					\
	if(cpt) m0=imin(m0,t),m1=imax(m1,t);				\
	else m0=m1=t;							\
	cpt++;								\
      }									\
    if(cpt==0) m0=m1=0;							\
    double avg=cpt?(double)sum/(double)cpt:0;				\
    printf("- minimum %s: %i%s\n",string,m0,cpt?"":" (undefined)");	\
    printf("- maximum %s: %i%s\n",string,m1,cpt?"":" (undefined)");	\
    printf("- average %s: %.2lf",string,avg);				\
    if(cpt) printf(" ± %.2lf (%li/%i)\n",				\
		   ECARTYPE(sum,sum2,cpt),sum,cpt);			\
    else printf(" (undefined)\n");					\
  }while(0)


/* calcule puis affiche la fréquence, le min et le max d'un tableau
   d'entiers */

// Ex:
//     FREQMINMAX(F,k,T,n,"size of T");
//       Calcule dans F la fréquence des éléments du tableau T
//       à n éléments, les valeurs de T étant des entiers de [0,k[.
//       Attention ! le tableau F est alloué par et initialisé par
//       la macro, et c'est la macro qui détermine k qui doit être
//       une variable. Le texte sert pour l'affichage.
//
#define FREQMINMAX(F,k,T,n,string)					\
  do{									\
     int u,v,i;								\
     F=SortInt(T,NULL,n,0,&k,SORT_FREQv);				\
     for(u=v=F[i=0];i<k;i++) u=imin(u,F[i]),v=imax(v,F[i]);		\
     printf("- balance ratio of the frequency %s: %.0lf%% (ceil{n/k}=%i max=%i min=%i)\n", \
	      string,100.0*u/(double)v,iceil(n,k),v,u);			\
     }while(0)

void ruling(char *s,int n){
/*
  Affiche n fois la chaîne de caractère s
*/
  int i;
  for(i=0;i<n;i++) printf("%s",s);
}
#define BARRE do{ruling("―",LEN_BARRE);printf("\n");}while(0)
#define RULING(x) ruling(STR_DISTRIB,imin((int)(LEN_DISTRIB*x),LEN_DISTRIB-5))


void PrintDistribution(const int *Z,const int n,int r,const char *message){
/*
  Affiche la distribution des valeurs entières contenues dans le
  tableau Z (de taille n>0) selon r>0 intervalles et avec le message
  m. Les intervalles consécutifs nuls sont regroupés. Donc le nombre
  d'intervalles affichés peut être < r. Si r<0, alors on affiche
  toutes les valeurs (de fréquence non nulle) plutôt que des
  intervalles. En plus de la distribution, on affiche la valeur min,
  max et la moyenne avec l'écart type.

  Exemple avec r=10 et m="routing table size":

  - routing table size distribution: 10 ranges
    [339,381[ 	01% ▪ [× 168] 
    [381,423[ 	08% ▪▪▪▪ [× 698] 
    [423,465[ 	15% ▪▪▪▪▪▪▪▪▪ [× 1282] 
    [465,507[ 	57% ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ [× 4848] 
    [507,549[ 	11% ▪▪▪▪▪▪ [× 975] 
    [549,591[ 	02% ▪ [× 230] 
    [591,633[ 	02% ▪ [× 182] 
    [633,675[ 	00%  [× 18] 
    [675,717[ 	00%  [× 4] 
    [717,752[ 	00%  [× 1] 
  - minimum routing table size: 339
  - maximum routing table size: 751
  - average routing table size: 481.48 ± 13.84 (4047360/8406)
*/
  if((n<=0)||(Z==NULL)){ printf("- empty distribution\n"); return; }
  int i,j,k,m0,m1; /* calcule m0=min_i{Z[i]} et m1=max_i{Z[i]} */
  const int b=(r<0);
  if(b) r=n;
  m0=m1=Z[0]; 
  for(i=1;i<n;i++) m0=imin(m0,Z[i]),m1=imax(m1,Z[i]);
  int d=(++m1)-m0; /* d=|[m0,m1[| = plage des valeurs de Z, d>=1 */
  r=imin(r,d); /* r=nombre d'intervalles maximum, ne peut dépasser d */
  d=iceil(d,r); /* d=ceil(d/r)=taille des intervalles sauf le dernier
		   qui est plus court. NB: d>=1 car d>=r>=1 */
  NALLOCZ(int,F,r,0); /* initialise F */
  for(i=0;i<n;i++) F[(Z[i]-m0)/d]++; /* compte les valeurs de chaque intervalle */
  for(i=k=0;i<r;k++,i=j){ /* k=nombre d'intervalles (ou valeurs) qui vont être affichés */
    j=i+1;
    if(F[i]==0){
      while((j<r)&&(F[j]==0)) j++; /* j=prochain intervalle */
      if(j==r) break;
    }
  }
  printf("- %s distribution: %i %s%s\n",message,k,b?"value":"range",PLURIEL(k));

  /* affichage des intervalles (ou valeurs), éventuellement en
     fusionnant les intervalles vides consécutifs */

  double x;
  for(i=0;i<r;){
    j=i+1;
    if(F[i]==0){
      while((j<r)&&(F[j]==0)) j++; /* j=prochain intervalle */
      if(j==r) break;
    }
    x=(double)F[i]/(double)n;
    if(b){ if(F[i]) printf("  %i: \t",m0+i*d); }
    else printf("  [%i,%i[ \t",m0+i*d,imin(m0+j*d,m1));
    printf("%02i%% ",(int)(100*x));
    RULING(x);
    printf(" [× %i]\n",F[i]);
    i=j; /* ici j>=i+1 */
  }
  free(F); /* ne sert plus à rien */

  /* affichage min/max/moy écart type pour Z */
  MINMAXMOY(Z[_i],n,1,message);
  return;
}


char *millier(const long i){
/*
  Renvoie une chaîne de caractères correspondant à l'entier long i
  écrit avec le séparateur des milliers ','. Ne pas faire de free()
  sur le pointeur renvoyé.

  Ex: > printf("n=%s\n",millier(-123456789));
      > n=-123,456,789
*/  
  static char r[64];
  char *p,*s=r; /* r=résultat final, s=r ou r+1 si i<0 */
  int n;
  
  snprintf(r,sizeof(r),"%li",i); /* écrit i dans r */
  if(i<0) s++; /* si i<0, r="-123..." */
  n=strlen(s); /* n=longueur(r) sans le signe */
  p=s+n-3;
  
  while(p>s){
    memmove(p+1,p,n);
    *p=',';
    n+=4; /* 4 chiffres avec la virgule */
    p-=3; /* paquet de 3 chiffres */
  }
  
  return r;
}


int lg(unsigned n){
/*
  Renvoie le plus petit entier k tel que n < 2^k. Autrement dit k est
  tel que 2^{k-1} <= n < 2^k. D'où:

            / 1+floor{log_2(n)} si n>0     n  | 0 1 2 3 4 5 6 7 8 9 ...
   lg(n) := |                            -----o------------------------
            \ 0 sinon                    lg(n)| 0 1 2 2 3 3 3 3 4 4 ...

  Lorsque n>0, c'est aussi le nombre de bits dans l'écriture binaire
  de n, et c'est également un de plus que la position du bit de poids
  fort.

  NB: on peut calculer lg(lg(lg(n))) sans erreur, car lg(n)≥0 est
  défini pour tout entier n≥0.
*/
  int k=0;
  while(n>0) n>>=1,k++;
  return k;
}


/*
  Calcule rapide d'un hash sur 32-bits, d'après Bob Jenkins (2006),
  disponible ici: http://burtleburtle.net/bob/c/lookup3.c
*/
#define rot(x,k) (((x)<<(k)) | ((x)>>(32-(k))))
#define mix(a,b,c) \
{ \
  c ^= b; c -= rot(b,14); \
  a ^= c; a -= rot(c,11); \
  b ^= a; b -= rot(a,25); \
  c ^= b; c -= rot(b,16); \
  a ^= c; a -= rot(c,4);  \
  b ^= a; b -= rot(a,14); \
  c ^= b; c -= rot(b,24); \
}


int hash_mix(const int u){
/*
  Renvoie un entier h(u) positif (sur 31 bits donc) basé sur le
  mélange de deux constantes entières aléatoires A et B de 32 bits. Si
  u<0, alors les constantes A,B sont initialisées et on renvoie -1.
*/
  static unsigned A,B;

  if(u<0){
    A=random();
    B=random();
    return -1;
  }

  unsigned a=A,b=B,c=u; // NB: a et b sont modifiées
  mix(a,b,c);

  return c>>1; // valeurs sur 31 bits, toujours positive
}


int hash_prime(const int u){
/*
  Renvoie un entier h(u) de [0,p[ où p=2^31-1=0x7FFFFFFF est un nombre
  premier (on pourrait aussi utiliser 2^61-1 qui est premier). Il faut
  0<=u<=p. Si u<0, alors les constantes A,B (qui accélèrent les
  calculs) sont initialisées et on renvoie -1.  On prend la fonction
  de hashage de Carter & Wegman (1978) avec h(u)=(A*u+B) mod p (puis
  h(u)%k) où A,B sont dans [0,p[ et A impair. Il est important de
  faire les calculs en "long unsigned" car la valeur A*x+B dépasse 32
  bits en général. Sinon, le calcul de (A*x+B)%p est erroné (s'il est
  réalisé en 32 bits).
*/
  static const long unsigned p=0x7FFFFFFF;
  static long unsigned A,B;
  const long unsigned x=u;

  if(u<0){
    A=(1+random())%p;
    A|=1; // force A à être impair
    B=random()%p;
    return -1;
  }

  return (A*x+B)%p; // en principe, le modulo ne sert à rien
}


int hash_shuffle(int u,const int n){
/*
  Calcule un entier p(u) de [0,n[ qui est une permutation
  p:[0,n[->[0,n[ basée sur deux constantes aléatoires R1,R2 de
  [0,n[. Le temps de calcul est d'environ (logn)/2. Si u<0, alors les
  constantes R1,R2,K,N0,N1 (qui accélèrent les calculs) sont
  initialisées et on renvoie -1. Cette permutation est basée sur deux
  "shuffles" qui sont les permutations P0 et P1 suivantes:

    P0(u)=(u>>1)+((u&1)==0)*N0; // avec N0=floor(n/2)
    P1(u)=(u>>1)+((u&1)==1)*N1; // avec N1=ceil(n/2)

  Pour n=10, la permutation P0(u) donne: 0-5-1-6-2-7-3-8-4-9. Donc on
  écrit 0,1,2,... aux positions paires 0,2,4,6,...

  Pour n=10, la permutation P1(u) donne: 5-0-6-1-7-2-8-3-9-4. Donc on
  écrit 0,1,2,... aux positions impaires 1,3,5,6,...

  On applique la permutation P0 ou P1 suivant les bits de poids faible
  de R1, en répétant cette opération K fois, où K=(logn)/2 est environ
  la moitié des bits nécessaires pour écrire n en binaire. En effet,
  on retrouve la valeur de départ après avoir effectué logn fois la
  permutation P0 ou logn fois la permutation P1.

  Enfin, on inverse le résultat (u->n-u) et on le décale de R2
  position modulo n (u->(u+R2)%n). Ainsi la permutation finale de u,
  p(u), peut valoir n'importe quelle valeur entre [0,n[.
*/

  /* variables globales / locales */
  static int R1=0;
  static int R2=0;
  static int N0=0;
  static int N1=0;
  static int K=0;

  if(u<0){
    R1=random()%n; /* R1 contient au moins k bits aléatoires */
    R2=random()%n; /* pour le décalage final */
    N0=N1=(n>>1);  /* N0=floor(n/2) */
    N1 += (n&1);  /* N1=ceil(n/2) */
    K=lg(n)/2;     /* K=la moitié de floor(logn) */
    return -1;
  }
  
  /* calcule u=p(u) */
  int i,b,r=R1;
  for(i=0;i<K;i++){
    b=(r&1); /* lit un bit b de R1 et applique P0 ou P1 */
    u=(u>>1)+((u&1)==b)*(b?N1:N0); /* u/2, u/2+N0 ou u/2+N1 */
    r>>=1; /* supprime le bit lu */
  }

  /* inverse et décalage aléatoire */
  return (n+R2-u)%n;
}


int *MakeHash(int *H,const int n,const int k,const int M){
/*
  Calcule un hash h:[0,n[->[0,k[ avec k<=n selon la méthode M (voir la
  variable HASH). Remplit et renvoie le tableau H de taille n (ou bien
  un nouveau tableau alloué si H=NULL). Si M n'est pas une valeur
  reconnue, on fait comme si M=H_MOD. On renvoie NULL si n<=0. En
  général, il faut éviter de faire un modulo (%) qui est jusqu'à 200
  fois plus lent qu'une addition entière.

  Il y a d'autres fonctions de hash très intéressantes à voir (comme
  http://burtleburtle.net/bob/hash/doobs.html). On y trouve en
  particulier, lookup3.c avec mix(a,b,c) de Bob Jenkins 2006
  http://burtleburtle.net/bob/c/lookup3.c. Il y a aussi MurmurHash2
  (2010).
*/

  if(n<=0) return NULL;
  if(H==NULL) ALLOC(H,n);
  
  int u,a;
  
  switch(M){
    
  case H_PRIME:
    hash_prime(-1); /* initialisation */
    for(u=0;u<n;u++) H[u]=hash_prime(u)%k;
    break;
    
  case H_SHUFFLE:
    hash_shuffle(-1,n); /* initialisation */
    for(u=0;u<n;u++) H[u]=hash_shuffle(u,n)%k;
    break;
    
  case H_MIX:
    hash_mix(-1); /* initialisation */
    for(u=0;u<n;u++) H[u]=hash_mix(u)%k;
    break;
    
  case H_MOD:
  default:
    a=random()%k; // attention que u+a reste >=0
    for(u=0;u<n;u++) H[u]=(u+a)%k;
    break;
  }
  
  return H;
}


rs_cluster_tables *rs_cluster(graph *G,int k){
/*
  Calcule pour le graphe G le routing scheme cluster de paramètre
  k>0. On renvoie les tables calculées. Le stretch est toujours <=
  5. Cependant lorsque k=1 le stretch devient <= 3.

  Si bit-0 de VARIANT=1, alors les tables B sont vides.
  Si bit-1 de VARIANT=1, alors les tables W sont vides.
  Si bit-2 de VARIANT=1, alors les tables B des voisins du cluster sont vides.
  Si bit-3 de VARIANT=1, alors les tables B des voisins du cluster n'ont pas de sommets du cluster.
*/
  
  const int n=G->n;
  param_bfs *X0,*X;
  int u,v,i,j,r,t,center;
  table **S,**B,**W,**R;
  int *Z,*F;

  printf("\nCLUSTER\n");
  BARRE;
  
  TopChrono(1); /* reset du chrono tmp */
  TopChrono(2); /* reset du chrono total */

  /* trouve un sommet center de degré max, v=deg(center) */
  for(u=v=center=0;u<n;u++)
    if(G->d[u]>v) v=G->d[center=u];

  if((VARIANT&2)&&(k<v+1)) Erreur(6); // paramètres incohérents
  k=imin(k,v+1); /* k=taille effective du cluster=#colors */
  printf("- cluster size: %i",k);
  printf(" (ceil{√n}=%i)\n",(int)ceil(sqrt(n)));
  printf("- degree of the center: %i (id:%i)\n",v,center);

  /* construit le cluster C de taille k à partir des k-1 sommets de
     plus haut degré voisins du center */

  /* construit une table H des degrés des voisins de center */
  /* on en profite pour calculer le degré min pour optimiser SortInt() */
  u=v; /* u=futur degré min, v=deg(center) */
  NALLOCZ(int,H,v,G->d[G->L[center][_i]]); /* H[i]=degré du i-ème voisin du center */
  for(i=0;i<v;i++) u=imin(u,H[i]); /* u=voisin de center de deg min */

  /* trie selon les listes des degrés H: Z[i]=indice dans H du degré de rank i */
  r=v-u+1; /* H[i] valeur dans [u,v]=[u,u+r[ */ 
  Z=SortInt(H,NULL,v,u,&r,SORT_INDEXi); /* H[Z[v-1]]=+grand degré parmi les v voisins du center */
  free(H); /* H ne sert plus à rien, seul Z sert encore */

  /* remplit C selon les degrés décroissant */
  NALLOC(int,C,k);
  C[0]=center; /* met le center dans C, forcément de plus haut degré */
  for(i=1;i<k;i++) C[i]=G->L[center][Z[--v]]; /* les k-1 voisins de + haut degré */
  free(Z); /* Z ne sert plus à rien */
  
  /* affiche un aperçu des degrés du cluster */
  printf("- cluster degree: ");
  APERCU(G->d[C[_i]],k,10,2);

  /* on trie C pour déterminer la couleur de u dans C rapidement. La
     couleur de u dans C est alors i=SetSort(C,k,1). Attention ! la
     couleur du center n'est pas forcément 0. */
  QSORT(C,k,fcmp_int);
  
  X0=bfs(G,center,NULL); /* calcule un BFS depuis center */
  printf("- eccentricity of the center: %i\n",X0->radius);
  printf("- time to construct the cluster with its BFS: %s\n",TopChrono(1));

  /*
    Calcule, pour chaque sommet u de C, le nombre de sommets qui ne
    sont pas dans C et qui ont u comme plus proche ancêtre. On calcule
    ensuite le maximum de ces nombres.

    On utilise un simple tableau H indiquant si un sommet a déjà été
    visité ou pas. De plus, si c'est un sommet de C, H[u] indique
    combien de sommets lui sont directement descendant. Au départ
    H[u]=1 si u est dans C, et H[u]=-1 sinon. Puis, pour chaque sommet
    u qui n'est pas déjà visité ou dans C, on marque u et on remonte
    dans l'arbre u=parent(u) jusqu'à trouver soit un sommet de l'arbre
    déjà visité ou bien un sommet de C. Pour marqué un sommet u que
    l'on visite, on pose H[u]=0. Si on tombe sur un sommet du cluster,
    on incrémente H[u]. Donc H[u]>0 si u est dans C, H[u]<0 si u n'a
    pas été visité, et H[u]==0 sinon (pas dans C mais déjà
    visité). Cet algo prend un temps O(n) car chaque arête de l'arbre
    n'est visitée qu'une seule fois.
   */
  ALLOCZ(H,n,-1);
  for(i=0;i<k;H[C[i++]]=1);
  for(u=0;u<n;u++){
    v=u;
    while(H[v]<0){ /* v n'est pas marqué */
      H[v]=0; /* on marque v */
      v=X0->P[v]; /* v=parent(v) */
    }
    if(H[v]>0) H[v]++;
    }
  j=0; /* cherche u de C avec H[u] maximum */
  for(i=1;i<k;i++) if(H[C[i]]>H[C[j]]) j=i;
  printf("- maximum number of cluster direct descendants: %i (id: %i)\n",H[C[j]],C[j]);
  free(H);

  /* efface les pères des sommets de C pour test d'appartenance rapide */
  for(i=0;i<k;i++) X0->P[C[i]]=-2; /* si -2 alors dans C */

  DEBUG(PRINTT(C,k););

  /****************************/
  /* calcule des stats / hash */
  /****************************/

  /* calcule le nombre d'arêtes intra et extra C, de voisins de C */
  TopChrono(1);
  ALLOCZ(H,n,0);
  t=r=0;
  for(i=0;i<k;i++){
    u=C[i]; // u=i-ème sommet de C
    for(j=0;j<G->d[u];j++){ /* pour chaque voisin de u */
      v=G->L[u][j]; /* v=j-ème voisin de u */
      if(X0->P[v]==-2) t++; // arête dans C
      else{ r++; H[v]=1; } // arête hors C; marque les sommets
    }
  }
  printf("- #edges inside the cluster: %i\n",t/2);
  printf("- #edges outgoing the cluster: %i\n",r);
  for(u=t=0;u<n;u++) t+=H[u];
  printf("- #neighbors of the cluster: %i\n",t);
  printf("- their average #neighbors in the cluster: %.2lf (%i/%i)\n",(double)r/(double)t,r,t);
  
  /* construit le hash H[u]=0..k-1 de chaque sommet u=0..n-1 */
  /* le sommet C[i] est en charge des sommets de hash i=0..k-1 */

  MakeHash(H,n,k,HASH);
  FREQMINMAX(F,k,H,n,"the hash");
  DEBUG(PRINTT(H,n);PRINTT(F,k););
  printf("- time to compute stats and hash: %s\n",TopChrono(1));
  BARRE;
  
  /************************/
  /* calcule les tables B */
  /************************/
  
  /*
    Tables seulement définies pour les sommets u qui ne sont pas dans
    C. Il s'agit dans un premier temps des sommets à distance < r de u
    où r est la distance du plus proche sommet de C. Cependant, si le
    plus proche sommet est center (cela ne peut pas se produire si C
    contient le centre et tous ces voisins, en particulier si k=n), il
    convient alors de vérifier si à distance r-1 de u il n'y a pas un
    sommet de C. Si oui, il faut alors prendre r-1 et non
    r. Eventuellement tree(center) peut-être modifié. Le sommet u est
    supprimé de B[u]. Dans un second temps on ajoute tous les voisins
    de u à B[u] si |B[u]| était vide, c'est-à-dire si u était voisin
    de C. Le choix par défaut (->vpd) est le voisin de u allant vers
    le sommet de B[u] le plus proche du center (dans tree(center)).
  */

  TopChrono(1);
  ALLOCZ(B,n,NULL); /* tableau de n tables B (vides au départ) */
  X=new_param_bfs(); /* pour ne pas le faire à chaque sommet u */
  X->clean=1; /* initialisation complète de X->D, puis partielle */

  int bmax,bsum; /* stat pour les "vraie" boules */
  bmax=bsum=0;
  
  for(u=0;u<n;u++){ /* pour chaque sommet u du graphe */

    /* Cherche d'abord le rayon r de B[u] pour calculer ensuite plus
       rapidement B[u] avec un BFS de profondeur r. Pour cela on
       remonte dans tree(center) jusqu'à trouvé un sommet de C.
       Cependant, ce n'est pas forcément le plus proche de u. Après le
       BFS depuis u, il faut vérifier si dans la dernière couche du
       BFS on a pas un sommet de C. Si oui, on décrémente le rayon et
       on modifie tree(center) de sorte à trouver un sommet de C
       plutôt et éviter de futurs problèmes. Dans la variante "boules
       vides" on ne touche pas à l'arbre. */

    if(VARIANT&1){ // variante "boules vides"
      if(X0->P[u]==-2) continue; // on ne fait rien si u est dans C
      B[u]=new_table(0); /* on crée une table B[u] vide, NB: B[u]->n=0 */
      B[u]->vpd=X0->P[u]; // voisin par défaut est le père de u
      continue;
    }

    t=v=u; r=-1;
    while(X0->P[v]!=-2){ /* si v pas dans C */
      v=X0->P[t=v]; /* t=v et v=parent(t) dans l'arbre X0->P */
      r++;
    }
    /* Ici:
       v=1er sommet ancêtre dans tree(center) de u qui est dans C
       t=dernier sommet de B[u] avant C (=fils de v dans tree(center))
       r=dist(u,t)=dist(u,v)-1
    */
    if(r<0) continue; /* rien à faire, on était parti d'un sommet de C */
    
    B[u]=new_table(0); /* on crée la table B[u] vide, NB: B[u]->n=0 */
      
    if(r>0){ /* ne rien faire ici si r=0 (=si u voisin de C) */
      X->hmax=r; bfs(G,u,X); /* calcule un BFS depuis u de profondeur r */
      i=B[u]->n=(X->n)-1; /* taille de B[u] sans u, si r est correct */
      if(v==center){ /* il faut vérifier si u n'a pas de sommet v dans C à distance r */
	/* on part des sommets i les plus éloignés de u, donc à distance r */
	/* NB: il y a toujours dans X->D au moins un sommet à distance < r, c'est u */
	while((X->D[X->file[i]]==r)&&(X0->P[X->file[i]]!=-2)) i--;
	if(X0->P[X->file[i]]==-2){ /* on a trouvé un sommet dans C à distance r */
	  v=X->file[i]; /* v=sommet dans C le plus proche de u */
	  t=X->P[v]; /* t=dernier sommet de B[u] avant C (=père de v dans le BFS de u) */
	  X0->P[t]=v; /* modifie le père de t dans tree(center) pour accélérér les tests suivants */
	  while(X->D[X->file[i]]==r) i--; /* on peut effacer les sommets de i+1 à X->n-1 */
	  B[u]->n=i; /* la taille de B[u] est i+1 avec u, donc i sans u */
	  r--; /* r=dist(u,v)-1=dist(u,t) */
	}
      }
    }
    
    /* Ici:
       v=le plus proche sommet de u qui est dans C
       t=dernier sommet de B[u] avant C (=fils de v dans tree(center))
       r=dist(u,t)=dist(u,v)-1>=0
       B[u]->n=(taille de boule de rayon r sans u)-1
     */

    /* pour les stats avant ajoût des voisins et optimisation */
    bmax=imax(bmax,1+B[u]->n);
    bsum += (1+B[u]->n);

    if(r==0){ /* si u voisin de C, alors B[u] contient certains voisins de u mais pas u */
      B[u]->vpd=v; /* v=voisin de u dans C */
      if(VARIANT&4){ /* alors B[u] = boule vide (sans u) */
	B[u]->n=B[u]->radius=0; // NB: B[u]->node et B[u]->dist sont à NULL par new_table(0)
	continue;
      }
      if(VARIANT&8){ /* alors B[u] = voisins de u pas dans C */
        ALLOC(B[u]->node,G->d[u]); /* taille deg(u) pour avoir suffisamment de place */
        B[u]->n=0; // sert d'indice pour B[u]->node
        for(i=0;i<G->d[u];i++){
          v=G->L[u][i]; // v=voisin de u
          if(X0->P[v]!=-2) B[u]->node[B[u]->n++]=v; // si v pas dans C, on l'ajoute à B[u]
        }
        B[u]->radius=(B[u]->n>0); // =1 ssi au moins un sommet dans B[u]
        REALLOC(B[u]->node,B[u]->n); // redimensionne B[u]
      }
      else{
	B[u]->n=G->d[u]; /* |B[u]|=deg(u) */
        B[u]->radius=1; /* rayon de B[u], car deg(u)>0 */ 
        ALLOCZ(B[u]->node,B[u]->n,G->L[u][_i]); /* copie tous les voisins de u dans B[u]->node */
      }
      ALLOCZ(B[u]->dist,B[u]->n,1); /* les sommets de B[u] sont tous à distance 1 */
      continue; /* on passe au sommet u suivant */
    }

    /* ici r>0, et donc B[u]->n>0, et t<>u */
    /* on cherche le voisin par défaut de u, cad le fils de u dans
       BFS(u) qui mène à t. Là, on pourrait tester le père dans X0
       pour éventuellement le corriger et éviter les problèmes de
       correction de r en r-1 comme plus haut. */
    while(X->P[t]!=u) t=X->P[t];
    B[u]->vpd=t; /* voisin de u par défaut, ne peut pas être dans C */
    B[u]->radius=r; /* rayon de la boule */
    
    /* copie les sommets (et leur distance) du BFS couvrant B[u] en supprimant u */
    ALLOCZ(B[u]->node,B[u]->n,X->file[_i+1]); /* B[u]->node[i]=i-ème sommet de B[u] */
    ALLOCZ(B[u]->dist,B[u]->n,X->D[B[u]->node[_i]]); /* B[u]->dist[i]=distance du i-ème à u */
  }
  printf("- time to construct tables B: %s\n",TopChrono(1));
  printf("- maximum non-cluster ball size: %i\n",bmax);
  printf("- average non-cluster ball size: %.2lf (%i/%i)\n",
	 (double)bsum/(double)n,bsum,n);
  MINMAXMOY(B[_i]->n,n,B[_i],"table B size");

  /*************************/
  /* optimise les tables B */
  /*************************/
  
  /*
    Dans un premier temps on met des -1 aux sommets à supprimer.
    Ensuite on réorganise les tables B en gardant l'ordre des sommets
    restant. NB: on peut avoir B[u]<>NULL et B[u]->n=0. B[u]->radius
    n'est pas mis à jour et représente toujours la distance entre u et
    le sommet de B[u] le plus loin avant "clean". Il est important
    qu'avant l'appel, les sommets de B[u] soient rangés dans l'ordre
    du bfs(), donc non triées.
  */

  X->clean=1; /* initialisation complète de X->D, puis partielle */
  for(u=0;u<n;u++){
    if(B[u]){ /* il faut que la table B[u] existe */
      /* si deg[u]<=1 */
      if(G->d[u]<2) t=0;
      else{
	/* si deg(u)>1 */
	v=B[u]->vpd; /* v=voisin par défaut */
	X->hmax=(B[u]->radius)-1;
	bfs(G,v,X); /* calcule un BFS depuis v de profondeur rayon de B[u] -1 */
	for(i=0;i<X->n;i++){ /* passe en revue les sommets du bfs(v...) */
	  t=X->file[i]; /* t=sommet du bfs(v,...) */
	  r=SetSearch(t,B[u]->node,B[u]->n,0); /* r=indice tq B[u]->node[r]=t */
	  if((r>=0)&&(B[u]->dist[r]==1+X->D[t])) B[u]->node[r]=-1; /* supprime t */
	}
	/* supprime les -1 de B[u] en décalant B[u]->node et B[u]->dist */
	for(i=t=0;i<B[u]->n;i++)
	  if(B[u]->node[i]>=0){
	    B[u]->node[t]=B[u]->node[i];
	    B[u]->dist[t]=B[u]->dist[i];
	    t++;
	  }
      }
      /* redimensionne les tables à t = nouvelle taille de B[u] */
      B[u]->n=t;
      REALLOC(B[u]->node,t);
      REALLOC(B[u]->dist,t);
    }
  }
  printf("- time to clean tables B: %s\n",TopChrono(1));
  MINMAXMOY(B[_i]->n,n,B[_i],"table B size");

  DEBUG(for(u=0;u<n;u++)
	  if(B[u]){
	    printf("u=%i B[u]->vpd=%i B[u]->radius=%i ",u,B[u]->vpd,B[u]->radius);
	    if(B[u]) PRINTT(B[u]->node,B[u]->n);
	  });
  
  /************************/
  /* calcule les tables S */
  /************************/

  /*
    Tables seulement définies pour les sommets u de C sauf center (car
    le centre est la racine de l'arbre et donc déjà sur le plus court
    chemin).  S[u] est la liste des sommets de hash i où u=C[i] et qui
    sont à distance au plus r=min(2,logn/loglogn) de u. On ne met ni u
    ni center dans S[u]. NB: Il est possible d'avoir dans S[u] des
    sommets de C (et même tous sauf u et center).
  */
  
  r=imax(2,lg(n));
  r=ceil((double)r/(double)lg(r));
  const int depthS=imax(2,r); /* peut-être > hauteur(tree(center)) */

  ALLOCZ(S,n,NULL); /* tableau de n tables S (vides au départ) */
  X->hmax=depthS; /* profondeur max du BFS */
  X->clean=1; /* initialisation complète de X->D, puis partielle */

  for(i=0;i<k;i++){ /* pour chaque sommet u du cluster sauf center */
    u=C[i]; if(u==center) continue;
    bfs(G,u,X); /* calcule un BFS depuis u<>center de profondeur depthS */
    t=imin(X->n-1,F[i]-(H[u]==i)-(H[center]==i)); /* t=nombre max de sommets dans S[u] */
    S[u]=new_table(t); /* table S pour u */
    S[u]->n=0; /* on se sert de S[u]->n comme indice */
    /* on va parcourir les sommets du bfs() de profondeur <= depthS et
       ne garder que ceux de hash = i (ici i = hash(u) = hash(C[i]) */
    S[u]->vpd=-1;
    /* S[u]->vpd est l'indice (dans S[u]->node) du 1er sommet v de
       S[u] qui est à une distance >= depthS-1 de u.  Dans S[u], les
       sommets sont rangés selon le bfs(), donc en fonction de leur
       distance à u. On s'en sert plus tard pour accélérer la
       suppression dans R[u] des sommets v déjà dans S[u] (lire
       commentaires sur le remplissage de R[u]). NB: si aucun sommet
       n'est de hash = i (mauvaise fonction de hashage !), alors
       S[u]->vpd<0. */
    for(j=1;(j<X->n)&&(S[u]->n<t);j++){ /* pour chaque v du bfs(), v<>u */
      /* lorsque S[u] contient t sommets, on peut s'arrêter. */
      v=X->file[j]; /* v=j-ème sommet de la file, v<>u */
      if((H[v]==i)&&(v!=center)){ /* v a le bon hash et v<>center */
	S[u]->node[S[u]->n]=v;
	S[u]->dist[S[u]->n]=X->D[v];
	if((X->D[v]>=depthS-1)&&(S[u]->vpd<0))
	  S[u]->vpd=S[u]->n; /* 1er sommet de hash i à distance >= depthS-1 */
	S[u]->n++; /* un sommet de plus dans S[u] */
      }
    }
    if(S[u]->n){ /* redimensionne les tables de S[u] pour pas gaspiller */
      REALLOC(S[u]->node,S[u]->n);
      REALLOC(S[u]->dist,S[u]->n);
    }else{ /* si aucun sommet dans S[u], on efface S[u] */
      free_table(S[u]);
      S[u]=NULL;
    }
  }
  free_param_bfs(X); /* plus de BFS à faire, X ne sert plus à rien */
  printf("- time to construct tables S: %s (radius %i)\n",TopChrono(1),r);
  
  DEBUG(
	for(i=0;i<k;i++){
	  u=C[i];
	  printf("u=%i hash=%i |S[u]|=%i depth=%i\n",u,i,S[u]?S[u]->n:0,depthS);
	  if(S[u]){
	    PRINTT(S[u]->node,S[u]->n);   
	    PRINTT(S[u]->dist,S[u]->n);
	  }
	}
	);
  
  /************************/
  /* calcule les tables R */
  /************************/
  
  /*
    Tables seulement définies pour les sommets u de C avec u=C[i].
    R[u] est la liste des sommets dont le hash est i et qui ne sont
    pas déjà dans S.  On ne met ni u ni center dans R[u]. C'est donc
    comme S[u] mais sans la restriction de distance. Pour les sommets
    de R[u], on appliquera le routage dans tree(center).
  
    Les sommets de hash=couleur(center) sont tous dans R[center] (sauf
    center lui-même si hash(center)=couleur(center)). Si u<>center,
    alors pour qu'un sommet v (de hash i) soit dans R[u] il faut qu'il
    soit à une distance >= depthS de center. S'il est à distance =
    depthS-1 (ou inférieure), alors en passant par center on obtient
    une route de u à v de longueur <= (depthS-1)+1=depthS ce qui n'est
    pas mieux que d'utiliser S[u]. Le sommet v sera donc déjà dans
    S[u]. Maintenant, les sommets candidats v à distance >= depthS de
    center ne peuvent être qu'à distance au moins depthS-1 de u. S'ils
    étaient à distance < depthS-1, alors ils seraient à distance <
    depthS de center. Or, ils sont déjà à distance >= depthS. Donc
    pour chercher si v est déjà dans S[u], on peut se contenter de
    vérifier les sommets de S[u] à distance >= depthS-1 de u. Ces
    sommets candidats sont stockés à partir de l'indice S[u]->vpd dans
    S[u]->node (car S n'est pas encore retriée). Dans R[u]->dist on
    stocke la distance réalisée dans l'arbre pour faire u->v.
  */
  
  ALLOCZ(R,n,NULL); /* tableau de tables R (vides au départ) */
  for(i=0;i<k;i++){ /* seulement pour les k sommets u de C */
    u=C[i]; /* u=sommet de C */
    t=F[i]-(H[u]==i); /* t=#sommets à mettre dans R[u] (moins éventuellement u) */
    if((u!=center)&&(H[center]==i)) t--; /* enlève u et aussi center si H[center]=i */
    if(S[u]) t-=S[u]->n; /* on enlève ceux déjà dans S[u] */
    if(t>0){ /* si table pas vide */
      R[u]=new_table(t); /* table R pour u */
      R[u]->n=0; /* on se sert de R[u]->n comme indice */
    }
  }
  free(F); /* ne sert plus à rien */
  /* on remplit les tables R[u]->node et R[u]->dist */
  /* R[u]->dist[i]=dist. dans tree(center) entre u et le i-ème sommet de R[u] */
  for(v=0;v<n;v++){ /* écrit chaque sommet v dans la bonne table R[u] */
    u=C[H[v]]; /* u=sommet de C en charge du hash de v (donc de couleur H[v]) */
    if((u==v)||(v==center)) continue; /* on ne met ni u ni center dans R[u] */
    if((X0->D[v]<depthS)&&(u!=center)) continue; /* seulement pour v à dist > depthS-1 de center */
    /* On cherche si v n'est pas déjà dans S[u]. Il suffit de chercher
       des sommets v à distance >= depthS-1 de u */ 
    if(S[u]){ /* cherche v dans S[u], si elle n'est pas vide bien sûr */
      j=S[u]->vpd; /* 1er sommet à distance >= depthS-1 de u */
      if(j<0) j=S[u]->n; /* ne rien faire si pas de 1er sommet (est-ce nécessaire ?) */
      for(;j<S[u]->n;j++) if(S[u]->node[j]==v) break; /* v est dans S[u] */
    }
    if((S[u]==NULL)||(j==S[u]->n)){ /* v n'a pas été trouvé dans S[u] */
      R[u]->node[R[u]->n]=v; /* ajoute v à R[u] */
      // Calcule R[u]->dist: on remonte de v jusqu'au cluster et on
      // voit si l'on passe par u ou pas. La distance d de v à u dans
      // tree(center) est la suivante:
      //   si u=center, alors          d=X0->D[v]
      //   si u<>center et t=u, alors  d=X0->D[v]-1
      //   si u<>center et t<>u, alors d=X0->D[v]+1
      if(u==center) /* remontée inutile si u est le center */
	R[u]->dist[R[u]->n]=X0->D[v];
      else{ /* on remonte à partir de v jusqu'à t=sommet dans C */
	t=v; while(X0->P[t]!=-2) t=X0->P[t]; /* si t pas dans C, t=parent(t) */
	R[u]->dist[R[u]->n]=X0->D[v]+((t==u)?-1:1);
      }
      R[u]->n++; /* un sommet de plus dans R[u] */
    }
  }
  printf("- time to construct tables R: %s\n",TopChrono(1));
  free_param_bfs(X0); /* ne sert plus à rien */

  DEBUG(
	for(i=0;i<k;i++){
	  u=C[i];
	  printf("u=%i hash=%i |R[u]|=%i\n",u,i,R[u]?R[u]->n:0);
	  if(R[u]) PRINTT(R[u]->node,R[u]->n);   
	  if(R[u]) PRINTT(R[u]->dist,R[u]->n);
	}
	);
    
  /************************/
  /* calcule les tables W */
  /************************/

  /*
    Tables seulement définies pour tous les sommets u de C, W[u] donne
    la liste des couleurs (=indice dans C) des voisins de u dans C. La
    couleur du center n'est jamais ajouté à W[u], car c'est le port
    par défaut. NB: Le nom du sommet de couleur i est C[i]. On ne se
    sert pas de W[u]->dist. Si W[u] est vide, alors on ne la supprime
    pas de façon à avoir le port par défaut. Dans le cas bit-1 de
    VARIANT=1, toutes les tables sont vides (routage dans l'étoile
    sans table pour le center et les feuilles).
  */

  ALLOCZ(W,n,NULL); /* tableau des tables W (vides au départ) */
  for(i=0;i<k;i++){ /* seulement pour les sommets u de C, y compris center */
    u=C[i]; /* u=sommet de C de couleur i */
    W[u]=new_table(0); /* table W pour u, NB: W[u]->n=0 */
    ALLOC(W[u]->node,k-1); /* au plus k-1 voisins */
    W[u]->vpd=center; /* voisin par défaut = center */
    if(VARIANT&2) continue; /* toutes les tables W seront vides */
    for(j=0;j<G->d[u];j++){ /* pour chaque voisin de u */
      v=G->L[u][j]; /* v=j-ème voisin de u */
      if(B[v]||(v==center)) continue; /* on veut v dans C et v<>center */
      W[u]->node[W[u]->n++]=SetSearch(v,C,k,1); /* ajoute la couleur de v à W[u]->node */
    }
    if(W[u]->n) /* réajuste la taille */
      REALLOC(W[u]->node,W[u]->n);
  }
  printf("- time to construct tables W: %s\n",TopChrono(1));

  DEBUG(
	for(i=0;i<k;i++){
	  u=C[i];
	  printf("u=%i hash=%i |W[u]|=%i\n",u,i,W[u]?W[u]->n:0);
	  if(W[u]) PRINTT(W[u]->node,W[u]->n);   
	}
	);

  /*******************/
  /* trie les tables */
  /*******************/
  for(u=0;u<n;u++){
    if(B[u]) QSORT2(B[u]->node,B[u]->dist,B[u]->n,fcmp_int);
    if(R[u]) QSORT2(R[u]->node,R[u]->dist,R[u]->n,fcmp_int);
    if(S[u]) QSORT2(S[u]->node,S[u]->dist,S[u]->n,fcmp_int);
    if(W[u]) QSORT(W[u]->node,W[u]->n,fcmp_int);
  }
  printf("- time to sort tables B,R,S,W: %s\n",TopChrono(1));
  BARRE;
  
  /* calcule de la taille Z[u] des tables de chaque sommet u */
  ALLOCZ(Z,n,1); /* taille=1 au moins pour chaque sommet u */
  for(u=0;u<n;u++){ /* pour chaque sommet u */
    if(B[u]) Z[u] += B[u]->n;
    if(S[u]) Z[u] += S[u]->n;
    if(R[u]) Z[u] += R[u]->n;
    if(W[u]) Z[u] += W[u]->n;
  }

  /* affiche taille min/max et moyenne/écart type des différentes tables */
  MINMAXMOY(B[_i]->n,n,B[_i],"table B size");
  MINMAXMOY(S[C[_i]]->n,k,S[C[_i]],"table S size");
  MINMAXMOY(R[C[_i]]->n,k,R[C[_i]],"table R size");
  MINMAXMOY(W[C[_i]]->n,k,W[C[_i]],"table W size");

  /* affiche la distribution des tailles de table */
  PrintDistribution(Z,n,10,"routing table size");
  free(Z); /* ne sert plus à rien */
  BARRE;

  /* assemble les tables en une seule pour le retour */
  NALLOC(rs_cluster_tables,RT,1);
  RT->B=B;
  RT->S=S;
  RT->R=R;
  RT->W=W;
  RT->C=C; // besoin pour le routage avec W
  RT->H=H; // besoin du hash des sommets
  RT->n=n; // besoin pour libérer les n tables
  RT->center=center; // besoin pour le routage en étoile (tables W vides)

  printf("total time: %s\n",TopChrono(2));
  return RT;
}


int rs_cluster_length(int u,int v,rs_cluster_tables *X){
/*
  Renvoie le nombre de sauts du routage selon les tables générées par
  rs_cluster() pour router un message de u à v, ou bien -1 si la route
  n'a pu être déterminée. Dans X on a toutes les tables nécessaire au
  schéma, notamment les tables B,S,R,W. Si u<0 alors on teste la
  validité des tables (et on renvoie une valeur non-nulle en cas
  d'erreur). L'algorithme dépend du bit-1 de VARIANT.

  Amélioration possible: Lorsque u a des voisins dans C, alors, plutôt
  que d'aller vers ->vdp, choisir un landmark de C parmi un ensemble
  prédéterminé de p = ceil(2m/n) = O(1) sommets (le degré moyen)
  d'éventuellement la bonne couleur h(v). Aller vers ->vdp seulement
  si cette recherche à échoué. Pour le schéma théorique, on peut
  découper la table B en deux, B1 et B2, où B2 serait les sommets
  voisins de u dans C, et doubler B2 (double accès par sommet et par
  couleur) de façon à garantir un temps constant dans tous les cas, et
  pas seulement un temps égale au degré moyen (garanti seulement avec
  grande proba dans les RPLG.
*/

  if(u<0){
    if(X==NULL) return 1;
    if(X->B==NULL) return 1;
    if(X->S==NULL) return 1;
    if(X->R==NULL) return 1;
    if(X->W==NULL) return 1;
    if(X->H==NULL) return 1;
    if(X->C==NULL) return 1;
    return 0;
  }
  
  DEBUG(printf("  \nu=%i v=%i: ",u,v););

  // on est arrivé
  if(u==v){ DEBUG(printf("u=v\n");); return 0; }
  int i;

  // si u n'est pas dans C
  
  if(X->B[u]){
    // v dans B[u] ?
    i=SetSearch(v,X->B[u]->node,X->B[u]->n,1);
    // si v dans B[u]:u -> v
    if(i>=0){ DEBUG(printf("in B[u]\n");); return X->B[u]->dist[i];}
    // si v pas dans B[u]: u -> B[u]->vpd -> v
    DEBUG(printf("via B[u]->vpd=%i\n",X->B[u]->vpd););
    return 1 + rs_cluster_length(X->B[u]->vpd,v,X);
  }
  
  // si u est dans C

  if(X->S[u]){ // il faut S non vide, donc u<>center
    // v dans S[u] ?
    i=SetSearch(v,X->S[u]->node,X->S[u]->n,1);
    // si v est dans S
    if(i>=0){ DEBUG(printf("in S[u]\n");); return X->S[u]->dist[i];}
  }

  // si v n'est pas dans S
  // v dans R[u] ?

  if(X->R[u]){ // il faut R non vide
    i=SetSearch(v,X->R[u]->node,X->R[u]->n,1);
    // si v est dans R, on cherche si v descendant de u dans tree(center)
    if(i>=0){ DEBUG(printf("in R[u]\n");); return X->R[u]->dist[i]; }
  }
  
  // si v n'est pas ni dans R ni dans S
  // H[v] dans W[u] ?

  if(X->W[u]){ /* ne peut pas être vide */
    i=SetSearch(X->H[v],X->W[u]->node,X->W[u]->n,1);
    // si H[v] est dans W[u]: u -> C[H[v]] -> v
    // si H[v] pas dans W[u]: u -> center -> v
    if((i>=0)||((VARIANT&2)&&(u==X->center))) u=X->C[X->H[v]];
    else u=X->W[u]->vpd;
    DEBUG(
	  if(u==X->C[X->H[v]])
	    printf("%s: u -> C[H[v]]=%i\n",(VARIANT&2)?"sans tables W":"H[v] in W[u]",u);
	  else
	    printf("%s: u -> W[u]->vpd=%i\n",(VARIANT&2)?"sans tables W":"H[v] not in W[u]",u);
	  );
    return 1 + rs_cluster_length(u,v,X);
  }

  // y'a un problème
  DEBUG(printf("fail: W[u] does not exist\n"););
  return FAIL_ROUTING;
}


rs_dcr_tables *rs_dcr(graph *G,int k){
/*
  Calcule pour le graphe G le routing scheme DCR ou AGMNT de paramètre
  k>0 correspondant au nombre de couleurs. On renvoie les tables ainsi
  calculées. Le stretch est toujours <= 5 (DCR) ou <=3 pour AGMNT. La
  distinction entre DCR et AGMNT se fait par le bit-2 de VARIANT. Les
  champs ->CONT et ->F qui sont également NULL pour DCR. Le bit-0 de
  VARIANT est réservé au choix des landmarks de plus haut degré,
  variante valable aussi bien pour DCR que pour AGMNT.

  Si une couleur n'existe pas, alors cela marche quand même: les
  boules couvrent tout le graphe et le stretch est alors 1.

  Spécificitées pour AGMNT: Pour le routage de u -> v, où
  hash(v)=C(u), il faut prendre la meilleure des options parmi:

  1) routage via un landmark s et w=nca(u,v) dans T(s). Dans ce cas,
     on a besoin de stocker s, éventuellement w pour éviter de le
     recalculer.

  2) routage via les boules contiguës via une arête x-y. Plus
     précisément, on prend la plus courte des routes, si elles
     existent, de la forme u->s->x-y->v tel que B(s) contient u et x,
     B(v) contient y, et x,y sont voisins. NB: u=s et y=v sont
     possibles.

  La route u->s->x-y->v, si elle existe, a la propriété que s->x-y->v
  est un plus court chemin. De plus le plus petit ancêtre commun dans
  T(s) entre u et x doit être s. Si c'était un sommet w<>s, alors w
  serait un meilleur candidat vérifiant toutes les propriétés. Donc le
  sommet x à chercher ne peut pas être un descendant de u' dans T(s)
  où u' est le voisin de s menant à u. Donc la longueur du routage est
  alors d(u,s)+d(s,x)+1+d(y,v).

  La recherche des routes selon les boules contiguës peut se faire
  ainsi (le calcul réellement réalisé est cependant un peu différent):

  pour tout sommet u:
    pour tout sommet s de B^-1(u):
      pour tout x de B(s) qui n'est pas descendant de u':
        pour tout voisin y de x:
	  pour tout v de B^-1(y) tq hash(v)=C(u):
	    calculer d(u,s)+d(s,x)+1+d(y,v).

  Il faut calculer B^-1(u) de chaque sommet, puis trier leurs sommets
  selon leur hash de sorte qu'on puisse trouver très rapidement tous
  les sommets de B^-1(y) ayant un hash donné (=C(u)). La complexité
  est grosso-modo O(m*n*ln(n)). Les distances entres les centres des
  boules n'ont pas à être calculées.

  Pour savoir si x est descendant ou pas de u', il faut deux choses:
  1) pour chaque sommet x de B(s), avoir un numéro dfs(x) dans l'arbre
  T(s) couvrant B(s).  2) pour chaque sommet s de B^-1(u), avoir
  l'intervalle [a,b] des dfs(x) pour les sommets x de B(s) descendant
  de u'. Il faut donc stocker 1 entier (dfs) pour chaque sommet des
  boules, et deux entiers (a,b) pour chaque sommet d'une boule
  inverse. En pratique cela nécessite des calculs et de la mémoire
  supplémentaires, et il n'est pas clair que cela soit plus efficace.
*/

  const int n=G->n;
  const int agmnt=VARIANT&4; // agmnt vrai ssi c'est le schéma AGMNT, sinon c'est DCR
  int u,v,w,i,j,t,c,nc,nl;

  printf("\n%s\n",agmnt? "AGMNT" : "DCR");
  BARRE;

  TopChrono(1); /* reset du chrono tmp */
  TopChrono(2); /* reset du chrono total */

  printf("- wanted number of colors: %i\n",k);

  /**********************************/
  /* calcule H & C, hash et couleur */
  /**********************************/

  // H[u]=0..k-1, hash du sommet u
  // C[u]=0..k-1, couleur du sommet u
  // u est landmark ssi C[u]=0

  DEBUG(PRINT(HASH););
  int *H=MakeHash(NULL,n,k,HASH);
  NALLOC(int,C,n);

  if((VARIANT&1)&&(k>1)){
    // variante où les landmarks sont les sommets de plus haut
    // degré. On va les colorier 0 puis colorier les autres sommets
    // avec une couleur de 1 à k-1. NB: il faut au moins 2 couleurs, k>1.
    for(u=0;u<n;u++) C[u]=1; // au départ les sommets n'ont pas de couleur
    t=iceil(n,k); // t=nombre de landmarks, estimé à la moyenne
    int *D=SortInt(G->d,NULL,n,1,NULL,SORT_INDEXi); // trie le degré des sommets
    for(i=0;i<t;i++) C[D[n-i-1]]=0; // le i-ème sommet de plus haut degré est un landmark
    for(u=0;u<n;u++) if(C[u]) C[u]=1+random()%(k-1); // colorie dans [1,k[ les sommets non landmark
  }else
    // variante par défaut: les sommets ont une couleur aléatoire dans [0,k[
    for(u=0;u<n;u++) C[u]=random()%k;

  DEBUG(PRINTT(C,n);PRINTT(H,n););

  /* affichage de stats sur H & C */
  /* calcule:
     nl=nombre de landmarks
     nc=nombre de couleurs
     F[c]=nombre de sommets de hash c
  */
  int *F;
  FREQMINMAX(F,k,C,n,"the colors"); /* fréquence des couleurs (ne sert qu'à l'affichage) */
  nl=F[0]; /* nl=nombre de landmarks, ceux de couleurs 0 */
  for(i=nc=0;i<k;i++) nc += (F[i]>0); /* nc=nombre de couleurs différentes, nc<=k */
  free(F); /* nécessaire car réalloué par le prochain FREQMINMAX */
  FREQMINMAX(F,k,H,n,"the hash"); /* fréquence des hashs (sert pour taille des tables) */
  for(i=u=0;i<k;i++) u += (F[i]>0); /* u=nombre de hashs différents, u<=k */

  printf("- real number of colors: %i\n",nc);
  printf("- real number of hashs: %i\n",u);
  printf("- number of landmarks: %i\n",nl);
  printf("- time to compute hash, color & stats: %s\n",TopChrono(1));
  BARRE;

  /*********************************************/
  /* calcule S, la liste des bfs des landmarks */
  /*********************************************/

  // S[u]=bfs(u,G) pour u landmark, NULL si u non landmark
  
  NALLOC(param_bfs*,S,n);
  for(u=0;u<n;u++)
    if(C[u]) S[u]=NULL; // C[u]=0 ssi u landmark
    else{
      S[u]=bfs(G,u,NULL);
      free(S[u]->file); /* libère les pointeurs inutilisés */
      S[u]->file=NULL;  /* important pour le free_param_bfs() plus tard */
    }
  printf("- time to construct landmark bfs (array S): %s\n",TopChrono(1));

  DEBUG(
	for(u=0;u<n;u++)
	  if(C[u]==0){
	    PRINT(u);
	    PRINTT(S[u]->P,n);
	    PRINTT(S[u]->D,n);
	    printf("\n");
	  }
	);
  
  /*****************************/
  /* calcule les tables B et W */
  /*****************************/

  /*
    B[u]=boule de voisinage de u, la plus "petite" contenant toutes
         les couleurs (composé de la liste des sommets et leur
         distance à u), la dernière couronne étant ordonnée selon les
         identifiants des sommets.

	 Il est important que les boules contiennent k couleurs, de
         façon à toujours pouvoir être capable de router vers un hash
         donné. Si une couleur n'apparaît pas dans le graphe, alors
         chaque boule devra contenir tous les sommets.

	 B[u] contient toujours au moins u, bien qu'on ne se sert
         jamais de cette entrée pour le routage. Par contre, pour la
         construction, c'est plus simple d'avoir u dans B[u]. On
         enlèvera donc cette entrée pour le routage (pour les stats et
         établir la taille des tables), après avoir construit toutes
         les tables.

    B[u]->node[i]=i-ème sommet de la boule de u
    B[u]->dist[i]=distance entre u et B[u]->node[i]
    B[u]->vpd=+proche landmark de u (+proche sommet de couleur 0 dans B[u])
    W[u]->node[c]=next-hop vers le 1er sommet de couleur c de B[u]->node
    (la table W, indexée par les couleurs, sert pour avoir le temps constant)

    Pour AGMNT:

    K[u]=taille de la boule inverse de B[u], ie K[u]=|{v : u in B[v]}|
    B[u]->radius=rayon de B[u]

    Algorithme:

      On fait un bfs(G,u) par couche (avec cont=1). On trie chaque
      couche par identifiant (si on a visité au moins assez de
      sommets), et on regarde à chaque fois le premier moment où l'on
      voit une couleur donnée (remplissage de W). Au passage on
      détermine, lorsque la couleur est 0, le +proche landmark t de
      u. On s'arrête quand on a visité toutes les couleurs ou bien que
      tous les sommets du graphe ont été visités (cela peut arriver si
      une couleur n'est portée par aucun sommet).
   */

  NALLOCZ(table*,B,n,NULL);  // B=tableau de n tables B (vides au départ)
  NALLOCZ(table*,W,n,NULL);  // W=tableau de n tables W (vides au départ)
  param_bfs *X=new_param_bfs(); // X=le résultat du bfs() depuis u
  NALLOC(int,T,n); /* T=tableau pour la dernière couche de sommets */
  int *K; if(agmnt) ALLOCZ(K,n,0); /* K=taille des boules inverses */
  X->clean=1; /* initialisation complète des distances, puis partielle */

  for(u=0;u<n;u++){ /* calcule B[u] & W[u] pour chaque u */
    c=0; /* c=nombre de couleurs (différentes) déjà rencontrées */
    X->cont=1; /* évite de recalculer le début de l'arbre */
    X->hmax=-1; /* au départ on visite seulement u */
    W[u]=new_table(0); /* W[u]->node[c]=-1 si la couleur c n'a pas été visitée */
    ALLOCZ(W[u]->node,k,-1); /* NB: on n'utilise pas ->dist */
    W[u]->n=k; /* taille de W[u] */
    do{
      
      X->hmax++; /* couche suivante (au départ hmax=0) */
      bfs(G,u,X); /* on parcoure tous les sommets jusqu'à distance (incluse) hmax de u */
      for(j=0,i=X->tf;i<X->n;i++) T[j++]=X->file[i]; /* copie la dernière couche dans T */
      if((c+j>=k)||(X->n>=n)) QSORT(T,j,fcmp_int); /* trie la dernière couche (=T), obligé à cause de W[u] */
      /* on ne trie T que si on a le potentiel pour avoir terminé:
	 avoir toutes les couleurs ou avoir visité tout le graphe.  Si
	 le nombre de sommets de la dernière couche (=j) + le nombre
	 de couleurs déjà rencontrées (=c) est au moins k alors on a
	 potentiellement atteint la dernière couche. */

      DEBUG(
	    //PRINT(u);PRINT(X->n);PRINT(X->tf);
	    //PRINTT(T,j);PRINT(c);PRINT(j);
	    //printf("\n");
	    );
      
      /* remplit W en parcourant les sommets de la dernière couche */
      /* le +proche landmark de u, lorsque rencontré, est stocké dans t */
      for(i=0;(i<j)&&(c<k);i++){
	v=T[i]; /* v=i-ème sommet de la dernière couche */
	if(W[u]->node[C[v]]>=0) continue; /* couleur déjà rencontrée */
	/* ici la couleur C[v] n'a jamais été rencontrée */
	w=v; /* NB: si v=u, alors il faut effacer le -1 dans W[u] */
	if(v!=u) while(X->P[w]!=u) w=X->P[w]; /* cherche le next-hop w depuis u vers v */
	W[u]->node[C[v]]=w; /* met le next-hop dans W[u], ou bien v si u=v */
	if(C[v]==0) t=v; /* si v est un landmark, alors c'est le +proche de u */
	c++; /* une couleur de plus */
      }

      /* On sort de la boucle si: soit on a toutes les couleurs (=k)
	 dans W[u]->node, ou bien on a visité tout le graphe, ce qui
	 est possible si toutes les couleurs n'étaient pas
	 représentées.  NB: ici, dans tous les cas, i est nombre de
	 sommets de la dernière couche à recopier partiellement. */
      
    }while((c<k)&&(X->n<n));

    /* on construit B[u] à partir de X->file et de T */
    B[u]=new_table(X->tf+i); // table de B[u]
    B[u]->vpd=t; // t=+proche landmark de u=W[u]->node[0]
    B[u]->radius=X->hmax; // rayon de B[u], ne sert à rien pour DCR

    /* construit B[u] */
    /* j=indice pour B[u]->node[] */
    /* t=indice pour X->file[] */
    for(j=t=0;t<X->tf;j++){ /* copie X->file sauf la dernière couche */
      v=X->file[t++];
      B[u]->node[j]=v; if(agmnt) K[v]++; // car v dans B[u]
      B[u]->dist[j]=X->D[v];
    }
    for(t=0;t<i;j++){ /* copie la partie traitée de la dernière couche */
      v=T[t++];
      B[u]->node[j]=v; if(agmnt) K[v]++; // car v dans B[u]
      B[u]->dist[j]=X->D[v];
    }
  }
  printf("- time to construct tables B & W: %s\n",TopChrono(1));
  free_param_bfs(X);
  free(T); /* ne sert plus à rien */

  /* tri des tables B */
  /* c'est important pour le routage et pour optimiser AGMNT */
  /* NB: B[u]->n>0, car u est toujours dans B[u] */
  for(u=0;u<n;u++) QSORT2(B[u]->node,B[u]->dist,B[u]->n,fcmp_int);
  printf("- time to sort tables B: %s\n",TopChrono(1));
  BARRE;


  /*****************************/
  /* Partie spécifique à AGMNT */
  /*****************************/

  contigue **CONT=NULL;

  if(agmnt){
    /*********************************/
    /* calcule les boules inverses I */
    /*********************************/

    // I[u][i]=i-ème sommets de B^-1[u] trié selon la distance à u
    // K[u]=nombre de sommets dans I[u] (déjà calculé)
    //
    // Les boules inverses ainsi triées servent à accélérer le calcul
    // des routes via les boules contiguës.

    // alloue les boules inverses (I) et les distances au centre (D)
    NALLOC(int*,I,n);
    NALLOC(int*,D,n); // D ne sert que pour le tri des boules inverses
    for(u=0;u<n;u++){
      ALLOC(I[u],K[u]); // allocation des boules inverses
      ALLOC(D[u],K[u]); // allocation des distances au centre
    }

    // remplit les boules inverses et les distances au centre
    for(u=0;u<n;u++) K[u]=0; // K[u]=nombre de sommets déjà mis dans I[v]
    for(u=0;u<n;u++) // pour toutes les boules
      for(i=0;i<B[u]->n;i++){ // pour tous les sommets v de B[u]
	v=B[u]->node[i];
	I[v][K[v]]=u; // ajoute u à I[v]
	D[v][K[v]]=B[u]->dist[i]; // dist(u,v)
	K[v]++;
      }
    // ici K[u]=taille de I[u] pour tout u
    printf("- time to construct inverse balls (array I): %s\n",TopChrono(1));
    DEBUG(
	  PRINTT(K,n);
	  for(u=0;u<n;u++){
	    PRINT(u);
	    PRINTT(I[u],K[u]);
	    PRINTT(D[u],K[u]);
	    printf("\n");
	  }
	  );
    /* trie les boules inverses selon la distance à u */
    /* NB: B[u] contient u, alors I[u] est de taille K[u]>0 */
    for(u=0;u<n;u++) QSORT2(D[u],I[u],K[u],fcmp_int); // trie I[u] selon D[u]
    FREE2(D,n); // ne sert plus à rien
    printf("- time to sort inverse balls: %s\n",TopChrono(1));

    /****************************************************/
    /* calcule le meilleur chemin de u à v de hash C[u] */
    /****************************************************/

    // CONT[u][i] = meilleur chemin pour aller de u à v, le i-ème
    // sommet dont le hash vaut C[u] ordonné selon i
    //
    // UC[i]=i-ème sommet trié u selon sa couleur C[u]
    // UH[i]=i-ème sommet trié v selon son hash H[v] puis selon v
    //
    // Ces deux tableaux servent à lister en temps n*k*H(k) ~ n^1.5
    // les paires (u,v) avec hash de v = C[u], les sommets de même
    // hash se retrouvant consécutifs dans UH. On pourrait construire
    // UC et UH en temps O(n) (au lieu de O(nlogn) comme on le fait)
    // en réutilisant les tableaux de fréquence des hashs et des
    // couleurs. Mais c'est plus complexe, et le gain en temps sera au
    // final négligeable.
    
    // construit la liste UC
    NALLOCZ(int,UC,n,_i); // UC[i]=i au départ
    fcmp_tabint(NULL,C); // tri selon C
    QSORT(UC,n,fcmp_tabint); // tri UC selon C

    // construit la liste UH
    NALLOCZ(int,UH,n,_i); // UH[i]=i au départ
    fcmp_tabinteq(NULL,H); // tri selon H
    QSORT(UH,n,fcmp_tabinteq); // tri UH selon H et v

    printf("- time to sort nodes according to hash and colors: %s\n",TopChrono(1));

    DEBUG(
	  PRINTT(UC,n);
	  PRINTT(UH,n);
	  );
    
    int l,d,d1,d2,s,h,p,q,r,x,y;
    ALLOC(CONT,n); // CONT[u] va être alloué pour tous les sommets u

    // Calcule la distance via les landmarks ou les boules contiguës.
    // On balaye les sommets u dans l'ordre croissant des couleurs.
    // Les sommets v, de hash C[u], sont consécutifs dans UH.

    // Pour les landmarks u, le routage vers v (de hash 0 donc)
    // s'effectue comme un routage via le landmark u (c'est ici un
    // meileur choix que le landmark de v), ce qui va produire une
    // route de plus court chemin dans T_u.

    for(j=t=0;j<n;j++){ // pour tous les sommets
      u=UC[j];   // u=sommet courant
      h=F[C[u]]; // h=nombre de sommets de hash C[u], h=0 est possible
      ALLOC(CONT[u],h); // alloue pour h sommets (ne fait rien si h=0)
      
      // calcule l'indice t de UH[t..t+h[ où sont rangés les v de hash C[u]
      q=j? UC[j-1] : u; // q=sommet juste avant u dans UC, q=u au départ (si j=0)
      if(C[u]>C[q]) t+=F[C[q]]; // changement de couleur ?

      // balaye tous les sommets v de hash C[u], de plus par ordre
      // croissant (important pour faire une recherche binaire au
      // moment du routage vers v)

      for(i=0;i<h;i++){
	v=UH[t+i]; // v=sommet de hash C[u]
	CONT[u][i].v=v; // stocke v dans tous les cas, important pour
			// la recherche binaire lors du routage
	if((C[v]==0)||(u==v)) continue; // ne rien à faire si v est landmark ou si u=v
	if(C[u]==0){ // NB: on peut avoir t=0 et C[u]<>0 si F[0]=0 par exemple
	  // si u est landmark, on provoque un routage de plus court
	  // en codant une solution via boule-contiguë de distance
	  // optimale dist(u,v) dans T_u. NB: CONT[u][i].w n'est pas défini
	  CONT[u][i].s=-1; // pour dire via boule-contiguë et forcer un +cc
	  CONT[u][i].d=S[u]->D[v]; // dist(u,v)
	  continue; // v suivant
	}
	if(SetSearch(v,B[u]->node,B[u]->n,1)>=0) continue; // si v est dans B[u], on a rien à faire
	// ici v n'est pas un landmark et pas dans B[u] => dist(u,v)>=B[u]->radius

	q=SetSearch(u,B[v]->node,B[v]->n,1); // est-ce que v est dans B^-1[u] ?
	if(q>=0){ // alors routage +cc
	  CONT[u][i].s=-1; // pour dire via boule-contiguë et forcer un +cc
	  CONT[u][i].d=B[v]->dist[q]; // dist(u,v)
	  continue; // v suivant
	}

	// route via landmark ? On vise une route du type u->w->v,
	// avec w=nca(u,v,T_s) pour un certain landmark s. On commence
	// par le +proche landmark s de v (ce qui suffit pour garantir
	// un stretch <= 3), si bien qu'on ne change de landmark que si
	// on fait strictement mieux.

	CONT[u][i].d=INT_MAX; // important car le champs d pourrait ne pas être initialisé (si nl=0)
	for(l=-1;l<nl;l++){ // balaye tous les landmarks et aussi B[v]->vpd
	  s=(l<0)? B[v]->vpd : UC[l]; // s=landmark=un sommet de couleur 0 ou bien B[v]->vpd
	  w=nca_bfs(u,v,S[s]); // w=nca(u,v,T_z)
	  d=dist_nca(u,v,w,S[s]->D); // dist(u,v) dans T_s
	  if(d<CONT[u][i].d){ // on a trouvé une meilleure route via w dans T_s
	    CONT[u][i].s=s;
	    CONT[u][i].w=w;
	    CONT[u][i].d=d;
	  }
	}
	DEBUG(
	      if(u==28 && v==36){
		printf("via best landmark:\n");
		PRINT(CONT[u][i].s);
		PRINT(CONT[u][i].w);
		PRINT(CONT[u][i].d);
	      }
	      );
	
	// route via boule-contiguë ? On vise une route du type
	// u->s->x-y->v.  La plus courte de ces routes fait que
	// nécessairement s=nca(u,x,T_s) et que s->x-y->v est un +cc.
	// NB: grâce au tri des B^-1[u], les sommets s sont parcourus
	// par distance croissante à u. Cela permet ainsi de stoper
	// plus rapidement la recherche d'un bon s.

	for(l=0;l<K[u];l++){ // balaye tous les sommets s de B^-1[u] = I[u]
	  s=I[u][l]; // s=sommet contenant u dans B[s], NB: s<>v car v pas dans B^-1[u] qui contient s
	  d=B[s]->dist[SetSearch(u,B[s]->node,B[s]->n,1)]; // d=|u->s|
	  q=SetSearch(v,B[s]->node,B[s]->n,1);
	  DEBUG(
		if(u==28 && v==36){
		  PRINT(s);
		  PRINT(CONT[u][i].d);
		  PRINT(B[s]->radius);
		  PRINT(d);
		  PRINT(q);
		  printf("\n");
		}
		);
	  if(q>=0){ // si v est dans B[s]
	    d += B[s]->dist[q]; // d=|u->s->v|
	    if(d<CONT[u][i].d){ // on a trouvé une route meilleure pour u->v
	      CONT[u][i].d=d;   // d=distance trouvée
	      CONT[u][i].s=-1;  // pour dire via boule-contiguë
	    }
	    continue; // on peut passer au s suivant
	  }
	  // ici v n'est pas dans B[s], donc dist(s,v) >= B[s]->radius >= 1 (s<>v)
	  if(d+1>=CONT[u][i].d) break; // aucun s ne pourra faire mieux -> v suivant
	  if(d+B[s]->radius>=CONT[u][i].d) continue; // ce s n'est pas assez bon -> s suivant
	  for(r=0;r<B[s]->n;r++){ // balaye tous les sommets x de B[s], x=s compris
	    // on suppose que x (=B[s]->node[r]) est le plus loin
	    // possible de s (ou le +proche de v). Il doit alors être
	    // à distance rayon de B[s] ou rayon de B[s] - 1 de
	    // s. C'est sans perte de généralité car v n'est pas dans
	    // B[s] et x sur un +cc entre s et v.
	    if(B[s]->dist[r]<B[s]->radius-1) continue; // x n'est pas le plus loin possible
	    d1 = d+B[s]->dist[r]+1; // d1=|u->s->x-y|
	    if(d1>=CONT[u][i].d) continue; // la route ne fera pas mieux
	    x=B[s]->node[r]; // x=le sommet de B[s]
	    for(p=0;p<G->d[x];p++){ // balaye les y voisins de x
	      y=G->L[x][p]; // y=voisin de x
	      q=SetSearch(y,B[v]->node,B[v]->n,1); // on cherche y dans B[v]
	      DEBUG(
		    if(u==28 && v==36){
		      PRINT(d1);
		      PRINT(x);
		      PRINT(y);
		      PRINT(q);
		    }
		    );
	      if(q<0) continue; // il faut y dans B[v]
	      // ici on a trouvé une route du bon type
	      d2 = d1+B[v]->dist[q]; // d2=d1+|y->v|=|u->s->x-y->v|
	      if(d2<CONT[u][i].d){   // on a trouvé une route meilleure pour u->v
		CONT[u][i].d=d2;     // d2=distance trouvée
		CONT[u][i].s=-1;     // pour dire via boule-contiguë
	      }
	    }
	  }
	}
      }
    }
    
    free(UH);
    free(UC);
    FREE2(I,n);
    printf("- time to compute best routes via landmark or contigue-ball: %s\n",TopChrono(1));
    BARRE;
  }// fin du "if(agmnt){ ..."
  
  DEBUG(
	for(u=0;u<n;u++){
	  PRINT(u);
	  PRINT(B[u]->vpd);
	  PRINT(B[u]->radius);
	  PRINTT(B[u]->node,B[u]->n);
	  PRINTT(W[u]->node,W[u]->n);
	  printf("\n");
	}
	if(agmnt){
	  for(u=0;u<n;u++){
	    printf("\nCONT[%i].v = ",u);for(i=0;i<F[C[u]];i++) printf("%i ",CONT[u][i].v);
	    printf("\nCONT[%i].s = ",u);for(i=0;i<F[C[u]];i++) printf("%i ",CONT[u][i].s);
	    printf("\nCONT[%i].d = ",u);for(i=0;i<F[C[u]];i++) printf("%i ",CONT[u][i].d);
	    printf("\nCONT[%i].w = ",u);for(i=0;i<F[C[u]];i++)
					  if(CONT[u][i].s<0) printf("- ");
					  else printf("%i ",CONT[u][i].w);
	    printf("\n");
	  }
	}
	);

  /****************************/
  /* taille totale des tables */
  /****************************/

  /* pour chaque sommet, il faut en moyenne:
     
     T1: une boule -> k*H(k)
     T2: une table de landmark -> n/k
     T3: une table des hashs -> k-1 (à cause du temps constant, sinon -> 0)
     T4: une table des sommets de même hash que C[u] -> n/k
   
     A cause de l'algorithme de routage, on peut optimiser les tables
     ainsi. Si v est dans B[u] et est un landmark, on peut l'enlever
     de T2, ce qui en moyenne revient à enlever H(k) entrées à
     T2. Donc |T2| = n/k - H(k) en moyenne. De même, on peut enlever
     de T4 les entrées des sommets qui sont dans T1 ou T2, mais le
     gain est très léger.
   
     La fonction donnant le nombre d'entrées est alors
     
              f(k,n) = 2n/k + (k-1)*(H(k)+1)

     voir la fonction func1(). Une façon de le voir l'optimisation de
     T2 est que les entrées de B[u] pourraient être découpées en deux
     tables: B1[u] et B2[u] où B1[u] seraient les entrées pour les
     sommets non-landmarks et B2[u] pour les landmarks. Pour B2[u] la
     structure seraient un peu différente, elle aurait à la fois celle
     de B1[u] et celle des landmarks.
  */

  NALLOCZ(int,Z,n,1); /* taille de la table de u, au moins 1 pour chaque u */
  NALLOCZ(int,Z1,n,nl); /* taille de la table des landmarks qui ne sont pas dans B[u] */
  NALLOCZ(int,Z2,n,F[C[_i]]); /* taille de la table des hash = C[u] hors de B[u] */

  for(u=0;u<n;u++){ /* pour chaque sommet u */
    for(i=0;i<B[u]->n;i++){ // calcule Z1[u] et Z2[u]
      v=B[u]->node[i]; // v dans B[u]
      Z1[u] -= (C[v]==0); /* enlève les landmarks v qui sont dans B[u] */
      Z2[u] -= (H[v]==C[u]); /* enlève les hashs C[u] qui sont dans B[u] */
    }
    Z[u] += B[u]->n-1; /* taille de B[u] moins le sommet u */
    Z[u] += Z1[u];     /* table des landmarks hors B[u] */
    Z[u] += k-1;       /* taille de W (rien pour la couleur 0) */
    Z[u] += Z2[u];     /* table des sommets dont le hash est la couleur de u hors de B[u] */
  }

  /* affiche taille min/max et moyenne/écart type des différentes tables */
  MINMAXMOY(B[_i]->n-1,n,1,"table B size"); /* on retire u de B[u] pour le routage */
  if(agmnt) MINMAXMOY(K[_i],n,1,"inverse table B size");
  MINMAXMOY(Z1[_i],n,1,"landmark table size");
  MINMAXMOY(k-1,n,1,"table W size");
  MINMAXMOY(Z2[_i],n,1,"own color/hash table size");

  /* pointeurs qui ne servent plus à rien */
  free(Z1);
  free(Z2);
  if(agmnt) free(K);
  else{ // ne sert plus pour DCR
    free(F);
    F=NULL;
  }

  /* affiche la distribution des tailles de table */
  PrintDistribution(Z,n,10,"routing table size");
  free(Z); /* ne sert plus à rien */
  printf("- theoretical average: 2√(n*ln(n*ln(n))) = %i\n",(int)(2*sqrt(n*log(n*log(n)))));
  BARRE;
  
  /* assemble les tables en une seule pour le retour */  
  NALLOC(rs_dcr_tables,RT,1);
  RT->B=B;
  RT->W=W;
  RT->S=S;
  RT->H=H;
  RT->C=C;
  RT->n=n;
  RT->CONT=CONT; // pour AGMNT seulement, NULL sinon
  RT->F=F; // pour AGMNT seulement, NULL sinon
  ALLOCZ(RT->dist,n,S[_i]?S[_i]->D:NULL); /* distances partielles */

  printf("total time: %s\n",TopChrono(2));
  return RT;
}


int rs_dcr_length_rec(int u,int v,int w,int a,rs_dcr_tables *X){
/*
  Fonction récursive donnant la longueur de la route de u à v avec les
  tables X initialisées par dcr, et l'en-tête (w,a):
  - w=landmark intermédiaire (s ou t, w=s par défaut)
  - a=nca(u,v) dans l'arbre du landmark w (a<0 par défaut)
*/
  
  DEBUG(printf("  \nu=%i v=%i: ",u,v););

  // on est arrivé
  if(u==v){ DEBUG(printf("u=v\n");); return 0; }
  int i;

  // v dans B[u] ? (NB: B[u] existe toujours)
  
  i=SetSearch(v,X->B[u]->node,X->B[u]->n,1);
  // si v dans B[u]:u -> v
  if(i>=0){ DEBUG(printf("in B[u]\n");); return X->B[u]->dist[i];}
  
  // v n'est pas dans B[u]
  // est-ce que v est dans la table des landmarks ?

  if(X->C[v]==0){ // 0=couleur des landmarks
    DEBUG(printf("%i is a landmark\n",v););
    return X->S[v]->D[u]; // dist(u,v), v=landmark
  }

  // v n'est pas landmark
  // est-on arrivé au responsable de v ?
  
  if((a<0)&&(X->C[u]!=X->H[v])){ // si pas encore au responsable
    u=X->W[u]->node[X->H[v]]; // next-hop de u vers le responsable de v
    DEBUG(printf("go to %i the closest node of color %i, the hash of %i (a=%i and w=%i)",u,X->H[v],v,a,w););
    return 1 + rs_dcr_length_rec(u,v,w,a,X);
  }

  // on est passé par (ou on est sur) le responsable de v ?
  
  if(a<0){ // on est arrivé au responsable de v
    int w1=X->B[w]->vpd; // landmark le +proche de la source (=w si a<0)
    int w2=X->B[v]->vpd; // landmark le +proche de cible (=v)
    param_bfs *Y1=X->S[w1]; // l'arbre de w1
    param_bfs *Y2=X->S[w2]; // l'arbre de w2
    int a1=nca_bfs(u,v,Y1); // ancêtre commun entre u et v dans Y1
    int a2=nca_bfs(u,v,Y2); // ancêtre commun entre u et v dans Y2
    int d1=dist_nca(u,v,a1,Y1->D); // d1=dist(u,v) dans Y1
    int d2=dist_nca(u,v,a2,Y2->D); // d2=dist(u,v) dans Y2
    if(d1<d2) w=w1,a=a1; else w=w2,a=a2;
    
    /* la ligne suivante est une optimisation: si d1=d2, alors on a
       intérêt de choisir l'ancêtre a1 ou a2 le plus loin de u cela
       laisse plus de chances à l'algo de court-circuiter le routage
       dans l'arbre. */
    if((d1==d2)&&(Y1->D[u]-Y1->D[a1]>=Y2->D[u]-Y2->D[a2])) w=w1,a=a1;
    
    DEBUG(
	  printf("node in charge of %i reached\n",v);
	  printf("go to a=%i, the nca of %i and %i in the tree of landmark w=%i\n",a,u,v,w);
	  );
  }

  // ici on a:
  //  w=landmark intermédiaire (landmark de s ou de t)
  //  a=ancêtre commun entre u et v dans l'arbre de w
  if(u==a){ // si u est arrivé au bon ancêtre
    DEBUG(printf("ancestor a=%i reached\n",a););
    return X->S[w]->D[v] - X->S[w]->D[a];
  }

  // ici il faudrait ajouter les racourcis, si un des ancêtres de v
  // dans l'arbre de w (stockés dans l'étiquette de v) est dans B[u]
  
  u=X->S[w]->P[u]; // on remonte dans l'arbre de w
  DEBUG(printf("go up to a=%i in the tree of landmark w=%i",a,w););
  return 1 + rs_dcr_length_rec(u,v,w,a,X);
}


int rs_dcr_length(int u,int v,rs_dcr_tables *X){
/*
  Renvoie le nombre de sauts du routage selon les tables générées par
  rs_dcr() pour router un message de u à v, ou bien -1 si la route n'a
  pu être déterminée. Dans X on a toutes les tables nécessaires au
  schéma. Si u<0, on réalise quelques tests de bases sur les tables
  X. La fonction fait essentiellement appel à une fonction recursive
  où la gestion d'un en-tête est nécessaire.
*/

  if(u<0){
    if(X==NULL) return 1;
    if(X->B==NULL) return 1;
    if(X->W==NULL) return 1;
    if(X->S==NULL) return 1;
    if(X->H==NULL) return 1;
    if(X->C==NULL) return 1;
    return 0;
  }

  return rs_dcr_length_rec(u,v,u,-1,X);
}


int rs_agmnt_length_rec(int u,int v,int s,int w,rs_dcr_tables *X){
/*
  Fonction récursive donnant la longueur de la route de u à v avec les
  tables X initialisées par agmnt, et l'en-tête (s,w):
  - s=landmark pour un routage via landmark (<0 par défaut)
  - w=nca(u,v,T_s) si s>=0

  L'algorithme est le suivant:
  1. si v est dans la boule de u, alors on route directement vers v
  2. si v est un landmark, alors on route directement vers v
  3. sinon, on route vers le plus proche sommet de la boule u dont
     la couleur = hash(v). Une fois en ce sommet on route via la
     meilleure des deux options:
     - via un des landmarks (via le nca w)
     - via une boule contiguë
  Le point 3 nécessite la gestion d'un en-tête.
*/
  
  DEBUG(printf("  \nu=%i v=%i: ",u,v););

  // on est arrivé
  if(u==v){ DEBUG(printf("u=v\n");); return 0; }
  int i;

  // v dans B[u] ? (NB: B[u] existe toujours)
  
  i=SetSearch(v,X->B[u]->node,X->B[u]->n,1);
  // si v dans B[u]: u -> v
  if(i>=0){ DEBUG(printf("v is in B[u]\n");); return X->B[u]->dist[i]; }
  
  // v n'est pas dans B[u]
  // est-ce que v est un landmark ? (dans la table des landmarks)

  if(X->C[v]==0){ // 0=couleur des landmarks
    DEBUG(printf("%i is a landmark\n",v););
    return X->S[v]->D[u]; // dist(u,v), v=landmark
  }

  // v n'est pas landmark
  // est-on arrivé au responsable de v ?

  if((s<0)&&(X->C[u]!=X->H[v])){ // si pas encore au responsable
    u=X->W[u]->node[X->H[v]]; // next-hop de u vers le responsable de v
    DEBUG(
	  printf("go to %i the next-hop to the closest node of color %i, the hash of %i",
		 u,X->H[v],v);
	  );
    return 1 + rs_agmnt_length_rec(u,v,s,w,X);
  }

  // on est passé par (ou on est sur) le responsable de v
  // est-ce qu'on vient d'arrivé sur le responsable ?
  
  if(s<0){ // on arrive au responsable de v
    // ici v n'est pas landmark ni dans la boule de u
    // cherche i tq v=CONT[u][i].v (cherche binaire, il doit y être)
    contigue cont; cont.v=v;
    contigue *p=bsearch(&cont,X->CONT[u],X->F[X->C[u]],sizeof(contigue),fcmp_contigue);
    DEBUG(
	  if(p==NULL) return FAIL_ROUTING; // ne devrait pas arriver car v est dans CONT[u]
	  );
    cont=*p; // ici cont.v=v
    if(cont.s<0){
      DEBUG(printf("routage via contigue-ball, d=%i\n",cont.d););
      return cont.d; // routage via boule contiguë
    }
    // routage via landmark s et nca w
    DEBUG(printf("routage via landmark %i and nca %i",cont.s,cont.w););
    return rs_agmnt_length_rec(u,v,cont.s,cont.w,X);
  }

  // on est passé par le responsable de v
  // on monte vers le père du landmark, sauf si on arrive au nca w
  if(u==w){ // w=s est possible
    DEBUG(printf("arrive at nca %i of landmark tree %i\n",w,s););
    return X->S[s]->D[v] - X->S[s]->D[w]; // distance entre w et v dans T_s
  }
  u=X->S[s]->P[u]; // u=père(u), NB: u<>root car u<>w=nca
  DEBUG(printf("routing to the parent of u in the landmark tree %i",s););
  return 1 + rs_agmnt_length_rec(u,v,s,w,X);
}


int rs_agmnt_length(int u,int v,rs_dcr_tables *X){
/*
  Renvoie le nombre de sauts du routage selon les tables générées par
  rs_dcr() pour AGMNT pour router un message de u à v, ou bien -1 si
  la route n'a pu être déterminée. Dans X on a toutes les tables
  nécessaires au schéma. Si u<0, on réalise quelques tests de bases
  sur les tables X. La fonction fait essentiellement appel à une
  fonction recursive où la gestion d'un en-tête est nécessaire.
*/

  if(u<0){
    if(rs_dcr_length(-1,0,X)==1) return 1;
    if(X->CONT==NULL) return 1;
    if(X->F==NULL) return 1;
    return 0;
  }

  return rs_agmnt_length_rec(u,v,-1,-1,X);
}


rs_tzrplg_tables *rs_tzrplg(graph *G,double t){
/*
  Calcule les tables de routage selon le schéma tz_rplg pour le graphe
  G, une adaptation du schéma de Thorup et Zwick avec les landmarks
  sur les sommets de plus haut degré.

  Cet algorithme est spécialisé pour les graphes de type RPLG(n,t). Le
  paramètre t attendu (power-law exponent) est celui du graphe G
  fournit. Les performances de ce schéma sont meilleures si le
  paramètre t est le bon, mais marche quel que soit t>1.5. Les valeurs
  de t et de VARIANT permettent de calculer le nombre de landmarks.
  L'ensemble des landmarks est aussi nommé "core" du graphe par Chung
  & Lu.

  On calcule deux tables: B et L. La table B[u], pour chaque sommet u,
  contient tous les sommets strictement plus proche que son plus
  proche landmark. La table L[u], pour chaque sommet u, contient le
  next-hop vers de u vers le landmark le plus proche de u selon
  l'arbre de +cc enraciné dans le landmark.

  Rem: dans [CSTW12], pour les BC graphs, pour les moyennes ils
  prennent probablement n=10K et non n=+grande composante connexe
  (7K), ce qui change les choses.
*/

  printf("\nTZ RPLG\n");
  BARRE;

  TopChrono(1); /* reset du chrono tmp */
  TopChrono(2); /* reset du chrono total */

  const int n=G->n;
  double gamma; /* gamma = paramètre fonction de t */
  int core_size; /* core_size = nombre de landmarks */
  int i,u;

  /*** Calcul des landmarks ***/
  /* C[i]=liste des landmarks, i=0..core_size-1 */

  // calcule D = liste triée (ordre croissant) des degrés des sommets
  // ou permutation aléatoire des sommets (si VARIANT=2)
  int *D;
  if(VARIANT<2) D=SortInt(G->d,NULL,n,1,NULL,SORT_INDEXi);
  if(VARIANT==2){ // permutation aléatoire
    ALLOCZ(D,n,_i);
    Permute(D,n);
  }

  // calcule core_size, dépend de VARIANT et de t
  if((VARIANT==0)&&(t>1.5)){
    gamma=(double)(t-2.0)/(double)(2.0*t-3.0);
    core_size=ceil(pow(n,gamma));
  }
  if((VARIANT==1)&&(t>1.5)){
    // C = { u : deg(u)>(n^gamma')/4 }, gamma'=(1-gamma)/(t-1)=1/(2t-3) 
    gamma=1.0/(double)(2.0*t-3.0);
    u=ceil(pow(n,gamma));
    // cherche le 1er sommet de degré <= u
    i=n-1; // part de la fin (haut degré)
    while((i>=0)&&(G->d[D[i--]]>u));
    core_size=i/4;
  }
  if((VARIANT==2)&&(t>0)) core_size=(int)t;
  if(t<0) core_size=(int)(-t);
  if(t==0) core_size=ceil(sqrt(n));

  core_size=imax(core_size,0);    // core_size >= 0
  core_size=imin(core_size,n); // core_size <= n

  // C=liste des landmarks
  NALLOCZ(int,C,core_size,D[n-_i-1]); // de +haut degré (ou sommets aléatoires)
  
  // affiche C et sa densité
  for(i=u=0;i<core_size;i++) u += G->d[C[i]];
  printf("- core degree: ");
  APERCU(G->d[D[n-_i-1]],core_size,10,2);
  free(D); /* ne sert plus à rien */
  printf("- core size: %i",core_size);
  if((VARIANT==0)&&(t>1.5)) printf(" (n^%g)",gamma);
  if((VARIANT==1)&&(t>1.5)) printf(" (deg>n^%g)",gamma);
  if(t==0) printf(" (sqrt(n))");
  printf("\n");
  printf("- sum of cluster's degrees: %i",u);
  if(u*100>n) printf(" (%.2lfn)",(double)u/(double)n);
  printf("\n- time to construct core: %s\n",TopChrono(1));
  
  /*** Construction des tables L ***/

  // L[u]->node[i] = identifiant du i-ème landmark
  // L[u]->dist[i] = distance entre u et L[u]->node[i]
  // lmin[u] = distance entre u et son plus proche landmark
  // label[u] = indice i du landmark le plus proche de u

  NALLOCZ(table*,L,n,new_table(core_size)); // n tables L de taille core_size
  NALLOCZ(int,lmin,n,n); // par défaut lmin[u]=n
  NALLOC(int,label,n);
  
  int l;
  param_bfs *X;
  
  for(i=0;i<core_size;i++){ // on fait un bfs() pour chaque landmark
    l=C[i]; // l=i-ème landmark
    X=bfs(G,l,NULL); // bfs() depuis l
    for(u=0;u<n;u++){
      L[u]->node[i]=l;
      L[u]->dist[i]=X->D[u];  // dist(u,l)
      if(X->D[u]<lmin[u]){ // landmark plus proche ?
	lmin[u]=X->D[u]; // NB: le rayon de B[u] sera lmin[u]-1
	label[u]=i; /* NB: L[u]->node[label[u]]=landmark le plus proche de u */ 
      }
    }
    free_param_bfs(X);
  }
  free(C);
  printf("- time to construct tables L: %s\n",TopChrono(1));
  
  /*** Construction des tables B ***/
  
  NALLOCZ(table*,B,n,NULL); // tableau de n tables B (vides au départ)
  X=new_param_bfs(); // pour bfs() depuis u
  X->clean=1; // initialisation complète des distances, puis partielle
  
  /* construit la boule B[u] de rayon lmin[u]-1 avec bfs(u,.) */

  for(u=0;u<n;u++){
    X->hmax=imax(0,lmin[u]-1); // la boule contient les sommets
			      // strictement plus proche que le coeur
    B[u]=new_table(0);
    // si B n'est pas vide, alors on fait un BFS pour calculer B
    if(X->hmax){
      bfs(G,u,X);
      B[u]->n=(X->n)-1; // taille de B[u] (sans u)
      /* copie les sommets (et leur distance) du bfs() en supprimant u */
      ALLOCZ(B[u]->node,B[u]->n,X->file[_i+1]); /* =i-ème sommet de B[u] */
      ALLOCZ(B[u]->dist,B[u]->n,X->D[B[u]->node[_i]]); /* =distance du i-ème à u */
    }
    // si B est vide alors:
    else B[u]->n=0;
  }
  free_param_bfs(X);
  free(lmin); /* ne sert plus à rien */
  
  printf("- time to construct tables B: %s\n",TopChrono(1));
  
  /* tri des tables B */
  for(u=0;u<n;u++) if(B[u]) QSORT2(B[u]->node,B[u]->dist,B[u]->n,fcmp_int);
  printf("- time to sort tables B: %s\n",TopChrono(1));
  BARRE;

  /* taille totale des tables */
  NALLOCZ(int,Z,n,1); /* taille=1 au moins pour chaque sommet u */
  for(u=0;u<n;u++){ /* pour chaque sommet u */
    if(B[u]) Z[u] += B[u]->n;
    if(L[u]) Z[u] += L[u]->n;
  }
  
  /* affiche taille min/max et moyenne/écart type des différentes tables */
  MINMAXMOY(B[_i]->n,n,B[_i],"table B size");
  MINMAXMOY(L[_i]->n,n,L[_i],"table L size");

  /* affiche la distribution des tailles de table */
  PrintDistribution(Z,n,10,"routing table size");
  free(Z); /* ne sert plus à rien */
  BARRE;
  
  /* assemble les tables en une seule pour le retour */
  NALLOC(rs_tzrplg_tables,RT,1);
  RT->B=B;
  RT->L=L;
  RT->label=label;
  RT->n=n; // besoin pour libérér les n tables
  
  printf("total time: %s\n",TopChrono(2));
  return RT;
}


int rs_tzrplg_length(int u,int v,rs_tzrplg_tables *X){
/*
  Renvoie le nombre de sauts du routage selon les tables générées par
  rs_tzrplg() pour router un message de u à v, ou bien -1 si la route
  n'a pu être déterminée. Dans X on a toutes les tables nécessaire au
  schéma, notamment les tables B et L. Si u<0 alors on teste la
  validité des tables (et on renvoie une valeur non-nulle en cas
  d'erreur).
*/

  if(u<0){
    if(X==NULL) return 1;
    if(X->B==NULL) return 1;
    if(X->L==NULL) return 1;
    return 0;
  }

  DEBUG(printf("  \nu=%i v=%i: ",u,v););

  // on est arrivé
  if(u==v){ DEBUG(printf("u=v\n");); return 0; }

  // routage dans la boule de u
  // v dans B[u] ?

  if(X->B[u]){
    int i=SetSearch(v,X->B[u]->node,X->B[u]->n,1);
    // si v dans B[u]: u -> v
    if(i>=0){
      DEBUG(
	    printf("in B[u], distance %i: (", X->B[u]->dist[i]);
	    int _v;
	    for(_v=0;_v<X->B[u]->n;_v++)
	      printf("%i,", X->B[u]->node[_v]);
	    printf(")\n");
	    );
      return X->B[u]->dist[i];
    }
  }

  // routage via le plus proche landmark de v
  // route: u -> lv -> v (aucun raccourci, voir Algo. 1 dans CSTW12)

  if(X->L[u]){
    // On doit récupérer 2 entrées: L[u][lv] et L[v][lv]. Pour cela on
    // doit d'abord trouver l'identité de lv qui se trouve dans
    // l'étiquette de v:
    int lv = X->label[v];
    if(X->L[v])
        return X->L[u]->dist[lv] + X->L[v]->dist[lv];
  }

  // y'a un problème, pas d'entrée pour v
  DEBUG(printf("fail: no route found when routing from %i to %i \n",u,v););
  return FAIL_ROUTING;
}


rs_bc_tables *rs_bc(graph *G, int k){
/*
  Schéma de routage selon Brady-Cowen 2006.  Le stretch est additif
  est <= 2k. En particulier, si k=0, il s'agit d'un routage de plus
  court chemin.

  Principe: On construit un arbre BFS (=T) enraciné dans le sommet de
  plus haut degré (=center). Le coeur (=C) est la boule de rayon k
  depuis la racine de T. Dans l'article d'origine, k=d/2 avec d pair.

  On construit une liste (=L) de BFS couvrant G ainsi qu'une forêt
  (=H) de BFS de G comme suit. Au départ, L={T}, et H est la forêt
  T\C. Puis, pour chaque arête {u,v} de G\C\T, on vérifie si l'ajoût
  de {u,v} à H crée un cycle ou pas. Si on ne crée pas de cycle, on
  met à jour la forêt H en lui ajoutant {u,v}. Si on crée un cycle, on
  calcule un BFS de G de racine u (ou v) qu'on ajoute à L (on favorise
  le sommet de plus grand degré sans les arêtes de T).

  Une fois toutes les arêtes {u,v} ainsi balayées, on cacule pour
  chaque composantes connexes de H un BFS. On obtient une forêt
  couvrante qu'on ajoute à L.

  L'algorithme de routage de u à v consiste simplement à router dans
  l'arbre A de L contenant u et v et qui minimise dist_A(u,v).
*/

  const int n=G->n;
  int u,v,center,x,y,d,i;

  printf("\nBRADY-COWEN\n");
  BARRE;

  TopChrono(1); /* reset du chrono tmp */
  TopChrono(2); /* reset du chrono total */

  /* trouve un sommet center de degré max, v=deg(center) */
  for(u=v=center=0;u<n;u++) if(G->d[u]>v) v=G->d[center=u];
  printf("- degree of the center: %i (id:%i)\n",v,center);

  /* calcule l'arbre T (bfs) */

  // on construit l'arbre en deux temps, d'abord on récupère le coeur
  // de rayon k, puis on poursuit en couvrant tout le graphe

  param_bfs* T_bfs=new_param_bfs(); // arbre de racine center
  T_bfs->clean=1;
  T_bfs->hmax=k;
  bfs(G,center,T_bfs);
  int core_size=T_bfs->n; // taille du coeur
  NALLOCZ(int,C,core_size,T_bfs->file[_i]); // C=liste des sommets du coeur
  printf("- core: "); APERCU(C[_i],T_bfs->n,10,2);
  printf("- core size: %i\n",core_size);

  // on finit la construction de T
  T_bfs->cont=1; 
  T_bfs->hmax=-1; // le bfs() va jusqu'au bout
  bfs(G,center,T_bfs); // T_bfs=bfs() de tout G
  printf("- excentricity of the center: %i\n",T_bfs->radius);
  printf("- time to construct the core and T: %s\n",TopChrono(1));

  // S[u]=-2 ssi u est dans le coeur, -1 sinon
  NALLOCZ(int,S,n,-1);
  for(i=0;i<core_size;i++) S[C[i]]=-2;
  free(C);

  /* calcule le graphe H */

  // au début H=T\C, ajoute à H les arête de T
  graph* H=new_subgraph(G); // on créer un graphe vide H de la taille de G
  for(u=0;u<n;u++){
    x=T_bfs->P[u]; // x=père de u dans T
    // ajoute {u,x} à H si u et son père x ne sont pas dans le coeur
    if((x>=0)&&(S[x]!=-2)&&(S[u]!=-2)) ADD_EDGE(H,u,x);
  }

  /* calcule les composantes connexes de H, via un DFS */
  param_dfs* H_dfs=new_param_dfs(n);
  ALLOCZ(H_dfs->C,n,S[_i]); // copie S dans ->C pour les sommets interdits
  H_dfs=dfs(H,0,H_dfs); // DFS depuis un sommet arbitraire
  const int nc=H_dfs->nc; // sauvegarde le nombre de composantes de H
  printf("- time to construct and to traverse (dfs) T\\C: %s\n",TopChrono(1));
  
  /* convertit les composantes en "représentant" pour FindSet() */
  // au départ couleur[u] est un indice
  // à la fin couleur[u] est un sommet source de dfs sur H

  int *couleur=H_dfs->C; // couleur[u]=componsante du sommet u de H
  for(u=0;u<n;u++){
    if(S[u]==-2) continue;
    couleur[u]=H_dfs->R[couleur[u]];
  }

  /*** calcule la liste L des BFS ***/

  // pour toutes arêtes uv de G, ajouter uv à H si cela ne crée pas de
  // cycle, sinon on calcule un BFS enraciné en u ou v (deg max dans
  // G\H) qu'on ajoute à L

  NALLOCZ(int,rang,n,0); // pour l'heuristic de FindSet()
  NALLOC(param_bfs*,L,n); // L=liste des BFS, L[i]=i-ème BFS
  NALLOCZ(param_bfs*,L_bfs,n,NULL); // L_bfs[u]=BFS pour le sommet u
  L[0]=L_bfs[center]=T_bfs; // on met T_bfs dans L et L_bfs
  int nbfs=1; // nombre de BFS dans la liste L
  int no_cycle=0; // nombre d'arêtes ne créeant pas de cycle

  for(u=0;u<n;u++){ // pour tous les sommets u de G
    if(S[u]==-2) continue; // ne rien faire si u est dans le coeur
    d=G->d[u]; // d=degré(u) dans G
    for(i=0;i<d;i++){
      v=G->L[u][i]; // v=i-ème voisin de u
      if(S[v]==-2) continue; // // ne rien faire si v est dans le coeur
      if((T_bfs->P[u]==v)||(T_bfs->P[v]==u)) continue; // ne rien faire si {u,v} dans E(T)
      if(u>=v) continue;
      // ici u<v
      couleur[u]=x=FindSet(u,couleur); // x=représentant de u
      couleur[v]=y=FindSet(v,couleur); // y=représentant de v
      if(x==y){ // si même représentant alors on fait un BFS depuis u ou v
	if((L_bfs[u])||(L_bfs[v])) continue; // BFS déjà calculé
	// ici, on a jamais fait de BFS ni depuis u ni depuis v
	// on enracine le BFS en u ou v selon le degré de G-H
	x=(G->d[u]-H->d[u] > G->d[v]-H->d[v])? u : v;
	L[nbfs++]=L_bfs[x]=bfs(G,x,NULL); // met à jour L et L_bfs
      }
      else{ // si pas même couleur alors pas de cycle
	ADD_EDGE(H,u,v);
	no_cycle++;
	UNION(x,y,couleur,rang); // fusion pour FindSet()
      }
    }
  }

  free(rang);
  x=no_cycle+nbfs-1; // x=nb d'arêtes dans G\C\T
  printf("- #edges in C: %i\n",NbEdges(G)-x-n+core_size);
  printf("- #edges in G\\C\\T: %i\n",x);
  printf("- #edges added to T\\C to make H: %i\n",no_cycle);
  printf("- #bfs trees computed from edges not in H: %i\n",nbfs-1);
  printf("- time to construct H and these bfs trees: %s\n",TopChrono(1));

  /*** calcule un BFS pour chaque composante de H ***/

  param_bfs* H_bfs=new_param_bfs();
  H_bfs->clean=0; // très important pour faire l'union de BFS (ne pas initialiser H_bfs->D)
  H_bfs->D=S; // pour sélectionner les sommets de G\S (surtout ne pas faire free(S) ...)
  
  // lance un BFS depuis chaque composante de H, à partir des racines
  // R[i] (du dfs de H) qui sont leur propre représentant
  for(i=x=0;i<nc;i++){ // x=nombre de BFS calculés dans H
    u=H_dfs->R[i];
    if(u==FindSet(u,couleur)){
      bfs(H,u,H_bfs); // superpose les BFS de H
      x++;
    }
  }
  printf("- #components of T\\C: %i\n",nc);
  printf("- #components of H: %i\n",x);
  printf("- time to construct the bfs forest for H: %s\n",TopChrono(1));
  free_param_dfs(H_dfs);
  free_graph(H);
  BARRE;
  
  /* taille totale des tables */
  NALLOC(int,Z,n); /* Z[u]=taille de la table de u, doit être au moins 1 */
  for(u=0;u<n;u++){ /* pour chaque sommet u */
    Z[u]=nbfs+1; // doit être au moins 1
    if(S[u]!=-2) Z[u]++; // si pas dans le coeur (à cause de H)
  }
  
  /* affiche la distribution des tailles de table */
  PrintDistribution(Z,n,10,"routing table size");
  free(Z); /* ne sert plus à rien */
  BARRE;
  
  /* transforme H en véritable forêt couvrante (les champs ->P et ->D
     des sommets du coeur n'étant pas forcément cohérent), puis
     l'ajoute à la liste L. Il faut bien sûr que la forêt ne soit pas
     vide. */

  if(nc){ // si H possède au moins un sommet, alors nc>0
    // corrige H_bfs pour qu'il soit une forêt couvrante, sinon
    // dist_bfs() dans _length() peuvent ne pas marcher
    // puis ajoute H_bfs à L
    for(u=0;u<n;u++)
      if(S[u]==-2){
	H_bfs->P[u]=-1;
	H_bfs->D[u]=0;
      }
    L[nbfs++]=H_bfs;
  }
  REALLOC(L,nbfs); // réduit L

  /* assemble les tables en une seule pour le retour */
  NALLOC(rs_bc_tables,RT,1);
  RT->L=L;
  RT->Lu=L_bfs;
  RT->nbfs=nbfs;
  ALLOCZ(RT->dist,n,L_bfs[_i]?L_bfs[_i]->D:NULL); // distances partielles

  printf("total time: %s\n",TopChrono(2));
  return RT;
}


int rs_bc_length(int u,int v,rs_bc_tables *X){
/*
  Renvoie le nombre de sauts du routage selon les tables générées par
  rs_bc() pour router un message de u à v, ou bien -1 si la route n'a
  pu être déterminée. Dans X on a toutes les tables nécessaire au
  schéma, notamment les tables ... Si u<0 alors on teste la validité
  des tables (et on renvoie une valeur non-nulle en cas d'erreur).

  Rem: on pourrait obtenir des routes plus courtes en faisant du
  pas-à-pas plutôt que de router dans l'arbre.
*/
  if(u<0){
    if(X==NULL) return 1;
    if(X->L==NULL) return 1;
    if(X->Lu==NULL) return 1;
    return 0;
  }

  DEBUG(printf("Routing from %i to %i\n",u,v););

  // 1) on est arrivé
  if(u==v) return 0;

  // 2) est-ce que u ou v est la racine d'un des BFS de L ?
  if(X->Lu[u]) return X->Lu[u]->D[v];
  if(X->Lu[v]) return X->Lu[v]->D[u];
  DEBUG(
	printf("ni u ni v ne sont racines d'un BFS\n");
	PRINT(X->nbfs);
	);

  // 3) calcule la distance dans entre u et v dans le meilleur arbre de L
  int d,i,d0=INT_MAX; // d0=distance min recherchée
  for(i=0;i<X->nbfs;i++){ // parcours les arbres de L
    d=dist_bfs(u,v,X->L[i]); // d=dist(u,v) dans le i-ème arbre de L
    if(d>=0) d0=imin(d0,d); // met à jour la distance min
  }
  
  if(d0<INT_MAX) return d0;
  return FAIL_ROUTING;
}


rs_hdlbr_tables *rs_hdlbr(graph *G,int k){
/*
  Calcule les tables de routage selon le schéma HDLBR pour le graphe
  G, une adaptation du schéma tz_rplg en name-independant.

  Cet algorithme est spécialisé pour les graphes de type RPLG(n,t).
  L'algorithme HDLBR, à l'instar de TZ_RPLG, est paramétré par défaut
  à un nombre de landmarks fixé à k=n^x=n^1/2:
  
    " We set x=1/2 because this setting minimizes the storage overhead, 
    making routing table size of landmarks and average routing table 
    size of non-landmark nodes both bounded by Õ(n^1/2) bits. "
*/

  printf("\nHDLBR\n");
  BARRE;

  TopChrono(1); /* reset du chrono tmp */
  TopChrono(2); /* reset du chrono total */

  const int n=G->n;
  int i,u,v;

  /*** Calcul des landmarks ***/
  // Core=liste des landmarks, les k sommets de plus haut degré

  // récupère la liste des identifiants des noeuds triée par degré
  int *D=SortInt(G->d,NULL,n,1,NULL,SORT_INDEXi); /* D=liste triée par degrés croissant des noeuds */

  // Core=liste des landmarks=liste des sommets de plus haut degré
  NALLOCZ(int,Core,k,D[n-_i-1]);

  // Trie le coeur par identifiants (utile pour setSearch lors du routage)
  QSORT(Core,k,fcmp_int);
  printf("- time to construct the Core: %s\n",TopChrono(1));
  BARRE;

  free(D);

  printf("- #landmarks: %i",k);
  printf(" (ceil{√n}=%i)\n",(int)ceil(sqrt(n)));
  printf("- landmark list: ");
  APERCU(Core[_i],k,10,2);
  printf("- landmark degrees: ");
  APERCU(G->d[Core[_i]],k,10,2);
  for(i=u=0;i<k;i++) u += G->d[Core[i]]; // somme des degrés des landmarks
  printf("- average landmark's degree: %.2lf\n",(double)u/(double)n);
  printf("- time to construct landmarks: %s\n",TopChrono(1));
  
  /*** Construction des tables L ***/
  // L[u] est définie pour tout les sommets u
  // L[u]->node[i]=parent du sommet u dans l'arbre BFS du i-ème landmark
  // L[u]->vpd[i]=indice de l(u), le plus proche landmark de u
  
  NALLOCZ(table*,L,n,new_table(k)); // n tables de taille k>0
  // par défaut L[u]->vpd=-1, il faut l'initialiser à 0
  for(u=0;u<n;u++) {
    L[u]->vpd=0; // L[u]->vpd=indice du landmark le plus proche de u
    L[u]->n=k;
    ALLOCZ(L[u]->node,k,-1); 
    ALLOCZ(L[u]->dist,k,-1); 
  }
  param_bfs *X=new_param_bfs(); // structure utilisée pour le résultat du BFS depuis u
  X->clean=1; // initialisation complète des distances, puis partielle

  for(i=0;i<k;i++){
    v=Core[i]; // v=landmark numéro i
    bfs(G,v,X);
    for(u=0;u<n;u++){
      L[u]->node[i]=X->P[u]; // père de u dans l'arbre de racine v   
      L[u]->dist[i]=X->D[u]; // dist(u,v)
      if(X->D[u]<L[u]->dist[L[u]->vpd]) L[u]->vpd=i; // met à jour L[u]->vpd
    }
  }
  free_param_bfs(X);
  printf("- time to construct tables L: %s\n",TopChrono(1));

  // calcule la distance inter-landmark
  for(i=u=0;u<k;u++) // i=distance inter-landmark
    for(v=u+1;v<k;v++)
      i=imax(i,L[Core[u]]->dist[v]);
  printf("- inter-landmark distance: %i\n",i);
  printf("- time to compute this distance: %s\n",TopChrono(1));
    

  /*** Construction des tables B ***/
  // B[u] est définie pour tout sommet u
  // B[u]->node=liste des sommets à distance < dist(u,l(u))
  // B[u]->dist=distance de ces sommets à u
  // B[u] contient en plus l(u) qui est toujours B[u]->node[0]
  // Si u est un landmark, alors B[u]={u}, sinon u n'est pas
  // stocké dans B[u].
  
  NALLOCZ(table*,B,n,new_table(0)); // n tables vides au départ
  X=new_param_bfs(); // structure utilisée pour le résultat du bfs depuis u
  X->clean=1; // initialisation complète des distances, puis partielle
  
  for(u=0;u<n;u++){
    X->hmax=imax(0,L[u]->dist[L[u]->vpd]-1); // hmax=dist(u,l(u))-1
    bfs(G,u,X);
    B[u]->n=X->n; // NB: ici X->node[0]={u}
    // si u est landmark alors B[u]={u}, sinon on remplace dans B[u]
    // le sommet u par l(u). Dans tous les cas, |B[u]| = X->n.
    ALLOCZ(B[u]->node,B[u]->n,X->file[_i]); // B[u]->node[i]=i-ème sommet de B[u]
    ALLOCZ(B[u]->dist,B[u]->n,X->D[B[u]->node[_i]]); // B[u]->dist[i]=distance du i-ème à u
    i=L[u]->vpd; // i=indice de l(u)
    if(Core[i]!=u){ // si u n'est pas un landmark
      B[u]->node[0]=Core[i];
      B[u]->dist[0]=L[u]->dist[i];
    }
    //free(L[u]->dist); // ne sert plus à rien
    //L[u]->dist=NULL;
  }
  free(X);
  printf("- time to construct tables B: %s\n",TopChrono(1));
  
  /*** Calcul de la taille des boules inverse ***/
  // D[u]=taille de la boule inverse de u, c'est-à-dire de B[u] mais
  // sans l(u).

  ALLOCZ(D,n,0); // compteurs à 0
  for(u=0;u<n;u++)
    for(i=1;i<B[u]->n;++i) // i>0, ne compte pas B[u]->node[0]=l(u)
      D[B[u]->node[i]]++;

  printf("- time to compute inverse ball sizes: %s\n",TopChrono(1));

  /* H[u]=0..k-1, hash du sommet u */
  int *H=MakeHash(NULL,n,k,HASH);

  /*
    On ne construit pas la table couleur car elle n'est pas nécessaire
    pour le routage, en distribué elle sert à faire la "traduction" du
    schéma étiquetté vers le schéma name-independant. Ici cette table
    est un entier qui compte le nombre d'entrées qu'elle aurait dû
    avoir.
  */

  int *F;
  FREQMINMAX(F,k,H,n,"the hash"); /* fréquence des hashs (sert pour taille des tables) */
  printf("- time to compute hash values: %s\n",TopChrono(1));

  /* tri des tables B */
  for(u=0;u<n;u++)
    if(B[u]) QSORT2(B[u]->node,B[u]->dist,B[u]->n,fcmp_int);
  printf("- time to sort tables B: %s\n",TopChrono(1));
  BARRE;

  /* taille totale des tables */  
  NALLOCZ(int,Z,n,1); /* taille=1 au moins pour chaque sommet u */
  for(u=0;u<n;u++){ /* pour chaque sommet u */
    if(B[u]) Z[u] += B[u]->n-1; /* -1 car le landmark l(u) est toujours dans B[u] */
    if(L[u]) Z[u] += L[u]->n;
    Z[u] += D[u]; /* vraie taille des boules inverses */ 
    i=L[u]->vpd;
    if(Core[i]==u) Z[u] += F[i]; // si u est un landmark
  }
  
  /* affiche taille min/max et moyenne/écart type des différentes tables */
  MINMAXMOY(B[_i]->n,n,B[_i],"table B size");
  MINMAXMOY(D[_i],n,1,"inverse ball size");
  MINMAXMOY(L[_i]->n,n,L[_i],"table L size");
  MINMAXMOY(F[_i],k,F[_i],"table Color size");
  free(D);
  free(F);

  /* affiche la distribution des tailles de table */
  PrintDistribution(Z,n,10,"routing table size");
  free(Z); /* ne sert plus à rien */
  BARRE;
  
  /* assemble les tables en une seule pour le retour */
  NALLOC(rs_hdlbr_tables,RT,1);
  RT->B=B;
  RT->L=L;
  RT->Core=Core;
  RT->core_size = k;
  RT->H=H;
  RT->n=n; // besoin pour libérér les n tables
  
  printf("total time: %s\n",TopChrono(2));
  return RT;
}


int rs_hdlbr_length_rec(int u,int v,rs_hdlbr_tables *X, int lv){
/*
  Renvoie le nombre de sauts du routage selon les tables générées par
  rs_hdlbr() pour router un message de u à v, ou bien -1 si la route
  n'a pu être déterminée. Dans X on a toutes les tables nécessaire au
  schéma, notamment les tables B, C et L. Si u<0 alors on teste la
  validité des tables (et on renvoie une valeur non-nulle en cas
  d'erreur).

  Le paramètre lv est l'indice du plus proche landmark de, ou -1 s'il
  n'a pas pu être déterminé.

  L'algorithme de routage HDLBR de u à v est le suivant: (comme on est 
  en name-independant on doit parfois trouver l'étiquette de v qui se 
  trouve sur le landmark de couleur h(v)).
    1) si u=v alors on est arrivé
    2) si v appartient à B[u], u appartient à B[v] ou v est un landmark, router vers v en plus court chemin
    3) si lv<0, alors router vers le landmark de couleur h(v)
    4) si u=lv router vers v avec la boule inverse étendue de u
    5) sinon, router vers le landmark de v (lv)
*/

  DEBUG(printf("Routing from %i to %i \n", u,v);)
  // 1) on est arrivé
  if(u==v) return 0;

  /** Routage dans la boule/boule inverse de u **/
  // 2) v dans B[u] ?
  if(X->B[u]){
    int i=SetSearch(v,X->B[u]->node,X->B[u]->n,1);
    // si v dans B[u]: u -> v
    if(i>=0) {
      DEBUG(printf("\t - Direct routing %i in B[%i] (distance %i) \n", v, u,X->B[u]->dist[i]);)
      return X->B[u]->dist[i];
    }
  }
  // 2) u dans B[v] (ie., v dans B^-1[u]) ?
  if(X->B[v]){
    int i=SetSearch(u,X->B[v]->node,X->B[v]->n,1);
    // si u dans B[v]: u -> v
    if(i>=0) {
      DEBUG(printf("\t - Direct routing %i in C[%i] (distance %i) \n", v, u,X->B[v]->dist[i]);)
      return X->B[v]->dist[i];
    }
  }
  // 2) v est un landmark ?
  if(X->L[u]){
    int i=SetSearch(v,X->Core,X->core_size,1);
    if (i>=0) {
      DEBUG(printf("\t - Direct routing %i is a landmark (distance %i)  \n", v, X->L[u]->dist[i]);)
      return X->L[u]->dist[i];
    }
  }

  // 3) Si le header est vide (lv<0), alors on est pas encore passé
  // par h(v) et on veut router un pas vers h(v)
  if(lv<0){
    int hash = X->H[v]; // indice du landmark dans le tableau Core
    int hv = X->Core[hash]; // landmark responsable de v
    // si on est arrivé en h(v), routage vers le plus proche landmark de v
    if(u==hv) {
      DEBUG(printf("\t - Manager hv=%i reached, now going toward lv=%i \n", hv, X->Core[X->L[v]->vpd]););
      return rs_hdlbr_length_rec(u, v, X, X->L[v]->vpd);
    }

    DEBUG(printf("\t - Routing toward v's manager (hv=%i) \n", hv););
    // sinon on monte vers le père de u dans l'arbre de racine h(v)
    return 1 + rs_hdlbr_length_rec(X->L[u]->node[hash],v,X,-1);
  }

  // 4) routage via le plus proche landmark de v
  if(lv>=0){
    DEBUG(printf("\t - Routing toward lv=%i \n", X->Core[lv]););
    return 1 + rs_hdlbr_length_rec(X->L[u]->node[lv],v,X,lv);
  }
  // y'a un problème, pas d'entrée pour v
  DEBUG(printf("fail: no route found when routing from %i to %i \n",u,v););
  return FAIL_ROUTING;
}


int rs_hdlbr_length(int u,int v,rs_hdlbr_tables *X){

  if(u<0){
    if(X==NULL) return 1;
    if(X->B==NULL) return 1;
    if(X->L==NULL) return 1;
    if(X->Core==NULL) return 1;
    if(X->H==NULL) return 1;
    return 0;
  }

  return rs_hdlbr_length_rec(u,v,X,-1); // header vide au départ
}


enum{
  SC_NONE,
  SC_ALL,
  SC_ONE,
  SC_NPAIRS,
  SC_PAIR,
  SC_EDGES,
  SC_UV,
  SC_UNTIL,
};


static inline long route_uv(const graph *G,
			    const int u,const int v,
			    const int h,const int hmax,
			    int **dist,long **stat,int *L,int *M,param_bfs *X){
/*
  Remplit la table de statistiques stat[][] et aussi de distance
  dist[][] des sommets testés avec h=longueur du routage de u à v
  (dépend aussi de SCENARIO.dist). Dans le scenario SC_EDGES les
  distances ne sont pas mise à jour (la distance étant toujours
  1). Renvoie une valeur <0 si une erreur est survenue ou bien la
  distance entre u et v sinon. Le tableau n'est pas remplit si une
  erreur est survenue, en particulier si h>hmax.

  On pourrait calculer ici h=length(u,v,X), plutôt que de le passer
  comme paramètre. On pourrait ainsi calculer avant dist(u,v), et
  ensuite la passer comme paramètre à length(u,v,X,d) qui pourrait
  éventuellement l'exploiter.
*/

  if((h<0)||(h>hmax)){
    DEBUG(printf("FAIL !!!  %i->%i: #hop=%i\n",u,v,h););
    return -1; /* y'a un problème, sommet suivant */
  }

  int j,k,l,t;

  /* calcule k=dist(u,v) */
  if(SCENARIO.mode==SC_EDGES) k=1; /* pas de bfs() dans ce cas */
  else{ /* est-ce que dist(u,v) est déjà connue ? */
    if((dist[u]==NULL)&&(dist[v]==NULL)){
      if(SCENARIO.dist){ /* on calcule et stocke dist[u] */
	bfs(G,u,X);
	ALLOCZ(dist[u],G->n,X->D[_i]); /* alloue dist[u] et y copie X->D[] */
      }else{ /* ici on ne stocke pas dist[u] */
	X->hmax=h; /* borne sup sur la distance de u à v */
	bfs(G,u,X);
      }
      k=X->D[v]; /* k=dist(u,v) */
    }
    else k=dist[u]? dist[u][v] : dist[v][u]; /* k=dist(u,v) */
  }
  
  if(h<k){
    DEBUG(printf("FAIL !!!  %i->%i: #hop=%i dist(u,v)=%i\n",u,v,h,k););
    return -1; /* il y'a un problème, trop court */
  }

  /* ici il faut faire +1 dans stat[k][j] où j=h-k */
  j=h-k;
  if(j>=L[k]){ /* tableau stat[k] n'est pas assez grand */
    l=L[k]; /* sauvegarde L[k] */
    L[k]=imax(j+1,l)<<1; /* on double la taille courante ou de h-k+1 */
    if(stat[k]==NULL) ALLOCZ(stat[k],L[k],0L); /* première fois */
    else{
      REALLOC(stat[k],L[k]); /* aggrandit stat[k] */
      for(t=l;t<L[k];t++) stat[k][t]=0L; /* initialise la fin du tableau */
    }
  }
  stat[k][j]++; /* une route de longueur h=k+j de plus */
  M[k]=imax(M[k],j); /* détour max pour la distance k */

  return k; /* tout c'est bien passé */
}


int pgcd(int a,int b){
/*
  Renvoie le plus grand commun diviseur de a et b. NB: a et b premier
  entre eux ssi |pgcd(a,b)|=1.
*/
  int c;
  while(a){ c=a; a=b%a; b=c; }
  return b;
}


int routing_test(graph *G,void *T,rt_length length,int hmax,int **distp){
/*
  Teste un scenario de routage pour le graphe G avec les tables de
  routage T et la fonction de longueur length(), puis affiche les
  distributions des longueurs de route, de distance et de stretch. Le
  scenario est décrit par la variable globale SCENARIO. La matrice
  distp[][] est une matrice partielle de distances du graphe
  (éventuellement calculées lors de la construction de T), permettant
  d'accélérer les tests. En effet, pour calculer le stretch on doit
  calculer la distance.

  Plus précisément, pour tout sommet u de G, distp[u], si non NULL,
  doit être la distance de u vers tous les sommets de G. Attention !
  un vecteur partiel de distance n'est pas autorisé. Il ne faut pas
  que distp[u] soit construit à l'aide de bfs() partiel par
  exemple. Il est cependant possible d'avoir distp[u]=NULL (car
  matrice partielle <> vecteur partiel). Si distp=NULL, alors aucune
  distance n'est pré-calculées. Le tableau distp n'est pas modifié.
  Seules les vecteurs de distance calculés par routing_test() sont
  libérés à la fin, charge à l'appelant de libérer les vecteurs de
  distp.

  La fonction renvoie une valeur non-nulle si une erreur s'est
  produite, et 0 sinon. La valeur de hmax indique la longueur maximum
  de routage autorisée. Cela permet de détecter des routes trop
  longues pour être correctes (et éviter une boucle infinie). Si
  hmax<0, alors on prend hmax=2n comme valeur par défaut.

  Pour trouver un petit graphe dont le stretch est > 3.00, on peut
  faire le script suivant:

  n=10; a=0; i=0; while [[ $(echo "$a <= 3.00" | bc) -eq 1 ]]; do a=$(./gengraph gabriel $n -seed $i -check routing scenario all agmnt -1 | grep -e "- maximum stretch:" | awk '{print $4}'); echo "n=$n, stretch=$a, seed=$i"; i=$((i+1)); done
*/
  if(SCENARIO.mode==SC_NONE) return 0;

  TopChrono(1);
  printf("\nROUTING TEST\n");
  BARRE;

  /* vérification de la table T */
  if(length(-1,0,T)) return 1;

  const int n=G->n;
  int u,v,i,j,k,t,h;
  int dmax; /* plus longue distance */
  int lmax; /* plus longue route */
  double x; /* pour les affichage de distributions */

  /* type "long" pour le comptage */
  long p,s,s2;
  long err=0L; /* nb de routage erronés */

  param_bfs *X=new_param_bfs(); /* pour le calcul de la distance */
  X->clean=1; /* pour appels multiples */

  /* dist[u][v]=distance entre u et v=0..n-1, n valeurs différentes au plus */
  /* on a dist[u]=NULL si on a pas encore calculé de bfs(u,...) */
  NALLOCZ(int*,dist,n,distp?distp[_i]:NULL); /* dist[u]=distp[u] ou NULL */

  /* stat[k][j]=#routage entre sommets à distance k et de longueur de k+j (j=détour) */
  NALLOCZ(long*,stat,n,NULL); /* stat[k]=NULL si pas encore eut de stat pour dist k */
  /* L[k]=taille du tableau stat[k] */
  NALLOCZ(int,L,n,0); /* taille nulle par défaut */
  /* M[k]=détour maximum (=longueur-k) rencontrée pour la distance k, M[k]<= L[k] */
  NALLOCZ(int,M,n,-1); /* par défaut M[k]<0 */
  if(hmax<0) hmax=n<<1;

  switch(SCENARIO.mode){
    
  case SC_ALL:
    p=(long)n*(long)(n-1); /* cast "(long)" important */
    if(2*lg(n)>(int)(8*sizeof(p)-1)) Erreur(33); /* dépassement arithmétique, p est trop grand */
    printf("- all-to-all pairs: %s pairs ...",millier(p));
    fflush(stdout);
    for(u=0;u<n;u++)
      for(v=0;v<n;v++)
	if(u!=v) err += (route_uv(G,u,v,length(u,v,T),hmax,dist,stat,L,M,X)<0);
    break;

  case SC_UNTIL:
    p=s=(long)n*(long)(n-1); /* cast "(long)" important */
    if(2*lg(n)>(int)(8*sizeof(p)-1)) Erreur(33); /* dépassement arithmétique, p est trop grand */
    printf("- all pairs until stretch ≥ %g ...",SCENARIO.stretch);
    fflush(stdout);
    t=random()%n;
    for(i=0;i<n;i++)
      for(j=0;j<n;j++){
	if(i==j) continue;
	u=(i+t)%n,v=(j+t)%n;
	h=length(u,v,T);
	k=route_uv(G,u,v,h,hmax,dist,stat,L,M,X); /* k=dist(u,v) */
	if(k<0){ err++; continue; }
	if(h>=(k*SCENARIO.stretch)){ /* stop: stretch max atteint */
	  s=(long)i*(long)(n-1)+(long)(j+(j<i)); /* s=#paires testées */
	  i=j=n; /* pour arrêter les deux boucles. NB: après i>n */
	}
      }
    break;

  case SC_PAIR:
  case SC_NPAIRS:
    x=(SCENARIO.mode==SC_NPAIRS);
    p=x?n:SCENARIO.u;
    printf("- %srandom pairs: %s pair%s ...",x?"n ":"",millier(p),PLURIEL(p));
    fflush(stdout);
    for(i=0;i<p;i++){
      u=random()%n;
      v=random()%(n-1); // NB: n>=2 car au moins une arête
      if(v>=u) v++; // u,v aléatoires dans [0,n[ avec u<>v
      err += (route_uv(G,u,v,length(u,v,T),hmax,dist,stat,L,M,X)<0);
    }
    break;
  
  case SC_EDGES:
    p=NbEdges(G)<<1;
    printf("- all neighbor pairs: %s pair%s ...",millier(p),PLURIEL(p));
    fflush(stdout);
    for(u=0;u<n;u++){
      p=G->d[u];
      for(i=0;i<p;i++){
	v=G->L[u][i];
        err += (route_uv(G,u,v,length(u,v,T),hmax,dist,stat,L,M,X)<0);
      }
    }
    break;
  
  case SC_ONE:
    u=SCENARIO.u;
    if(u>=n) Erreur(17);
    p=n-1;
    if(u<0) u=random()%n;
    printf("- one-to-all pairs from %i: %s pair%s ...",u,millier(p),PLURIEL(p));
    fflush(stdout);
    for(v=0;v<n;v++)
      if(u!=v) err += (route_uv(G,u,v,length(u,v,T),hmax,dist,stat,L,M,X)<0);
    break;
  
  case SC_UV:
    u=SCENARIO.u;
    v=SCENARIO.v;
    if((u>=n)||(v>=n)) Erreur(17);
    if(u<0) u=random()%n;
    if(v<0) v=random()%n;
    printf("- routing from %i to %i: 1 pair ...",u,v);
    fflush(stdout);
    err += (route_uv(G,u,v,length(u,v,T),hmax,dist,stat,L,M,X)<0);
    break;

  default:
    Erreur(23);
  }

  printf(" %s (%s)\n",((SCENARIO.mode==SC_UNTIL)&&(i==n))?"not found":"Ok",TopChrono(1));
  if(SCENARIO.mode==SC_UNTIL){
    printf("- #tested routings: %s",millier(s));
    if(s==p) printf(" (all the pairs)");
    printf("\n");
    p=s; /* maintenant p=nombre de paires testées */
  }
  
  /* libère les tableaux devenus inutiles: X et L */
  free_param_bfs(X);
  free(L);
  
  /* on ne libère que les distances calculées par routing_test() */
  for(u=0;u<n;u++) if((distp==NULL)||(distp[u]==NULL)) free(dist[u]);
  free(dist);

  /*
    Affiche la distribution:
    - des longueurs de route,
    - des distances,
    - des stretch.
  */

  /* détermine dmax (=distance max) et lmax (=longueur routage max) */
  dmax=lmax=0;
  for(k=0;k<n;k++)
    if(stat[k]){
      dmax=imax(dmax,k);
      lmax=imax(lmax,k+M[k]);
  }

  /* pour un petit gain mémoire: réajuste les tableaux dépendant de la distance k */
  /* NB: le tableau L n'existe plus */
  REALLOC(stat,dmax+1); /* stat[0..dmax] */
  REALLOC(M,dmax+1);    /*    M[0..dmax] */

  /* calcule le nombre p de routages corrects (non erronés) */
  for(k=0,p=0L;k<=dmax;k++)
    if(stat[k])
      for(j=0;j<=M[k];j++)
	p += stat[k][j];

  printf("- #failed routings: %s",millier(err));
  if(err) printf(" (%i%%)",p?(int)ceil((100.0*err)/(double)(p+err)):100);
  printf("\n");
  if(p==0L){
    printf("  all routings failed!\n");
    goto fin_routing_test;
  }
  /* NB: ici p>0 */
  
  BARRE;
  printf("- route length distribution:\n");

  /* calcule F[i]=nombre de routages de longueur i */
  NALLOCZ(long,F,lmax+1,0L);
  for(k=0;k<=dmax;k++)
    if(stat[k])
      for(j=0;j<=M[k];j++)
	F[k+j] += stat[k][j];
  
  /* s=somme totale des longueurs, s2=somme du carré des longueurs */
  for(i=0,s=s2=0L;i<=lmax;i++) s+=((long)i)*F[i],s2+=((long)i)*((long)i)*F[i];

  for(i=0;i<=lmax;i++){
    if(F[i]==0L) continue;
    x=(double)F[i]/(double)p;
    printf("    %i \t%02i%% ",i,(int)(100.0*x));
    RULING(x);
    printf(" [× %li] \n",F[i]);
  }
  printf("- average route length: %.02lf ± %.02lf (%li/%li)\n",s/(double)p,ECARTYPE(s,s2,p),s,p);
  printf("- maximum route length: %i\n",lmax);
  free(F);

  BARRE;
  printf("- distance distribution:\n");
  ALLOCZ(F,dmax+1,0L); /* F[k]=nombre de distances de longueur k */
  for(k=0;k<=dmax;k++)
    if(stat[k])
      for(j=0;j<=M[k];j++)
	F[k] += stat[k][j];

  /* s=somme totale des distances, s2=somme du carré des distances */
  for(i=0,s=s2=0L;i<=dmax;i++) s+=((long)i)*F[i],s2+=((long)i)*((long)i)*F[i];

  for(i=0;i<=dmax;i++){
    if(F[i]==0) continue;
    x=F[i]/(double)p;
    printf("    %i \t%02i%% ",i,(int)(100.0*x));
    RULING(x);
    printf(" [× %li] \n",F[i]);
  }
  printf("- average distance: %.02lf ± %.02lf (%li/%li)\n",s/(double)p,ECARTYPE(s,s2,p),s,p);
  printf("- maximum distance: %i\n",dmax);
  free(F);

  BARRE;
  printf("- stretch distribution:\n");

  /* t=borne sup sur le nombre de stretch différents */
  for(k=t=0;k<=dmax;k++)
    if(stat[k]) t += M[k]+1;

  triplet *ptr;
  triplet e={0,0,0L}; /* triplet nul */
  NALLOCZ(triplet,P,t,e); /* P=tableau de triplets (j,k,stat[k][j]) */

  /* on construit P */
  /* t=nombre de triplets ajoutés à P */

  /* TODO: t pourrait être très grand ≃ n^2 et donc il faudrait qu'il
     soit de type long. Ensuite, si t est très grand, il ne faut pas
     tous les afficher, mais plutôt une distribution (10 ranges par
     exemples). La même remarque s'applique à la distribution des
     distances (cependant, il ne peut avoir que n-1 distances
     différentes). */

  for(k=t=0;k<=dmax;k++)
    if(stat[k])
      for(j=0;j<=M[k];j++)
	if(stat[k][j]){
	  if(k){ u=pgcd(k,j); e.x=j/u; e.y=k/u; }
	  else{ e.x=0; e.y=1; } /* pour dist=0, le stretch est 1+0 */
	  e.z=stat[k][j];
	  ptr=bsearch(&e,P,t,sizeof(triplet),fcmp_stretch);
	  if(ptr) P[(int)(ptr-P)].z += e.z;  /* e est déjà dans P, ajoute e.z */
	  else{
	    P[t++]=e; /* on ajoute une nouvelle entrée e à P */
	    QSORT(P,t,fcmp_stretch); /* trie le nouveau tableau */
	  }
	}

  /* affiche les stretchs contenus dans P */
  for(i=0;i<t;i++){
    x=(double)(P[i].z)/(double)p;
    printf("  %0.3lf (%i/%i)",1.0+(double)P[i].x/(double)P[i].y,P[i].x+P[i].y,P[i].y);
    printf("\t%02i%% ",(int)(100.0*x));
    RULING(x);
    printf(" [× %li] \n",P[i].z);
  }
  free(P);

  /* stretch moyen, écart type et stretch max */
  /* on pourrait l'avoir directement avec P[], mais l'avantage avec
     stat[][] est qu'on peut aussi avoir la distance max qui atteint
     le stretch max. On perd l'info dans P[] à cause du pgcd). */

  double r=0,r2=0; /* r=somme total des stretch, r2=somme du carré des stretch */
  double smax=1; /* smax=stretch max */

  u=v=1; /* v=distance et u=déviation pour le stretch max */
  for(k=0;k<=dmax;k++)
    if(stat[k])
      for(j=0;j<=M[k];j++){
	x=1; if(k) x += (double)j/(double)k; /* x=stretch, x=1.0 si k=0 */
	r += (double)stat[k][j]*x; r2 += (double)stat[k][j]*x*x;
	if(x>=smax) smax=x,u=j,v=k;
      }
  printf("- average stretch: %.03lf ± %.03lf (%.02lf/%li)\n",(double)r/(double)p,ECARTYPE(r,r2,p),r,p);
  printf("- maximum stretch: %.02lf (%i/%i)\n",1.0+(double)u/(double)v,u+v,v);

  if(err){
    BARRE;
    printf("- Warning! there were %s failed routings\n",millier(err));
  }

 fin_routing_test:
  BARRE;
  free(M);
  FREE2(stat,dmax+1);
  return 0;
}


/***********************************

       FIN ROUTINES POUR LES
          ROUTING SCHEMES

***********************************/


void RemoveVertex(graph *G,const int u){
/*
  Supprime toutes les occurences du sommet u dans le graphe G non
  null, sans renuméroter les sommets et sans réallouées les
  listes. Après cet appel, les listes de G ne sont plus forcément
  triées, mais la liste de G->L[u] existe toujours. Simplement G->d[u]
  est à zéro.

  Effet de bord: modifie G->sort, G->d[.], G->L[.].
  Attention ! G->n n'est pas modifié.
*/
  const int du=G->d[u];
  int i,j,v,dv;
  G->d[u]=0;
  for(i=0;i<du;i++){
    v=G->L[u][i];
    dv=G->d[v];
    for(j=0;j<dv;j++)
      if(G->L[v][j]==u)	G->L[v][j]=G->L[v][--dv];
    G->d[v]=dv;
  }
  G->sort=0;
  return;
}


int AdjGraph(const graph *G,const int u,const int v){
/*
  Renvoie 1 ssi dans le graphe G le sommet v apparaît dans la liste
  d'adjacence de u. Si les listes du graphe G sont est triées (G->sort
  est vrai), une dichotomie est effectuée (donc en log(deg(u))) sinon
  c'est un parcours séquentiel (en deg(u)).
*/

  /* recherche dichotomique si le graphe est trié */
  if(G->sort)
    return (bsearch(&v,G->L[u],G->d[u],sizeof(int),fcmp_int)!=NULL);

  /* sinon, recherche séquentielle */
  int i;
  for(i=G->d[u];i>0;)
    if(G->L[u][--i]==v) return 1; /* on a trouvé v dans la liste de u */
  return 0;
}


int Treewidth(graph *H,const int code)
/*
  Calcul approché ou excate de la treewidth de H. Si code=0, on
  calcule une borne supérieure (heuristique degmin). Si code=1 on
  calcule la valeur exacte. H n'est pas modifié. Dans H->int1 on
  renvoie le nb de tests effectués.
*/
{
  if(H==NULL) return 0;
  int n,j,t,tw,l,v,n1;
  int d,r,i,u,k;
  int d0,r0,i0,u0;
  graph *G;

  H->int1=0;
  n=H->n;
  n1=n-1;
  tw=(code==1)? Treewidth(H,0) : imin(n1,NbEdges(H));

  /* tw=upper bound sur tw. On a intérêt d'avoir la valeur la plus
     faible possible, car on va alors éliminer les mauvais ordres
     d'éliminations plus rapidement.

    tw max en fonction du nb m d'arêtes:
    m=0 => tw=0
    m=1 => tw=1
    m=2 => tw=1
    m=3 => tw<=2
    m=4 => tw<=2
    m=5 => tw<=2?
    m=6 => tw<=3
    ...
   */

  NALLOCZ(int,P,n,_i); /* permutation initiales des sommets */
  NALLOCZ(int,D,n,n1);

  do{
    G=ExtractSubgraph(H,NULL,0,0); /* on copie H */
    GraphRealloc(G,D); /* plonge G dans le graphe complet */
    k=0; /* tw par défaut si G sans d'arêtes */

    for(j=0;j<n1;j++){ /* n-1 fois */

      H->int1++;
      u0=P[i0=j];
      d0=r0=G->d[u0];
      if(code==0){
	for(i=j;i<n1;i++){
	  u=P[i]; d=G->d[u];
	  if(d<=d0){
	    /* calcule r=nb d'arêtes dans N(u) */
	    for(r=l=0;l<d;l++)
	      for(t=l+1;t<d;t++)
		r += AdjGraph(G,G->L[u][l],G->L[u][t]);
	    if((d<d0)||(r<r0)){ d0=d;r0=r;i0=i;u0=u; }
	  }
	} /* fin for(i=...) */
	P[i0]=P[j]; /* décale P[i], que si code=0 */
      }
      k=imax(k,d0); /* met à jour tw */
      if(k>=tw) goto nextP; /* on ne fera pas moins */
      RemoveVertex(G,u0); /* supprime u */
      /* remplace N(u) par une clique */
      for(i=0;i<d0;i++){
	u=G->L[u0][i];
	for(t=i+1;t<d0;t++){
	  v=G->L[u0][t];
	  if(!AdjGraph(G,u,v)) ADD_EDGE(G,u,v);
	}
      }
      
    } /* fin for(j=...) */

    tw=imin(tw,k);
  nextP:
    free_graph(G);
    if((code==0)||(tw<2)) break; /* si tw(G)=0 ou 1, alors on ne
				    trouvera pas plus petit */
  }while(NextPermutation(P,n,NULL));

  free(D);
  free(P);
  return tw;
}


int *Prune(const graph *G,int *z)
/*
  Algorithme permettant de supprimer les sommets du graphe G par
  dégrés minimum. Un ordre d'élémination des sommets est renvoyé sous
  la forme d'un tableau de taille n. Plus précisément, on renvoie un
  tableau de sommets T de sorte que u=T[i] est un sommet de degré
  minimum du graphe G\(T[0]...T[i-1]). La complexité est linéaire,
  i.e. O(n+m). Si *z<>NULL, on retourne dans z le degré maximum
  rencontré donc la dégénérescence du graphe.
*/
{
  int i,j,u,v,d,k,p,r,t,c,n1;
  const int n=G->n;

  /*
    1. On construit les tableaux suivants:

    T[0..n[ = tableau de sommets triés par degré croissant (le résultat)
    R[0..n[ = tableau de positions dans T des sommets
    D[0..n[ = tableau des degrés des sommets
    P[0..n[ = tableau des positions dans T du premier sommet de degré donné

    2. On traite les sommets dans l'ordre T[0],T[1],... Supposons que
    u=T[i]. Alors, pour chaque voisin v dans G, encore existant (donc
    situé après i dans le tableau T), on met à jour la structure
    T,R,D,P.
   */

  /* initialise T,R,D,P */
  if(n<1) return NULL;
  NALLOCZ(int,D,n,G->d[_i]);
  int *R=SortInt(D,NULL,n,0,NULL,SORT_INDEXe);
  NALLOC(int,T,n); for(u=0;u<n;u++) T[R[u]]=u;
  NALLOCZ(int,P,n,-1);
  for(i=n1=n-1;i>=0;i--) P[D[T[i]]]=i; /* i=n-1,...,0 */

  for(c=i=0;i<n1;i++){ /* pour chaque sommet, pas besoin de faire le dernier */
    u=T[i]; /* u = sommet de degré (D[u]) minimum */
    if(D[u]>c) c=D[u]; /* mémorise le plus grand degré rencontré */
    d=G->d[u]; /* d = deg(u) */
    for(j=0;j<d;j++){ /* pour chaque voisin */
      v=G->L[u][j]; /* v=voisin de u */
      r=R[v]; /* v=T[r]; */
      if(r>i){ /* mettre à jour que si v existe encore */
	k=D[v]--; /* le degré de v diminue, k=D[v] juste avant */
	p=P[k--]++; /* on va échanger v et T[p]. Les sommets de degré k commencent juste après */
	if(P[k]<0) P[k]=p; /* y'a un sommet de degré k en plus */
	if(p>i){ /* si p est avant i, ne rien faire. Cel arrive que si k=0 */
	  t=T[p]; /* t = 1er sommet de de degré k */
	  T[p]=v; /* dans T, on avance v en position p (à la place de t donc) */
	  T[r]=t; /* on met t à la place de v */
	  R[v]=p; /* on met à jour la nouvelle position de v dans T */
	  R[t]=r; /* on met à jour la nouvelle position de t dans T */
	}
      }
    }
  }

  if(z!=NULL) *z=c; /* retourne la dégénérescence */
  free(R);
  free(D);
  free(P);
  return T;
}


int *GreedyColor(graph *G,int *R)
/*
  Algorithme permettant de colorier de manière gloutonne un graphe G à
  n sommets. La complexité en temps est linéaire, i.e. O(n+m). Il faut
  G non NULL. On renvoie un tableau C[0..n[ où C[u] est la couleur du
  sommet u, un entier entre 0 et n-1. Le tableau R[0..n[ donne un
  ordre (inverse) de visite des sommets: On colorie le sommet u avec
  la plus petite couleur possible dans le graphe induit par les
  sommets situé après u dans R (on commence donc avec R[n-1]). Si
  R=NULL, alors l'ordre R[i]=i est utilisé. On renvoie dans G->int1 la
  couleur maximum utilisée (donc le nb de couleurs-1). Cette valeur
  vaut -1 si n<1.  */
{
  int i,j,u,d,l,c;
  const int n=G->n;

  /*
    On utilise 3 tableaux:

    C[0..n[ = tableau final des couleurs, au départ =-1
    L[0..n-1[ = liste de couleurs utilisées par le sommet courant
    M[0..n-1[, M[i]=1 ssi la couleur i est utilisée par un voisin du
    sommet courant, au départ =0. On met une sentiennelle à M si bien
    que toujours on M[n-1]=0.

  */

  G->int1=-1;
  if(n<1) return NULL;
  NALLOCZ(int,C,n,-1);
  NALLOCZ(int,M,n,0);
  i=n-1;
  NALLOC(int,L,i);

  for(;i>=0;i--){ /* en partant de la fin */
    u=(R==NULL)? i : R[i];
    d=G->d[u]; /* d = deg(u) */

    /* on liste dans L et M les couleurs rencontrées */
    for(j=l=0;j<d;j++){ /* pour chaque voisin v de u */
      c=C[G->L[u][j]]; /* c=couleur du voisin v */
      if(c>=0){ /* voisin déjà colorié ? */
	L[l++]=c; /* on ajoute la couleur c à la liste */
	M[c]=1; /* la couleur c n'est pas à prendre */
      }
    }
    
    /* on cherche la 1ère couleur à 1 (=non-rencontrée) */
    j=0; while(M[j]) j++; /* s'arrête toujours à cause de la sentiennelle */
    C[u]=j; /* j est la couleur recherchée pour u */
    G->int1=imax(G->int1,j); /* couleur maximum rencontrée */

    /* il faut ré-initialiser rapidement M */
    for(j=0;j<l;j++) M[L[j]]=0; /* on efface les 1 qu'on à mis dans M */
  }
  
  free(L);
  free(M);
  return C;
}


void HalfGraph(graph *G,int code){
/*
  Transforme G en un graphe asymétrique (avec G->sym=1) en graphe
  symétrique en ne gardant que les arcs u->v tels que v<u si code=0 ou
  bien tels que u<v si code=1. En particulier les boucles sont
  supprimés. Les listes d'adjacence sont aussi triées par ordre
  croissant et les tables raccourcies (realloc). Cela permet de faire
  des parcours de graphe deux fois plus rapidement, par exemple, pour
  vérifier qu'une coloration est propre. Le champs ->sym est mis à
  jour. L'exécution est plus rapide si code=0.
*/
  if(G==NULL) return;
  if(!G->sym) return;

  const int n=G->n;
  int i,u,d;
  NALLOC(int,D,n); // tableau des nouveaux degrés

  SortGraph(G,0); // tri par ordre croissant les listes
  
  for(u=0;u<n;u++){
    d=G->d[u];
    for(i=0;i<d;i++)
      if(G->L[u][i]>=u) break;
    if(code){ // il faut v>u
      while((i<d)&&(G->L[u][i]==u)) i++; // les boucles
      memmove(G->L[u],G->L[u]+i,(G->d[u]-i)*sizeof(int));
      i=G->d[u]-i;
    }
    D[u]=i; // coupe la liste à i
  }
  
  GraphRealloc(G,D); // modifie toutes les listes d'adjacence
  free(D);
  G->sym=0; // G est désormais asymétrique
  return;
}


void kColorSat(graph *G,const int k){
/*
  Affiche les contraintes SAT pour la k-coloration de G au format
  Dimacs. Attention ! G est modifié (manque la moitié de ses arcs).

  Le nombre de variables est n*k.
  Le nombre de clause est n+m*k.
*/

  if(G==NULL) return;
  const int n=G->n;
  const int m=NbEdges(G);
  int c,u,i,d;
  
  printf("p cnf %i %i\n",n*k,n+m*k);

  /*
    Variables: chaque sommet i a k variables numérotés x(i,c) avec
    i=0..n-1 et c=0..k-1. Il faut x(i,c)>0. On pose x(i,c)=1+k*i+c.

    Contraintes sommets: il faut x(i,0) v ... v x(i,k-1), et la
    conjonction pour tous les sommets.

    Containtes arêtes: pour chaque arête {i,j} et couleur c il ne faut
    pas que x(i,c)=x(j,c). Autrement dit il faut -x(i,c) v -x(j,c). Il
    faut la conjonction de toutes les couleurs c et toutes les arêtes
    {i,j}.
  */

  /* liste chaque sommet */
  for(u=0;u<n;u++){
    for(c=0;c<k;c++) printf("%i ",1+u*k+c);
    printf("0\n");
  }

  HalfGraph(G,0); /* enlève la moitié des arcs */

  /* liste chaque arc */
  for(u=0;u<n;u++){
    d=G->d[u];
    for(i=0;i<d;i++)
      for(c=0;c<k;c++)
	printf("-%i -%i 0\n",1+u*k+c,1+G->L[u][i]*k+c);
  }
  
  return;
}


void kIndepSat(graph *G,const int k){
/*
  Affiche les contraintes SAT pour ensemble indépendant de taille k
  pour le graphe G au format Dimacs. Attention ! G est modifié (manque
  la moitié de ses arcs).

  Variables:

    Pour chaque i=0..n-1:
    X(i)=1 ssi le sommet i est dans l'ensemble indépendant.
    Si i-j est une arête alors -X(i) v -X(j).
    NB: pour Vertex Cover, il faudrait suffit d'ajouter X(i) v X(j).

    Pour chaque t=0..n et b=0..k:
    Y(t,b)=1 ssi la somme des t variables X_0+...+X(t-1) est au
    moins b. On a:

       Y(n,k) = 1
       Y(t,0)=1 pour tout t>=0
       Y(t,b)=0 pour tout 0<=t<b et b>0
       Y(t,b) => Y(t-1,b) v (Y(t-1,b-1) ^ X(t-1))

       <=> -Y(t,b) v Y(t-1,b) v Y(t-1,b-1) ET
           -Y(t,b) v Y(t-1,b) v X(t-1)

    #variables X(i): n
    #variables Y(t,b): (n+1)*(k+1)
    #variables totales: n+(n+1)*(k+1)

*/

  if(G==NULL) return;

  const int n=G->n;
  const int m=NbEdges(G);
  int i,j,t,d,b;  

  // il faut entrer le nombre exactes de variables et de clauses
  printf("p cnf %i %i\n",n+(n+1)*(k+1),m+n+2+k*(k+1)/2+k*(2*n+1-k));

  /* numéros des variables: attention de ne surtout pas utiliser la
     variable numéro 0 qui signifie "fin de ligne" */

#define X(i)   ((i)+1)             // numéro de la variable X(i)
#define Y(t,b) (n+1+(t)*(k+1)+(b)) // numéro de la variable Y(t,b)

  HalfGraph(G,0); /* rend le graphe simple et asymétrique */

  // pour chaque arêtes i-j: -X(i) v -X(j)
  // #clauses: m
  for(i=0;i<n;i++){
    d=G->d[i];
    for(j=0;j<d;j++)
      printf("-%i -%i 0\n",X(i),X(G->L[i][j]));
  }

  // Y(n,k)=1
  // #clause: 1
  printf("%i 0\n",Y(n,k));

  // cas b=0 et t>=0: Y(t,0)=1
  // #clauses: n+1
  for(t=0;t<=n;t++) printf("%i 0\n",Y(t,0));

  // cas b>=1 et 0<=t<b: Y(t,b)=0
  // #clauses: 1+2+3+...+k = k*(k+1)/2
  for(b=1;b<=k;b++)
    for(t=0;t<b;t++)
      printf("-%i 0\n",Y(t,b));

  // cas b>=1 et t>=b: récurrence
  // #clauses: 2*∑_{b=1}^k (n-b+1) = 2*(n+1)*k - k(k+1) = k*(2*n+1-k)
  for(b=1;b<=k;b++)
    for(t=b;t<=n;t++){
      printf("-%i %i %i 0\n",Y(t,b),Y(t-1,b),Y(t-1,b-1));
      printf("-%i %i %i 0\n",Y(t,b),Y(t-1,b),X(t-1));
    }

  return;
}


int *kColor(graph *G,const int k){
/*
  Algorithme permettant de colorier en au plus k couleurs un graphe G,
  si c'est possible. La complexité est (n+m)*k^{n-1} dans le pire des
  cas.  On renvoie un tableau C[0..n[ où C[u] est la couleur du sommet
  u, un entier entre 0 et k-1. On renvoie NULL s'il n'est pas possible
  de colorier G en k couleurs, si G=NULL ou k<1. On renvoie dans
  G->int1 la couleur maximum utilisée (qui peut être < k-1). On
  utilise toutes les couleurs de [0,G->int1]. On symétrise le graphe,
  qui est donc modifié, afin d'enlever la moitié des arêtes à
  vérifier.

  La stratégie est la suivante. On part d'une coloration initiale des
  sommets C=[0,...,0,k-1] où la couleur du dernier sommet u=n-1 est
  fixée par C[n-1]=k-1. Puis on vérifie si C est propre ou non. Si
  c'est non, on incrémente C comme un compteur, et on recommence avec
  la coloration suivante. Ainsi toutes les colorations possibles de G
  sont passées en revue. On teste toujours en priorité la dernière
  arête qui a posé problème avant de vérifier tout le graphe.

  Optimisations possibles à faire:
  
  1. Réduction de données. On peut supprimer récursivement les sommets
     de degré < k. On pourra toujours les ajouter à la fin si la
     coloration a réussie.

  2. On peut décomposer le graphe en composantes connexes, il y a
     ainsi moins de colorations possibles à tester.

  3. On peut renuméroter les sommets selon un parcours BFS. De cette
     façon la vérification pourrait être accélérée lorsqu'on change
     une seule couleur.

*/

  int c,i,d,u,v,b,*T;
  if((G==NULL)||(k<1)) return NULL;
  HalfGraph(G,0); /* enlève la moitié des arêtes */
  const int n=G->n;

  NALLOCZ(int,C,n,0); /* C[u]=couleur du sommet u */
  C[n-1]=k-1; /* on peut supposer que le sommet n-1 a une couleur fixée */
  if(n<2) goto fin_kcolor; /* s'il y a un seul sommet */
  b=1; /* b=vrai ssi la dernière arête coloriée est propre */

  do{
    /* vérifie si C est une coloration propre */
    if(b){ /* on le fait que si b est vrai */
      for(u=0;u<n;u++){ /* pour chaque sommet */
	c=C[u]; d=G->d[u]; /* degré et couleur de u */
	for(i=0;i<d;i++) /* pour chaque voisin de i */
	  if(c==C[G->L[u][i]]) break; /* coloration pas propre */
	if(i<d) break; /* coloration pas propre */
      }
      if(u==n) goto fin_kcolor; /* la coloration est propre */
    }
    /* ici l'arête (u,i) n'est pas correctement coloriée */
    
    /* on change C en incrémentant le compteur */
    v=0;
  loop_kcolor:
    C[v]++;
    if(C[v]==k){
      C[v++]=0;
      if(v==n) break;
      goto loop_kcolor;
    }
    
    /* est-ce que l'arête (u,i) est toujours mal coloriée ?
       si oui, pas la peine de vérifier tout le graphe */
    b=(C[u]!=C[G->L[u][i]]); /* b=vrai ssi (u,i) est propre */

  }while(v<n);

  /* aucune coloration trouvée */
  free(C);
  return NULL;

 fin_kcolor:
  /* on a trouvé une coloration propre */
  /* on réduit l'espace des couleurs utilisées.  NB: ici on a encore
     C[n-1]=k-1 */

  ALLOCZ(T,k,-1);
  /* si T[c]=i, alors la couleur c se renumérote en i. Si c n'a jamais
     été rencontrée, alors i=-1 */

  for(i=u=0;u<n;u++){
    if(T[C[u]]<0) T[C[u]]=i++; /* la couleur C[u] n'a pas jamais été vue */ 
    C[u]=T[C[u]]; /* ici la couleur C[u] a déjà été vue */
  }

  free(T); /* plus besoin de T */
  G->int1=i-1; /* couleur max utilisée */
  return C;
}


int *power_law_seq(const int n,const double t,int *T){
/*
  Ecrit dans le tableau T une distribution de degré pour un graphe à
  n>0 sommets selon une lois en puissance d'exposant t>0. On renvoie
  NULL si les paramètres n et t ne sont pas corrects. La distribution
  est codée par une suite (2k,n_1,d_1,...,n_k,d_k) de paires (n_i,d_i)
  qui signifie n_i sommets de degré d_i. Si T=NULL, alors T est alloué
  et renvoyé, sinon il doit être assez grand pour recevoir la
  distribution. La taille de T, sans jamais dépasser 2n+1, vaut en
  théorie |T| = O(n^{1/t}). En pratique on a:

    t    n     |T|        t    n     |T|        t    n     |T|
   ----------------      ----------------      ----------------
    2.0  10^3  50         2.5  10^3  28         3.0  10^3  18
    2.0  10^4  156        2.5  10^4  70         3.0  10^4  40
    2.0  10^5  494        2.5  10^5  176        3.0  10^5  86
    2.0  10^6  1560       2.5  10^6  446        3.0  10^6  188
    2.0  10^7  4932       2.5  10^7  1122       3.0  10^7  404
   ----------------      ----------------      ----------------

  La distribution est la suivante:
  - d_1=1, n_1 = ⌊exp(a)⌋ + n-s(a)
  - d_i=i, n_i = ⌊exp(a)/i^t⌋ pour i dans [2,p(a)]
  - a est un réel minimisant |n-s(a)| avec
    s(a) := ∑_{i=1}^p(a) ⌊exp(a)/i^t⌋
    p(a) := ⌊exp(a/t)⌋

  Attention ! Dans l'article original de [Lu01], il y a une erreur
  dans le signe du terme correctif r=n-s(a) pour les sommets de degré
  1. Il faut faire +r et non -r comme c'est écrit.

  Le nombre de paires (n_i,d_i) dans T est exactement p(a). Si n>0,
  alors T contient au moins 1 paire (et au plus n), car les d_i sont
  différents. On a s(0)=1 car p(0)=1. Aussi si a>=ln(n)+1, alors s(a)
  ≥ floor{exp(a)} ≥ n.

  En choisissant a0=0 et a1=ln(n)+1, on a alors s(a0) <= n <=
  s(a1). La fonction s(a) étant croissante, pour minimiser n-s(a) on
  réalise une recherche binaire pour a dans [a0,a1]. Le nombre
  d'itérations est au plus 64 (car sizeof(double)*8=64).
*/
  if((t<=0)||(n<=0)) return NULL;

  /* calcule la valeur de 'a' optimale */
  
  double a0=0,a1=log(n)+1.0; // intervalle pour a
  double a,b,e; // a=milieu de [a0,a1], b=meilleur a
  int p,i,j,s,cont=1; // cont=1 ssi on continue le calcul
  int r=0; // r=valeur minimum (=0 au départ)
  
  do{
    a=(a0<a1)? (a0+a1)/2 : a0; // si a0>=a1, on calcule s puis on s'arrête
    if((a==a0)||(a==a1)) cont=0; // intervalle trop faible, on calcule s puis on s'arrête
    e=exp(a); p=(int)exp(a/t); // p=nombre de paires, NB: p>=1
    for(s=0,i=1;i<=p;i++) s += (int)(e/pow(i,t)); // calcule s=s(a)
    if(s==n) cont=0; // valeur optimale, on va avoir b=a
    if((r==0)||(abs(n-s)<abs(r))) b=a,r=n-s; // NB: si s=n, alors b=a
    if(s<n) a0=a; // on est avant n
    if(s>n) a1=a; // on est après n
  }while(cont);

  /* ici on a calculé b, le meilleur a */

  e=exp(b);
  p=(int)exp(b/t); // nombre de paires
  if(T==NULL) ALLOC(T,2*p+1);
  T[0]=2*p; // taille du tableau

  /* écrit la distribution dans T */
  
  for(j=i=1;i<=p;i++){
    T[j++]=(int)(e/pow(i,t)); // NB: si i=1, T[1]=floor(exp(b))
    T[j++]=i;
  }
  DEBUG(PRINT(j););

  T[1]+=r; // correction pour les sommets de degré 1
  
  return T;
}


/***********************************

         ROUTINES POUR LES
       FONCTIONS D'ADJACENCE

***********************************/

// Pour avoir la liste de toutes les fonctions d'adjacence:
// grep '^int.*(query\* const Q){' gengraph.c

/* code souvent utilisé dans les fonctions d'adjacence */

#define TEST_n       do{ if(Q->n<=0){ Q->n=0; return 0; }}while(0)
#define SET_n(X)     do{ Q->n=(X); TEST_n; }while(0)
#define RET_n(X)     do{ Q->n=(X); if(Q->n<=0) Q->n=0; return 0; }while(0)
#define RET_a(X)     do{ Q->a=(X); return 0; }while(0)
#define RET_error(X) do{ Q->error=(X); return 1; }while(0)


int free_rep(query* const Q){
/*
  Routine permettant de libérer les tableaux Q->rep. Utilise Q->n.

  La fonction devrait toujours être exécutée dans une fonction
  d'adjacence avec Q->code=QUERY_END. Renvoie 0 si tout c'est bien
  passé, 1 sinon.
*/
  if((Q->n<=0)||(Q->code!=QUERY_END)) return 1; // fin anormale
  FREE2(Q->rep,Q->n);
  return 0; // fin normale
}


int free_pos(query* const Q){
/*
  Routine permettant de libérer les tableaux XPOS, YPOS ...

  La fonction devrait toujours être exécutée dans une fonction
  d'adjacence avec Q->code=QUERY_END. Renvoie 0 si tout c'est bien
  passé, 1 sinon.
*/
  if(Q->code!=QUERY_END) return 1; // fin anormale
  free(XPOS),free(YPOS);
  free(XSEED),free(YSEED);
  XPOS=YPOS=XSEED=YSEED=NULL;
  return 0;
}


/***********************************

         GRAPHES DE BASE

***********************************/


/*
  Une fonction d'adjacence de graphe, disons adj(Q), effectue un
  calcul en fonction de la requête Q->code et des divers paramètres
  contenus dans Q. Le code typique ressemble donc à un switch(Q->code)
  comme ci-dessous.

  int adj(query* const Q){
    switch(Q->code){
    case QUERY_INIT: Q->n=...; return 0; // ou RET_n(...);
    case QUERY_ADJ:  Q->a=...; return 0; // ou RET_a(...);
    }
    return 1;
  }

  Toutes les fonctionalités ne sont pas forcément implémentées.
  Cependant, a minima, QUERY_INIT (pour déterminer Q->n) et QUERY_ADJ
  (pour déterminer Q->a) doivent être implémentées. La fonction doit
  renvoyer 0 si tout c'est bien passée, et 1 sinon. Dans ce cas, il
  est possible de renseigner Q->error.

  Rem1: Plutôt que d'utiliser "return 0" à la fin de chaque "case", on
  pourrait mettre un "break" puis un "return 0" en dehors du "switch"
  avec aussi un "default: return 1". Cependant un "return" direct est
  plus rapide, ce qui peut être important pour QUERY_ADJ et QUERY_NAME
  par exemple.

  Rem2: Attention aux subtilités de l'instruction "switch". Les
  différents "case" correspondent à des étiquettes comme un "goto". En
  particulier, les déclarations avec initialisation ne sont pas
  valides juste après un "case". Une solution est de déplacer (plus
  loin) la déclaration ou d'inserer une instruction (éventuellement
  vide ";") entre le "case" et la déclaration, ou encore de créer un
  block. Voici un résumé:

    (incorrect)        (correct)          (correct)           (correct)
  
  case QUERY_ADJ:    case QUERY_ADJ:    case QUERY_ADJ:;    case QUERY_ADJ:
    int k=42;          int k;             int k=42;           { int k=42;
    ...;               k=42;              ...;                  ...;
                       ...;                                     }
   
  Rem3: Il faut eviter de traiter l'erreur directement dans QUERY_INIT
  avec Erreur(...) car les fonctions peuvent s'appeler en cascade, et
  on peut ne plus comprendre d'où vient l'erreur. C'est à l'appelant
  de gérer l'erreur. Dans la plupart des cas on peut se contenter de
  renvoyer le graphe vide avec RET_n(0) si les paramètres sont
  incorrects. Si l'on souhaite renvoyer une erreur particulière on
  peut (devrait) faire:

    case QUERY_INIT:
      ...;
      if(...) RET_error(42);
      ...;

  Les différentes fonctionnalités sont:

    QUERY_INIT: détermine Q->n, initialise la fonction
    QUERY_END: termine la fonction
    QUERY_ADJ: détermine Q->a en fonction de Q->i, Q->j et des paramètres
    QUERY_NAME: détermine Q->name, le nom du sommet Q->i

    [et à finaliser ...]

    QUERY_LIST: détermine Q->L, la liste des voisins de Q->i
    QUERY_DOT: détermine le dessin d'une arête au format dot
    QUERY_DEG: détermine le degré du sommet Q->i
    QUERY_LNAME: pour la liste des noms des sommets

  Certaines fonctions d'adjacence suppose Q->i<Q->j. Cela signifie que
  le résultat peut ne pas être correcte, voir produire une erreur
  grave, si l'option -DIRECT est présente.

  C'est une mauvaise idée que d'utiliser des variables statiques dans
  les fonctions d'ajacence car elles peuvent s'appeler entre
  elles. Par exemple arboricity(Q) est utilisée par plusieurs autres
  fonctions. Si les pré-calculs sont les bienvenus dans QUERY_INIT
  pour optimiser QUERY_ADJ, il faut impérativement stocker leurs
  résultats dans les champs de la variable Q.

  On devrait systématiquement, lors de l'initialisation avec
  QUERY_INIT, avoir le test "if(Q->n<=0){ Q->n=0; return 0; }" ou
  encore la macro "TEST_n;" (le "=0" dans le "if(Q->n<=0)" est
  important) permettant de générer le graphe vide sans rien faire
  d'autre surtout pour les graphes allouant de la mémoire (Q->rep,
  Q->xpos, ...) qui nécessitent que Q->n>0.

  Les graphes utilisant le tableau Q->rep[u] (représentation
  implicite) pour chaque sommet u devraient se terminer par la
  libération des tableaux comme ceci:

    case QUERY_END:
      return free_rep(Q);

  De même, les graphes géométriques qui utilisent les tableaux XPOS et
  YPOS devraient se terminer par leur libération comme ceci:
  
    case QUERY_END:
      return free_pos(Q);
*/


int load(query* const Q){
/*
  Graphe défini par un fichier (ou l'entrée standard). A
  l'initialisation, il est chargé en mémoire dans la variable Q->G de
  type "graph". Suivant la valeur de Q->G->sort, le test d'adjacence
  est linéaire ou logarithmique en min{deg(i),deg(j)}.
*/
  switch(Q->code){
    
  case QUERY_END:
    free_graph(Q->G);
    Q->G=NULL;
    return 0;

  case QUERY_INIT:
    Q->G=File2Graph(Q->sparam,2); /* remplit Q->G à partir du fichier */
    if(Q->G->f>0){ /* si c'est une famille, on sélectionne le premier graphe */
      graph *G=ExtractSubgraph(Q->G->G[0],NULL,0,0); /* copie le premier graphe */
      free_graph(Q->G); /* libère complètement la famille Q->G */
      Q->G=G; /* Q->G=premier graphe */
    }
    if(!Q->G->sym) DIRECTED=1; /* si le graphe est asymétrique, affichage DIRECTED */
    RET_n(Q->G->n);

  case QUERY_ADJ:
    /* pour avoir du min{deg(i),deg(j)} en non-orienté */
    if((!DIRECTED)&&(Q->G->d[Q->i]>Q->G->d[Q->j])) Q->a=AdjGraph(Q->G,Q->j,Q->i);
    else Q->a=AdjGraph(Q->G,Q->i,Q->j);
    return 0;

  }

  return 1;
}


int prime(query* const Q){
  switch(Q->code){
  case QUERY_INIT: RET_n(Q->param[0]);
  case QUERY_ADJ: RET_a((Q->i>1)? ((Q->j%Q->i)==0) : 0);
  }
  return 1;
}


int paley(query* const Q){
/*
  Le résidu est r=|i-j|. Pour savoir si r est un carré, on teste s'il
  existe un entier k<=(n-1)/2 tel que (k*k)%n=r. Le nombre de carrés
  si dans Z/nZ, si n est premier vaut (n+1)/2.
*/
  const int n=Q->param[0];
  switch(Q->code){
  case QUERY_INIT: RET_n(n);
  case QUERY_ADJ:;
    int k;
    const int q=n/2;
    const int r=abs(Q->i-Q->j);
    for(k=1;k<=q;k++) if(r==((k*k)%n)) RET_a(1);
    RET_a(0);
  }
  return 1;
}


int comb(query* const Q){
/*
  Utilise i<j.
  On pourrait utiliser "grid 1 n -apex -1".
*/
  const int n=Q->param[0];
  switch(Q->code){
  case QUERY_INIT: RET_n(2*n);
  case QUERY_ADJ: RET_a( ((Q->j<n)&&(Q->j==Q->i+1))||(Q->j==Q->i+n) );
  }
  return 1;
}


int sunlet(query* const Q){
/*
  Utilise i<j.
  On pourrait utiliser: "grid 1 -n -apex -1" 
*/
  if((Q->code==QUERY_ADJ)&&(Q->i==0)&&(Q->j==Q->param[0]-1)) RET_a(1);
  return comb(Q);
}


/* types pour alkane() */
enum{
  ALK_NOR,
  ALK_CYC,
  ALK_ISO,
  ALK_NEO,
  ALK_SEC,
  ALK_TER,
};


int alkane(query* const Q){
/*
  Utilise i<j.

  Les sommets de C ont des numéros < n.
  Les sommets de H ont des numéros >= n.
  Pour C, dans l'ordre on trouve:
  - les sommets de degré 1
  - les sommets de degré 3
  - les sommets de degré 4
  - les sommets de degré 2
  Cette numérotation permet de gérer de manière uniforme et simple
  les adjacences C-H. Les adjacences C-C se font aux cas par cas.

  Q->wrap[0]: borne inférieure pour n
  Q->wrap[1]: nombre de sommets C de degré 1
  Q->wrap[2]: nombre de sommets C de degré 3
  Q->wrap[3]: nombre de sommets C de degré 4

  Vidéos sur les noms des alkanes:
  https://fr.khanacademy.org/science/organic-chemistry/bond-line-structures-alkanes-cycloalkanes#naming-alkanes
*/
  const int type=Q->param[0];
  const int n=Q->param[1]; // n=nombre d'atomes de C
  static int P_NOR[]={1,2,0,0};
  static int P_CYC[]={3,0,0,0};
  static int P_ISO[]={4,3,1,0};
  static int P_NEO[]={5,4,0,1};
  static int P_SEC[]={6,3,1,0};
  static int P_TER[]={7,4,0,1};
  int i=Q->i;
  int j=Q->j;
  int a,b;

  switch(Q->code){

  case QUERY_END:
    Q->wrap=NULL; // ne pas faire free() car variable static
    return 0;

  case QUERY_INIT:
    switch(type){
    case ALK_NOR: Q->wrap=P_NOR; break;
    case ALK_CYC: Q->wrap=P_CYC; break;
    case ALK_ISO: Q->wrap=P_ISO; break;
    case ALK_NEO: Q->wrap=P_NEO; break;
    case ALK_SEC: Q->wrap=P_SEC; break;
    case ALK_TER: Q->wrap=P_TER; break;
    }
    if(n<Q->wrap[0]) RET_n(0);
    RET_n(3*n+2*(type!=ALK_CYC));

  case QUERY_NAME:
    strcpy(Q->name,(i<n)?"C":"H");
    return 0;
    
  case QUERY_ADJ:
    if((i<n)&&(n<=j)){ // i dans C et j dans H
      j-=n;
      if(j<3*Q->wrap[1]) RET_a(j/3==i); // deg(j)=1
      j-=3*Q->wrap[1],i-=Q->wrap[1]; 
      if(j<1*Q->wrap[2]) RET_a(j/1==i); // deg(j)=3
      j-=1*Q->wrap[2],i-=Q->wrap[2]+Q->wrap[3];
      RET_a(j/2==i);                    // deg(j)=2
    }
    if(j>=n) RET_a(0); // i et j >=n => pas d'arête

    switch(type){ // i et j sont dans C

    case ALK_NOR:
      // 0-2-...-1
      RET_a((j-i==n-1)||((j==i+1)&&(i>=1)));

    case ALK_CYC:
      // 0-1-2-...-0
      RET_a((j-i==n-1)||(j==i+1));

    case ALK_ISO:
      // 0
      //  \
      // 1-3-4-...-2
      RET_a( ((i<2)&&(j==3))||
	     ((i==2)&&(j==n-1))||
	     ((j==i+1)&&(i>=3))
	     );

    case ALK_NEO:
      // 0
      //  \
      // 1-4-5-...-3
      //  /
      // 2
      RET_a( ((i<3)&&(j==4))||
	     ((i==3)&&(j==n-1))||
	     ((j==i+1)&&(i>=4))
	     );
      
    case ALK_SEC:
      //     4-...-1
      //    /
      // 0-3-a-...-2
      a=(n-4)/2+4;
      RET_a( ((i==0)&&(j==3))||
	     ((i==1)&&(j==a-1))||
	     ((i==2)&&(j==n-1))||
	     ((i==3)&&(j==a))||
	     ((j==i+1)&&(i>=3)&&(j<a))||
	     ((j==i+1)&&(i>=a))
	     );

    case ALK_TER:
      //     5-...-1
      //    /
      // 0-4-a-...-2
      //    \
      //     b-...-3
      a=(n-5)/3+5;
      b=a+(n-a)/2;
      RET_a( ((i==0)&&(j==4))||
	     ((i==1)&&(j==a-1))||
	     ((i==2)&&(j==b-1))||
	     ((i==3)&&(j==n-1))||
	     ((i==4)&&(j==a))||
	     ((i==4)&&(j==b))||
	     ((j==i+1)&&(i>=4)&&(j<a))||
	     ((j==i+1)&&(i>=a)&&(j<b))||
	     ((j==i+1)&&(i>=b))
	     );
    }

  }
  return 1;
}


int mycielski(query* const Q){
/*
  Utilise i<j.
  Modifie Q->i et Q->j.
*/
  int ki,kj,b,k,t;
  switch(Q->code){

  case QUERY_INIT:
    k=Q->param[0];
    RET_n((k<2)? 0 : 3*(1<<(k-2))-1);

  case QUERY_ADJ:
    ki=ceil(log2((double)(Q->i+2)/3.));
    kj=ceil(log2((double)(Q->j+2)/3.));
    k=3*(1<<kj)-2; /* rem: k est pair */
    b=(Q->j==k);
    if(ki==kj) RET_a(b);
    if(b) RET_a(0);
    t=Q->j-(k/2);
    if(t==Q->i) RET_a(0);
    if(Q->i<t){ Q->j=t; return mycielski(Q); }
    Q->j=Q->i;Q->i=t;
    return mycielski(Q);

  }
  return 1;
}


int windmill(query* const Q){
  switch(Q->code){
  case QUERY_INIT: RET_n(2*Q->param[0]+1);
  case QUERY_ADJ: RET_a((Q->i==0)||((Q->i&1)&&(Q->j==Q->i+1)));
  }
  return 1;
}


int matching(query* const Q){
/*
  Utilise i<j.
  On pourrait utiliser aussi: "fdrg 2n 1 ."
*/
  switch(Q->code){
  case QUERY_INIT: RET_n(2*Q->param[0]);
  case QUERY_ADJ: RET_a((Q->i==Q->j-1)&&(Q->j&1));
  }
  return 1;
}


int ring(query* const Q){
  const int n=Q->param[0];
  switch(Q->code){
  case QUERY_INIT: RET_n(n);
  case QUERY_ADJ:;
    int k;
    const int t=Q->param[1]+2;
    for(k=2;k<t;k++)
      if((Q->j==((Q->i+Q->param[k])%n))||
	 (Q->i==((Q->j+Q->param[k])%n)))
	RET_a(1);
    RET_a(0);
  }
  return 1;
}


int cage(query* const Q){
  const int n=Q->param[0];
  const int k=Q->param[1];
  switch(Q->code){
  case QUERY_INIT:
    if(k<1) RET_n(0);
    RET_n(n);
  case QUERY_ADJ:
    if(Q->j==(Q->i+1)%n) RET_a(1);
    if(Q->i==(Q->j+1)%n) RET_a(1);
    if(Q->j==(Q->i+Q->param[(Q->i%k)+2])%n) RET_a(1);
    if(Q->i==(Q->j+Q->param[(Q->j%k)+2])%n) RET_a(1);
    RET_a(0);
  }
  return 1;
}


int crown(query* const Q){
/*
  Utilise i<j.
*/
  int k=Q->param[0];
  switch(Q->code){
  case QUERY_INIT: RET_n(2*k);
  case QUERY_ADJ: RET_a((Q->i<k)&&(Q->j>=k)&&(Q->i!=(Q->j-k)));
  }
  return 1;
}


int fan(query* const Q){
/*
  Utilise i<j.
  On pourrait utiliser aussi: "grid 1 p -apex q"
*/
  const int p=Q->param[0];
  const int q=Q->param[1];
  switch(Q->code){
  case QUERY_INIT: RET_n(p+q);
  case QUERY_ADJ: RET_a(((Q->j==Q->i+1)&&(Q->j<p))||((Q->i<p)&&(Q->j>=p)));
  }
  return 1;
}


int split(query* const Q){
/*
  Utilise i<j.
  On pourrait utiliser aussi: "ring p 0 -not -apex q"
*/
  switch(Q->code){
  case QUERY_INIT: RET_n(Q->param[0]);
  case QUERY_ADJ: RET_a((Q->i<Q->param[1])||(Q->j<Q->param[1]));
  }
  return 1;
}


int chess(query* const Q){
  const int p=Q->param[0];
  const int q=Q->param[1];
  
  switch(Q->code){

  case QUERY_NAME:
    name_base(Q->name,Q->i,Q->param[0],2,",","()",1);
    return 0;
    
  case QUERY_INIT: RET_n(p*q);
    
  case QUERY_ADJ:;
    const int x=Q->param[2];
    const int y=Q->param[3];
    const int xi=Q->i%p;
    const int yi=Q->i/p;
    const int xj=Q->j%p;
    const int yj=Q->j/p;  
    RET_a(((abs(xi-xj)==x)&&(abs(yi-yj)==y))||((abs(xi-xj)==y)&&(abs(yi-yj)==x)));
  }
  
  return 1;
}


int grid(query* const Q){
/*
  Modifie Q->i et Q->j.
*/
  int x,y,k,z,p,h,b,d=Q->param[0];

  switch(Q->code){

    case QUERY_NAME:
      z=1; /* z=vrai ssi toutes les dimensions sont < 10 */
      int R[DIMAX];
      for(k=0;k<d;k++){
	b=Q->param[k+1];
	R[k]=Q->i%b;
	Q->i /= b;
	z &= (b<11);
      }
      name_vector(Q->name,R,d,(z?"":","),(z?"":"()"),0,"%i");
      return 0;

  case QUERY_INIT:
    free(Q->wrap);ALLOC(Q->wrap,d);
    for(Q->n=1,k=0;k<d;k++){
      p=Q->param[k+1];
      Q->wrap[k]=(p<0);
      p=abs(p);
      Q->param[k+1]=p;
      Q->n *= p;
    }
    return 0;

  case QUERY_END:
    free(Q->wrap);Q->wrap=NULL;
    return 0;

  case QUERY_ADJ:
    z=h=k=b=0;
    while((k<d)&&(b<2)&&(h<2)){
      p=Q->param[k+1];
      x=Q->i%p;
      y=Q->j%p;
      h=abs(x-y);
      if(h==0) z++;
      if((h>1)&&(Q->wrap[k])&&(((x==0)&&(y==p-1))||((y==0)&&(x==p-1)))) h=1;
      if(h==1) b++;
      Q->i /= p;
      Q->j /= p;
      k++;
    }
    RET_a((b==1)&&(z==(d-1)));
  }

  return 1;
}


int clebsch(query* const Q){
  if((Q->code==QUERY_ADJ)&&((Q->i|Q->j)==Q->n-1)&&((Q->i&Q->j)==0)) RET_a(1);
  /* sommets opposés */
  return grid(Q);
}


int rplg(query* const Q){
/*
  Q->param[0]=n
  Q->dparam[0]=t
  Q->dparam[1]=∑_i w_i
  Q->drep[i][0]=degré moyen du sommet i = (n/(i+1))^(1/(t-1))
*/
  switch(Q->code){

  case QUERY_END:
    if(Q->n>0) FREE2(Q->drep,Q->n);
    return 0;

  case QUERY_ADJ:
    RET_a( (RAND01 < (Q->drep[Q->i][0]*Q->drep[Q->j][0]/Q->dparam[1])) );
    
  case QUERY_INIT:
    SET_n(Q->param[0]);
    const double n=Q->n;
    const double c=1/(Q->dparam[0]-1.0);
    double s=0;
    int k;

    ALLOCMAT(Q->drep,Q->n,1);
    for(k=0;k<Q->n;k++) s += (Q->drep[k][0]=pow(n/((double)k+1.0),c));
    Q->dparam[1]=s;

    return 0;
  }
  return 1;
}


int butterfly(query* const Q){
/*
  Modifie Q->i et Q->j.
*/
  int d=Q->param[0]+1; /* d=dim+1 */
  switch(Q->code){

  case QUERY_INIT:
    d--;         /* d=dim */
    Q->n=(1<<d); /* n=2^dim */
    Q->n *= d+1; /* n=(dim+1)*2^dim */
    return 0;

  case QUERY_ADJ:;
    int x=Q->i/d;Q->i%=d; /* i -> (x,i) = (mot binaire,niveau) */
    int y=Q->j/d;Q->j%=d; /* j -> (y,j) = (mot binaire,niveau) */
  
    if(Q->j==Q->i+1) RET_a((x==y)||((x^y)==(1<<Q->i)));
    if(Q->i==Q->j+1) RET_a((x==y)||((x^y)==(1<<Q->j)));
    RET_a(0);
  }
  
  return 1;
}


int debruijn(query* const Q){
  const int d=Q->param[0];
  const int b=Q->param[1];
  int k;

  switch(Q->code){

  case QUERY_INIT:
    if((d<0)||(b<1)) RET_n(0); // il faut d>=0 et b>=1
    for(Q->n=1,k=0;k<d;k++) Q->n *= b;
    return 0;

  case QUERY_ADJ:;
    const int x=Q->j-(Q->i*b)%Q->n;
    const int y=Q->i-(Q->j*b)%Q->n;
    RET_a(((0<=x)&&(x<b))||((0<=y)&&(y<b)));
  }

  return 1;
}


int barbell(query* const Q){
/*
  Utilise i<j.
*/
  const int n1=abs(Q->param[0]);
  const int n2=abs(Q->param[1]);
  const int p=Q->param[2];
  switch(Q->code){

  case QUERY_INIT:
    RET_n(n1+n2+p-1);
  
  case QUERY_ADJ:
    if(Q->j<n1){
      if(Q->param[0]<0) RET_a((Q->j==Q->i+1)||((Q->i==0)&&(Q->j==n1-1))); /* cycle 1 */
      RET_a(1); /* clique 1 */
    }
    if(Q->i>=n1-1+p){
      if(Q->param[1]<0) RET_a((Q->j==Q->i+1)||((Q->i==n1-1+p)&&(Q->j==n1+n2+p-2))); /* cycle 2 */
      RET_a(1); /* clique 2 */
    }
    if((n1-1<=Q->i)&&(Q->j<n1+p)) RET_a((Q->j-Q->i==1)); /* chemin */
    RET_a(0);
  }

  return 1;
}


int shuffle(query* const Q){
  int n=Q->param[0];
  switch(Q->code){

  case QUERY_NAME:
    name_base(Q->name,Q->i,2,n,"","",-1);
    return 0;

  case QUERY_INIT:
    RET_n(1<<n);
  
  case QUERY_ADJ:
    if((Q->i>>1)==(Q->j>>1)) RET_a(1);
    n=Q->n/2;
    if(Q->j==((Q->i<<1)%Q->n+(Q->i>=n))) RET_a(1);
    if(Q->j==((Q->i>>1)+((Q->i&1)?n:0))) RET_a(1);
    RET_a(0);
  }

  return 1;
}


int kautz(query* const Q){
/*
  A chaque sommet i, qui est entier de [0,b*(b-1)^(d-1)[, on associe
  un mot (x_1,...,x_d) codée sous la forme d'un entier rep[i][0] de
  [0,b^d[ (d lettres sur un alphabet à b lettres). Alors i et j sont
  adjacents dans le Kautz ssi rep[i][0] et rep[j][0] sont adjacents
  dans le De Bruijn.

  param[] = d b x y où
   x = #sommets du De Buijn = b^d
   y = #sommets du Kautz    = b*(b-1)^(d-1)
*/
  int d,b,k,u,q,r,s,t,x;
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);
      
  case QUERY_INIT:
    d=Q->param[0];
    b=Q->param[1];
    if((d<1)||(b<2)) RET_n(0); // il faut d>=1 et b>=2
    Q->n=x=b;
    t=b-1;
    for(k=1;k<d;k++){
      Q->n *= t;
      x *= b;
    }
    TEST_n; // fin si n<0
    Q->param[2]=x;    /* #sommets de De Bruijn */
    Q->param[3]=Q->n; /* #sommets de Kautz */
    ALLOCMAT(Q->rep,Q->n,1);

    for(u=0;u<Q->n;u++){ /* pour tous les sommets faire .. */
      /* On voit un sommet u comme (r,s2...sd) avec r dans [0,b[ et
	 s_i de [0,b-1[. On le converti en (x1,...,xd) avec xd=r et
	 x_(d-1)=s2 ou s2+1 suivant si s2<r ou pas, etc. */
      r=x=u%b;
      q=u/b;
      for(k=1;k<d;k++){
	s=q%t;
	s+=(s>=r);
	r=s;
	q=q/t;
	x=x*b+r;
      }
      Q->rep[u][0]=x;
    }
    return 0;

  case QUERY_ADJ:
    Q->n=Q->param[2]; /* modifie le nb de sommets */
    Q->i=Q->rep[Q->i][0];
    Q->j=Q->rep[Q->j][0];
    debruijn(Q); /* adjacence écrite dans Q->a */
    Q->n=Q->param[3]; /* rétablit le nb de sommets */
    return 0;
  }
  
  return 1;
}


int ggosset(query* const Q){
/*
  ggosset p k d_1 v_1 ... d_k v_k, mais
  Q->param[] = p k d d_1 v_1 ... d_k v_k

  Sommets: les permutations du vecteur (v_1...v_1, ..., v_k...v_k) (ou
  de son opposé) telles que le nombre de valeurs entières v_t est
  d_t. On pose Q->rep[i] = vecteur du sommet i qui est de taille
  d=d_1+...+d_k. On a l'arête i-j ssi le produit scalaire entre
  Q->rep[i] et Q->rep[j] vaut p.

  Pour calculer tous les vecteurs possibles de taille d (les sommets)
  à partir de (v_1...v_1,...,v_k...v_k) on procède comme suit (codage
  d'un multi-ensemble à k valeurs):

  On choisit les positions de v_1 (il y en a Binom(d,d_1) possibles,
  puis on choisit les positions de v_2 parmi les d-d_1 restantes (il y
  en a Binom(d-d_1,d_2)), puis les positions de v_3 parmi les
  d-d_1-d_2 retantes, etc. Pour chaque t, les positions des d_t
  valeurs v_t sont codées par un sous-ensemble S[t] de [0,d[ de taille
  d_t.
*/
  const int d=Q->param[2];
  const int k=Q->param[1];
  int t,p,m,u,v,c,z,**S;

  switch(Q->code){
    
  case QUERY_END:
    return free_rep(Q);

  case QUERY_NAME:
    name_vector(Q->name,Q->rep[Q->i],d,",","[]",1,"%i");
    return 0;

  case QUERY_ADJ:
    /* Calcule le produit scalaire de Q->rep[i] et Q->rep[j] */
    for(p=t=0;t<d;t++) p += Q->rep[Q->i][t]*Q->rep[Q->j][t]; 
    RET_a((p==Q->param[0]));
    
  case QUERY_INIT:

    /* Calcule Q->n */
    if(d<=1) Erreur(6); /* paramètre incorrect */
    m=d; Q->n=2; 
    for(t=3;m>0;t += 2){
      p=Q->param[t]; /* p=d_i */
      if(p<0) Erreur(6); /* paramètre incorrect */
      Q->n *= Binom(m,p);
      m -= p;
    }
    TEST_n; // fin si graphe vide
    ALLOCMAT(Q->rep,Q->n,d); /* vecteur de taille d représentant les sommets */
    NALLOC(int,P,d); /* tableau intermédiaire */
    ALLOCMAT(S,k,d); /* on réserve k tableaux (sous-ensembles) de taille <= d */
    for(t=0;t<k;t++) NextSet(S[t],-1,d); /* initialise les sous-ensembles */
    /* NB: taille |S[t]|=d_{t-1}=Q->param[2t+3] et v_{t-1}=Q->param[2t+4] */

    for(u=0;u<Q->n;u+=2){
      /* Pour chaque sommet u on fait:

	 1. on remplit rep[u] et rep[u+1] à partir des sous-ensembles S[0]...S[k-1]
	 2. on incrémente les sous-ensembles S[0]...S[k-1]
	 
	 Pour l'étape 1, il faut passer par un tableau intermédiaire P
	 puis remplir rep[u] et rep[u+1] à partir de P. Supposons d=5,
	 k=3, et S[0]={1,3} S[1]={1}, S[2]={0,1}.  On met dans P les
	 indices t des S[t] aux bons emplacements comme suit:

           - au départ P est vide: P={-,-,-,-,-}
	   - puis on ajoute S[0]:  P={-,0,-,0,-}
	   - puis on ajoute S[1]:  P={-,0,1,0,-}
	   - puis on ajoute S[2]:  P={2,0,1,0,2}
      */

      /* Calcule P */
      for(t=0;t<d;t++) P[t]=-1; /* tableau vide au départ */
      for(t=0;t<k;t++){ /* pour chaque sous-ensemble S[t] */
	m=-1;
	z=Q->param[2*t+3];
	for(p=v=0;p<z;p++){ /* on parcoure S[t] */
	  /* mettre t dans la S[t][p]-ème case vide de P à partir de l'indice v */
	  /* utilise le fait que dans S[t] les nombres sont croissant */
	  /* v=position courante de P où l'on va essayer d'écrire t */
	  /* c=combien de cases vides de P à partir de v faut-il encore sauter ? */
	  /* si P[v]=-1 et c=0 alors on écrit P[v]=t */
	  c=S[t][p]-m;
	  m=S[t][p]; /* mémorise l'ancienne valeur de S[t][p] */
	  while(c>0){
	    if(P[v]<0){
	      c--;
	      if(c==0) P[v]=t; /* écrit t et passe à la case suivante */
	    }
	    v++;
	  }
	}
      }

      /* Remplit rep[u] et rep[u+1] grâce au tableau P (et incrémenter u) */
      for(t=0;t<d;t++){
	v=Q->param[2*P[t]+4]; /* valeur v_t à écrire */
	Q->rep[u][t]=v;
	Q->rep[u+1][t]=-v;
      }

      /* Incrémente S[0]...S[k-1] grâce à NextSet() */
      t=0; /* on commence avec S[0] */
      m=d; /* S[t] dans {0...d-2} */
      v=3; /* Q->param[v] = taille de S[t] */
      while((t<k)&&(!NextSet(S[t],m,Q->param[v]))){
	t++; /* si S[t] fini on passe à S[t+1] */
	m -= Q->param[v];
	v += 2;
      }
      /* si t=k alors on a atteint le dernier sommet */
    }

    free(P);
    FREE2(S,k);
    return 0;

  }
  return 1;
}


int schlafli(query* const Q){
/*
  Sous-graphe induit par le voisinage d'un sommet du graphe de gosset
  qui a deux fois plus de sommets. Les voisins T du sommets 0 sont:

       T = 2 4 6 8 11 12 14 17 19 20 22 25 27 29 30 32 35 37 39 41 42 44 47 49 51 53 55
    2u+2 = 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54
       D = 0 0 0 0 +1 0  0  +1 +1 0  0  +1 +1 +1 0  0  +1 +1 +1 +1 0  0  +1 +1 +1 +1 +1
         = 00  001    0011        00111          001111            0011111

  En fait on pourrait calculer T[u] en fonction de u comme suit.  La
  position des début de chaque bloque dans D = T - (2u-2) vaut
  0,2,5,9,14,20 soit (b+1)(b/2+2) pour b=0,1,...,4. On a donc
  T[u]=2u+2+D[u]. Si u=0...26 est le numéro de voisin, alors on
  cherche le plus grand entier k possible avec (k+1)(k/2+2) <= u <
  (k+2)(k/2+5/2). On pose k=ceil(sqrt(2u))-1, si x=(k+1)(k/2+2) > u,
  alors on fait k--. Du coup D[u]=(u-x>1).

  Modifie Q->param.
*/
  switch(Q->code){

  case QUERY_END:
  case QUERY_ADJ:
    return ggosset(Q);

  case QUERY_INIT:
    /* paramètres (modifiés) identiques à gosset */
    Q->param[0]=8; Q->param[1]=2; Q->param[2]=8;
    Q->param[3]=2; Q->param[4]=3; Q->param[5]=6; Q->param[6]=-1;
    ggosset(Q);

    int u;
    /* T[u]=u-ème voisin du sommet 0, cf. gengraph gosset -format vertex 0 */
    const int T[]={2,4,6,8,11,12,14,17,19,20,22,25,27,29,
		   30,32,35,37,39,41,42,44,47,49,51,53,55};
    for(u=0;u<27;u++){
      free(Q->rep[u]); // NB: u<T[u]
      Q->rep[u]=Q->rep[T[u]]; // ici Q->rep[T[u]] ne peut pas être NULL ou libéré
      Q->rep[T[u]]=NULL;
    }
    for(;u<56;u++) free(Q->rep[u]);
    // ici il n'y a que Q->rep[0]...Q->rep[26] qui sont encore aloués, et
    // qui reste à libérer à la fin de la fonction (ADJ_END)
    Q->n=27;
    return 0;

  }
  return 1;
}


int linial(query* const Q){
/*
  Q->rep[i][0...t-1] = représentation du nom du sommet i (sa vue).
*/
  const int m=Q->param[0];
  const int t=Q->param[1];
  int k,u;
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);
    
  case QUERY_NAME:
    if(m<10) name_vector(Q->name,Q->rep[Q->i],t,"","",1,"%i");
    else name_vector(Q->name,Q->rep[Q->i],t,",","()",1,"%i");
    return 0;

  case QUERY_INIT:;
    if((m<t)||(t<1)) RET_n(0); /* graphe vide si n<t */
    for(Q->n=1,k=m-t+1;k<=m;k++) Q->n *= k; /* calcul Q->n */
    ALLOCMAT(Q->rep,Q->n,t);
    NALLOC(int,S,t);
    NALLOC(int,P,t);
    NextArrangement(S,P,-1,t); /* initialisation de S et P */
    for(u=0;u<Q->n;u++){ /* génère tous les arrangements */
      for(k=0;k<t;k++) Q->rep[u][k]=S[P[k]];
      NextArrangement(S,P,m,t);
    }
    free(P);
    free(S);
    return 0;
    
  case QUERY_ADJ:
    if((Q->rep[Q->i][0]!=Q->rep[Q->j][t-1])||(Q->param[0]==t)){
      for(u=1,k=0;u<t;u++,k++) if(Q->rep[Q->i][u]!=Q->rep[Q->j][k]) break;
      if(u==t) RET_a(1);
    }
    if(Q->i<Q->j){ SWAP(Q->i,Q->j,k); return linial(Q); }
    RET_a(0);
  }

  return 1;
}


int linialc(query* const Q){
  if(Q->code==QUERY_INIT){
    int u,k,v,m1,x,y;
    int m=Q->param[0];
    int t=Q->param[1];
    Q->n=m; m1=m-1;
    for(Q->n=m,u=m1,k=1;k<t;k++) Q->n *= u; /* calcule Q->n=m*(m-1)^t */
    ALLOCMAT(Q->rep,Q->n,t);
    /* on transforme u en (x0,x1,...x_t) avec x0 in [0,m[ et x_i in [0,m-1[ */
    for(u=0;u<Q->n;u++){
      x=Q->rep[u][0]=(u%m);
      for(v=u/m,k=1;k<t;k++){
	y = v%m1; /* y in [0,m-1[ */
	v /= m1;
	x=Q->rep[u][k]=y+(y>=x); /* si x=y, on incrémente y */
      }
    }
    return 0;
  }
  return linial(Q);
}


int gpstar(query* const Q){
/*
  Q->rep[i][0...n-1] = représentation de la permutation du sommet i.
*/
  int k,u;
  const int n=Q->param[0];
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);

  case QUERY_NAME:
    if(n<10) name_vector(Q->name,Q->rep[Q->i],n,"","",1,"%i");
    else name_vector(Q->name,Q->rep[Q->i],n,",","()",1,"%i");
    return 0;

  case QUERY_INIT:
    if(n<=0) RET_n(0);
    for(Q->n=1,k=2;k<=n;k++) Q->n *= k;
    ALLOCMAT(Q->rep,Q->n,n); /* ici Q->n>0 */
    NALLOCZ(int,P,n,_i); /* initialise P */
    /* génère toutes les permutations */
    for(u=0;u<Q->n;u++){
      for(k=0;k<n;k++) Q->rep[u][k]=P[k]+1; /* copie P dans Q->rep */
      NextPermutation(P,n,NULL);
    }
    free(P);
    return 0;

  case QUERY_ADJ:
    /* distance de Hamming: on compte les différences entre les
       tableaux rep[i] et rep[j] */
    for(k=u=0;k<n;k++) u += (Q->rep[Q->i][k]!=Q->rep[Q->j][k]);
    RET_a((u==Q->param[1]));
  }
  
return 1;
}


int pancake(query* const Q){
/*
  Utilise i<>j.
  Q->rep[i][0...n-1] = représentation de la permutation du sommet i.
  On utilise Q->param[1] pour le signe de la permutation:
   s=+1 pour pancake()
   s=-1 pour bpancake()
*/
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);

  case QUERY_INIT:
    gpstar(Q); // génère les n! permutations
    Q->param[1]=1; // signe>0
    return 0;

  case QUERY_ADJ:
    /* test d'adjacence à partir de rep[i] et rep[j] */
    if(Q->i==Q->j) RET_a(0); // il est important d'avoir i<> j */
    const int s=Q->param[1]; // signe pour le reversal
    int k=Q->param[0],z=0;
    do k--; while(Q->rep[Q->i][k]==Q->rep[Q->j][k]); /* i<>j, donc on s'arrête toujours */
    while((k>=0)&&(Q->rep[Q->i][k]==s*Q->rep[Q->j][z])) k--,z++; /* teste le "reversal" signé */
    RET_a((k<0)); /* adjacent si préfixe=reversal */
  }
  
  return 1;
}


int bpancake(query* const Q){
/*
  rep[i][0...n-1] = représentation de la permutation signée du sommet
  i. Il s'agit d'une valeur de {1,2,...,n,-1,-2,...,-n} (on évite
  soigneusement 0, car -0=+0.
*/
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);

  case QUERY_ADJ:
    return pancake(Q);

  case QUERY_NAME:
    name_vector(Q->name,Q->rep[Q->i],Q->param[0],"","",1,"%+i");
    return 0;

  case QUERY_INIT:
    Q->param[1]=-1; // signe<0
    const int n=Q->param[0];
    int k,p,q,t,c,u;
    if(n<=0) RET_n(0);
    for(t=1,k=2;k<=n;k++) t *= k; /* calcule t=n! */
    p=(1<<n); /* p=2^n */
    SET_n(p*t); /* Q->n=#sommets, forcément >0 */
    ALLOCMAT(Q->rep,Q->n,n); /* permutations signées représentant les sommets */
    NALLOCZ(int,P,n,_i); /* initialise une permutation P (non signée) */

    /* Génère toutes les permutations signées. On représente les
       signes par les bits de l'entier q=0...2^n-1. Si bit à 1 -> +1,
       sinon -> -1 */

    for(c=u=0;c<t;c++){ /* répète n! fois */
      for(q=0;q<p;q++,u++) /* répète 2^n fois */
	for(k=0;k<n;k++){ /* (1<<k)=mask=bit-0,bit-1,bit-2...,bit-(n-1) */
	  Q->rep[u][k]=P[k]+1; /* copie P dans Q->rep avec le signe +/-1 */
	  if(q&(1<<k)) Q->rep[u][k] *= -1; /* copie P dans Q->rep avec le signe +/-1 */
	}
      NextPermutation(P,n,NULL);
    }

    free(P);
    return 0;
  }

  return 1;
}


int pstar(query* const Q){
/*
  rep[i][0...n-1] = représentation de la permutation du sommet i.
*/

  switch(Q->code){
    
  case QUERY_END:
  case QUERY_INIT:
    return gpstar(Q);

  case QUERY_ADJ:
    /* Il faut deux différences dont le premier chiffre */
    Q->a=(Q->rep[Q->i][0]!=Q->rep[Q->j][0]);
    if(Q->a) return gpstar(Q);
    return 0; // ici Q->a=0
  }

  return 1;
}
  

int gabriel(query* const Q){
/*
  L'initialisation et la terminaison sont communes à beaucoup de
  graphe géométriques. L'adjacence est en O(n). Il est important de
  tester les n sommets (même ceux supprimés par -delv), sinon le
  résultat n'est pas un sous-graphe.
*/
  const int n=Q->param[0];

  switch(Q->code){

  case QUERY_END:
    return free_pos(Q);
 
  case QUERY_INIT:
    SET_n(n);
    InitXY(Q);
    return 0;

  case QUERY_ADJ:;
    /* c=(xc,yc)=milieu du segment i-j */
    const double xc=(XPOS[Q->i]+XPOS[Q->j])/2;
    const double yc=(YPOS[Q->i]+YPOS[Q->j])/2;
    const double r=fmin(Norme_dxy(fabs(xc-XPOS[Q->i]),fabs(yc-YPOS[Q->i])),
			Norme_dxy(fabs(xc-XPOS[Q->j]),fabs(yc-YPOS[Q->j])));
    /*
      r=dist(c,i)=rayon du disque centré au milieu de i-j. Attention !
      en général r<>dist(i,j)/2 car ce n'est pas forcément la norme
      L2. A cause des arrondis, il est possible que r soit surestimé,
      et donc que dist(c,i)<>dist(c,j).  C'est pour cela qu'on prend
      r=min(dist(c,i),dist(c,j)).
    */

    int k;
    for(k=0;k<n;k++)
      if(Norme_dxy(fabs(xc-XPOS[k]),fabs(yc-YPOS[k]))<r) RET_a(0);
    RET_a(1);
  }

  return 1;
}
  

/*
  Une surface S de genre g orientée ou non, avec ou sans bord est
  modélisée par un polygone convexe régulier P à p=4g cotés numérotés
  de 0 à p-1 de manière alternée comme suit (on parle de polygone
  fondamental https://en.wikipedia.org/wiki/Fundamental_polygon). Le
  coté 0 est positionné verticalement et le plus à droit. Une fois le
  coté i positionné, le coté i+1 est, si i est pair, situé deux cotés
  après le coté i en tournant dans le sens direct (le sens contraire
  des aiguilles d'une montre), et si i est impair, situé deux cotés
  avant. Voir l'exemple suivant pour g=2.

                                   1
                                _______
                               /       \
                            3 /         \ 2
                             /           \
                             |           |
                           4 |     S     | 0
                             |           |
                             \           /
                            6 \         / 7  
                               \_______/
                                   5    
 
  Les cotés sont appariés pour former les coutures de S. La surface S
  possède donc 2g coutures. C'est la taille du tableau XYsurface
  codant la signature de S, c'est-à-dire la façon dont les coutures
  sont réalisées. Si l'on découpe S le long de ses coutures, alors on
  obtient tout simplement P. Le coté i est apparié au coté i+1 si i
  est pair, et au coté i-1 sinon. Dit autrement, le coté i est apparié
  avec le coté i^1 (i xor 1 = inversion du bit de poids faible de
  i). Les coutures sont numérotées de 0 à 2g-1, le numéro de la
  couture du coté i étant simplement i>>1 (soit floor(i/2)). Les cotés
  appariés doivent avoir la même longueur. Pour simplifier, on
  supposera que tous les cotés de P sont de même longueur, soit
  2sin(π/p) en supposant que P est inscrit dans un cercle de rayon
  unité.

  Chaque couture de numéro c correspond donc à l'appariement des cotés
  2c et 2c+1. Elle possède une signature XYsurface[c] qui peut prendre
  trois valeurs (+1,-1,0) codant les trois façons de réaliser une
  couture. La signature >0 signifie que les cotés 2c et 2c+1 sont
  recollés sans inversion (handle) comme dans une surface
  orientable. La signature <0 signifie que les cotés 2c et 2c+1 sont
  recollés avec une inversion (crosscap) comme dans une surface
  non-orientable. Enfin, la signature =0 signifie que les cotés 2c et
  2c+1 forment un bord (ou trou) de la surface (border).

  On peut assembler des copies de P pour former un graphe de
  polygones. Ce graphe est a priori infini. Les sommets sont des
  copies de P recollés par les cotés appariés (les coutures) après
  rotation, translation voir renversement. Un polygone Q a autant de
  voisins qu'il a de coté i qui ne correspond pas à un bord (il faut
  donc que XYsurface[i>>1] soit non nul). On va voir ci-après qu'on
  peut en fait limiter le graphe des polygones à un arbre d'environ
  g^O(g) polygones. Un chemin sur la surface correspond donc à un
  chemin dans le graphe des polygones, passant d'un polygone à un
  polygone voisin chaque fois qu'il tranverse une couture.

  Pour obtenir le voisin R du coté i du polygone Q, on procède comme
  suit (qui reste valide pour le disque hyperbolique):

  1) on fait une copie de Q, appelée R, qu'on supperpose à Q;

  2) on fait une rotation de R correspondant à deux cotés du polygone
     (4π/p) dans le sens direct si i est impair, et inverse si i est
     pair, si bien que les cotés i de Q et i^1 de R sont confondus;

  3) on effectue une inversion de R autour du coté i de Q, soit une
     symétrie par rapport au coté i. (Dans le cas hyperbolique, c'est
     une inversion par rapport à le cercle représentant le coté i);

  4) si la signature de la couture i>>1 est négative il faut de plus
     faire une symétrie axiale de R selon la perpendiculaire passant
     par le coté i^1 (=renversement).

  Sans avoir de preuve, on va supposer qu'une géodésique d'un point u
  à un point v sur la surface S correspond sur P (une fois S découpée)
  à une suite de segments de droite, chaque segment ayant comme
  extrémité u, v ou un point du périmètre de P. De plus, ces segments
  peuvent être juxtaposés pour former une droite allant de u à une
  copie de v dans l'un des polygones du graphe des polygones. On va
  aussi supposer qu'une géodésique ne coupe jamais deux fois la même
  couture. Ainsi on peut restreindre le graphe des polygones à un
  arbre dont la racine est une copie de P et de profondeur <= 2g, avec
  la restriction que polygone n'a pas de fils ni avec les cotés
  correspondant avec une couture déjà utilisée par un de ces ancêtres
  (en particulier par son père si ce n'est pas la racine).

  Ainsi, chaque point u de P se retrouve copié dans chacun des
  polygones de l'arbre. On repère chaque copie par la suite des
  numéros des cotés permettant de passer d'un polygone à son voisin,
  soit un chemin dans l'arbre des polygones. On supposera que les
  chemins sont valides, c'est-à-dire qu'il ne contient pas de coté
  correspondant à un bord et que deux cotés appariés ne se suivent
  jamais.

  double surface_geodesic(point u,point v,int *C,point *D)

  // construit la géodésique allant du point u au point v(C). Renvoie
  // la longueur de la géodésique ou -1 si elle n'existe pas. Pour
  // qu'elle existe il faut que la droite entre u et v(C) traverse
  // effectivement dans l'ordre les cotés indiqués par C. On renvoie
  // également dans le tableau D (qui sera supposé de taille assez
  // grande) la liste des points d'intersections de u à v(C).

  Ainsi pour trouver une géodésique entre u et v, il faut calculer une
  géodésique vers v(C) pour toutes les copies possibles de v, et de
  prendre la plus courte. Notons qu'elle existe bien toujours car la
  géodésique entre u et v(C) avec C={-1} existe toujours.

  Pour la construction d'une droite sur le disque poincaré, voir:
  https://en.wikipedia.org/wiki/Poincar%C3%A9_disk_model#Compass_and_straightedge_construction
*/


int surface_next(int *C){
/*
  Détermine le chemin valide (sans bord ni deux fois la même couture)
  immédiatement après le chemin C, vu comme un compteur qu'on essaye
  d'incrémenter. Renvoie 1 si on a pu le faire, 0 sinon. Le chemin
  doit être valide en entrée. Le chemin C identifie la copie d'un
  polygone en précisant la suite (terminée par -1) des cotés qu'il
  faut traverser depuis le polygone racine pour y arriver.
  
  On utilise la fonction comme ceci:

    C[0]=-1; // initialisation du chemin
    do{
    ...; // traitement du chemin C
    }while(surface_next(C)); // chemin suivant

  Ex:  signature: bb      signature: hb      signature: hh
       1 chemin:          3 chemin:          13 chemins:
         C = -1             C = -1             C = -1
                            C = 0 -1           C = 0 -1
                            C = 1 -1           C = 1 -1
			                       C = 2 -1
			                       C = 3 -1
			                       C = 2 0 -1
					       C = 3 0 -1
					       C = 2 1 -1
					       C = 3 1 -1
					       C = 0 2 -1
					       C = 1 2 -1
					       C = 0 3 -1
					       C = 1 3 -1

  Comme on le voit, les chemins sont parcourus par longueur
  croissante, puis par valeur (comme un entier écrit en base p=#coté
  du polygone) lu de droit à gauche. Certains chemins générés peuvent
  mener au même polygone, comme [0 2] et [2 0] dans le cas d'un
  tore. Il y a donc 13 chemins valides, mais seulement 9 polygones
  différents. L'orientation de la surface n'a pas d'impact sur la
  sortie de la fonction. Seuls les bords en ont un.

  Pour améliorer la complexité, on utilise des variables statiques en
  supposant que les appels à la fonction se suivent, c'est-à-dire
  qu'on ne change pas à la main C[] (sinon en faisant C[0]=-1) entre
  deux appels à la fonction. La complexité amortie est en 2^O(g), où
  g=XYsurface, car pour visiter tous les environ g! chemins valides on
  énumère environ g^g chemins (compteurs en base g de taille g). En
  moyenne, cela fait donc g^g/g! ≃ e^g pour chaque appel, à des
  polynômes en g près. Comme il s'agit d'une permutation contrainte,
  on pourrait s'inspirer de l'algo de NextPermutation() qui serait
  bien plus efficace (en g! au lieu de g^g) et qui est le suivant: (il
  faudrait en fait raisonner sur les coutures plutôt que les cotés)

  1. Trouver le plus grand index i tel que C[i] < C[i+1].
     S'il n'existe pas, la dernière permutation est atteinte.
  2. Trouver le plus grand indice j tel que C[i] < C[j].
  3. Echanger C[i] avec C[j].
  4. Renverser la suite de C[i+1] jusqu'au dernier élément.

  Autres points d'améliorations:

  1. Certains chemins mêmes au même polygone, et donc il n'ont pas
     besoin d'être tous testés. Par exemple, pour le tore, [0 2] et [2
     0]. Il faudrait arriver à n'en prendre qu'un seul ou le rendre
     canonique.

  2. Certaines successions de cotés ne peuvent correspondre à un plus
     court chemin. La distance la plus grande dans un polygone inscrit
     dans un cercle de rayon 1 est 2. Cela interdit a priori des
     successions de cotés, comme la succession de trois cotés
     diamétralement opposés.
  
*/
  
  static int F[SURFACEMAX]; // F[c]=fréquence des coutures de C[]
  static int ok=0; // ok=vrai ssi le chemin courant est valide
  static int p; // p=nombre de cotés du polygone, p>=4

  int i,j,c,z;

  if((C[0]<0)||(!ok)){ // initialisation de F[]
    for(c=i=0;c<XYsurfacesize;F[c++]=0); // met tout à 0
    while(C[i]>=0) F[C[i++]/2]=1; // met des 1 pour chaque couture (suppose C valide)
    ok=1; p=2*XYsurfacesize;
  }

  for(i=0;;){

    // ici on suppose qu'on a pas réussit à incrémenter aucune des
    // valeurs d'indice < i, que F[] et ok sont à jour. On essaye
    // d'incrémenter C[i] de sorte que le nouveau coté C[i] ne soit ni
    // un bord ni corresponde à une couture déjà utilisée.
    
    if(C[i]<0){ // on arrive à la fin du chemin courant
      if(i==XYsurfacesize) return ok=0; // on a tout exploré
      C[i+1]=-1; // il reste à mettre à jour C[i] qui vaut ici -1
    }
    
    // ici on cherche le plus petit coté après C[i] (modulo p) qui
    // n'est pas un bord. On passe en revue tous les cotés possibles.
    // Il faut répéter le for(j=...) p fois si C[i]=-1, sinon p-1 fois
    // pour ne pas retomber sur C[i]. NB: on entre toujours dans le
    // for(j=...) car p>1, et à la fin de la boucle, on a toujours
    // C[i]<>c.
    
    c=C[i]; // c=C[i] initial
    z=1; // z=vrai si le coté trouvé est > c (sinon on est passé par 0)
    for(j=(C[i]>=0);j<p;j++){ // répète p ou p-1 fois
      C[i]++;
      if(C[i]==p) C[i]=z=0; // il faudra incrémenter C[i+1] (retenue)
      if(XYsurface[C[i]/2]) break; // on a trouvé un coté qui n'est pas un bord 
    }
    if(j==p) return ok=0; // on a pas trouvé de coté valide

    // on met à jour F[] et ok car on a modifié C[i]
    if(c>=0){ // enlève la couture initiale, si elle existe
      F[c/2]--; // NB: au départ F[c/2]>0, forcément
      // met à jour ok qui peut rester constant ou passer de faux à vrai
      if(!ok){ ok=1; j=0; while(C[j]>=0) ok &= (F[C[j++]/2]<=1); }
    }
    F[C[i]/2]++; // ajoute la nouvelle couture
    ok &= (F[C[i]/2]<=1); // met à jour ok qui peut rester constant ou passer à faux

    if(z){ // on a trouvé le coté sans boucler (incrément sans retenue)
      if(ok) return 1; // on a terminé
      i=0; // on aurait du finir ici, mais ok=0, donc on recommence
    }else i++; // on a trouvé un coté mais il y a eut une retenue (-> i+1)
  }
  
}


double surface_geodesic(point u,point v,int *C){
  u.x=v.x=C[0]; // pour éviter le Warning à la compilation
  return -1;
}


point surface_image(point u,int *C){
/*
  Donne les coordonnées du point u(C) correspondant à la copie de u en
  suivant le chemin C dans l'arbre des polygones. C est une suite
  d'entiers terminée par -1 et supposé valide. En particulier elle ne
  contient aucun coté qui est un bord. L'algorithme est en complexité
  O(|C|).

  Algorithme: On se déplace de polygone en polygone suivant le chemin
  C en construisant le chemin décrit par les centres. Le polygone
  courant est repéré par trois éléments:

    (1) son centre (x,y);
    (2) son décalage d (d=numéro de coté du polygone le + à droit);
    (3) son orientation s (s=1 pour le sens direct, s=-1 sinon).

  Au départ, x=y=0, d=0, et s=1. L'objectif est donc de mettre à jour
  le centre, le décalage et l'orientation pour le polygone voisin en
  traversant le coté i. Une fois le polygone final calculé (centre,
  décalage, orientation), on calcule les coordonnées de v dans ce
  polygone là. Pour cela on effectue une rotation dont l'angle dépend
  du décalage d et dont le sens (direct ou non) dépend de
  l'orientation s.

  Pour trouver de numéro de secteur j correspondant au coté de numéro
  i, c'est-à-dire l'entier j de [0,2p[ tel que le coté i est la corde
  du cône d'angle compris entre j*2π/p et (j+1)*2π/p, il suffit
  d'échanger les deux derniers bits de i:

                         coté    T   secteur

                         ...00  <->  ...00
			 ...01  <->  ...10
			 ...10  <->  ...01
			 ...11  <->  ...11

  Et bien sûr la même transformation, notée T, permet de passer du
  numéro de secteur j au numéro de coté i. On a donc i=T(T(i)).
  L'échange des deux derniers bits de la variable i peut se faire
  ainsi:

     Méthode 1: b=i&3; j=i+(b==1)-(b==2);
     Méthode 2: const int T[]={0,-1,1,0}; j=i+T[i&3];

  Pour calculer le nouveau décalage d' du polygone P' voisin par le
  coté i d'un polygone P de décalage courant d et d'orientation s, on
  peut faire ainsi. On imagine qu'on se déplasse du centre de P au
  centre de P' en traversant le coté i de P. On arrive alors par le
  coté i'=i^1 de P'. Si j est un numéro de secteur, on note opp(j) le
  numéro de secteur opposé à j. Il ne dépend pas de l'orientation.
  Bien sûr on a opp(j)=(j+p/2)%p où p est le nombre de cotés des
  polygones. Le secteur correspondant à d' vaut alors:

          ( opp(T(i')) + T(d)-T(i) + p/2 )%p  si s'=s
	  ( opp(T(i')) + T(i)-T(d) + p/2 )%p  sinon

  Pour les décalages (d), on a donc a priori intérêt à travailler avec
  les numéros de secteur plutôt que les numéros de coté.
*/
  static const int T[]={0,-1,1,0}; // transformation coté <-> secteur
  const int p=2*XYsurfacesize; // p=nombre de coté du polygone
  const double a0=M_2PI/p; // a0=angle du segment reliant deux centres de polygones voisins
  const double r0=2*cos(a0/2); // r0=distance entre deux centres de polygones voisins

  double a;
  int t,j,i;

  // configuration du polygone courant
  double x=0,y=0; // x,y=position du centre
  int d=0; // d=numéro de coté vertical le plus à droite
  int s=1; // s=orientation

  // calcule la position du centre du polygone déterminé par le chemin C
  for(t=0;C[t]>=0;t++){ // pour chaque coté C[t]
    i=C[t]; // i=coté courant
    j=i+T[i&3]; // j=numéro de secteur
    a=s*a0*(j-d); // a=angle du segment reliant le centre courant au prochain
    x += r0*cos(a);
    y += r0*sin(a);
    s *= -XYsurface[t];
    // A FINIR: mise à jour de d ? (voir plus haut)
  }

  // A FINIR: il reste à ajouter le point v au polygone final
  point v={x,y};
  u.x += x;
  u.y += y;
  return v;
}


int sgabriel(query* const Q){
/*
  Graphe de Gabriel défini sur une surface, elle-même définie par
  l'option -xy surface. En partie codée par Louis Démoulins
  (stagiaire de Licence 2 en juin/juillet 2016).
*/
  if(XYsurfacesize==0) return gabriel(Q);
  // ici genre>0
  
  switch(Q->code){

  case QUERY_END:
  case QUERY_INIT:
    return gabriel(Q);

  case QUERY_DOT:
    USERDOT.adj=sgabriel;
    return 0;

  case QUERY_ADJ:;
    const int n=Q->param[0];
    const double dmax=Norme_dxy(2,2); // distance max dans la surface [-1,-1]x[+1,+1]
    const point pi={XPOS[Q->i],YPOS[Q->i]}; // pi=point i
    const point pj={XPOS[Q->j],YPOS[Q->j]}; // pj=point j
  
    //USERDOT.i=i;
    //USERDOT.j=j;

    NALLOC(int,cj,XYsurfacesize+1);
    NALLOC(int,cz,XYsurfacesize+1);

    double dz,d; // distances
    double r; // r=distance entre milieu m de pi-pj et pi
    point m,pz;
    int z;
    cj[0]=-1; // cj=chemin pour la copie de pj

    do{ // pour toutes les copies cj de j faire:
      z=0; // il faut définir z à cause du "continue"
      if(surface_geodesic(pi,pj,cj)<0) continue; // la géodésique n'existe pas
      m=surface_image(pj,cj); // m=pj(cj)
      m.x=(m.x+pi.x)/2, m.y=(m.y+pi.y)/2; // m=milieu entre pi et pj(cj)
      r=Norme_dxy(fabs(m.x-XPOS[Q->i]),fabs(m.y-YPOS[Q->i])); // r=distance(pi,m)
      // il faut vérifier qu'aucun sommet est dans le disque de rayon r autour de m
      for(;z<n;z++){ // NB: ici z=0
	if(z==Q->j) continue;
	// cherche la copie de z la plus proche de m
	pz.x=XPOS[z], pz.y=YPOS[z]; // pz=point d'indice z, cz=chemin de z
	cz[0]=-1; // cz=chemin pour le copie de pz
	dz=dmax; // dz=distance entre le meilleur z et m
	do{
	  d=surface_geodesic(m,pz,cz);
	  if(d>=0) dz=fmin(dz,d); // met à jour la distance
	}while(surface_next(cz));
	if(dz<r) break; // pz est dans le disque -> on sort du for(z=...) avec z<n
      }
    }while((z<n)&&(surface_next(cj))); // si z<n, il faut essayer une autre copie de pj
    
    // ici:
    // z=n (-> 1, car il y a aucun point dans le disque de rayon r)
    // z<n (-> 0, un des points est dans le disque)
    
    free(cz), free(cj);
    RET_a((z==n));
  }

  return 1;
}
  
  
int pat(query* const Q){
/*
  Graphe issu du jeu de Pat Morin. Un sommet correspond à un sommet
  (x,y) d'une des k grilles. On suppose que i<j, ce qui revient à dire
  que la grille de i est placée avant ou est égale à celle de j.

  Exemple: p=q=3 et k=4

  06 07 08  15 16 17  24 25 26  33 34 35
  03 04 05  12 13 14  21 22 23  30 31 32
  00 01 02  09 10 11  18 19 20  27 28 29

  Grille 0  Grille 1  Grille 2  Grille 3

*/
  const int p=Q->param[0]; // p=colonne
  const int q=Q->param[1]; // q=ligne
  const int r=Q->param[2]; // r=round
  const int pq=p*q;

  switch(Q->code){

  case QUERY_END:
    return free_pos(Q);

  case QUERY_INIT:
    SET_n(pq*r);
    if((p<=0)||(q<=0)||(r<=0)) RET_n(0);
    
    /* Détermine les coordonnées des points pour en controler le
       dessin. Les grilles sont placées à des hauteurs de plus en plus
       grandes pour "voir" les arêtes. Au départ (z=0), la première
       grille est mise à une hauteur 0, puis à une hauteur 1, puis à
       une hauteur 3, puis à une hauteur 6, etc. La hauteur de la
       grille pour z quelconque est z(z+1)/2. En fait, on utilise
       z(z+1)/2.5 pour une meilleure lisibilité. */

    ALLOC(XPOS,Q->n);
    ALLOC(YPOS,Q->n);
    int u,x,y,z;
    double h;

    for(u=0;u<Q->n;u++){// pour tous les sommets, faire:
      z=u/(pq);      // z=grille
      x=u%p;         // x=colonne
      y=(u%(pq))/p;  // y=ligne
      h=z*(z+1)/2.5; // h=décalage vers le haut de la grille z
      XPOS[u]=(double)(x+z*p)/(double)imax(p,q);
      YPOS[u]=(double)(y+h*q)/(double)imax(p,q);
    }

    XYtype=XY_USER; /* coordonnées fixées par l'utilisateur */
    InitXY(Q); // pour les options -xy noise/scale ... */
    return 0;

    case QUERY_ADJ:;
      // z=grille
      const int zi = Q->i/pq;
      const int zj = Q->j/pq;

      /* ici: zi<=zj car i<j */

      // x=colonne
      const int xi = Q->i%p;
      const int xj = Q->j%p;

      // y=ligne
      const int yi = (Q->i%pq)/p;
      const int yj = (Q->j%pq)/p;
      
      if(zj==zi) RET_a(((xj>=xi)&&(yj<=yi)) || ((xj<=xi)&&(yj>=yi)));
      /* ici: zj>zi */

      RET_a(((xj>=xi)&&(yj==yi)) || ((xj==xi)&&(yj>=yi)));
  }

  return 1;
}


int uno(query* const Q){
  switch(Q->code){

  case QUERY_END:
    return free_pos(Q);

  case QUERY_NAME:
    sprintf(Q->name,"(%i,%i)",(int)XPOS[Q->i],(int)YPOS[Q->i]);
    return 0;

  case QUERY_INIT:
    SET_n(Q->param[0]);
    const int p=Q->param[1];
    const int q=Q->param[2];
    if((p<=0)||(q<=0)) RET_n(0);
    ALLOC(XPOS,Q->n);
    ALLOC(YPOS,Q->n);
    int u;
    
    for(u=0;u<Q->n;u++){
      XPOS[u]=(int)(RAND01*p);
      YPOS[u]=(int)(RAND01*q);
    }
    
    XYtype=XY_USER; /* coordonnées fixées par l'utilisateur */
    InitXY(Q); // pour les options -xy noise/scale ... */
    return 0;
    
  case QUERY_ADJ:;
    RET_a((XPOS[Q->i]==XPOS[Q->j])||(YPOS[Q->i]==YPOS[Q->j]));
  }

  return 1;
}
  
  
int unok(query* const Q){
  if(Q->code==QUERY_INIT){
    SET_n(Q->param[0]);
    const int p=Q->param[1];
    const int q=Q->param[2];
    const int kp=Q->param[3];
    const int kq=Q->param[4];
    if((p<1)||(q<1)||(kp<1)||(kq<1)||(Q->n>imin(p*kp,q*kq))) RET_n(0);

    ALLOC(XPOS,Q->n);
    ALLOC(YPOS,Q->n);

    NALLOCMATZ(int,G,p,q,0); // G[i][j]=1 ssi on a mis un sommet en (i,j)
    NALLOCZ(int,A,p,0); // A[i]=nombre de sommets de G placés dans la ligne i
    NALLOCZ(int,B,q,0); // B[j]=nombre de sommets de G placés dans la colonne j
    int u,z,c,i,j,t;

    /* Principe: 

       Soit z le nombre de cases libres où l'on peut placer un point,
       c'est-à-dire telles que G[i][j]=0, A[i]<kp et B[j]<kq. On
       choisit une des z cases libres avec c=random()%z, puis on
       parcoure la grille, et quand la case libre numéro c arrive on y
       met un nouveau sommet (i,j), et on met à jour z. Pour la mise à
       jour on vérifie lorsque la ligne ou la colonne devient saturée
       en comptant les cases libres sur la croix de centre (i,j). La
       complexité en temps en environ npq/4 (en moyenne le point
       aléatoire se situe en p*q/2 et la somme (n-i)*p*q/2 donne du
       n*p*q/4).

       On optimise cet algorithme en commençant d'abord par un
       exécuter un algorithme à rejets, qui ne nécessite pas le
       parcours de toute la grille. Il est rapide au début puis va
       ralentir en fonction de la densité des cases libres de la
       grille et des rejets qui se produisent. On estime alors le
       temps dépensé par cet algorithme depuis le départ (t1) et celui
       restant si l'on continuait par l'algorithme classique (t2).
       Tant que t1<t2 on applique l'algorithme à rejets. Ensuite on
       change pour l'algorithme classique.
    */
    
    // partie commune aux deux algorithmes: ajoût du point(i,j) dans
    // la grille avec la mise à jour des positions XPOS, YPOS, de z,
    // des vecteurs A et B.  le point (i,j). I1 et I2 servent à la
    // mise à jour de l'estimation du temps pour l'algorithme à rejet
    // (non utilisé pour l'algorithme classique).
#define UPDATE(I1,I2)							\
    do{									\
      G[i][j]=1,z--;							\
      A[i]++; if(A[i]==kp){ I1; for(t=0;t<q;t++) z -= (B[t]<kq)&&(G[i][t]==0); } \
      B[j]++; if(B[j]==kq){ I2; for(t=0;t<p;t++) z -= (A[t]<kp)&&(G[t][j]==0); } \
      XPOS[u]=i,YPOS[u]=j;						\
    }while(0)
    
    u=0;           // u=nombre de sommets déjà tirés
    z=c=p*q;       // z=nombre de cases encore libres
    c/=4;          // c=temps moyen/sommet du temps de l'algo classique
    int t1=0;      // t1=temps dépensé par l'algo à rejet
    int t2=Q->n*c; // t2=temps estimé restant pour l'algo classique

    // Algorithme à rejets
    do{
      i=random()%kp;
      j=random()%kq;
      if((A[i]<kp)&&(B[j]<kq)&&(G[i][j]==0)){
	UPDATE(t1+=q,t1+=p);
	u++; if(u==Q->n) break;
	t2 -= c; // le temps de l'algo classique diminue
      }
    }while(t1++<t2);

    // Algorithme classique
    for(;u<Q->n;u++){
      c=random()%z; // c=numéro de case libre aléatoire uniforme parmi les libres
      for(i=0;i<p;i++){
	if(A[i]<kp) // sinon aucune case libre dans cette ligne
	  for(j=0;j<q;j++)
	    if((B[j]<kq)&&(G[i][j]==0)){
	      if(c==0){ // on a trouvé la case libre aléatoire
		UPDATE(,);
		i=j=p+q; // termine les deux for()
	      }else c--; // attend de trouver la case libre c
	    }
      }
    }
    
    FREE2(G,p);
    free(A);
    free(B);
    XYtype=XY_USER; /* coordonnées fixées par l'utilisateur */
    InitXY(Q); // pour les options -xy noise/scale ... */
    return 0;
  }
  return uno(Q);
}


int ngon(query* const Q){
/*
  Utilise i<j.

  codage de x:
  - bit-0 (x&1): AC connecté à A
  - bit-1 (x&2): BC connecté à B
  - bit-2 (x&4): carré
  - bit-3 (x&8): asymétrie
*/
  const int p=Q->param[0];
  int c=Q->param[1];
  const int x=Q->param[2];

  switch(Q->code){

  case QUERY_END:
    return free_pos(Q);
 
  case QUERY_INIT:
    if(x==-1){
      if(p<3) RET_n(0);
      SET_n(p); // n=p
    }
    if(x==-2){
      if(p<1) RET_n(0);
      SET_n(3*p); // n=3p
    }
    if(x>=0){
      // il faut c dans [0,p/2] et p>=1
      if((c<0)||(c>p/2)||(p<1)) RET_n(0);
      if((x<0)||(x>15)) RET_n(0);
      if(x&4) SET_n(4*p); // n=4p
      else SET_n(3*p); // n=3p
    }

    // donne des coordonnées aux points
    ALLOC(XPOS,Q->n);
    ALLOC(YPOS,Q->n);
    const double t=M_2PI/Q->n;
    double a=0;
    for(int i=0;i<Q->n;i++){
      XPOS[i]=cos(a);
      YPOS[i]=sin(a);
      a += t;
    }
    return 0;

  case QUERY_ADJ:;
    const int n=Q->n;
    int u=Q->i;
    int v=Q->j;
    if((u==0)&&(v==10))
      u=u+0;

    // teste le cycle
    if((v==u+1)||((v==n-1)&&(u==0))) RET_a(1);

    // triangulation à trois fans
    if(x==-1){
      if((v==c)||(u==n-c)) RET_a(1);
      if((u==0)&&(v==n-c)) RET_a(1);
      if((u==0)&&(c<v)&&(v<n-c)) RET_a(1);
      RET_a(0); // u,v non adjacents
    }

    const int p2=2*p; // p2=2p
    const int p3=3*p; // p3=3p
    int s=0; // vrai si asymétrique

    // diagonale du carré central
    if((u==0)&&(x&4)&&(v==p2)) RET_a(1);

    // on se ramène à u<v dans l'intervalle [0,p]
    for(;;){
      if((u<p)&&(v<=p)) break;                               // intervalle [0,p]
      if((p<=u)&&(u<p2)&&(v<=p2)){ u-=p; v-=p; s=1; break; } // intervalle [p,2p]
      if((p2<=u)&&(u<p3)&&(v<=p3)){ u-=p2; v-=p2; break; }   // intervalle [2p,3p]
      if(p3<=u){ u-=p3; v-=p3; s=1; break; }                 // intervalle [3p,n]
      if(u==0){ // cas particulier avec u=0
	if((x&4)&&(x>=0)){
	  if(v>=p3){ u=v-p3; v=p; s=1; break; }              // intervalle [3p,n]
	}else if(v>=p2){ u=v-p2; v=p; break; }               // intervalle [2p,n]
      }
      RET_a(0); // u,v non adjacents
    }
    // ici u<v dans [0,p]

    // triangulation récursive
    if(x==-2){
      int k=1;
      if(v==p){ // si v=p -> v = puissance de 2 >=v
	while(v>k) k<<=1;
	v=k;
      }
      k=1;
      const int w=v-u;
      while(w>k) k<<=1;
      if((w==k)&&(v%k==0)) RET_a(1);
      RET_a(0); // u,v non adjacents
    }

    if((u==0)&&(v==p)) RET_a(1); // grand triangle ou carré
    if((x&8)&&s) c=p-c; // asymétrie

    const int A=x&1; // vrai si étoile depuis A
    const int B=x&2; // vrai si étoile depuis B

    if(v<=c){ // si u,v dans [A,C]
      if(A) RET_a(u==0); // étoile depuis A
      RET_a(v==c);       // étoile depuis C
    }
    if(u>=c){ // si u,v dans [C,B]
      if(B) RET_a(v==p); // étoile depuis B
      RET_a(u==c);       // étoile depuis C
    }
    RET_a(0);
  }

  return 1; // fin anormale
}


int behrend(query* const Q){
/*
  Utilise i<j.
  Adjacence en O(log(p)).
*/
  const int p=Q->param[0];
  const int k=Q->param[1];
  int a,b,w,s,i;

  switch(Q->code){

  case QUERY_END:
    free(Q->wrap),Q->wrap=NULL;
    return 0;
    
  case QUERY_INIT:
    if((p<2)||(k<2)) RET_n(0);
    SET_n(p*k);
    /*
      On construit un ensemble X d'éléments de [0,p[ à partir de
      toutes les permutations P de [0,a[ interprétés comme entier en
      base b. Donc |X|=a! et |X|<p. Les entiers a et b doivent
      vérifier:
      
      1) a<b/k
      2) max{X}<p/(k-1) où max{X}=sum_{i=0}^{a-1} i*b^i
      3) a est maximum
      
      Principe pour trouver a et b: on part de a=2, et on essaye de
      voir s'il existe un b. Pour cela on choisit b=a*k+1, la plus
      petite valeur possible vérifiant 1), on calcule la somme
      (max{X}) et on vérifie qu'elle ne dépasse pas p/(k-1). Si la
      somme est correcte, on passe à a+1, sinon on a trouvé a et b.
     */

    const int s0=iceil(p,k-1); // NB: s<p/(k-1) => s<Ceil(p/(k-1))

    /* cherche a et b */
    i=a=1,s=0,b=k+1; // NB: b=a*k+1, w=b^i, s=sum i*b^i=max{X}
    while(s<s0)
      if(i==a) a++, b+=k, w=1, i=s=0;
      else s+=i*w, w*=b, i++;

    a--, b-=k; // ici a et b vérifient 1),2),3)

    /* crée l'ensemble X */
    NALLOCZ(int,P,a,_i); // P = permutation sur [0,a[
    ALLOC(Q->wrap,p); // Q->wrap = ensemble X de taille <= p

    s=1; // s=indice courant de Q->wrap, i=0 est réservé pour (a!)
    do{
      Q->wrap[s]=0, w=1;
      for(i=0;i<a;i++) Q->wrap[s]+=P[i]*w, w*=b;
      s++;
    }
    while(NextPermutation(P,a,NULL));
    free(P);
    REALLOC(Q->wrap,s); /* ici s=a!+1, s<=p */
    Q->wrap[0]=--s; // wrap[0]=s=|X|, wrap[1..s]=X

    /* trie X pour test d'adjacence plus rapide, en log(s) */
    QSORT(Q->wrap+1,s,fcmp_int); // ordre croissant
    return 0;

  case QUERY_ADJ:;
    int si=Q->i/p; // stable de i
    int sj=Q->j/p; // stable de j
    if(si==sj) RET_a(0); // pas adjacent
    a=Q->i%p; // rang de i dans son stable
    b=Q->j%p; // rang de j dans son stable
    if((si==0)&&(sj==k-1)){ SWAP(a,b,s); sj=1; } // NB: sj-si=1
    if(sj-si>1) RET_a(0); // pas adjacent
    // ici i-j ssi b=a+s mod p
    s=(b-a+p)%p;
    RET_a(bsearch(&s,Q->wrap+1,Q->wrap[0],sizeof(int),fcmp_int)!=NULL);
  }
  return 1;
}

int thetagone(query* const Q){
/*
  L'adjacence est en O(k*N).
*/
  const int p=Q->param[1];
  const int k=Q->param[2];

  switch(Q->code){

  case QUERY_END:
    return free_pos(Q);
      
  case QUERY_INIT:
    if(p<3) Q->param[1]=-1; /* p infini */
    if(k<1) RET_n(0);
    return gabriel(Q);

  case QUERY_ADJ:;
    /*
      Adjacence: pour tous les axes t, on cacule P_t(i,j)=distgone(),
      puis on détermine s'il existe un autre sommet z avec P_t(i,z)
      plus petit. Si c'est non (et que P_t(i,j) est finie), alors i
      est adjacent à j, sinon ils ne le sont pas (mais j peut être
      adjacent à i !).
    */
    int t,z;
    double d;

    const int n=Q->param[0];
    const double w=Q->dparam[0];

    for(t=0;t<k;t++){ /* pour tous les axes t */
      d=distgone(Q->i,Q->j,t,p,k,w); /* calcule P_t(i,j) */
      if(d<DBL_MAX){ /* si la distance d est finie */
	for(z=0;z<n;z++) /* pour tous les autres sommets z, même supprimés ! */
	  if((z!=Q->i)&&(distgone(Q->i,z,t,p,k,w)<d)) z=n; /* z plus proche ? */
	/* si oui, on arrête: P_t(i,j) contient z */
	if(z==n) RET_a(1); /* on a pas trouvé de sommet z plus proche que j */
      }
    }

    /*
      A priori ici il n'y a pas d'arête entre i et j. Il faut
      cependant tester aussi entre j et i car la distance P_t(i,j)
      n'est pas symétrique.
    */
    if(Q->i>Q->j) RET_a(0); /* pas d'arêtes i-j ni j-i */
    SWAP(Q->i,Q->j,t);
    return thetagone(Q);

  }
  return 1;
}
  

int squashed(query* const Q){
/*
  rep[i][0..k-1]: mot du sommet i sur {0,1,2} où la lettre '2'
  remplace le symbole '*'.

  Valeur expérimentale de connectivité pour p=1/3:

   n | 50 | 100 | 200 | 300 | 1000 | 10000 |
  ---o----o-----o-----o-----o------o-------o
   k |  4 |  15 |  18 |  19 |   22 |    28 |

  La probabilité que la distance de Hamming entre deux lettres soit
  non-nulle est q=2((1-q)/2)^2. C'est la probabilité d'avoir 0-1 ou
  1-0. Donc, la probabilité d'arête P_e que deux mots donnés soient
  voisins est

            P_e = 2kq (1-q)^{k-1} = 2q/(1-q) k*(1-q)^k

  C'est en effet deux fois celle de l'arc u->v, c'est-à-dire k*q (la
  distance est 1 à l'une des k positions) et (1-q) pour les k-1
  autres. Par analogie avec les graphes d'Erdos-Reny, lorsque P_e >
  ln(n)/n, alors le graphe est connexe. On a:

  P_e > ln(n)/n  <=>  k*(1-q)^k > ((1-q)/(2q))*ln(n)/n
                 <=>  k*ln(1-q)*e^(k*ln(1-q)) < ln(1-q)*((1-q)/(2q))*ln(n)/n
		 <=>  k*ln(1-q) < W( ln(1-q)*((1-q)/(2q))*ln(n)/n )
		 <=>  k < W(...)/ln(1-q)

  où W(x) est la fonction W de Lambert. Elle vérifie:

            W(x)=y  <=> x=y*e^y

  Elle est croissante pour x>=-1/e. W(-1/e)=-1, W(0)=0, W(e)=1. Si
  x>=e, alors A(x) + (1/2)*B(x) <= W(x) <= A(x) + (e/(e-1))*B(x) où
  A(x)=ln(x)-ln(ln(x)) et B(x)=ln(ln(x))/ln(x). En fait, pour x>=e,

            W(x) ~ ln(x) - ln(ln(x)) + ln(ln(x))/ln(x).

  Si 0<x<1, alors W(1)*x < W(x) < x avec
  W(1) = 0.5671432904097838729999686622... la constante Omega.
*/

  int k=Q->param[1];
  int u,t;

  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);

  case QUERY_NAME:
    if(k>NAMEMAX+1) Erreur(17);
    VIDE(Q->name);
    t=0;
    while(k--) 
      if(Q->rep[Q->i][k]==2) t += sprintf(Q->name+t,"%c",'*');
      else t += sprintf(Q->name+t,"%i",Q->rep[Q->i][k]);
    return 0;

  case QUERY_INIT:
    SET_n(Q->param[0]);
    double p=Q->dparam[0]; /* proba d'avoir '*' */
    if(k<0){ // valeur par défaut pour k
      if((p==0)||(p==1)) RET_n(0); // il faut 0<p<1
      double q=(1-p)/2; q=2*q*q; // on veut: q = 2*((1-p)/2)^2;
      double x=log(1-q)*((1-q)/(2*q))*log(Q->n)/Q->n;
      // on veut: k > W(x)/ln(1-q)
      //printf("x=%lf\n",x);
      x=fmax(-1/M_E,x); // on force x>=-1/e
      // calcule une borne sup sur W(x)
      // si x<e, W(x)>W(1)*x, sinon W(x) >= ln(x) - ln(ln(x)) + 0.5*ln(ln(x))/ln(x)
      if(x<M_E) x*=0.5671;
      else{
	x=log(x)-log(log(x))+0.5*log(log(x))/log(x);
	x /= log(1-q); // x=W(x)/ln(1-q)
      }
      Q->param[1]=k=floor(1+x); // pour avoir k > x
      //printf("k=%i p=%lf q=%lf\n",k,p,q);
    }
    if((k<1)||(k>=Q->n)||(Q->n<2)) RET_n(0);
    ALLOCMAT(Q->rep,Q->n,k);
    for(u=0;u<Q->n;u++) /* pour tous les sommets u du graphe */
      for(t=0;t<k;t++) /* pour toutes les lettres */
	if(RAND01<p) Q->rep[u][t]=2; /* lettre '*' avec proba p */
	else Q->rep[u][t]=(RAND01<0.5); /* sinon lettre 0 ou 1 */
    return 0;
    
  case QUERY_ADJ:;
    /* est-ce que la distance de Hamming entre rep[i] et rep[j] vaut 1 ? */
    Q->a=0;
    while(k--)
      if(Q->rep[Q->i][k]==1-Q->rep[Q->j][k]){
	Q->a++;
	if(Q->a>1) RET_a(0); /* pas adjacent */
      }
    return 0; /* ici Q->a=0 ou 1 */
  }

  return 1;
}


int udg(query* const Q){
/*
  NB: deux points peuvent avoir les même coordonnées et la norme n'est
  pas forcément symétrique en i,j.
*/
  switch(Q->code){

  case QUERY_END:
  case QUERY_INIT: 
    return gabriel(Q);

  case QUERY_ADJ:
    if(dist_ij(Q->i,Q->j)<=Q->dparam[0]) RET_a(1);
    if(Q->i>Q->j) RET_a(0);
    int t; SWAP(Q->i,Q->j,t);
    return udg(Q);
  }

  return 1;
}


int rng(query* const Q){
/*
  Adjacence en O(N).
*/
  switch(Q->code){

  case QUERY_END:
  case QUERY_INIT:
    return gabriel(Q);

  case QUERY_ADJ:;
    int k;
    const int n=Q->param[0];
    const double r=dist_ij(Q->i,Q->j); // !!! pas forcément symétrique en i,j

    for(k=0;k<n;k++) /* teste même les sommets supprimés */
      if(fmax(dist_ij(Q->i,k),dist_ij(Q->j,k))<r){
	if(Q->i>Q->j) RET_a(0); /* norme pas forcément symétrique */
	SWAP(Q->i,Q->j,k);
	return rng(Q);
      }
    
    RET_a(1);
  }

  return 1;
}


int nng(query* const Q){
/*
  Adjacence en O(N).
*/
  switch(Q->code){

  case QUERY_INIT: 
  case QUERY_END:
    return gabriel(Q);

  case QUERY_ADJ:;
    int k;
    const int n=Q->param[0];
    const double r=dist_ij(Q->i,Q->j); // !!! pas forcément symétrique en i,j

    for(k=0;k<n;k++) /* teste même les sommets supprimés */
      if((dist_ij(Q->i,k)<r)&&(k!=Q->i)) k=n; /* alors d(k,i)<d(i,j) */
    
    if(k==n) RET_a(1);

    /* Avant de dire que i a pour voisin j, il faut tester si j a pour
       voisin i car le test n'est pas symétrique. */

    if(Q->i>Q->j) RET_a(0);
    SWAP(Q->i,Q->j,k);
    return nng(Q);
  }
  
  return 1;
}


int hexagon(query* const Q){
/*
  Utilise i<j.
    
  On voit le graphe comme une grille de p+1 lignes de 2q+2 colonnes
  allant du coin en bas à droite (0,0) au coin (p,2q+1) (en notation
  (ligne,colonne)), et dans laquelle deux sommets ont été supprimés: le
  coin (0,2q+1) et le coin (p,2q+1) si p est impair, le coin (p,0) si
  p est pair. Les numéros sont consécutifs sur une ligne, de haut en
  bas.

  Ex:

  hexagon 3 2   hexagon 2 3

  o-o-o-o-o x   o-o-o-o-o-o-o x
  |   |   |     |   |   |   |
  o-o-o-o-o-o   o-o-o-o-o-o-o-o
    |   |   |     |   |   |   |
  o-o-o-o-o-o   x o-o-o-o-o-o-o
  |   |   |
  o-o-o-o-o x

*/
  int li,lj,ci,cj;

  const int p=Q->param[0];
  const int q=Q->param[1];
  const int t=(q<<1)+2; /* longueur d'une ligne */

  switch(Q->code){

  case QUERY_INIT:
    RET_n((p+1)*t-2);

  case QUERY_ADJ:;
    // pour éviter de modifier Q
    int i=Q->i;
    int j=Q->j;

    if(i>=t-1) i++; /* on insère virtuellement le coin (0,2q+1) */
    if(j>=t-1) j++;
    if((p&1)==0){ /* si p est impair, on a rien à faire */
      if(i>=p*t) i++; /* on insère virtuellement le coin (p,0) si p est pair */
      if(j>=p*t) j++;
    }

    /* on calcule les coordonnées de i et j placés sur cette
       grille (avec les coins manquant) */

    li=i/t;ci=i%t;
    lj=j/t;cj=j%t;
  
    /* utilise le fait que i<j: dans le dernier cas lj=li+1 */
    RET_a(((li==lj)&&(abs(ci-cj)==1)) || ((ci==cj)&&(lj==(li+1))&&((ci&1)==(li&1))));
  }

  return 1;
}


int whexagon(query* const Q){
/*
  Utilise i<j et hexagon(Q).
  Modifie Q->i et Q->j.
*/
  int li,ci,lj,cj;
  const int p=Q->param[0];
  const int q=Q->param[1];
  switch(Q->code){

  case QUERY_INIT:
    hexagon(Q);
    RET_n(Q->n+p*q);

  case QUERY_ADJ:;
    int t=Q->n-p*q; /* t=nombre de sommets de l'hexagone */

    /* teste si i et j sont dans l'hexagone */
    if((Q->i<t)&&(Q->j<t)) return hexagon(Q);

    /* teste si i et j sont hors l'hexagone */
    if((Q->i>=t)&&(Q->j>=t)) RET_a(0);

    /* on a i dans l'hexagone et j en dehors car i<j */
    lj=(Q->j-t)/q;cj=(Q->j-t)%q; /* j est le centre de l'hexagone (lj,cj) */
    t=(q<<1)+2; /* t=longueur d'une ligne de l'hexagone */
    
    /* on calcule les coordonnées de i dans l'hexagone */
    if(Q->i>=t-1) Q->i++; /* on corrige */
    if(((p&1)==0)&&(Q->i>=p*t)) Q->i++; /* on recorrige */
    li=Q->i/t;ci=Q->i%t; /* (li,ci) coordonnées de i */
    
    RET_a( ((li==lj)||(li==(lj+1))) && (abs(2*cj+1+(lj&1)-ci)<2) );
  }

  return 1;
}


int hanoi(query* const Q){
/*
  Modifie Q->i et Q-j.

  Adjacence: on écrit i (et j) en base b, mot de n lettres. i et j
  sont adjacents ssi i=Puv...v et j=Pvu...u où P est un préfixe
  commun, et u,v des lettres qui se suivent (modulo b).
*/
  int ri,rj,u,v,k;
  const int n=Q->param[0];
  const int b=Q->param[1];
  switch(Q->code){

  case QUERY_NAME:
    if(b<10) name_base(Q->name,Q->i,b,n,"","",1);
    else name_base(Q->name,Q->i,b,n,",","",1);
    return 0;

  case QUERY_INIT:
    if(b<2) RET_n(0);
    for(Q->n=k=1;k<=n;k++) Q->n *= b; /* #nombre de sommets = b^n */
    return 0;

  case QUERY_ADJ: /* on égraine les chiffres en base b */
    
    for(ri=rj=k=0;k<n;k++){
      if(ri!=rj) break; /* on s'arrête dès qu'on diffère, on garde k */
      ri=Q->i%b; Q->i/=b; /* ri=dernier chiffre, i=i sans dernier chiffre */
      rj=Q->j%b; Q->j/=b; /* rj=dernier chiffre, j=j sans dernier chiffre */
    }
    if((((ri+1)%b)!=rj)&&(((rj+1)%b)!=ri)) RET_a(0); /* alors pas voisin */

    u=ri; v=rj; /* ici u et v sont consécutifs (mod b) */
    for(;k<n;k++){
      ri=Q->i%b; Q->i/=b;
      rj=Q->j%b; Q->j/=b;
      if(ri!=v) RET_a(0); /* pas bon */
      if(rj!=u) RET_a(0); /* pas bon */
    }
    RET_a(1);
  }

  return 1;
}


int sierpinski(query* const Q){
/*
  On utilise rep[i][k] pour représenter le sommet i. C'est un mot d'au
  plus n lettres de [0,b[. On pose rep[i][n]=L où L est la longueur du
  mot.

  Le mot du sommet i représente la suite des cycles auxquels il
  appartient, sauf s'il est l'un des b sommets initiaux. Prennons
  b=3. Les 3 sommets du triangle sont 0,1,2. Si n>1, les autres
  sommets commenceront tous par 0. On a alors 3 autres triangles,
  numérotés 0,1,2, le triangle i à pour sommet le sommet i. Les 3
  sommets internes sont 00,01,02. Le sommet 00 partage les triangles 0
  et 1, 01 les triangles 1 et 2, 02 les triangles 2 et 0. Si n>2, tous
  les autres sommets commenceront par 00, 01 ou 02 suivant le
  sous-triangle auquel ils appartiennent. Dans le sous-triangle 0, les
  3 sommets internes seront 000, 001, 002. Etc.

  Ex: n=2 b=3         0 
                     /  \
                   00 -- 02
                  /  \  /  \
                 1 -- 01 -- 2

  Adjacence:

  CAS 1: extrémité (|i|=1 et |j|=n)
  Si i=x et n>1, alors il est voisin de:
  - 0 x^{n-2} x
  - 0 x^{n-2} (x-1)
  Si i=x et n=1, alors c'est le CAS 2.
  
  Soit P=le plus préfixe commun entre i et j, et k=|P|

  CAS 2: sommet du triangle le plus interne (|i|=|j|=n)
  Si i=Px, k=n-1, et n>0, alors il est voisin de:
  - P (x+1)
  - P (x-1)

  CAS 3: sommet entre deux triangles (1<|i|<n et |j|=n)
  Si i=Px, alors:

  CAS 3.1: i=P=Qx est préfixe de j (k=|i|).
  Alors il est voisin de:
  - P (x+1)^{n-k-1} x
  - P (x+1)^{n-k-1} (x+1)

  CAS 3.2: i=Px.
  Alors il est voisin de:
  - P (x+1) x^{n-p-2} x
  - P (x+1) x^{n-p-2} (x-1)

*/
  int k,t,r,x,c,li,lj;
  const int n=Q->param[0];
  const int b=Q->param[1];
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);
    
  case QUERY_NAME:
    if(b<10) name_vector(Q->name,Q->rep[Q->i],Q->rep[Q->i][n],"","",1,"%i");
    else name_vector(Q->name,Q->rep[Q->i],Q->rep[Q->i][n],",","()",1,"%i");
    return 0;

  case QUERY_INIT:
    if((b<3)||(n<1)) RET_n(0); /* graphe vide, graphe non défini */
    for(Q->n=b,k=2;k<=n;k++) Q->n=b*Q->n-b;
    ALLOCMAT(Q->rep,Q->n,n+1); /* rep[t][k]=k-ème lettre du sommet t */
    for(t=0;t<Q->n;t++){ /* calcule les noms */
      x=t;
      k=r=0;
      if(x<b) Q->rep[t][0]=x; /* un des b sommets du 1er cycle ? */
      else{
	x -= b;	
	for(;;){
	  Q->rep[t][k++]=r;
	  if(x<b) { Q->rep[t][k]=x; break; }
	  x -= b; r=x%b; x /= b;
	}
      }
      Q->rep[t][n]=k+1-(n==0); /* longueur du mot, corrige si n=0 */
    }
    return 0;

  case QUERY_ADJ:
    li=Q->rep[Q->i][n]; /* longueur de i */
    lj=Q->rep[Q->j][n]; /* longueur de j */
    /* propriété: si i<j, alors li<=lj */
    if(lj<n) RET_a(0);
    
    /* CAS 1 */
    if((li==1)&&(n>1)&&(Q->rep[Q->j][0]==0)){
      x=Q->rep[Q->i][0];
      for(c=t=1;t<=n-2;t++) if(Q->rep[Q->j][t]!=x) { c=0; break; }
      if(c){ /* c=vrai ssi j=0x^(n-2) */
	if(Q->rep[Q->j][n-1]==x) RET_a(1);
	if(Q->rep[Q->j][n-1]==((x-1+b)%b)) RET_a(1);
      }
    }

    /* calcule k=longueur du préfixe commun */
    for(k=0;(k<li)&&(k<lj)&&(k<n);k++)
      if(Q->rep[Q->j][k]!=Q->rep[Q->i][k]) break;

    /* CAS 2 */
    if((li==n)&&(k==n-1)){
      x=Q->rep[Q->i][k];
      if(Q->rep[Q->j][k]==((x+1)%b)) RET_a(1);
      if(Q->rep[Q->j][k]==((x-1+b)%b)) RET_a(1);
    }

    /* CAS 3 */
    if((li==1)||(li==n)) RET_a(0);
    x=Q->rep[Q->i][li-1];
    /* ici on a 1<|i|<n, |j|=n, et x=dernière lettre de i */

  /* CAS 3.1 */
  if(k==li){
    for(t=k;t<=n-2;t++) if(Q->rep[Q->j][t]!=((x+1)%b)) break;
    if(t>n-2){
      if(Q->rep[Q->j][n-1]==x) RET_a(1);
      if(Q->rep[Q->j][n-1]==((x+1)%b)) RET_a(1);
    }
  }
  
  /* CAS 3.2 */
  if((k==li-1)&&(Q->rep[Q->j][k]==((x+1)%b))){
    for(t=k+1;t<=n-2;t++) if(Q->rep[Q->j][t]!=x) break;
    if(t>n-2){
      if(Q->rep[Q->j][n-1]==x) RET_a(1);
      if(Q->rep[Q->j][n-1]==((x-1+b)%b)) RET_a(1);
    }
  }
  
  RET_a(0); /* si i>j, alors on ne fait rien */
  }

  return 1;
}


int banana(query* const Q){
/*
  Utilise i<j.
  Modifie Q->param[1].

   232323
   \|\|\|
    0 0 0
    | | |
    1 1 1
    \ | /
     n-1
*/
  const int n=Q->param[0];
  const int k=Q->param[1];
  switch(Q->code){

  case QUERY_INIT:
    Q->param[1]++;
    if((n<=0)||(k<=0)) RET_n(0);
    RET_n(n*(k+1)+1);
    
  case QUERY_ADJ:
    if(Q->j==Q->n-1) RET_a(Q->i%k==1);
    RET_a((Q->i/k==Q->j/k)&&(Q->i%k==0));
  }
  return 1;
}


int rpartite(query* const Q){
/*
  Utilise i<j.

  Les sommets sont numérotés consécutivement dans chacune des
  parts. Utilise le tableau Q->wrap en interne: wrap[k]=a_1+...+a_k
  est la k-ième somme partielle, avec wrap[0]=0. Donc les sommets de
  la part i sont numérotés de wrap[i-1] à wrap[i] (exclu).
*/
  int k,s,r=Q->param[0];
  switch(Q->code){

  case QUERY_INIT:
    if(r<=0) RET_n(0);
    free(Q->wrap);ALLOC(Q->wrap,r);
    Q->n=Q->wrap[0]=0;
    for(k=1;k<=r;k++){
      if(Q->param[k]<=0) RET_n(0);
      Q->n += Q->param[k];
      Q->wrap[k]=Q->n;
    }
    return 0;
    
  case QUERY_END:
    free(Q->wrap);Q->wrap=NULL;
    return 0;

  case QUERY_ADJ:
    /*
      Pour calculer l'adjacence, avec i<j, on calcule les numéros de
      la part de i et de j: adjacent ssi part(i)<>part(j). Pour cela,
      on fait une recherche dichotomique de la part de i (c'est-à-dire
      d'un k tq wrap[k]<=i<wrap[k+1]). La complexité est O(log(r)),
      r=#parts.
    */
    
    /* Cherche la part de i dans [s,r[: au départ c'est dans [0,r[ */
    s=0;
    k=r/2;
    while(s<r){
      if(Q->i<Q->wrap[k]){
	if(Q->j>=Q->wrap[k]) RET_a(1); /* i et j sont dans des parts <> */
	r=k; /* ici i<j<wrap[k]: on cherche dans [s,k[ */
	k=(s+k)>>1;
      }else{ /* ici wrap[k]<=i<j */
	if(Q->j<Q->wrap[k+1]) RET_a(0); /* i et j sont dans la part k */
	s=k; /* ici wrap[k]<=i<j: on cherche dans [k,r] */
	k=(k+r)>>1; 
      }
    }
    RET_a((Q->j>=Q->wrap[k+1])); /* ici i est dans la part k. On vérifie si j aussi */
  }
  
  return 1;
}


int aqua(query* const Q){
/*
  On se sert de rep[i][0..n], n=param[0].
*/
  int n=Q->param[0]; /* n=nombre de paramètres */
  int *C=Q->param+1; /* C=tableau de contraintes */
  int k,x,y,t;
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);

  case QUERY_NAME:
    name_vector(Q->name,Q->rep[Q->i],n,",","",1,"%i");
    return 0;
    
  case QUERY_INIT:;
    x=*C; /* s=Q->param[1]=premier terme=la somme */
    if(n<0) n=0;
    int *S=NextPart(NULL,n,x,C);
    Q->n=0; /* calcule une fois Q->n pour ALLOCMAT() */
    do Q->n++; while(NextPart(S,n,x,C));
    ALLOCMAT(Q->rep,Q->n,n); /* ici Q->n>0 */
    Q->n=0; /* calcule les sommets */
    do{
      for(k=0;k<n;k++) Q->rep[Q->n][k]=S[k];
      Q->n++;
    }while(NextPart(S,n,x,C));
    free(S);
    return 0;

  case QUERY_ADJ:
    if(Q->i==Q->j) RET_a(0);

    /* compte et mémorise les différences entre rep[i] et rep[j]: il y en a au moins 2 */
    for(k=0,x=y=-1;k<n;k++)
      if(Q->rep[Q->i][k]!=Q->rep[Q->j][k]){
	if((x>=0)&&(y>=0)) RET_a(0); /* si plus de 2 différences, alors pas d'arc */
	if(x<0) x=k; else if(y<0) y=k;
      }
    
    /* soit on a versé x vers y (soit le contraire) */
    /* k=quantité que l'on peut verser de x à y */
    for(t=1;t<=2;t++){ /* les deux cas sont identiques en inversant x et y */
      k=imin(C[y],Q->rep[Q->i][x]+Q->rep[Q->i][y]) - Q->rep[Q->i][y];
      if((Q->rep[Q->j][y]==Q->rep[Q->i][y]+k)&&
	 (Q->rep[Q->j][x]==Q->rep[Q->i][x]-k)) RET_a(1);
      if(t) SWAP(x,y,k); /* on échange une fois, la 1ère fois */
    }
    RET_a(0);    
  }

  return 1;
}


int flower_snark(query* const Q){
  int i=Q->i,j=Q->j; // pour ne pas modifier Q->i et Q->j
  const int u=(j>>2)-(i>>2); /* i<j, donc u>=0 */
  const int k=Q->param[0];
  switch(Q->code){

  case QUERY_NAME:
    sprintf(Q->name,"%c%i",((i&3)==0)? 'c' : 't'+(i&3),i>>2);
    return 0;
    
  case QUERY_INIT:
    RET_n(k<<2);

  case QUERY_ADJ:
    i &= 3;
    j &= 3;
    if(u==0) RET_a((i==0));
    if((u==1)&&(i==j)) RET_a((i>0));
    if(u!=k-1) RET_a(0);
    i*=j;
    RET_a(((i==1)||(i==6)));
  }

  return 1;
}


int gear(query* const Q){
/*
  Utilise i<j.
*/
  int n=Q->param[0];
  switch(Q->code){

  case QUERY_INIT:
    n<<=1;  /* double le paramètre n */
    Q->param[0]=n; Q->param[1]=Q->param[2]=2; Q->param[3]=0;
    RET_n(n+1);
    
  case QUERY_ADJ:
    if((Q->j==n)&&((Q->i&1)==0)) RET_a(1);
    if(Q->j<n) return cage(Q);
    RET_a(0);
  }
  
  return 1;
}


int arboricity(query* const Q){
/*
  Utilise Q->rep pour la représentation implicite. Q->rep[i][0..k-1]
  sont les k pères du sommet i, où k=Q->param[1]. Si rep[i][j]<0,
  c'est que le père j de i n'existe pas (i est une racine de la forêt
  j par exemple).
*/
  const int k=Q->param[1];
  int t,v;
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);

  case QUERY_INIT:
    /* calcule Q->n et Q->rep[i][0..k-1] */
    SET_n(Q->param[0]);
    if(k<1) RET_n(0);       /* il faut k>0 */
    ALLOCMAT(Q->rep,Q->n,k);/* rep=représentation finale */
    NALLOCZ(int,P,Q->n,_i); /* P=permutation aléatoire */
    NALLOC(int,T,Q->n);     /* T=arbre aléatoire */
    
    for(t=0;t<k;){ /* pour chaque arbre */
      Dyck(T,Q->n-1,1,DYCK_TREE); /* calcule un arbre aléatoire T */
      T[0]=0; /* racine = 0 plutôt que -1, car P[-1] n'existe pas ! */
      for(v=0;v<Q->n;v++) Q->rep[P[v]][t]=P[T[v]]; /* copie et permute le père */
      Q->rep[P[0]][t]=-1; /* père de la racine = -1 */
      if(++t<k) Permute(P,Q->n); /* si pas fini, on permute aléatoirement P */
    }
    
    free(T);
    free(P);
    return 0;
    
  case QUERY_ADJ:
    for(t=0;t<k;t++) /* j est un père de i ou le contraire */
      if((Q->rep[Q->i][t]==Q->j)||
	 (Q->rep[Q->j][t]==Q->i)) RET_a(1);
    RET_a(0);
  }

  return 1;
}


int kpage(query* const Q){
/*
  Modifie Q->param.

  Utilise rep[] pour la représentation implicite. Chaque sommet i
  possède 2k pères: rep[i][2p] et rep[i][2p+1] sont les 2 pères du
  sommet i de la page p (p-ème outerplanar), p=0...k-1. Pour le
  générer, on fait l'union de k outerplanars connexes enracinés et
  plan ("outer-plan" en fait).  Chacun est numérotés selon un parcours
  de la face extérieure avec une permutation circulaire aléatoire,
  sauf le premier. Pour l'adjacence on fait comme pour un graphe
  d'arboricité 2k. Le pagenumber d'un ne peut dépasser ceil(n/2), mais
  des valeurs de k plus grandes sont quand même autorisées.

  Le tirage de chaque outerplanar est construit à partir d'un arbre
  bicolorié (couleur 0 ou 1). Le 2e parent de u existe si u est
  colorié (disons à 1). Il vaut alors le prochain sommet non-related
  de u. On ne tient pas compte de la dernière branche, ces sommets
  n'ayant jamais de 2e parent. Ce tirage est uniforme sur les
  "outer-plan" connexes.
*/
  switch(Q->code){

  case QUERY_ADJ:
    return arboricity(Q);

  case QUERY_END:
    Q->param[1]<<=1; // rétablit les paramètres
    return free_rep(Q);

  case QUERY_INIT:;
    int u,p,q,t,z,k,c;
    SET_n(Q->param[0]);
    k=Q->param[1]; /* k=nombre de pages <= ceil{n/2} */
    if(k<1) RET_n(0); /* il faut k>0 */
    Q->param[1]<<=1; /* double le paramètre pour arboricity(Q) */
    ALLOCMAT(Q->rep,Q->n,Q->param[1]);
    NALLOC(int,T,Q->n); /* T=arbre aléatoire */
    /*
      On veut, pour chaque page p:
      rep[u][2p+0]=père1 de u
      rep[u][2p+1]=père2 de u
    */

    for(p=q=0;p<k;p++,q++){ /* pour chaque page p=0...k-1. Attention! q=2*p */
      Dyck(T,Q->n-1,1,DYCK_TREE); /* calcule un arbre DFS aléatoire T */
      
      /* c=permutation circulaire aléatoire pour les noms de sommets */
      if(p) c=random()%Q->n; else c=0; /* aucune permutation pour p=0 */

      /* calcule le père1 des sommets, le père dans T */
      Q->rep[c][q]=-1; /* la racine n'a pas de père */ 
      for(t=1;t<Q->n;t++) /* parcoure les sommets t de T selon le DFS, sauf la racine */
	Q->rep[(t+c)%Q->n][q]=(T[t]+c)%Q->n;

      /* calcule le père2 des sommets */

      /* Principe: le sommet courant est t. Chaque fois qu'on démarre
	 une nouvelle branche (t-T[t]>1), alors on parcoure les
	 sommets u allant de t-1 à T[t] non compris (la branche donc)
	 et décide pour chaque u de le connecter ou pas vers t (t qui
	 est donc le prochain non-related de u). NB: On ne parcoure
	 qu'au plus deux fois chacun des sommets. */

      q++; /* pour le père2 */
      Q->rep[c][q]=-1; /* la racine n'a pas de père */
      for(t=1;t<Q->n;t++){ /* parcoure les sommets t de T selon le DFS, sauf la racine */
	u=(t+c)%Q->n; /* u=sommet du graphe correspondant à t */
	Q->rep[u][q]=-1; /* par défaut, pas de père2 */
	if(t-T[t]>1){ /* ici t démarre une nouvelle branche */
	  z=t-1; /* z=dernier sommet de la branche précédante */
	  while(z!=T[t]){ /* z parcoure la branche précédante */
	    if(random()&1) Q->rep[(z+c)%Q->n][q]=u; /* ajoute un voisin vers u si coloré */
	    z=T[z]; /* z descend le long de la branche */
	  }
	}
      }
    }

    DEBUG(
	  for(p=0;p<k;p++){
	    ruling("―",10);printf("\n");
	    printf("page %i:\n",p);
	    for(u=0;u<Q->n;u++)
	      printf("%i: %i %i\n",u,Q->rep[u][2*p],Q->rep[u][2*p+1]);
	  }
	  ruling("―",10);printf("\n");
	  );

    free(T);
    return 0;

  }
  return 1;
}


int cactus(query* const Q){
/*
  rep[i][0] et rep[i][1] sont les 2 pères du sommet i. L'algorithme
  est le même que pour un outerplanar (soit un k-page avec k=1) sauf
  qu'un sommet u démarrant une nouvelle branche ne peut avoir qu'au
  plus seul voisin z dans la branche précédante. Mais cela ne suffit
  pas. Il faut de plus qu'aucun sommet entre z et T[u] n'ait déjà de
  voisin dans sa branche précédante.
*/
  switch(Q->code){

  case QUERY_END:
  case QUERY_ADJ:
    return arboricity(Q);

  case QUERY_INIT:;
    int u,z,t;
    SET_n(Q->param[0]);
    Q->param[1]=2; /* pour arboricity(Q) */
    ALLOCMAT(Q->rep,Q->n,2);
    NALLOC(int,T,Q->n); /* T=arbre DFS aléatoire */
    Dyck(T,Q->n-1,1,DYCK_TREE); /* calcule un arbre DFS aléatoire T */
    for(u=0;u<Q->n;u++) Q->rep[u][0]=T[u],Q->rep[u][1]=-1; /* père0 = T[u], père1 = -1 au départ */
    for(u=1;u<Q->n;u++){ /* parcoure les sommets u de T selon le DFS, sauf la racine */
      if(u-T[u]>1){ /* ici u démarre une nouvelle branche */
	z=u-1; /* z=dernier sommet de la branche précédante */
	t=0; /* t=sommet de la branche vers lequel u va pointer */
	while(z!=T[u]){
	  if((t==0)&&(random()&1)) t=z; /* u veut choisir z */
	  if((t)&&(Q->rep[z][1]>0)) t=0; /* pas bon, z dessous t à déjà un voisin */
	  z=T[z]; /* sommet suivant de la branche */
	}
	if(t) Q->rep[u][1]=t; /* c'est bon, u peut choisir t */
      }
    }
    free(T);
    return 0;
  }

  return 1;
}


int planar(query* const Q){
/*
  Modifie Q->param.

  Q->rep[u][0..1] = représentation implicite de u, ces deux pères.
  Attention ! Q->param[1] est modifié (et perdu) à l'initialisation à
  cause de l'utilisation de arboricity().
*/
  switch(Q->code){
    
  case QUERY_END:
  case QUERY_ADJ:
    return arboricity(Q);
      
  case QUERY_INIT:;
    int n=Q->param[0]; // n=nombre de faces
    int f=Q->param[1]; // f=taille maximum des faces internes
    int d=Q->param[2]; // d=degré des sommets internes
    const int w=(f<0); // w=vrai ssi face de taille au plus r

    f=abs(f);
    if((f<3)||(n<=0)) RET_n(0);

    /* le nombre maximum de sommets du graphe est a priori
       f+(n-1)*(f-2) = n(f-2)+2: 1 cycle de taille f au départ et on
       crée, pour chaque autre face, au plus f-2 nouveaux sommets */

    int k=n*(f-2)+2; /* k=nombre max de sommets */
    if(d<0) d=k; /* pas de contraintes de degré */
    NALLOCZ(int,A,k,-1); /* A[u]=1er voisin de u, aucun par défaut */
    NALLOCZ(int,B,k,-1); /* B[u]=2e voisin de u, aucun par défaut */
    NALLOC(int,C,k); /* C[i]=u=i-ème sommets de la face extérieure */
    NALLOC(int,T,k); /* pour mettre à jour C */
    NALLOCZ(int,D,k,2); /* D[u]=degré du sommet u, 2 par défaut */

    /* initialisation de la 1ère face, la face extérieure */
    int fr=(w)? 3+(random()%(f-2)) : f; /* fr=taille de la 1ère face */
    for(Q->n=0;Q->n<fr;Q->n++){
      A[Q->n]=(Q->n+1)%fr; // père vers le prochain du cycle
      C[Q->n]=Q->n; // face extérieure
    }
    int c=Q->n; // c=|C|=nombre de sommets de la face extérieure */
    int a,b,t,s,p,u,*Z;
    
    /* ajoute toutes n-1 autres faces */
    /* ici Q->n est le nouveau sommet courant */

    while((--n)>0){ // tant qu'il reste une face à faire
      fr=(w)? 3+(random()%(f-2)) : f; // fr=taille de la nouvelle face
      k=random()%c; // C[k]=un sommet de C au hasard
      A[Q->n]=C[k]; // 1er père de C[k]
      D[C[k]]++; // un voisin de plus pour C[k]

      /* on va tester les voisins valides (degré interne >= d) autour
	 de C[k], puis en choisir un au hasard. Enfin on le connectera
	 par un chemin jusqu'au nouveau sommet courant, Q->n, et
	 mettra à jour la face extérieure. */

      p=imin(fr-2,c/2); // p=nombre maximum de sommets de part et
                        // d'autre de C[k] dont il faut tester le degré
      
      /* teste les successeurs de C[k] sur C: C[k+1]...C[k+p] */
      for(a=1;a<=p;a++)	if(D[C[(k+a)%c]]<d) break;
      if(a>p) a--;
      // ici on peut se connecter à n'importe quel sommet entre C[k+1]..C[k+a] 

      /* teste les prédécesseurs de C[k] sur C: C[k-1]...C[k-p] */
      for(b=1;b<=p;b++)	if(D[C[(k+c-b)%c]]<d) break;
      if(b>p) b--;
      // ici on peut se connecter à n'importe quel sommet entre C[k-1]..C[k-b] 

      t=random()%(a+b); // choisit 1 indice parmi [0,a+b[
      if(t<a) s=1; else{ s=-1; t-=a; } 
      p=fr-t-3; // p=nombre de sommets du chemin, Q->n non compris, p=0 possible

      // t random dans [0,a[ si s>0
      // t random dans [0,b[ si s<0
      // s>0 => chemin = C[k]-N-(N+1)- ... - (N+p)-C[k+t+1]
      // s<0 => chemin = C[k-t-1]-(N+p)- ... - (N+1)-N-C[k]
      
      b=(k+s*(t+1)+c)%c; // C[b]=C[k+t+1] ou C[k-t-1]=sommet à connecter
      D[C[b]]++; // un voisin de plus pour C[b]
      for(u=Q->n;u<Q->n+p;u++) B[u]=u+1; // adjacence du chemin
      B[u]=C[b]; // le dernier sommet du chemin pointe sur C[b], u=Q->n possible

      /* mise à jour de C: on supprime d'abord les sommets de C qui
	 vont passer à l'intérieur, puis on calcule dans T la nouvelle
	 face extérieur (en sautant les sommets effacés de C):
	 - si s>0, il faut supprimer C[k+1]...C[k+t]
	 - si s<0, il faut supprimer C[k-t]...C[k-1]
      */
      for(u=1,k+=c;u<=t;u++) C[(k+u*s)%c]=-1; // supprime des sommets de C

      /* créer la nouvelle face extérieure T, a=|T| */
      for(u=a=0;u<=p;u++) T[a++]=Q->n++; // commence par les p+1 nouveaux sommets
      for(u=0,b+=c;u<c;u++) // le reste de l'ancienne face extérieure
	if(C[(b+u*s)%c]>=0) T[a++]=C[(b+u*s)%c];

      SWAP(C,T,Z); // échange C et T
      c=a; // nouvelle taille de C
    }

    free(T);
    free(C);
    free(D);

    /* recopie A,B dans Q->rep */
    Q->param[1]=2;     // pour utiliser l'arboricité 2
    ALLOCMAT(Q->rep,Q->n,2);
    for(u=0;u<Q->n;u++){
      Q->rep[u][0]=A[u];
      Q->rep[u][1]=B[u];
    }
    free(A);
    free(B);
    return 0;
  }

  return 1;
}


int hyperbolic(query* const Q){
/*
  Modifie Q->param.

  Q->rep[u][0..2] = représentation implicite de u, ces au plus trois
  pères. Le 3e père ne peut se produire que si p=3. Attention !
  Q->param[1] est modifié (et perdu) à l'initialisation à cause de
  l'utilisation de la routine arboricity() pour le test d'adjacence.

  rep[u][0]=successeur(u) sur le dernier niveau (cycle)
  rep[u][1]=parent(u) vers le niveau précédant (-1 si non défini)
  rep[u][2]=parent(u) vers le succ du niveau précédant (-1 si non défini)

*/
  switch(Q->code){
    
  case QUERY_END:
  case QUERY_ADJ:
    return arboricity(Q);
      
  case QUERY_INIT:;
    int p=Q->param[0]; // p=taille des faces
    int k=Q->param[1]; // k=degré des sommets
    int h=Q->param[2]; // h=nombre de niveaux
    if((p<3)||(k<2)||(h<1)) RET_n(0);
    if(k==2) h=1;

    int t,c,u,i,m;
    int n2=p; // =nombre de sommets de degré 2
    int n3=0; // =nombre de sommets de degré 3
    int n4=0; // =nombre de sommets de degré 4
    Q->n=p;   // =nombre total de sommets au départ

    // calcule le nombre de sommets: pour chaque niveau, on compte le
    // nombre de sommets de degré 2, 3 ou 4. En général, les sommets
    // de degré 3 sont engendrés par chaque sommet du niveau
    // précédent. Et ceux de degré 2 sont engendrés par les chemins
    // entre les ces sommets de degré 3 consécutifs nouvellement
    // crées. Il faut ensuite tenir compte des cas particuliers. Pour
    // p=3, il faut fusionner le dernier et le premier sommets de
    // degré 3. Apparaît alors un sommet de degré 4. Pour k=3, il faut
    // tenir compte des sommets de degré 3 de la couche précédante qui
    // doivent être sautés. Le cas p=k=3 doit être géré différemment
    // encore.

    for(c=1;c<h;c++){ // pour chaque nouveau niveau
      t=n3+n2+n4; // total de la couche précédante
      u=n3;
      PRINT(n2);
      PRINT(n3);
      PRINT(n4);
      PRINT(t);
      printf("\n");
      if(t==0){ h=c; break; } // il n'y a plus de sommets
      n3=(k-2)*n2+(k-3)*n3+(k-4)*n4;
      if(p==3){ n3-=2*t,n4=t,n2=0; if(k==3) n3=1,n4=0; h=2; }
      else{ n2=n3*(p-3) - t; if(k==3) n2-=u; }
      Q->n += n2+n3+n4;
    }
    PRINT(h);
    PRINT(Q->n);

    // pour utiliser l'arboricité
    Q->param[1]=2+((p==3)&&(k>3));
    if(h==1) Q->param[1]=1; // c'est un cycle
    
    ALLOCMAT(Q->rep,Q->n,3);
    for(u=0;u<Q->n;u++) Q->rep[u][1]=Q->rep[u][2]=-1; // aucun parent par défaut

    // cas très particulier
    if((p==3)&&(k==3)){ // cycle ou K_4
      Q->rep[0][0]=1; Q->rep[1][0]=2; Q->rep[2][0]=0;
      if(h>1){
	Q->rep[3][0]=-1;
	Q->rep[0][1]=Q->rep[1][1]=Q->rep[2][1]=3;
      }
      return 0;
    } // ici p>3 ou k>3
    
    int f=0; // f=1er sommet du nouveau niveau
    int d=0; // d=1er sommet du niveau précédent
    t=p;     // t=nouveau sommet
    
    for(c=0;c<h;c++){ // pour chaque nouveau niveau
      for(u=d;u<f;u++){ // parcoure les sommets u du niveau précédent
	// on détermine le parent du nouveau sommet t (par défaut il
	// n'en a pas), le cycle est déterminé dans un 2e temps par un
	// parcours complet du niveau (=cycle) crée
	m=k-4+(Q->rep[u][1]<0)+(Q->rep[u][2]<0); // m=k-2 ou k-3 ou k-4 suivant deg(u)
	for(i=0;i<m;i++){ // pour chaque nouveau sommet de degré 3
	  Q->rep[t][1]=u; // t sommet de deg 3 et de parent u
	  t += p-2; // prochain sommet de deg 3 (saute les sommets de deg 2)
	}
	if(m) t--;
      }
      // parcoure le nouveau niveau pour créer le cycle, NB: ici u=f
      for(;u<t;u++) Q->rep[u][0]=u+1;
      Q->rep[u-1][0]=f; // termine le cycle (évite des calculs de modulo)
      d=f; f=t; // met à jour (d,f) <- (f,t)
    }
    
    return 0;
  }

  return 1;
}


int NextDyck(int *X,int n){
/*
  Calcule le prochain mot de Dyck de longueur 2n contenu dans X (qui
  doit donc être un tableau de taille au moins 2n). Renvoie 1 ssi le
  dernier mot de Dyck a été atteint. Si n<0, alors X est initialisé au
  mot X=(10)^n.  Les mots sont énumérés selon le nombre croissant de 1
  consécutifs à gauche.  L'algorithme est celui de
  https://github.com/cassioneri/Dyck

  Pour n=4:

  10101010 10101100 10110010 10110100 10111000 11001010 11001100
  11010010 11010100 11011000 11100010 11100100 11101000 11110000

*/
  const int m=2*n-1;
  int y=0;
  int x=0;
  int i;

  if(n<0){
    x=1;
    y=-(n<<1);
    for(i=0;i<y;i++){ X[i]=x; x=1-x; }
    return (y==2);
  }

  for(i=m;i>0;i--)
    if(X[i]){
      if(X[i-1]) x++;
      else{
	X[i-1]=1;
	X[i]=0;
	for(y=y-x;y!=0;y--) X[++i]=0;
	while(i<m){
	  X[++i]=1;
	  X[++i]=0;
	}
	return 0;
      }
    }else y++;

  return 1; /* dernier mot atteint */
}


int flip(query* const Q){
/*
  Q->rep[i] = mot de Dyck = [ 1 0 1 1 0 0 ]. Ici n est le nombre de un
  du mot = Q->param[0]-2 = 3.  Le mot représente un arbre binaire
  complet (chaque noeud interne à deux fils exactement). On obtient le
  codage de l'arbre en mot de Dyck par un parcours DFS et en écrivant
  1 si l'arête parcouru mène à un fils gauche et 0 si elle mène à un
  fils droit. On écrit rien lorsqu'on parcoure les arêtes vers le
  père.

  Rotation sur le noeud x:

               B   C                      A   B
                \ /                        \ /
             A   o           ->             o   C
              \ /                            \ /
               x                              x

  U=[... 1 A 0 1 B 0 C ...]  ->  V=[... 1 1 A 0 B 0 C ...]

  Algorithme d'adjacence entre les mots U et V:

  1. On cherche le plus grand suffixe commun. Soit i la position tel
     que U[i]<>V[i] et U[i+1...2n-1] = V[i+1...2n-1]. Pour cela on
     remonte de 2n-1 jusqu'à la première différence.

  2. On échange éventuellement U et V de sorte que U[i]=1 et V[i]=0.

  3. On essaye de lire [... 1 A 0 ...] dans V à partir de i et dans U
     à partir de i-1, toujours selon les indices décroissant. En même
     temps qu'on vérifie que V[i]=U[i-1], on calcule la hauteur h avec
     +1 si V[i]=0 et -1 si V[i]=1 (car on lit le mot à l'envers). On
     s'arête dès que h<0. Soit i l'indice dans V où h<0 pour la
     première fois.

  4. On vérifie alors que V[j-1]=1, puis que U[0..j-2]=V[0..j-2]. On
     conclut que U et V sont adjacents.

  https://fr.wikipedia.org/wiki/Nombre_de_Catalan#Chemins_sous-diagonaux_dans_le_carr.C3.A9
  https://en.wikipedia.org/wiki/Tree_rotation
*/
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);
    
  case QUERY_NAME:
    name_vector(Q->name,Q->rep[Q->i],2*(Q->param[0]-2),"","",1,"%i");
    return 0;
  
  case QUERY_INIT:;
    int n,t,u;
    /* calcule Q->n = Catalan(n) */
    n=Q->param[0]-2;
    if(n<=0) RET_n(0);
    t=(n<<1); /* t=2n */
    SET_n(Binom(t,n)/(n+1));
    ALLOCMAT(Q->rep,Q->n,t);
    NALLOC(int,X,t);
    NextDyck(X,-n); /* initialise le 1er mot */
    t *= sizeof(int);
    for(u=0;u<Q->n;u++){
      bcopy(X,Q->rep[u],t); /* copie le mot de Dyck courrant vers Q->rep[u] */
      NextDyck(X,n); /* calcule le mot suivant */
    }
    free(X);
    return 0;
  
  case QUERY_ADJ:;
    int* U=Q->rep[Q->i];
    int* V=Q->rep[Q->j];
    int* W;
    int h=1; /* hauteur de [ ... 1 A 0 ... ] */
    int k=2*(Q->param[0]-2)-1; /* k=dernière position de U ou V */

    /* calcule suffixe */
    while(U[k]==V[k]) k--;
    if(V[k]) SWAP(U,V,W); /* échange U et V */

    /* k=position du 0 dans V */
    while((V[k]==U[k-1])&&(h>0)) h += 1-2*V[--k];

    /* problème ? */
    if((V[--k]==0)||(h>0)) RET_a(0);

    /* préfixe */
    while((k>=0)&&(U[k]==V[k])) k--;
    RET_a((k<0));
  }

  return 1;
}


int linegraph(query* const Q){
/*
  Chaque sommet i possède 2 couleurs prises dans 1...k. Les sommets i
  et j sont adjacents si une couleur de l'un est une couleur de
  l'autre.  Utilise Q->rep[u][0..1] pour la représentation des
  couleurs du sommet u.
*/
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);

  case QUERY_INIT:;
    int u,k;
    SET_n(Q->param[0]);
    k=Q->param[1]; if(k<=0) RET_n(0);
    ALLOCMAT(Q->rep,Q->n,2);
    for(u=0;u<Q->n;u++){
      Q->rep[u][0]=random()%k;
      Q->rep[u][1]=random()%k;
    }
    return 0;
    
  case QUERY_ADJ:
    RET_a(
	    ((Q->rep[Q->i][0]==Q->rep[Q->j][0])||(Q->rep[Q->i][0]==Q->rep[Q->j][1])||
	     (Q->rep[Q->i][1]==Q->rep[Q->j][0])||(Q->rep[Q->i][1]==Q->rep[Q->j][1]))
	    );
  }
  
  return 1;
}


int ringarytree(query* const Q){
/*
  Le sommet i est un chemin P(i)=x_1,x_2,... allant de la racine (=0)
  à i, chaque lettre x_t est le numéro du fils, numéro dans [0,r[ pour
  la racine et dans [0,k[ pour les autres noeuds internes. P(0)={} est
  vide.

  Pour que i et j soient voisins, avec i<j, il faut:

  - si p=0 (seulement connexion dans l'arbre): P(j)=P(i),x contient
    une seule lettre supplémentaire

  - si p=1 (p=0 et chemin entre noeuds de même niveau):
    P(i)=C,x_1,...,x_k
    P(j)=C,y_1,...,y_k
    y_1=1+x_1
    y_t=0 et x_t=k-1 pour tout t>1

  - si p=2 (p=1 et cycle entre noeuds de même niveau):
    x_1=0 et y_1=r-1 ou k-1 (suivant si C={} ou pas)
    x_t=0 et y_t=k-1 pour tout t>1

    Principe: on calcule P(i) et P(j) en parallèle, lettre par lettre.
*/
  int h=Q->param[0];
  const int k=Q->param[1];
  const int r=Q->param[2];

  switch(Q->code){

  case QUERY_NAME:
    if(Q->i==0) strcpy(Q->name,"ε"); // pour la racine
    else{
      int d=r;    /* d=nombre de fils de la racine de T, d=r puis d=k */
      int t=Q->n; /* t=taille des fils */
      int p=0;    /* nombre de caractères écrit dans Q->name */
      const int v=(imax(k,r)>10); /* vrai ssi il faut une virgule */
      VIDE(Q->name);
      while(Q->i>0){
	Q->i--;t--;t/=d;
	p+=sprintf(Q->name+p,"%i",Q->i/t); /* écrit le numéro du fils f=i/t*/
	if(p>NAMEMAX) Erreur(17);
	Q->i%=t;
	d=k; /* maintenant d=k */
	if((Q->i>0)&&(v)) p+=sprintf(Q->name+p,","); /* ajoute une "," */
      }
    }
    return 0;
    
  case QUERY_INIT: /* calcule Q->n */
    if((h<0)||(k<0)||(r<0)){ Q->n=0; return 0; }
    if((h==0)||(r==0)){ Q->n=1; return 0; }
    if(k==0) Q->param[0]=h=1; /* si k=0 et h>0, alors la hauteur est 1 */
    if(k<=1) Q->n=1+r*h;
    else{ /* ici k>1, h>0 et r>0 */
      int t; /* taille sous-arbre = 1+r+r^2+...+r^{h-1} = (r^h-1)/(k-1) */
      for(t=0,Q->n=1;t<h;t++) Q->n *= k; /* après cette boucle, Q->n=r^h */
      Q->n=1+r*(Q->n-1)/(k-1); /* Q->n=racine + r x (taille sous-arbre) */
    }
    return 0;

  case QUERY_ADJ:;
    /* calcule le préfixe commun de P(i) et P(j) */
    const int p=Q->param[3]; /* p=0,1,2 */
    int x=Q->i,y=Q->j;/* copies de i et j, NB: x<y */
    int t=Q->n; /* t=taille de T, l'arbre où x et y sont */
    int d=r;    /* d=nombre de fils de la racine de T, d=r puis d=k */
    int fx,fy;
    
    /* ici x et y sont dans un arbre T de taille t ayant d fils. La
       taille des fils de T est t'=(t-1)/d. Pour trouver le fils fx de
       T contenant x il suffit de faire (x-1)/t'. Le suffixe de x dans
       ce nouveau sous-arbre est alors (x-1)%t'. */
    
    fx=fy=0;
    
    /* tant que x et y sont tout deux dans T (ils viennent du même
       fils fx=fy, et aucun d'eux n'est la racine de T: on calcule
       alors leurs fils fx et fy et met à jour la nouvelle taille de T
       ainsi que le suffixe de x et y dans ce nouveau T. */
    
    while((fx==fy)&&(x>0)&&(y>0)){
      x--,y--,t--; /* on enlève la racine de T */
      t/=d; /* t=taille des fils de T */
      fx=x/t,fy=y/t; /* fils des sous-arbres de x et de y */
      x%=t,y%=t; /* suppression du préfixe de x et de y */
      d=k; /* maintenant d=k */
    }
    
    /*
      Ici la situation est la suivante: x et y sont dans des
      sous-arbres isomorphes à T de taille t, dans les sous-arbres des
      fils fx et fy. Chacun de ces sous-arbres est à d fils. Si fx=fy
      c'est que x est la racine (x=0) car y=0 est impossible puisque
      x<y.
	
                       o
                   fx / \ fy
                     o   o
                    / \ / \
                     x   y

    */

    /* fx=fy: x et y sont dans le même arbre T */
    if(fx==fy) RET_a(((y-1)%((t-1)/d)==0)); /* y fils de x ? */
    if(p<=0) RET_a(0);
    
    /* fx<>fy: x et y sont dans des sous-arbres différents */
    
    /* si 1er et dernier sous-arbres (voisins dans le cycle), on échange x et y */
    if((p==2)&&(fx==0)){
      int b=(t==(Q->n-1)/r)? r:k; /* fx,fy fils de niveau 1 ou > 1 ? */
      if(fy==b-1){ SWAP(x,y,b); fy=1; }
    }
    
    if(fy-fx>1) RET_a(0); /* sous-arbres qui ne sont pas voisins */
    
    /* x et y sont dans des sous-arbres voisins, chacun de taille t */
    
    for(;;){
      if((x==0)&&(y==0)) RET_a(1); /* racines voisines */
      if((x==0)||(y==0)) RET_a(0); /* pas même niveau */
      x--,y--,t--,t/=k,fx=x/t,fy=y/t,x%=t,y%=t; /* met à jour x,y,t,fx,fy */
      if((fx!=0)&&(fy!=k-1)) RET_a(0); /* il faut x fils 0 et y fils d-1 */
    }
    return 0;
  }

  return 1;
}


int rarytree(query* const Q){
/*
  Utilise Q->rep. Les sommets sont numérotés selon un DFS modifié: on
  pose les fils avant la récursivité (voir Dyck()).
*/
  switch(Q->code){

  case QUERY_END:
  case QUERY_ADJ:
    return arboricity(Q);

  case QUERY_INIT:;
    int n=Q->param[0]; if(n<=0) RET_n(0);
    int b=Q->param[1]; if(b<2) RET_n(0);
    int z=Q->param[2]; if((z!=0)&&(z!=1)) RET_n(0); /* z=0 ou 1 */
    int *B=Dyck(NULL,n,b-1,DYCK_KTREE); /* B=arbre b-aire aléatoire avec n noeuds internes */
    SET_n(b*n+1+z);
    ALLOCMAT(Q->rep,Q->n,1); /* représentation implicite */
    for(b=0;b<Q->n-z;b++) Q->rep[b][0]=B[b]; /* copie l'arbre B dans Q->rep, |B|=N-z */
    free(B);
    if(z) Q->rep[Q->n-1][0]=0; /* le dernier sommet pointe vers la racine 0 */
    Q->param[1]=1; /* pour test d'aboricité k=1 */
    return 0;
  }
  
  return 1;
}


int ktree(query* const Q){
/*
  Q->rep[i][0...k-1] sont les k pères du sommet i. Cette fonction
  utilise un 3e paramètre caché, Q->param[2] qui vaut 0, 1 ou 2
  suivant s'il faut générer un arbre aléatoire (ktree), un chemin
  (kpath) ou un étoile (kstar).
*/
  switch(Q->code){

  case QUERY_END:
  case QUERY_ADJ:
    return arboricity(Q);

  case QUERY_INIT:; /* calcule Q->n et Q->rep[i][0..k-1] */
    const int n=Q->param[0];
    const int k=Q->param[1];
    if((k<0)||(n<=k)) RET_n(0);
    SET_n(n);

    int t,p,w,x,y,*T;
    if(Q->param[2]==0) T=Dyck(NULL,Q->n-k-1,1,DYCK_TREE); /* arbre de n-k noeuds */
    if(Q->param[2]==1) ALLOCZ(T,Q->n-k,_i-1); /* chemin de n-k noeuds */
    if(Q->param[2]==2) ALLOCZ(T,Q->n-k,0); /* une étoile à n-k noeuds */
    ALLOCMAT(Q->rep,Q->n,k); /* représentation implicite */

    /* Chacun des k+1 sommets de la racine (numéros de 0 à k) ont pour
       pères tous les autres sommets (=> clique) */

    for(t=0;t<=k;t++) /* pour les k+1 sommets */
      for(p=0;p<k;p++) /* pour les k pères 0..k-1 */
	Q->rep[t][p]=(t+p+1)%(k+1); /* il faut sauter t */

    /* On utilise le fait que les noeuds de T forment un DFS. En
       traitant les sommets dans l'ordre on est sûr que le père est
       déjà traité */
 
    for(t=k+1;t<Q->n;t++){ /* on démarre à k+1, les k+1 sommets de la racine sont
			      déjà traités, tout sommet a donc une racine */
      p=T[t-k]; /* p=noeud père du sommet t du graphe */
      w=random()%(k+1); /* indice d'un des sommets du noeud père qui ne sera pas choisi */
      p += k; /* p=nom réel du père dans le graphe */
      Q->rep[t][0]=p; /* remplit avec le père lui-même */
      x=0; /* x=indice des pères pour Q->rep[p], x=0..k-1 */
      y=(w>0); /* y=prochain indice des pères pour Q->rep[t]. Si w=0 on
		  saute le père */
      while(x<k){
	Q->rep[t][y++]=Q->rep[p][x++];
	if(w==x) y--;
      }
    }
    
    free(T);  /* libère l'arbre */
    return 0;
  }
  
  return 1;
}


int apollonian(query* const Q){
/*
  Utilise Q->rep pour la représentation implicite, Q->rep[i][0,1,2] sont les
  3 pères du sommet i. On pourrait factoriser avec polygon().
*/
  switch(Q->code){

  case QUERY_END:
  case QUERY_ADJ:
    return arboricity(Q);

  case QUERY_INIT: /* calcule Q->n et Q->rep[i][0..2] */
    Q->n=Q->param[0]; if(Q->n<4) RET_n(0); /* il faut au moins 4 sommets */
    Q->param[1]=3; /* pour test d'arboricité, 3 pères */
    const int n=Q->n-3; /* n=nombre de sommets internes */
    const int m=3*n+1; /* nombre de sommets de l'arbre ternaire */
    int *P=Dyck(NULL,n,2,DYCK_KTREE); /* arbre ternaire à n sommets internes */

    /*
      Principe de la construction. On part d'un arbre ternaire à n
      sommets internes (dont la racine), comme dans rarytree(). La
      racine correspond à un K_4 dont le centre est le sommet 3. Puis
      le k-ième noeud interne de l'arbre (donc qui n'est pas une
      feuille) correspond à un nouveau K_4 dont le centre est un
      nouveau sommet numéroté k et qui est connecté à un triangle
      parent. Il y en a trois possibles suivant que le numéro du fils
      où l'on est.

                             0          3    (=triangle:012, centre:3)
                            /|\        /|\
                           1 2 3      4 . .  (=triangle:301, centre:4)
                          /|\        /|\
                         4 5 6      . 5 .    (=triangle:401, centre:5)
                          /|\        /|\
                         7 8 9      . . .

	P = [-,0,0,0,1,1,1,5,5,5] (= père dans l'arbre)
	C = [3,4,-,-,-,5,-,-,-,-] (= centre des triangles des sommets internes de l'arbre)
	
    */

    ALLOCMAT(Q->rep,Q->n,3); /* représentation implicite */
    int u,p,c;

    /* calcule C[u]=numéro du centre du triangle ou 0 si feuille, pour
       tout noeud u de l'arbre */

    NALLOC(int,C,m);
    for(u=1;u<m;u++) C[u]=0,C[P[u]]=1; /* par défaut u est une feuille, et son père non */
    C[0]=c=3; /* racine = centre 3 */
    for(u=1;u<m;u++) if(C[u]) C[u]=++c; /* met le numéro de centre */
    
    for(u=0;u<4;u++)   /* pour le premier tirangle */
      for(c=0;c<3;c++) /* pour les 3 fils de chaque sommet du K_4 */
	Q->rep[u][c]=(u+c+1)%4;

    /* on calcule le triangle Q->rep[c], pour chaque centre c=3...N.
       Q->rep[c][0..2] représente le triangle dont le centre est c */

    for(u=1;u<m;u++){ /* on parcoure les noeuds de l'arbre */
      c=C[u];
      if(c){ /* si u est un noeud interne */
	p=C[P[u]]; /* p=centre du père de u dans l'arbre */
	Q->rep[c][0]=p; /* un sommet du triangle est le centre */
	Q->rep[c][1]=Q->rep[p][u%3]; /* on en prend deux autres parmis le triangle du père */
	Q->rep[c][2]=Q->rep[p][(u+1)%3]; /* en fonction du numéro du fils -> u%3 */
      }
    }

    free(C);
    free(P);
    return 0;
  }
  
  return 1;
}


int polygon(query* const Q){
/*
  Utilise Q->rep pour la représentation implicite. Chaque sommet u
  possède 2 pères: Q->rep[u][0] et Q->rep[u][1]. Utilise Q->param[1]=2 pour
  l'adjacence: arboricity(i,j).

  La construction est similaire à apollonian(), avec un arbre binaire
  au lieu de ternaire, les K_3 remplaçant les K_4. On pourrait
  factoriser polygon() et apollonian() en passant en paramètre un
  booléen, via Q->param[2] par exemple: 0 pour polygon() et 1 pour
  apollonian(). Pas clair qu'on puisse encore généraliser cette
  construction où apollonian() et polygon() seraient des instances.
*/
  switch(Q->code){

  case QUERY_END:
  case QUERY_ADJ:
    return arboricity(Q);

  case QUERY_INIT: /* calcule Q->n et Q->rep[i][0..1] */
    Q->n=Q->param[0]; if(Q->n<3) RET_n(0);
    ALLOCMAT(Q->rep,Q->n,2); /* chaque sommet a deux pères */
    Q->param[1]=2; /* pour arboricity(Q) */
    int u,c,p;
    const int n=Q->n-2; /* n=nombre de sommets internes de l'arbre binaire */
    const int m=2*n+1; /* m=nombre total de sommets de l'arbre */
    int *T=Dyck(NULL,n,1,DYCK_KTREE); /* calcule un arbre binaire aléatoire T */

    /* Principe: on parcours T selon un parcours en profondeur
       modifié, on pose les deux fils avant la récursion.

       Dans l'exemple ci-dessous, la racine correspond au triangle
       (0,a,b). Puis, le fils gauche interne (=1) correspond au
       triangle (1,0,a). Le fils gauche interne 4 correspond au
       triangle (4,1,0). Les autres sommets feuilles ne correspondent
       à aucun triangle. Dans le graphe, les numéros des sommets sont
       décalés de +2.

                             0           
			    / \            a---b      0---1
                           1   2          / \ /      / \ /
                          / \            1---0      3---2
                         3   4            \ /        \ /
                            / \            4          4
                           5   6

       Plus précisément, à chaque nouveau fils u>0 de T qui est un
       fils interne on associe un nouveau sommet c du graphe, un coin
       du triangle. Si u est un fils gauche de v=T[u] alors on
       connecte c aux sommets v et Q->rep[v][0] de G. Si u est un fils
       droit, on connecte c aux sommets v et Q->rep[v][1]. La parité d'un
       noeud détermine s'il s'agit d'un fils droit ou gauche.

	T = [-,0,0,1,1,4,4] (= père dans l'arbre)
	C = [2,3,-,-,-,4,-] (= derniers sommets des triangles des noeuds internes de l'arbre)
    */

    NALLOC(int,C,m);
    for(u=1;u<m;u++) C[u]=0,C[T[u]]=1; /* par défaut u est une feuille, et son père non */
    C[0]=c=2; /* dernier sommet du triangle de la racine */
    for(u=1;u<m;u++) if(C[u]) C[u]=++c; /* met le numéro de centre */
    
    Q->rep[0][0]=Q->rep[0][1]=0;    /* le sommet 0 n'a aucun père */
    Q->rep[1][0]=-1,Q->rep[1][1]=0; /* le sommet 1 a un seul père, 0 */
    Q->rep[2][0]=0,Q->rep[2][1]=1;  /* le sommet 2 a pour pères 0 et 1 */

    for(u=1;u<m;u++){ /* on parcoure les noeuds de l'arbre */
      c=C[u];
      if(c){ /* si u est un noeud interne */
	p=C[T[u]]; /* p=centre du père de u dans l'arbre */
	Q->rep[c][0]=p;
	Q->rep[c][1]=Q->rep[p][u%2];
      }
    }

    free(C);
    free(T);
    return 0;
  }
  
return 1;
}


int treep(query* const Q){
/*
  Utilise Q->wrap. Q->wrap[i] = père de i dans l'arbre.  Les feuilles
  de cet arbre sans sommet de degré deux sont les sommets de numéro <
  p.
*/
  int *P=Q->wrap; /* raccourci pour Q->wrap */
  switch(Q->code){

  case QUERY_END:
    if(Q->n>0) free(P);
    return 0;

  case QUERY_INIT:;
    const int p=Q->param[0]; /* p=#feuilles */
    if(p<3) RET_n(0);
    ALLOC(P,(p<<1)-1); /* wrap[] = 1 tableau d'au plus 2p-2 entiers */

    /*
      Principe du calcul d'un arbre aléatoire à p feuilles et sans
      sommets de degré deux.

      L'idée est de construire l'arbre à partir des feuilles en
      déterminant le père de ces sommets, tout en garantissant qu'un
      sommet interne soit bien de degré au moins deux (sans son père).

      Soit A la liste des sommets qui n'ont encore pas de père. Cette
      liste de taille au plus p contient initialement toutes les
      feuilles: A={0,1,...,p-1}. Tant que |A|>2 on répète la procédure
      suivante:

      1. On tire un tableau aléatoire R de taille |A| entiers >=0 dont
         la somme fait |A| et possédant au moins une valeur > 1. Pour
         cela on fixe une case aléatoire à 2 puis on répète |A|-2 fois
         R[random()%|A|]++.

      2. Puis on détermine les pères de certains sommets de A en
         fonction de R. L'idée est que si R[k]>1, alors les sommets
         A[u]...A[u+R[k]-1] vont recevoir le même père, un nouveau
         sommets crée et réinjecté à la liste A. Plus précisément, on
         parcoure séquentiellement R. Si R[k]=0, alors on passe
         simplement à k+1 sans rien faire d'autre. Si R[k]=1, on
         maintient A[u] dans la liste A on passe au prochain sommet
         u+1 de A. Le père du sommet A[u] n'est alors toujours pas
         fixé. Par contre, si R[k]>1, alors on crée un nouveau sommet
         v que l'on ajoute à la fin de la liste A, et les R[k] sommets
         A[u]...A[u+R[k]-1] ont alors tous pour père le sommet v.

      Lorsque |A|<=2 on ne rajoute plus aucun sommet, et le nombre de
      sommets N est alors déterminé. On fixe alors que le sommet A[0]
      est racine. Si |A|=2 alors le père de A[1] est fixé à A[0]. Il
      n'est pas possible d'avoir |A|=0 car tout regroupement de sommet
      crée au moins un nouveau sommet dans père (et donc ajouté à A).

      Une complexité de O(p^2) est possible car |A| diminue seulement
      d'une unité si à chaque étape R ne contient qu'une seule valeur
      > 1. Cependant, en moyenne O(log(p)) étapes suffisent, soit
      O(p*log(p)) en tout, car il est facile de voir que R contient
      une fraction de |A| valeurs > 2. (Dans la représentation unaire
      de R il y a n/2 blocks et la moitié sont de longueur > 2.) Et
      pour chacun de tels blocks tous sauf 1 seront enlevé de A.
    */

    int t;
    NALLOCZ(int,A,p,_i); /* tableau des sommets actifs */
    NALLOC(int,R,p); /* tableau des valeurs random */
    int u; /* u=indice dans A du prochain sommet actif */
    int q; /* q=indice dans A du prochain nouveau sommet actif, q<=u */
    int k; /* k=indice dans R, k>=u */
    int a; /* a=taille de A, a<=p */
    Q->n=a=p; /* Q->n=nb courant de sommets dans le graphe, Q->n>=p */

    while(a>2){
      for(k=0;k<a;R[k++]=0); /* tableau R à zéro */
      R[random()%a]=2; /* met un "2" quelque part */
      for(k=2;k<a;k++) R[random()%a]++; /* incrémente a-2 fois des positions de R */
      for(k=u=q=0;k<a;k++){ /* parcoure les valeurs de R */
	if(R[k]==0) continue;
	if(R[k]==1) { A[q++]=A[u++]; continue; }
	t=u+R[k]; /* ici t>=2 */
	for(;u<t;u++) P[A[u]]=Q->n; /* P[A[u]]=père de A[u]=nouveau sommet */
	A[q++]=Q->n++; /* un sommet de plus, et un nouveau actif de plus */
      }
      a=q; /* nouvelle taille de A=nb de nouveaux sommets actifs */
    }

    P[A[0]]=-1;
    if(a==2) P[A[1]]=A[0];
    
    free(A);
    free(R);
    REALLOC(P,Q->n); /* recalibrage du tableau */
    Q->wrap=P; /* met à jour Q->wrap */
    return 0;

  case QUERY_ADJ:
    RET_a(((P[Q->i]==Q->j)||(P[Q->j]==Q->i))); /* arbre */
    return 0;
  }

  return 1;
}


int halin(query* const Q){
/*
  Utilise treep().
*/
  const int p=Q->param[0]; /* p=#feuilles */
  if((Q->code==QUERY_ADJ)&&(Q->i<p)&&(Q->j<p)) /* cycle */
    RET_a(((Q->j==((Q->i+1)%p))||(Q->i==((Q->j+1)%p))));
  
  return treep(Q); /* arbre */
}


int permutation(query* const Q){
/*
  Utilise wrap, wrap[i] = permutation du sommet i.
*/
  switch(Q->code){

  case QUERY_END:
    if(Q->n>0) free(Q->wrap);
    return 0;
    
  case QUERY_NAME:
    sprintf(Q->name,"(%i,%i)",Q->i,Q->wrap[Q->i]);
    return 0;

  case QUERY_INIT: /* permutation aléatoire de [0,n[ dans Q->wrap */
    SET_n(Q->param[0]);
    ALLOCZ(Q->wrap,Q->n,_i);
    Permute(Q->wrap,Q->n);
    return 0;
    
  case QUERY_ADJ:
    RET_a(((Q->i-Q->wrap[Q->i])*(Q->j-Q->wrap[Q->j])<0));
  }

  return 1;
}


int interval(query* const Q){
/*
  A chaque sommet i correspond un intervalle [a,b] de [0,2n[, avec
  a=Q->rep[i][0] et b=Q->rep[i][1].
*/
  switch(Q->code){
    
  case QUERY_END:
    return free_rep(Q);

  case QUERY_NAME:
    sprintf(Q->name,"[%i,%i]",Q->rep[Q->i][0],Q->rep[Q->i][1]);
    return 0;

  case QUERY_INIT:
    SET_n(Q->param[0]);
    const int m=(Q->n<<1);
    int k,x;

    /* génère un intervalle Q->rep[k] pour k, [a,b] dans [0,2n[ avec a<=b */
    ALLOCMAT(Q->rep,Q->n,2);
    for(k=0;k<Q->n;k++){
      x=random()%m;
      Q->rep[k][0]=x;
      Q->rep[k][1]=x+random()%(m-x);
    }
    return 0;

  case QUERY_ADJ:
    RET_a(
	  ((Q->rep[Q->i][0]<=Q->rep[Q->j][0])&&(Q->rep[Q->j][0]<=Q->rep[Q->i][1])) ||
	  ((Q->rep[Q->i][0]<=Q->rep[Q->j][1])&&(Q->rep[Q->j][1]<=Q->rep[Q->i][1])) ||
	  ((Q->rep[Q->j][0]<=Q->rep[Q->i][0])&&(Q->rep[Q->i][1]<=Q->rep[Q->j][1]))
	  );
  }
  
  return 1;
}


int circle(query* const Q){
/*
  Graphe d'inclusion d'intervalle de [0,2n[, avec a=Q->rep[i][0] et
  b=Q->rep[i][1].
*/
  if(Q->code==QUERY_ADJ)
    RET_a(
	  ((Q->rep[Q->i][0]<=Q->rep[Q->j][0])&&(Q->rep[Q->j][1]<=Q->rep[Q->i][1])) ||
	  ((Q->rep[Q->j][0]<=Q->rep[Q->i][0])&&(Q->rep[Q->i][1]<=Q->rep[Q->j][1]))
	  );
  return interval(Q);
}


int sat(query* const Q){
/*
  Utilise: i<j.

  [0,2n[: les variables positives et négatives
  [2n+t*k,2n+(t+1)*k[: t-ème clique de taille k

  Chaque sommet-clause t est connecté à une variable (positive ou
  négative) aléatoire stockée dans Q->rep[t][0].
*/
  
  const int n=Q->param[0];
  const int m=Q->param[1];
  const int k=Q->param[2];
  const int n2=2*n;
  int t; 
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);
    
  case QUERY_INIT:
    if((n<1)||(m<1)||(k<1)) RET_n(0);
    SET_n(n2+m*k);
    ALLOCMAT(Q->rep,Q->n,1);    
    for(t=n2;t<Q->n;t++) Q->rep[t][0]=random()%n2;
    return 0;
    
  case QUERY_ADJ:
    if(Q->j<n2) RET_a(((Q->j==Q->i+1)&&(Q->j&1))); /* i-j et j impaire */
    if(Q->i>=n2) RET_a((Q->j-Q->i<=k)); /* i et j dans la même clique ? */
    RET_a((Q->rep[Q->j][0]==Q->i)); /* i dans le matching et j dans une clique */
  }
  
  return 1;
}


int gpetersen(query* const Q){
/*
  Utilise i<j.
  u_i dans [0,n[ et v_i dans [n,2n[, Q->n=2n.
*/
  const int n=Q->param[0];
  const int r=Q->param[1];
  switch(Q->code){

    case QUERY_NAME:
      sprintf(Q->name,"%c%i",(Q->i<n)?'u':'v',Q->i%n);
      return 0;

    case QUERY_INIT:
      RET_n(2*n);
      
    case QUERY_ADJ:
      /* u_i-v_i */
      if(Q->j==Q->i+n) RET_a(1);

      /* sinon, cas fréquent, pas d'arête entre u_i et v_j */
      if((Q->i<n)&&(n<=Q->j)) RET_a(0);

      /* u_i-u_{i+1 mod n}, ici i<j<n */
      if(Q->i<n) RET_a((Q->j==(Q->i+1)%n)||(Q->i==(Q->j+1)%n));
      
      /* v_i-v_{i+r mod n} */
      /* ici i,j<n mais j<i possible*/
      RET_a((Q->j-n==(Q->i+r)%n)||(Q->i-n==(Q->j+r)%n));
  }

  return 1;
}


int antiprism(query* const Q){
  if(Q->code==QUERY_ADJ)
    if(Q->j==Q->param[0]+(Q->i+1)%Q->param[0]) RET_a(1);
  return gpetersen(Q);
}


int deltohedron(query* const Q){
/*
  Utilise i<j.
  Les sommets de [0,n[ forme le cycle avec n=Q->param[0].
*/
  const int n=Q->param[0];
  switch(Q->code){

  case QUERY_INIT:
    if(n<=0) RET_n(0);
    RET_n(n+2);

  case QUERY_ADJ:
    RET_a(( (Q->j==(Q->i+1)%n) ||
	      ((Q->j==n-1)&&(Q->i==0)) ||
	      ((Q->j==n)&&((Q->i&1)==0)) ||
	      ((Q->j==n+1)&&((Q->i&1)==1)) )
	    );
  }

  return 1;
}


int helm(query* const Q){
/*
  Utilise i<j.
  On pourrait faire "wheel n -star -1".
*/
  int n=Q->param[0];
  switch(Q->code){

  case QUERY_INIT:
    if(n<3) RET_n(0);
    RET_n(2*n+1);

  case QUERY_ADJ:
    if(Q->j-Q->i==n) RET_a(1); // branche
    if(Q->j>n) RET_a(0);
    if(Q->i==0) RET_a(1); // roue à n rayons
    RET_a((Q->j-Q->i==1)||((Q->j==n)&&(Q->i==1)));
  }
  
  return 1;
}


int haar(query* const Q){
/*
  Utilise i<j.
  Utilise un paramètre auxiliaire Q->param[1].
*/
  const unsigned n=Q->param[0];
  const int k=Q->param[1];
  switch(Q->code){

  case QUERY_INIT:
    Q->param[1]=lg(n); // paramètre auxiliaire
    RET_n(2*Q->param[1]); // n=0 si k=0

  case QUERY_ADJ:
    if((Q->i>=k)||(Q->j<k)) RET_a(0); // 0 si i,j dans la même part
    const int j=(Q->i-(Q->j-k)+k)%k; // u_i adjacent à v_{i+j mod k} ?
    RET_a((n>>j)&1); // vrai ssi bit-j de n est à 1 ?
  }
  
  return 1;
}


int turan(query* const Q){
/*
  Utilise i<j.
*/
  const int n=Q->param[0];
  const int r=Q->param[1];
  switch(Q->code){

  case QUERY_INIT:
    if((n<0)||(r<=0)||(n<r)) RET_error(6); /* il faut n>=r>0 */
    RET_n(n);

  case QUERY_ADJ:
    RET_a((Q->j-Q->i)%r);
  }
  
  return 1;
}


int kneser(query* const Q){
/*
  Q->rep[i][0...k-1] sont les ensembles représentant les sommets.
*/
  int v,x,y;
  const int n=Q->param[0];
  const int k=Q->param[1];
  const int r=Q->param[2];
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);

  case QUERY_NAME:
    name_vector(Q->name,Q->rep[Q->i],k,",","{}",1,"%i");
    return 0;

  case QUERY_INIT:;
    if((n<0)||(k<0)||(k>n)) RET_n(0);
    SET_n(Binom(n,k)); /* on a besoin de Q->n pour allouer Q->rep */
    ALLOCMAT(Q->rep,Q->n,k); /* on connaît Q->n */
    if(Q->n==1) return 0; /* si Q->n=1, fin: graphe à 1 sommet */
    NALLOC(int,S,k);
    NextSet(S,-1,k); /* premier sous-ensemble */
    for(x=0;x<Q->n;x++){ /* pour tous les sommets x du graphe */
      for(y=0;y<k;y++) Q->rep[x][y]=S[y]; /* copie dans Q->rep[x] */
      NextSet(S,n,k); /* sous-ensemble suivant */
    }
    free(S);
    return 0;

  case QUERY_ADJ:
    /*
      Calcule si l'intersection possède au plus r éléments.
      L'algorithme ici est en O(k) en utilisant le fait que les
      éléments de Q->rep sont rangés dans l'ordre croissant.
    */
    v=x=y=0; /* indices pour i et j, v=nb d'élements commun */

    while((x<k)&&(y<k)&&(v<=r))
      if(Q->rep[Q->i][x]==Q->rep[Q->j][y]) v++,x++,y++;
      else if(Q->rep[Q->i][x]<Q->rep[Q->j][y]) x++; else y++;

    RET_a((v<=r));
  }

  return 1;
}


int rig(query* const Q){
/*
  Q->rep[i][0]=taille t_i de l'ensemble associé au sommet i
  Q->rep[i][1...t_i]=ensemble associé au sommet i
*/
  int x,y,k,t;
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);
    
  case QUERY_NAME:
    name_vector(Q->name,Q->rep[Q->i]+1,Q->rep[Q->i][0],",","{}",1,"%i");
    return 0;
    
  case QUERY_INIT:
    SET_n(Q->param[0]);
    const double p=Q->dparam[0];
    k=Q->param[1];
    NALLOC(int,S,k+1); /* ensemble S[1...k] temporaire pour un sommet */
    ALLOC(Q->rep,Q->n);
    for(x=0;x<Q->n;x++){ /* pour chaque sommet x */
      t=0; for(y=1;y<=k;y++) if(RAND01<p) S[++t]=y; /* t=taille de S */
      ALLOC(Q->rep[x],t+1); /* espace pour le sommet x */
      Q->rep[x][0]=t; /* écrit S dans Q->rep[x][1...t] */
      for(y=1;y<=t;y++) Q->rep[x][y]=S[y];
    }
    free(S);
    return 0;

  case QUERY_ADJ:
    /*
      Détermine si l'intersection de Q->rep[i][1...] et
      Q->rep[j][1...]  est vide ou pas.  L'algorithme utilise le fait
      que les éléments de Q->rep sont rangés dans un ordre croissant.
    */

    x=y=1; /* indices pour les ensemble de i et j */
    k=Q->rep[Q->i][0]; /* taille de l'ensemble de i */
    t=Q->rep[Q->i][0]; /* taille de l'ensemble de j */

    while((x<=k)&&(y<=t)){
      if(Q->rep[Q->i][x]==Q->rep[Q->j][y]) RET_a(1);
      if(Q->rep[Q->i][x]<Q->rep[Q->j][y]) x++; else y++;
    }
    RET_a(0);
  }

  return 1;
}


int bdrg(query* const Q){
/*
  param[]={2t,n_1,d_1, ... n_t,d_t}, donc param[0]=|param|-1.  Crée
  dans Q->G un graphe dont la distribution des degrés est donnés par
  (n_i,d_i). La construction est basée sur un matching aléatoire des
  demi-arêtes. Pour cela on construit un tableau T où chaque sommet u
  est répété deg(u) fois dans T (cela forme les demi-arêtes sortantes
  de u). Les arêtes du graphe sont alors les arêtes simples entre les
  sommets T[2*i] et T[2*i+1] (on supprime boucle et arête-multiple).
  La taille de T est égale à ∑d_i. Si elle n'est pas paire, on
  supprime une demi-arête à un sommet ayant un d_i>0. C'est toujours
  correct, car si d_i=0 pour tous les i, c'est que ∑d_i est paire.

  Considérons l'exemple suivant:

    - 1 sommet  de degré 1 (le sommet 0)
    - 2 sommets de degré 2 (les sommets 1 et 2)
    - 1 sommets de degré 3 (le sommet 3)

    => param[]={6,1,1,2,2,1,3}
    => T[]=[ 0 1 1 2 2 3 3 3 ] (on répète le sommet u deg(u) fois)
    => permute(T) = [ 2 3 2 0 3 1 3 3 ]
    => arêtes:        2-3,2-0,3-1,3-3
    => arêtes simples: 2-3, 2-0, 3-1
*/
  switch(Q->code){

  case QUERY_END:
  case QUERY_ADJ:
    return load(Q);

  case QUERY_INIT:;
    int c,u,u0,v,k,p;
    c=Q->param[0]; /* c=2*t=nombre de valeurs dans param[] sans le premier */
    if(c&1) RET_n(0); /* il faut un nombre paires de valeurs */
    for(Q->n=0,k=1;k<=c;k+=2) Q->n+=Q->param[k]; /* Q->n=nombre de sommets */
    TEST_n; // fin si graphe vide
    Q->G=new_graph(Q->n); /* crée le graphe (non vide) */
    for(k=u=0;k<c;k+=2){ /* alloue les listes */
      if(Q->param[k+2]>0) u0=u; /* mémorise un sommet avec un degré d_i>0 */
      for(p=0;p<Q->param[k+1];p++,u++){
	Q->G->d[u]=Q->param[k+2]; /* degré du sommet u */
	ALLOC(Q->G->L[u],Q->G->d[u]); /* liste d'adjacence de u */
      }
    }
    for(v=k=0;k<c;k+=2) v+=Q->param[k+1]*Q->param[k+2]; /* v=∑ d_i*n_i=2*|E| */
    if(v&1){ /* si la somme est impaire, on enlève une demi-arête au sommet u0 */
      Q->G->d[u0]--; /* u0 est nécessairement défini si la somme était impaire */
      REALLOC(Q->G->L[u0],Q->G->d[u0]); /* raccourcie la liste de u0 */
      v--; /* une demi-arête de moins */
    }
    c=v; /* c=nombre de demi-arêtes */
    NALLOC(int,T,c); /* T=tableau des demi-arêtes */
    for(u=k=0;u<Q->n;u++) /* parcoure tous les arcs et remplit T */
      for(p=0;p<Q->G->d[u];p++) T[k++]=u;
    Permute(T,c); /* mélange aléatoire les demi-arêtes */
    degres_zero(Q->G); /* pour ADD_EDGE */
    for(p=0;p<c;){ /* parcoure T et ajoute les arêtes simples */
      u=T[p++]; v=T[p++]; if(u==v) continue; /* pas de boucle */
      if(Q->G->d[u]>Q->G->d[v]) SWAP(u,v,k); /* pour une recherche plus rapide */
      if(SetSearch(v,Q->G->L[u],Q->G->d[u],0)<0) /* v dans Q->G->L[u] ? */
	ADD_EDGE(Q->G,u,v); /* non: ajoute u-v et met à jour la taille des listes */
    }
    free(T);
    GraphRealloc(Q->G,Q->G->d);
    return 0;
  }

  return 1;
}


// Deux routines qui servent plusieurs fois dans fdrg().

// Echange dans list[] les mini-sommets u1 et u2 tout en mettant à
// jour pos[]; t est une variable temporaire pour SWAP.
#define SWAP_MINI_SOMMET(u1,u2,t)		\
  SWAP(list[pos[u1]],list[pos[u2]],t);		\
  SWAP(pos[u1],pos[u2],t);


// Vérifie si les l'arêtes entre les sommets vi et vj formeraient une
// boucle ou une multi-arête. Si c'est le cas, on met t=1, sinon t=0.
// Pour savoir si vi et vj sont déjà voisins on cherche s'il existe un
// mini-sommet de vi qui a été pris (d'indice >= e) et qui a pour
// parent vj. Pour avoir un temps de recherche en
// O(min{deg(vi),deg(vj)}), cherche aussi dans les voisins de vj si
// deg[vj]<deg[vi].  NB: k^1 représente l'indice suivant k si k est
// pair, et précédant si k est impair.
#define	NOT_GOOD_PAIR(vi,vj,t)				\
  if(vi==vj) t=1;					\
  else{							\
    if(deg[vj]<deg[vi]) u=vj,v=vi; else u=vi,v=vj;	\
    int d=deg[u],f=first[u],k;				\
    for(t=0;t<d;t++){				      	\
      k=pos[f+t];					\
      if((k>=e)&&(parent[list[k^1]]==v)) t=d;		\
    }							\
    t=(t>d);						\
  }


int fdrg(query* const Q){
/*
  param[] = {2t, n_1,d_1,..,n_t,d_t}.  Crée dans Q->G un graphe dont
  la distribution des degrés est donnés par param[]. La construction
  est basée sur l'algorithme (avec rejet) de [BKS10] (section 6.1)
  basé sur l'implémentation de [SW97]. Une pré-version a été
  programmée par Aliénor Brabant, en stage de L2 en juin 2016.

  L'algorithme est le suivant (ici deg(i)=degré du sommet i souhaité,
  D[i]=degré restant du sommet i, M=sum_i deg(i)=2 fois le nombre
  total d'arêtes souhaitées):

  (1) E={}, D[i]=deg(i) pour i=0..n-1, V={0,...,n-1}

  (2) Choisir i,j de V avec proba p_ij = D[i]*D[j] *
      (1-deg(i)*deg(j)/2M) parmi toutes paires avec (i,j) tq i<>j et
      {i,j} pas dans E. Ajouter {i,j} à E puis décrémenter D[i] et
      D[j].

  (3) Répéter (2) jusqu'à ce qu'on ne puisse plus rajouter d'arête à E

  (4) Si |E|<M/2, recommencer sinon G=(V,E).

  L'algorithme termine à l'étape (4) la première fois avec probabilité
  1-o(1). Pour les graphes réguliers, on peut prendre plus simplement
  p_ij = D[i]*D[j]. Pour générer un biparti aléatoire (de séquence
  fixée) on peut remplacer (2) par p_ij = D[i]*D[j] *
  (1-deg(i)*deg(j)/M) et ajouter seulement les arêtes qui ne sont pas
  dans une même part.

  On implémente l'algorithme en gérant une liste de "mini-sommets",
  qui sont les brins incidents aux sommets et qu'on doit connectés.
  Chaque sommet i possède deg(i) mini-sommets. On choisit donc deux
  mini-sommets avec la bonne probabilité p_ij tels que si connectés,
  ils ne forment ni de boucles ni de multi-arêtes. Si c'est le cas, on
  sélectionne cette paire, pour former une nouvelle arête, et l'on
  déplace ces deux mini-sommets dans la 2e partie de la liste. A la
  fin, la liste représente toutes les paires de mini-sommets
  connectés, et donc la liste des arêtes.

  Pour obtenir une complexité en O(M*dmax+dmax^4) il faut découper
  l'étape (2) en 3 phases. Il faut aussi faire attention à l'étape
  (3).  Certes, on ne peut ajouter que O(M) arêtes, mais détecter si
  l'on peut ou non ajouter encore une arête n'est pas évident. En
  effet, il peut rester des paires de mini-sommets non connectés sans
  pourtant pouvoir ajouter une seule arêtes si toutes les paires de
  mini-sommets ont des sommets déjà connectés. Voici un exemple de
  blocage pour n=5 sommets et une séquence de degré 2,2,2,3,3.

                         sommets={0,1,...,4}
                    mini-sommets={a,b,...,l}
  
     0  1  2  3   4                            0  1  2  3   4
    /│ /│ /│ /│\ /│\                          /│ /│ /│ /│\ /│\  
    ab cd ef ghi jkl                          ab cd ef ghi jkl 
    ││ ││ ││ ││* ││*                          ││ ││ ││ │││ │││
    │└─┘│ ││ │└──┘│                           │└─┘└─┼┼─┘│└─┘││
    └───┼─┘│ │    │                           └─────┼┼──┘   ││
        └──┼─┘    │                                 └┼──────┼┘
           └──────┘                                  └──────┘

  choix: ae,bc,dg,fk,hj                    choix: ah,bc,dg,el,fk,ij
  dernier choix impossible

  Cependant, pour être bloqué il doit exister dans le graphe courant
  une clique de k sommets ayant chacun encore au moins un mini-sommet
  libre. Une clique de k sommets ayant encore un mini-sommet de libre
  implique k <= dmax: un sommet i de la clique doit être connecté à
  k-1 autre sommet et ce degré est <= dmax-1 puisque un mini-sommet
  est encore libre. Autrement dit, si le nombre d'arêtes restantes M -
  ∑ D[i] = M-e > dmax^2 on ne peut être bloqué puisqu'il reste plus de
  dmax sommets avec au plus dmax-1 mini-sommets de libres.
  
*/
  switch(Q->code){

  case QUERY_ADJ:
    return load(Q);

  case QUERY_INIT:;
    int c=Q->param[0]; // c=2*t=nombre de valeurs dans param[] sans le premier
    if((c&1)||(c<0)) RET_error(6); // il faut un nombre paires de valeurs
    Q->n=graphical(Q->param+1,c/2); // Q->n=nombre de sommets
    if(Q->n<0) RET_error(34); // erreur: séquence non graphique
    TEST_n; // fin si graphe vide (si n=0)
    int u,v,i,j,k,t,d,vi,vj,ui,uj,e,p,nv,nh;

    /* calcule le nombre M de mini-sommets, pour les ALLOC() suivants */
    for(t=0,i=1;i<=c;i+=2) t+=Q->param[i]*Q->param[i+1];
    const int M=t; // M=nombre total de mini-sommets=∑d_i

    /* Conventions:
         i = indice ou position, entier de [0,M[
         u = mini-sommet, entier de [0,M[
         v = sommet, entier de [0,Q->n[
    */
    
    NALLOC(int,list,M); // list[i]=u, liste des mini-sommets (au départ list[i]=i)
    NALLOC(int,pos,M); // pos[u]=i, position i du mini-sommets u dans list[]
    NALLOC(int,parent,M); // parent[u]=v, sommets correspondant au mini-sommets u
    NALLOC(int,listv,Q->n); // listv[i]=v, liste des sommets libres
    NALLOC(int,dfree,Q->n); // dfree[v]=nombre de mini-sommets libre pour le sommet v
    NALLOC(int,first,Q->n); // first[v]=1er mini-sommet du sommet v;
   // les mini-sommets de v sont donc first[v], first[v]+1,..., first[v]+deg[v]-1

    /* On parcoure les (n_i,d_i) et calcule:
       - le degré max
       - la taille des listes d'adjacence du graphe final Q->G
       - le tableau first[]
       - le tableau parent[]
    */

    Q->G=new_graph(Q->n); // crée le graphe, alloue Q->G->d
    t=u=v=0; // t=degré max, u=mini-sommet, v=sommet

    for(i=1;i<=c;i+=2){
      d=Q->param[i+1]; // degré courant
      t=imax(t,d); // t=max{d_i}
      for(j=0;j<Q->param[i];j++,v++){
	first[v]=u; // premier mini-sommet du sommet v
	for(k=0;k<d;k++) parent[u++]=v; // tous le même parent v
	Q->G->d[v]=d; // degré du sommet v
	ALLOC(Q->G->L[v],d); // liste d'adjacence de u
      }
    }

    const int M2=2*M; // constante pour les probas
    const int dmax=t; // dmax=degré max du graphe final
    const int ddmax=t*t; // constante dmax^2
    const int pair=(~1); // mask permettant (avec &) de rendre pair un entier 
    const int noreg=(c>2); // noreg=vrai ssi le graphe n'est pas régulier
    int * const deg=Q->G->d; // raccourci: deg[v]=degré de v dans le graphe final
    NALLOC(int,listh,4*ddmax); // listh[]=liste d'arêtes possible pour la PHASE 3

    for(;;){ // forever ... sauf si la PHASE 3 réussit. En principe
	     // elle réussit dès le 1er coup avec proba 1-o(1)

      e=M; // e=nombre de mini-sommets encore disponibles

      // PHASE 1
      //
      // On titre uniformément deux mini-sommets pris dans list[0..e[,
      // e étant le nombre de mini-sommets disponibles (non
      // connectés). La probabilité d'obtenir la paire de sommets
      // (vi,vj) est ainsi proportionnelle à D[vi]*D[vj]. On garde
      // cette paire avec proba 1-deg(vi)*deg(vj)/2M (ou proba 1 si
      // c'est un graphe régulier). Ainsi, la proba de succès (dans le
      // cas non-régulier) est proprotionnelle au produit (D[i]*D[j])
      // * (1-deg(vi)*deg(vj)/2M). On rejette si c'est une boucle ou
      // si vi et vj étaient déjà voisins, ce qui prend un temps de
      // O(dmax). Le nombre de répétitions est 2 en moyenne. La
      // complexité totale de la PHASE 1 est O(m*dmax) en moyenne.
      
      // initialise list[] et pos[]
      for(i=0;i<M;i++) list[i]=pos[i]=i;

      while(e>2*ddmax){ // on ne peut pas être bloqué avec ce test

	// choisit deux mini-sommets: list[i] et list[j]
	i=random()%e, vi=parent[list[i]];
	j=random()%e, vj=parent[list[j]];
	if(noreg) // NB: si le graphe est régulier, on accèpte toujours vi,vj
	  if(M2*RAND01>M2-deg[vi]*deg[vj]) continue; // rejet ?
	NOT_GOOD_PAIR(vi,vj,t); // vi=vj ou si vi-vj ?
	if(t) continue; // recommence si pas bon
	
	// on ajoute l'arête vi-vj en déplaçant les mini-sommets
	// list[i] et list[j] en position e-1 et e-2 dans list[] tout en
	// mettant à jour pos[] et e.
	ui=list[i],uj=list[j];
	k=list[--e]; SWAP_MINI_SOMMET(ui,k,t); // échange ui et k
	k=list[--e]; SWAP_MINI_SOMMET(uj,k,t); // échange uj et k

	// NB: Après le premier SWAP_MINI_SOMMET(), list[] est
	// modifiée, et donc il est possible que le mini-sommet uj ne
	// soit plus en position j dans list[] (si j=e-1 par exemple).
	// Il faut donc fixer ui=list[i] et uj=list[j] AVANT de faire
	// le premier SWAP_MINI_SOMMET().
      }
    
      // PHASE 2
      //
      // Cette phase est similaire à la PHASE 1, sauf qu'on considère
      // plutôt les sommets libres que les mini-sommets libres. On
      // reste dans cette phase tant que le nombre de sommets encore
      // libres (nv) est assez grand: nv > 2dmax. On choisit deux
      // sommets libres, et on répète jusqu'à en trouver deux qui ne
      // soient pas déjà voisins (2 répétitions en moyenne). On ne
      // peut pas être bloqué puisqu'on ne peut pas avoir de clique si
      // le nombre de sommets restant est > dmax. Puis on choisit
      // uniformément un mini-sommet pour chacun des sommets, et on
      // répète jusqu'à en avoir deux non connectés (O(dmax^2)
      // répétitions en moyenne). La complexité totale de la PHASE 2
      // est O(dmax^4) en moyenne.
      
      // initialise (1) listv[0..nv[, la liste des sommets libres,
      // c'est-à-dire ayant encore au moins un mini-sommet de libres;
      // (2) dfree[v], le nombre de mini-sommets de libres pour v si v
      // est dans listv[] (non défini sinon).
      for(v=nv=0;v<Q->n;v++){
	u=first[v]; // u=premier mini-sommet de v
	d=k=deg[v]; // d=dfree[v]=deg[v] au départ
	for(i=0;i<k;i++) if(pos[u+i]>=e) d--;
	if(d) dfree[v]=d,listv[nv++]=v; // on ajoute v à listv[] si d>0
      }

      while(nv>2*dmax){ // ne peut pas être bloqué si nv>dmax
	
	// choisit deux sommets de listv[] non voisins
	i=random()%nv, vi=listv[i];
	j=random()%nv, vj=listv[j];
	NOT_GOOD_PAIR(vi,vj,t); // vi=vj ou si vi-vj ?
	if(t) continue; // recommence si pas bon
	
	// met à jour dfree[], puis éventuellement listv[] et nv
	dfree[vi]--; if(dfree[vi]==0) listv[i]=listv[--nv]; // déplace en supprimant vi
	if(j==nv) j=i; // NB: dans ce cas vj a pu être déplacé en i
	dfree[vj]--; if(dfree[vj]==0) listv[j]=listv[--nv]; // déplace en supprimant vj

	// choisit un mini-sommet pour vi et pour vj et vérifie qu'ils
	// ne soient pas déjà connectés; on ne peut pas être bloqué
	// car vi et vj ne sont pas connectés et ils ont chacun un
	// mini-sommet de libre. Le nombre de répétitions est
	// O(dmax^2) en moyenne.
	do{
	  ui=first[vi]+(random()%deg[vi]);
	  uj=first[vj]+(random()%deg[vj]);
	}while((pos[ui]>=e)||(pos[uj]>=e));
	
	// on ajoute l'arête entre les mini-sommets ui et uj, et met à
	// jour list[], pos[] et e
	k=list[--e]; SWAP_MINI_SOMMET(ui,k,t); // échange ui et k
	k=list[--e]; SWAP_MINI_SOMMET(uj,k,t); // échange uj et k
      }
    
      // PHASE 3
      //
      // On considère le graphe H des arêtes encore possibles, donc
      // conctant les sommets de listv[]. Il possède nv <= 2dmax
      // sommets et <= nv*(nv-1)/2 < 2dmax^2 arêtes. On choisit une
      // arête aléatoire uniforme de H qu'on accepte ensuite avec une
      // probabilité D[i]*D[j]/dmax^2. Le nombre de répétitions est
      // donc O(dmax^2) en moyenne. On met à jour H et on répète
      // jusqu'à ce que H n'ait plus de sommets libres.  Pour mettre à
      // jour H, il faut supprimer l'arête i-j choisie mais aussi
      // supprimer toutes les arêtes de H où l'un des sommets i ou j
      // apparaît si jamais D[i] ou D[j] passe à 0. Si, à l'issue de
      // la PHASE 3 le nombre d'arêtes ajoutées au graphe final est <
      // M/2 (H devient vide trop tôt), il faut alors recommencer les
      // trois phases. La complexité totale de la PHASE 3 est
      // O(dmax^4) en moyenne.
      
      // On considère listh[0..nh[, la liste des arêtes possibles de H
      // (graphe qu'on ne construit pas en fait) où l'on stocke chaque
      // arête possible dans des cases consécutives de listh[]. La
      // taille de ce tableau doit être nv*(nv-1) < 4dmax^2. On
      // l'initialise en balayant toutes les paires possibles de
      // sommets de listv[] et en testant leur adjacence.  Aussi, on
      // va stocker dans list[0..p[ (au début donc) les arêtes
      // générées à la PHASE 3. On a la place car celles générées aux
      // phases précédantes ont été stockées dans list[e..M[ (à la fin
      // donc) et p<=e.

      for(i=nh=0;i<nv;i++)
	for(j=i+1;j<nv;j++){ // NB: i<j
	  vi=listv[i],vj=listv[j]; // vi,vj=sommets encore libre
	  NOT_GOOD_PAIR(vi,vj,t); // t=(vi=vj ou vi-vj)
	  if(t==0) listh[nh++]=vi,listh[nh++]=vj; // ajoute vi-vj
	}
      
      p=0; // p=position de la prochaine entrée libre dans list[]
      while(nh>0){ // tant qu'il reste une arête dans listh[]
      
	// choisit une arête i-j de listh[] uniformément
	t=(random()%nh)&pair; // position (paire) aléatoire
	vi=listh[t],vj=listh[t+1];
	if(RAND01*ddmax>dfree[vi]*dfree[vj]) continue; // rejet ?

	// ajoute l'arête vi-vj dans list[] et l'enlève de listh[]
	list[p++]=vi,list[p++]=vj;
	listh[t+1]=listh[--nh],listh[t]=listh[--nh];

	// met à jour dfree[]
	i=(--dfree[vi]==0); // i=vrai ssi il faut supprimer vi dans listh[]
	j=(--dfree[vj]==0); // j=vrai ssi il faut supprimer vj dans listh[]

	if(i||j){
	  // ici il faut supprimer vi ou vj (ou les deux) dans tout
	  // listh[]. On parcoure listh[] une 1ère fois et on met des
	  // paires (-1,-1) si vi ou vj apparaît, puis une 2e fois
	  // pour enlever les (-1,-1) en délaçant la dernière de
	  // listh[]. NB: si v est en position t dans listh[], alors
	  // son voisin est en position t^1. Cela prend un temps
	  // |listh|=O(dmax^2), le même temps pour accepter vi-vj en
	  // moyenne.
	  for(t=0;t<nh;t++)
	    if((i&&(listh[t]==vi))||(j&&(listh[t]==vj)))
	      listh[t]=listh[t^1]=-1; // efface l'arête contenant vi ou vj
	  t=0;
	  while(t<nh)
	    if(listh[t]<0) listh[t+1]=listh[--nh],listh[t]=listh[--nh];
	    else t+=2;
	}
      }
      
      if(p==e) break; // fin: list[] est complètement remplie
    }
    
    // on écrit les arêtes de list[] dans Q->G, le codage des arêtes
    // n'est pas le même avant et après la position e
    degres_zero(Q->G); // pour ADD_EDGE(Q->G,...)
    for(i=0;i<e;i+=2) ADD_EDGE(Q->G,list[i],list[i+1]);
    for(;i<M;i+=2) ADD_EDGE(Q->G,parent[list[i]],parent[list[i+1]]);
    
    free(list);
    free(pos);
    free(parent);
    free(first);
    free(listv);
    free(listh);
    free(dfree);

    return 0;
  }

  return 1;
}


int seg_intersection(point p1,point q1,point p2,point q2){
/*
  Renvoie vrai ssi le segment ]p1q1[ intersecte le segment
  ]p2q2[. Pour s'intersecter il faut que: (1) les points p2,q2 soient
  de part et d'autre de la droite portée par (p1q1); et (2) les points
  p1,q1 soient de part et d'autre de la droite portée par
  [p2q2]. Sinon, si les points sont colinéaires et il faut tester si
  l'un des points appartient au segment de l'autre.

  Pour déterminer si point X est au-dessus, en-dessous ou sur la
  droite portée par (AB) il suffit de calculer le signe de
  det(X-A,B-A).
*/
  const double dp2=det(q1.x-p1.x,q1.y-p1.y,p2.x-p1.x,p2.y-p1.y); // p2 au-dessus de (p1q1) ?
  const double dq2=det(q1.x-p1.x,q1.y-p1.y,q2.x-p1.x,q2.y-p1.y); // q2 au-dessus de (p1q1) ? 
  const double dp1=det(q2.x-p2.x,q2.y-p2.y,p1.x-p2.x,p1.y-p2.y); // p1 au-dessus de (p2q2) ?
  const double dq1=det(q2.x-p2.x,q2.y-p2.y,q1.x-p2.x,q1.y-p2.y); // q1 au-dessus de (p2q2) ? 
  
  if((dp2*dq2<0)&&(dp1*dq1<0)) return 1; // cas général

// est-ce que C est dans ]AB[ ?
#define IS_IN(C,A,B)				\
  ( (fmin(A.x,B.x)<C.x)&&(C.x<fmax(A.x,B.x))&&	\
    (fmin(A.y,B.y)<C.y)&&(C.y<fmax(A.y,B.y)) )
 
  if((dp2==0)&&(IS_IN(p2,p1,q1))) return 1; // p2 dans ]p1q1[
  if((dq2==0)&&(IS_IN(q2,p1,q1))) return 1; // q2 dans ]p1q1[
  if((dp1==0)&&(IS_IN(p1,p2,q2))) return 1; // p1 dans ]p2q2[
  if((dq1==0)&&(IS_IN(q1,p2,q2))) return 1; // q1 dans ]p2q2[
  return 0;
}
 

int rlt(query* const Q){
/*
  Pour simplifier, on prend V = [0,2p[ x [0,2q[ comme ensemble de
  points de la grille, c'est-à-dire les couples (i,j) avec i et j
  pairs. On en déduit que l'ensemble M des points milieu est [0,2p-1[
  x [0,2q-1[ \ V, c'est-à-dire des points qui ont au moins une
  coordonnées impaires.

  Pour chaque point milieu R, il faut prendre tous les points A
  possibles de la grille et vérifier que le point B symétrique de A
  selon R est dans la grille, que l'arête résultante (le segment [AB])
  soit de la bonne longueur (<=d) et n'intersectent aucune autre arête
  déjà en place dans E.
*/
  switch(Q->code){

  case QUERY_ADJ:
    return load(Q);

  case QUERY_END:
    return free_pos(Q);

  case QUERY_INIT:;
    const int p=Q->param[0];
    const int q=Q->param[1];
    if((p<=0)||(q<=0)) RET_n(0);
    Q->n=p*q;

    int d=Q->param[2];
    if(d<0) d=2*Q->n; // d=+∞
    d=2*d-1; // d=1->1, d=2->3, d=3->5, etc.
    
    const int m=(2*p-1)*(2*q-1)-p*q; // m=nombre d'arêtes=nombre de points dans M
    const int mod2=(~1); // mask (avec &) pour entier pair inférieur

    // t=nombre de points maximum que l'on peut mettre dans L
    int t=(p&1); t=(p==1)? 1 : p*(q-t)/2+t;

    NALLOC(int,M,2*m); // M=liste de points milieu de chaque arête
    NALLOC(int,E,4*m); // E=liste des segments [AB] qui sont des arêtes
    NALLOC(int,L,2*t); // L=liste des points A possibles pour un point R de M
    Q->G=new_graph(Q->n); // NB: les degrés sont nuls

    int i,j,k,u,v;
    int Ri,Rj,i0,i1,j0,j1,di,dj;
    point A,B,E0,E1;
    int nM=0;
    int nE=0;
    int nL;

    // initialise les points de M
    if(d) // si d=0, il faut nM=0
      for(i=0;i<2*p-1;i++)
	for(j=0;j<2*q-1;j++)
	  if((i&1)||(j&1)) M[nM++]=i,M[nM++]=j; // au moins une coord. impaire

    while(nM){ // tant que M n'est pas vide (et d<>0)
      t=(random()%nM)&mod2; // position (paire) aléatoire
      Ri=M[t],Rj=M[t+1]; // R=(Ri,Rj) points aléatoire de M
      M[t+1]=M[--nM],M[t]=M[--nM]; // supprime R de M

      // détermine la liste L des points A=(i,j) possibles pour R. Il
      // faut R milieu de [AB] donc B=2R-A, A et B dans [0,2p-1[ x
      // [0,2q-1[, et aussi que ]AB[ ne soit pas trop grand et ne
      // contienne aucun point entier à part R. Pour cela, il suffit
      // que les coordonnées de (R-A) soient premières entre elles. A
      // cause des symétries, seuls certains points A doivent être
      // testés dans le rectangle [i0,i1]x[j0,j1] comme suit:
      if(Rj<q) j0=0,j1=Rj&mod2; else j0=(Rj+1)&mod2,j1=2*q-2;
      if(Ri<p) i0=0,i1=2*Ri; else i0=2*Ri-2*p+2,i1=2*p-2;

      nL=0;
      for(i=i0;i<=i1;i+=2){
	A.y=i;B.y=2*Ri-i;di=abs(i-Ri);
	for(j=j0;j<=j1;j+=2){
	  dj=abs(j-Rj);
	  if((j==Rj)&&(i>Ri)) continue; // AR même colonne et A au-dessus
	  if((j==Rj)&&(di>1)) continue; // AR même colonne et A,R pas voisins
	  if((i==Ri)&&(dj>1)) continue; // AR même ligne et A,R pas voisins
	  if(imax(di,dj)>d) continue; // [AR] est trop long
	  if(pgcd(di,dj)!=1) continue; // ]AR[ contient un point entier
	  A.x=j;B.x=2*Rj-j;
	  // est-ce que ]AB[ intersecte une arête de E ?
	  for(k=0;k<nE;k+=4){ // pour toutes les arêtes e=(E0,E1) de E
	    E0.y=E[k+0],E0.x=E[k+1];
	    E1.y=E[k+2],E1.x=E[k+3];
	    if(seg_intersection(A,B,E0,E1)) break; // si intersection
	  }
	  if(k<nE) continue; // A n'est pas bon: suivant
	  L[nL++]=i,L[nL++]=j; // ajoute A=(i,j) à L
	}
      }

      // choisit un point A de L, met à jour le degré des sommets A et
      // B, et met l'arête correspondante dans E
      t=(random()%nL)&mod2; // position (paire) aléatoire
      u=(   L[t]/2)*q + (   L[t+1]/2); // u=numéro du sommet de A
      v=(Ri-L[t]/2)*q + (Rj-L[t+1]/2); // v=numéro du sommet de B
      Q->G->d[u]++;
      Q->G->d[v]++;
      E[nE++]=L[t],E[nE++]=L[t+1]; // A
      E[nE++]=2*Ri-L[t],E[nE++]=2*Rj-L[t+1]; // B
    }

    // copie les arêtes de E dans Q->G
    for(u=0;u<Q->G->n;u++) ALLOC(Q->G->L[u],Q->G->d[u]); // alloue les listes
    degres_zero(Q->G); // pour ADD_EDGE()
    for(k=0;k<nE;k+=4){
      u=(E[k+0]/2)*q + (E[k+1]/2); // u=numéro du sommet de A
      v=(E[k+2]/2)*q + (E[k+3]/2); // v=numéro du sommet de B
      ADD_EDGE(Q->G,u,v); // ajoute l'arête u-v
    }

    // fin
    free(M);
    free(E);
    free(L);
    return 0;
  }

  return 1;
}


int kout(query* const Q){
/*
  Utilise i<j.

  Q->rep[i]=tableau des voisins de i. Si Q->rep[i][j]<0, alors c'est
  la fin de la liste.  Attention ! Si i<j, alors Q->rep[i] ne peut pas
  contenir le sommet j (qui a été crée après i). Pour le test
  d'adjacence, il faut donc tester si i est dans Q->rep[j], et pas le
  contraire !  On a forcément k>0.
*/
  const int k=Q->param[1];
  int x,y;
  switch(Q->code){

  case QUERY_END:
    return free_rep(Q);

  case QUERY_INIT:
    SET_n(Q->param[0]);
    ALLOCMAT(Q->rep,Q->n,k);
    NALLOC(int,T,Q->n);
    int r,d,z;

    Q->rep[0][0]=-1;     /* le sommet 0 est seul !*/
    for(x=1;x<Q->n;x++){ /* x=prochain sommet à rajouter */
      r=imin(x,k);     /* le degré de x sera au plus r=min{x,k}>0 */
      d=1+random()%r; /* choisir un degré d pour x: d=[1,...,r] */
      for(y=0;y<x;y++) T[y]=y; /* tableau des voisins possibles */
      r=x;              /* r-1=index du dernier élément de T */
      for(y=0;y<d;y++){ /* choisir d voisins dans T */
	z=random()%r;   /* on choisit un voisin parmi ceux qui restent */
        Q->rep[x][y]=T[z]; /* le y-ème voisin de x est T[z] */
        T[z]=T[--r];    /* enlève T[z], met le dernier élément de T à sa place */
      }
      if(d<k) Q->rep[x][d]=-1; /* arrête la liste des voisins de x */
    }

    free(T);
    return 0;

  case QUERY_ADJ: /* i<j, donc j a été crée après i */
    Q->a=0;
    for(y=0;y<k;y++){
      if(Q->rep[Q->j][y]==Q->i) RET_a(1);
      if(Q->rep[Q->j][y]<0) return 0; // ici Q->a=0
    }
    return 0; // ici Q->a=0
    
  }
  return 1;
}


int expander(query* const Q){
/*
  Q->rep[i]=tableau des k>0 successeurs de i. Il est possible que le même
  voisin apparaisse plusieurs fois dans Q->rep[i].

  Algorithme:
   1. on part du tableau T[i]=i, pour i=0..n-1
   2. on forme le cycle T[0]-T[1]-...-T[n-1]-T[0]
   3. on permute seulement les n-1 premières valeurs de T
   4. on recommence en 2.

  Il faut répéter k fois la ligne 2, cependant on a pas besoin
  d'effectuer la dernière permutation de la ligne 3. Rem: il est
  inutile de permuter les n cases de T, seule les n-1 première cases
  suffisent car il s'agit de permutations circulaires.
*/
  switch(Q->code){

  case QUERY_END:
  case QUERY_ADJ:
    return arboricity(Q);

  case QUERY_INIT:;
    SET_n(Q->param[0]);
    int k=Q->param[1]; if(k<1) RET_n(0); /* k>0 */
    ALLOCMAT(Q->rep,Q->n,k);
    NALLOCZ(int,T,Q->n,_i); /* tableau pour une permutation aléatoire */

    int c,t;
    for(c=0;c<k;){ /* répète k fois, pour chaque cycle */
      for(t=0;t<Q->n;t++) Q->rep[T[t]][c]=T[(t+1)%Q->n]; /* suit et copie le cycle numéro c */
      if(++c<k) Permute(T,Q->n-1); /* permutation aléatoire (inutile le dernier coup) */
      /* seul Q->n-1 éléments ont besoin d'être permutés (permutation circulaire) */
    }
    
    free(T);
    return 0;
  }
  
  return 1;
}


/***********************************

           GRAPHES DEFINIS
        A PARTIR D'UN TABLEAU
          (DE TAILLE FIXEE)

***********************************/


/* codes<0 pour GraphFromArray() */
enum{
  GFA_END   =-1, // fin du graphe
  GFA_CUT   =-2, // fin d'une séquence
  GFA_PATH  =-3, // chemin de sommets consécutifs
  GFA_HAM   =-4, // cycle Hamiltonien
  GFA_STAR  =-5, // étoile
  GFA_WHEEL =-6, // camembert
};


int GraphFromArray(query* const Q,const int *T){
/*
  Fonction d'adjacence générique permettant pour des graphes de petite
  taille et sans paramètre. Plutôt que d'utiliser une matrice ou une
  liste, c'est un tableau T qui définit les adjacences du graphe.

  En gros, GraphFromArray(Q,T) fait Q->a=1 ssi Q->i et Q->j se suivent
  dans le tableau T, c'est-à-dire s'il existe un indice k tel que
  T[k]=Q->i et T[k+1]=Q->j ou le contraire.  Cette fonction suppose
  que Q->i < Q->j. Les arêtes du graphe sont ainsi couvertes par des
  chemins éventuellement non élémentaires. La fonction renvoie
  toujours 0 si tout c'est bien passé et 1 sinon, comme une fonction
  d'ajacence standard.

  Les valeurs négatives du tableau ont des significations
  particulières (voir ci-dessous la signification des codes GFA_xxx).
  Le principe de l'algorithme de décodage est le suivant: on lit le
  tableau T progressivement valeur après valeur, et soit x le dernier
  sommet mémorisé, c'est-à-dire une valeur entière >=0. Si on lit une
  valeur y>=0 alors y est un sommet et le graphe possède l'arête
  x-y. On remplace alors x par y et on continue la lecture de T. Si
  y<0, c'est un code GFA_xxx et le comportement est alors spécifique
  comme décrit ci-dessous:

  GFA_END: fin du tableau et donc du graphe.

  GFA_CUT: fin d'une séquence. Le prochain sommet lu ne sera pas
     connecté au dernier sommet mémorisé du tableau. Sans ce code, la
     prochaine valeur >=0 lue est toujours connectée au dernier sommet
     mémorisé. En fait, toute valeur <0 qui n'est pas un des codes
     reconnus à le même effet que GFA_CUT.

  GFA_PATH: chemin de sommets consécutifs, "a,GFA_PATH,n" représente
      un chemin de n arêtes suivant le sommet a. C'est équivalent à
      "a,a+1,...,a+n,GFA_CUT".

  GFA_HAM: cycle Hamiltonien, équivalent à
  "0,GFA_PATH,Q->n-1,Q->n-1,0,GFA_CUT".

  GFA_STAR: étoile. La séquence "a,GFA_STAR,a_1,...,a_n,GFA_CUT"
      représente une étoile de centre a et de feuilles a_1,...,a_n.
      C'est équivalent à "a,a_1,a,a_2,a,...,a,a_n,GFA_CUT". On peut
      remplacer le GFA_CUT par n'importe qu'elle valeur négative qui
      peut ainsi être enchaînée.

  GFA_WHEEL: comme GFA_STAR, sauf qu'en plus les sommets a_i,a_{i+1}
      sont adjacents. Attention ! a_n et a_1 ne sont pas adjacents.

  Une bonne stratégie pour trouver un code assez court pour un graphe
  (c'est-à-dire un tableau T de petite taille) est de l'éplucher par
  degré maximum jusqu'à obtenir des chemins ou cycles (non forcément
  simples). Les sommets ainsi épluchés peuvent être codés par GFA_STAR
  ou GFA_WHEEL, les chemins et de cycles par GFA_PATH.

    Ex:        5-6-7
               |  /      -> code = 0,GFA_PATH,10,GFA_CUT,4,8,GFA_END
       0-1-2-3-4-8-9-10

  Les sommets isolés n'ont pas à être codés. Cependant, pour
  déterminer Q->n il est nécessaire que le plus grand sommet du graphe
  figure dans T. Si cela n'est pas le cas (c'est-à-dire que le sommet
  n-1 se trouve être un sommet isolé), on peut ajouter la séquence
  "n-1,GFA_CUT," par exemple en début du tableau T.

*/
  const int i=Q->i;
  const int j=Q->j;
  int k,a,n;
  switch(Q->code){

  case QUERY_INIT:
    for(k=n=0;T[k]!=GFA_END;k++) n=imax(n,T[k]);
    RET_n(n+1);

  case QUERY_ADJ:
    for(k=0;;){
      switch(T[k]){
	  
      case GFA_END:
	RET_a(0); /* adjacence non trouvée */
	
      case GFA_PATH:
	a=T[k-1];
	n=T[++k];
	if((a<=i)&&(i<a+n)&&(j==i+1)) RET_a(1); /* suppose que i<j */
	k++;
	continue;
	
      case GFA_HAM:
	if(j==i+1) RET_a(1); /* suppose que i<j */
	if((i==0)&&(j==(Q->n-1))) RET_a(1);
	k++;
	continue;
	  
      case GFA_STAR:
	a=T[k-1];
	while(T[++k]>=0){
	  if((i==a)&&(j==T[k])) RET_a(1);
	  if((j==a)&&(i==T[k])) RET_a(1);
	}
	continue;
	  
      case GFA_WHEEL:
	a=n=T[k-1];
	while(T[++k]>=0){
	  if(((i==a)||(i==n))&&(j==T[k])) RET_a(1);
	  if(((j==a)||(j==n))&&(i==T[k])) RET_a(1);
	  n=T[k];
	}
	continue;
	
      default:
	/* ici T[k+1] existe car T[k]<>GFA_END. Si T[k]<0 alors cela
	   aura le même effet que de lire GFA_CUT, à savoir le lien
	   T[k]-T[k+1] est coupé puisque i et j sont >0. */
	if((i==T[k])&&(j==T[k+1])) RET_a(1);
	if((j==T[k])&&(i==T[k+1])) RET_a(1);
	k++;
	continue;
      }
    }
  }

  return 1;
}


int tutte(query* const Q){
  const int T[]={
    0,GFA_PATH,8,
    4,8,9,3,9,10,11,2,
    11,12,1,12,13,14,10,
    14,7,6,15,13,GFA_CUT,
    0,16,GFA_PATH,7,
    23,19,23,24,18,24,25,26,17,
    26,27,16,27,28,29,25,
    29,22,21,30,28,GFA_CUT,
    0,31,GFA_PATH,7,
    38,34,38,39,33,39,40,41,32,
    41,42,31,42,43,44,40,
    44,37,36,45,43,GFA_CUT,
    15,20,GFA_CUT,
    30,35,GFA_CUT,45,5,GFA_END};
  return GraphFromArray(Q,T);
}


int icosahedron(query* const Q){
  const int T[]={
    5,0,1,2,3,4,5,6,7,8,
    9,10,11,6,4,7,3,8,2,
    9,1,10,5,1,2,0,3,4,0,
    5,6,10,9,11,9,11,8,11,7,GFA_END};
  return GraphFromArray(Q,T);
}


int rdodecahedron(query* const Q){
  const int T[]={
    0,1,2,3,4,5,0,6,
    7,1,7,8,2,8,9,3,9,
    10,4,10,11,5,11,6,
    7,12,9,12,11,GFA_CUT,
    0,13,2,13,4,GFA_END};
  return GraphFromArray(Q,T);
}


int cuboctahedron(query* const Q){
  const int T[]={
    0,1,2,3,4,5,0,6,7,2,
    8,9,5,6,10,7,1,11,0,GFA_CUT,
    3,10,4,9,11,8,3,GFA_END};
  return GraphFromArray(Q,T);
}


int herschel(query* const Q){
  const int T[]={
    0,GFA_PATH,10,
    10,3,2,9,0,7,6,1,4,GFA_CUT,
    8,5,10,GFA_END};
  return GraphFromArray(Q,T);
}


int goldner_harary(query* const Q){
  const int T[]={
    0,1,2,0,3,4,1,5,6,7,2,8,3,
    1,6,2,3,9,6,10,1,9,2,10,GFA_CUT,
    4,9,5,GFA_CUT,
    8,9,7,GFA_END};
  return GraphFromArray(Q,T);
}


int triplex(query* const Q){
/* proche d'un flower_snark 3 */
  const int T[]={
    0,GFA_PATH,8,8,0,
    9,GFA_STAR,0,3,6,GFA_CUT,
    10,GFA_STAR,1,4,7,GFA_CUT,
    11,GFA_STAR,2,5,8,GFA_END};
  return GraphFromArray(Q,T);
}


int jaws(query* const Q){
  const int T[]={
    0,GFA_PATH,15,15,0,5,6,7,2,1,16,14,15,10,
    9,8,13,12,17,3,4,18,11,10,9,19,6,
    19,18,GFA_CUT,16,17,GFA_END};
  return GraphFromArray(Q,T);
}


int starfish(query* const Q){
/* proche d'un flower_snark 5 */
  const int T[]={
    0,1,2,3,4,0,GFA_CUT,
    5,6,7,8,9,5,GFA_CUT,
    10,11,12,13,14,10,GFA_CUT,
    15,GFA_STAR,0,5,10,GFA_CUT, 
    16,GFA_STAR,1,6,11,GFA_CUT, 
    17,GFA_STAR,2,7,12,GFA_CUT,
    18,GFA_STAR,3,8,13,GFA_CUT,
    19,GFA_STAR,4,9,14,GFA_END};
  return GraphFromArray(Q,T);
}


int fritsch(query* const Q){
  const int T[]={GFA_HAM,0,4,6,3,7,2,4,GFA_CUT,1,8,1,7,6,8,5,0,2,GFA_END};
  return GraphFromArray(Q,T);
}


int soifer(query* const Q){
  const int T[]={GFA_HAM,0,6,4,2,6,0,7,5,8,1,3,8,5,3,GFA_END};
  return GraphFromArray(Q,T);
}


int poussin(query* const Q){
  const int T[]={
    GFA_HAM,0,2,4,6,8,10,12,14,2,
    0,13,11,9,5,8,12,7,3,14,7,GFA_CUT,
    4,1,9,0,11,GFA_CUT,1,5,GFA_CUT,3,6,GFA_END};
  return GraphFromArray(Q,T);
}


int errera(query* const Q){
  const int T[]={
    0,GFA_STAR,10,9,5,15,11,16,GFA_CUT,
    1,GFA_STAR,6,16,12,13,4,7,GFA_CUT,
    2,GFA_WHEEL,6,7,8,9,10,GFA_CUT,
    3,GFA_WHEEL,11,12,13,14,15,11,GFA_CUT,
    4,GFA_STAR,7,8,5,14,13,GFA_CUT,
    5,GFA_STAR,8,9,14,15,GFA_CUT,
    16,GFA_STAR,6,11,12,10,GFA_END};
  return GraphFromArray(Q,T);
}


int kittell(query* const Q){
  const int T[]={
    0,GFA_PATH,22,
    10,12,14,5,13,11,3,10,2,9,1,6,0,
    5,16,6,17,7,1,8,22,20,18,7,GFA_CUT,
    5,15,21,14,22,9,12,GFA_CUT,
    9,14,GFA_CUT,18,8,20,GFA_CUT,
    19,GFA_STAR,15,16,17,21,GFA_CUT,
    2,0,3,0,4,13,4,11,GFA_END};
  return GraphFromArray(Q,T);
}


int frucht(query* const Q){
  const int T[]={GFA_HAM,0,4,5,3,2,7,6,8,9,11,10,1,GFA_END};
  return GraphFromArray(Q,T);
}


int moser(query* const Q){
  const int T[]={0,1,2,3,0,4,5,6,0,GFA_CUT,1,3,2,5,4,6,GFA_END};
  return GraphFromArray(Q,T);
}


int markstrom(query* const Q){
  const int T[]={
    8,0,GFA_PATH,12,5,13,GFA_PATH,7,
    3,21,22,23,0,23,1,2,21,22,18,GFA_CUT,
    15,20,14,13,4,GFA_CUT,
    6,12,7,GFA_CUT,17,19,GFA_CUT,
    9,11,10,16,GFA_END};
  return GraphFromArray(Q,T);
}


int robertson(query* const Q){
  const int T[]={GFA_HAM,0,4,8,13,17,2,6,10,14,0,
    1,9,16,5,12,1,GFA_CUT,3,11,18,7,15,3,GFA_END};
  return GraphFromArray(Q,T);
}


int headwood4(query* const Q){
  const int T[]={
    0,GFA_WHEEL,1,5,4,3,13,12,2,GFA_CUT,
    1,GFA_WHEEL,5,6,7,8,2,GFA_CUT,
    2,GFA_WHEEL,8,9,10,11,GFA_CUT,
    16,GFA_WHEEL,4,15,17,19,5,GFA_CUT,
    12,GFA_WHEEL,11,20,19,18,GFA_CUT,
    14,15,3,14,13,18,17,14,18,GFA_CUT,
    22,GFA_WHEEL,6,21,24,23,7,GFA_CUT,
    5,20,6,11,21,10,24,9,23,8,GFA_END};
  return GraphFromArray(Q,T);
}


int wiener_araya(query* const Q){
  const int T[]={
    0,GFA_STAR,1,4,12,15,GFA_CUT,
    1,GFA_PATH,39,
    41,GFA_STAR,20,23,36,GFA_CUT,
    18,36,35,16,17,1,2,19,GFA_CUT,
    3,21,22,5,4,8,7,26,27,9,10,29,28,39,25,24,6,GFA_CUT,
    8,12,11,31,32,13,14,34,33,37,38,23,GFA_CUT,
    30,40,33,GFA_END};
  return GraphFromArray(Q,T);
}


int zamfirescu(query* const Q){
  const int T[]={
    0,GFA_PATH,47,
    0,9,5,1,17,1,2,19,18,22,21,20,39,3,4,41,
    40,47,43,42,6,7,44,31,32,28,27,34,33,
    45,46,38,37,21,GFA_CUT,
    8,30,29,10,11,27,26,13,14,24,25,35,
    36,23,16,15,0,12,GFA_END};
  return GraphFromArray(Q,T);
}


int hatzel(query* const Q){
  const int T[]={
    0,GFA_PATH,55,
    46,30,31,11,12,8,7,26,27,9,10,29,28,
    47,48,44,45,38,39,43,56,52,51,23,22,
    5,6,24,25,50,49,56,GFA_CUT,
    8,4,0,12,13,32,33,37,36,40,41,42,53,
    54,20,21,3,2,19,18,55,41,GFA_CUT,
    14,34,35,16,17,1,0,15,GFA_END};
  return GraphFromArray(Q,T);
}


int harborth(query* const Q){
  const int T[]={
    0,GFA_PATH,19,GFA_CUT,19,
    0,20,1,21,2,22,3,23,4,24,5,
    25,6,26,7,27,8,28,9,29,10,
    30,11,31,12,32,13,33,14,34,15,
    35,16,36,17,37,18,38,19,39,0,
    39,38,40,41,37,36,35,42,41,43,
    44,45,21,20,45,43,40,39,GFA_CUT,
    46,44,22,23,24,46,25,26,27,
    47,48,28,29,48,49,47,46,GFA_CUT,
    50,49,51,30,31,51,50,32,33,34,
    42,50,GFA_END};
  return GraphFromArray(Q,T);
}


int doily(query* const Q){
  const int T[]={
    9,0,GFA_PATH,9,GFA_CUT,
    9,11,3,13,7,10,1,12,5,14,9,GFA_CUT,
    0,13,5,GFA_CUT,2,14,7,GFA_CUT,
    4,10,9,GFA_CUT,6,11,1,GFA_CUT,
    8,12,3,GFA_END};
  return GraphFromArray(Q,T);
}


int bidiakis(query* const Q){
  const int T[]={
    0,GFA_HAM,
    0,6,GFA_CUT,1,5,GFA_CUT,11,7,GFA_CUT,
    10,2,GFA_CUT,9,3,GFA_CUT,8,4,GFA_CUT,
    GFA_END};
  return GraphFromArray(Q,T);
}


int cricket(query* const Q){
  const int T[]={0,1,2,3,1,4,GFA_END};
  return GraphFromArray(Q,T);
}


int moth(query* const Q){
  const int T[]={0,1,2,3,4,1,5,1,3,GFA_END};
  return GraphFromArray(Q,T);
}


int suzuki(query* const Q){
  const int T[]={
    0,GFA_WHEEL,1,2,3,4,5,6,7,8,GFA_CUT,
    8,2,4,6,8,
    1,3,5,7,1,
    5,6,10,1,9,
    4,3,9,10,7,
    8,10,5,9,2,
    GFA_END};
  return GraphFromArray(Q,T);
}


int bull(query* const Q){
  const int T[]={0,1,2,3,4,3,1,GFA_END};
  return GraphFromArray(Q,T);
}


int hgraph(query* const Q){
  const int T[]={0,1,2,3,GFA_CUT,4,1,2,5,GFA_END};
  return GraphFromArray(Q,T);
}


int rgraph(query* const Q){
  const int T[]={0,1,2,3,4,1,5,GFA_END};
  return GraphFromArray(Q,T);
}


/********************

  OPERATEURS UNAIRES

  C'est un graphe comme les autres (défini par la fonction d'adjacence
  adj), sauf qu'il prend comme paramètre un autre graphe. En général,
  les opérateurs doivent mettre POS=0, sinon les positions XPOS/YPOS
  pourraient ne pas être définies pour certains indices.

  - apex(Q): ajoute des sommets universels
  - star(Q): ajoute des sommets de degré 1 à un graphe.

*********************/


int apex(query* const Q){
/*
  [NE MARCHE PLUS DEPUIS v4.5 -- A FINIR]

  Les sommets de 0..N0-1 sont ceux du graphe initial. Ceux de numéro
  ≥ N0 sont les apices (donc de degré N0). On se sert du fait que
  i<j. Désactive POS.
*/
  query* const P=Q->query;
  return P->adj(P);
}


int star(query* const Q){
/*
  [NE MARCHE PLUS DEPUIS v4.5 -- A FINIR]

  Les sommets de 0...N0-1 sont ceux du graphe initial. Ceux de numéro
  ≥ N0 sont ceux de degré 1.  Se sert du fait que i<j. Désactive POS.
*/
  query* const P=Q->query;
  return P->adj(P);
}


/***********************************

        FONCTIONS DU PROGRAMME
              PRINCIPAL

***********************************/

int InitVertex(const int n,double p){
/*
  Remplit le tableau V[i] donnant l'étiquette finale du sommet i et
  supprime certains des n sommets suivant la valeur p. Si p<0, alors
  on supprime exactement |p| sommets. Utilise aussi SHIFT.  Si PERMUTE
  est vrai V[] est remplit d'une permutation aléatoire de
  SHIFT+[0,n[. Si V[i]=-1 cela signifie que i a été supprimé (p).  La
  fonction retourne le nombre de sommets final du graphe, c'est-à-dire
  le nombre d'étiquettes >=0. Si k sommets ont été supprimés, alors
  les valeurs de V[] sont dans SHIFT+[0,n-k[.

  Initialise aussi le tableau VF[j], avec j=0...n-k, de sorte que
  VF[j]=i si VF[i] est le j-ème sommet non supprimé. Il est important
  que VF[] ait une taille de Q->n au départ. Un realloc() le
  redimensionne plus tard dans la fonction.
*/

  int i,j,k,r; /* r=n-(nb de sommets supprimés) */
  long seuil;

  /* supprime les sommets */
  if(p<0){ /* ici on en supprime exactement |p| sommets */
    for(i=0;i<n;i++) VF[i]=i; /* on va se servir temporairement de VF */
    r=-(int)p; /* les r premières valeurs de VF seront les sommets à supprimer */
    r=imin(r,n); /* r ne doit pas dépasser n */
    for(i=0;i<r;i++){
      j=i+(random()%(n-i));
      SWAP(VF[i],VF[j],k);
    }
    for(i=0;i<r;i++) V[VF[i]]=-1; /* on supprime ces r sommets */
    for(i=r=0;i<n;i++) /* on remplit V et VF */
      if(V[i]>=0) { VF[r]=i; V[i]=r++; }
  }
  else{ /* ici on supprime chacun des n sommets avec proba p */
    seuil=(double)p*(double)RAND_MAX;
    for(i=r=0;i<n;i++)
      if(random()<seuil) V[i]=-1;
      else { VF[r]=i; V[i]=r++; }
  } /* dans les deux cas, r=nb de sommets restant */

  /* réajuste le tableau VF à la taille minimum */
  REALLOC(VF,r);

  if(PERMUTE) Permute(V,n);

  /* ne rien faire si SHIFT=0 */
  if(SHIFT>0)
    for(i=0;i<r;i++)
      V[VF[i]] += SHIFT;

  return r;
}


void ScanINC(int *dmax,int *dmin,int *m){
/*
  Calcule, en fonction des tableaux INC[] et VF[], le degré max, le
  degré min et le nombre d'arêtes du graphe (final).
*/
  int d,i;
  *m=*dmax=0;
  *dmin=INT_MAX;
  if(NF<=0) *dmin=0;
  for(i=0;i<NF;i++){
    d=INC[VF[i]]; /* d=degré du i-ème sommet existant */
    *m += d;
    if(d>*dmax) *dmax=d;
    if(d<*dmin) *dmin=d;
  }
  *m >>= 1;
  return;
}


char *MakeCMD(char *s,const int deb,const int fin){
/*
  Routine permettant de recomposer la ligne de commande, en
  particulier pour obtenir le nom du graphe généré (avec les
  options). On ajoute à la fin de la chaîne s les arguments ARGV[i]
  pour i∈[deb,fin[. Si s=NULL alors un pointeur statique sur la chaîne
  est renvoyée, pointeur qui n'a donc pas besoin d'être libéré.

  Si un argument comprend un espace, il est alors parenthésé par '...'
  de façon à être interprété comme un seul argument. Les arguments
  sont séparés par un espace. Le dernier argument est toujours suivi
  d'un espace, éventuellement inutile.
*/
  static char r[CMDMAX];
  if(s==NULL){ s=r; VIDE(s); }
  
  int i;
  for(i=deb;i<fin;i++)
    if(index(ARGV[i],' ')) /* si l'argument est en plusieurs mots */
      strcat(strcat(strcat(s,"'"),ARGV[i]),"' ");
    else strcat(strcat(s,ARGV[i])," ");

  return s;
}


/* pour DateHeure() */
#define DATE_FORMAT "%d/%m/%Y - %Hh%M'%S"
#define SIZE_DATE 22 // avec le 0 final


char *DateHeure(void){
/*
  Renvoie la date et l'heure courante. Il n'est pas (et il ne faut pas
  !) faire de free() sur le pointeur retourné.
*/
  static char date[SIZE_DATE];
  time_t t=time(NULL);
  struct tm *tm=localtime(&t);
  strftime(date,SIZE_DATE*sizeof(char),DATE_FORMAT,tm);

  date[SIZE_DATE-1]=0;
  return date;
}


void Header(const int c){
/*
  Affiche et calcule le préambule (date, nom du graphe, nb de
  sommets). Suivant la valeur de c, le nombre d'arêtes est donnée.  Si
  bit-0 de c=1, alors l'entête de base est affichée.  Si bit-1 de c=1,
  alors on affiche le nombre d'arêtes, le degré min et max.

*/

  /* affichage de la date, nom du graphe, ligne de commande et de n */
  if(c&1){
    printf("//\n");
    printf("// %s - seed=%u\n",DateHeure(),SEED);
    printf("// %s\n",MakeCMD(NULL,0,ARGC)); // ne pas libérer ce pointeur
    printf("// n=%i",NF);
  }

  /* affichage du nombre d'arêtes, maxdeg et mindeg */
  if(c&2){
    int maxdeg,mindeg,nbedges;
    ScanINC(&maxdeg,&mindeg,&nbedges);
    if(!(c&1)) printf("//\n//");
    printf(" m=%i",nbedges);
    printf(" maxdeg=%i",maxdeg);
    printf(" mindeg=%i",mindeg);
  }

  if(c) printf("\n//\n");
  return;
}


static inline char *ComputeName(query* const Q){
/*
  Routine servant plusieurs fois dans Out(Q) pour l'affichage à la
  volée du sommet Q->i. On renvoie le nom du sommet à afficher dans
  Q->name. Par défaut c'est V[Q->i]. Dans le cas du format standard on
  affiche le nom donné par Q->adj(Q) et QUERY_NAME, ou bien les
  coordonnées suivant le cas (|LABEL|==3).

  NB: Dans le cas du format dot, il faut laisser les indices (donc
  V[Q->i]) car les "label" des sommets sont codés seulement à la fin
  dans l'attribut "label".
*/
  sprintf(Q->name,"%i",(V==NULL)?Q->i:V[Q->i]); /* par défaut name="V[i]" */
  if(FORMAT==F_standard){
    if(abs(LABEL)==1){ Q->code=QUERY_NAME; Q->adj(Q); }
    if(POS && abs(LABEL)==3) sprintf(Q->name,"(%g,%g)",XPOS[Q->i],YPOS[Q->i]);
  }
  return Q->name;
}


void Out(query* const Q){
/*
  Affiche l'arête i-j suivant le format FORMAT.

  Si HEADER=1, alors Out() doit se charger d'afficher l'en-tête.

  Si CHECK>0 alors Out() doit se charger de créer et de remplir la
  liste d'adjacence du graphe GF et de déterminer son nombre de
  sommets NF. Pour cela une liste (de type "list") est créee et
  progressivement rempli. A la fin, on calcule GF avec List2Graph().

  Variables globales modifiées:
  - N, GF, NF, VF, NAME
  - CAPTION, NPAL, PALETTE
  - VERTEX0, VERTEXN
  - USERDOT

  Autres variables globales utilisées:
  - CHECK, FORMAT, ROUND, WIDTH
  - HEADER, DIRECTED, VCOLOR
  - XMIN, XMAX, YMIN, YMAX
  - VSIZEK, VSIZESTD, VSIZEXY
  - POS, LABEL, LEN, XPOS, YPOS
  - PARAM_PAL, CPARAM
  - COLORCHAR, COLORBASE
*/
  int x,y,z;
  double w;
  
  /* variables qui conservent leurs valeurs d'un appel à l'autre */
  static list *L0,*L1,*L2; /* pour -check: tête, dernier, avant-dernier */
  static int cpt;  /* compteur d'arêtes ou de sommets isolés par ligne */
  static int last; /* extrémité de la dernière arête affichée */
  static char *sep1; /* séparateur d'arête: "" (pour standard) ou ";" (dot) */
  static char *sep2; /* autre séparateur */
  static char *sepe; /* caractère pour une arête: "-" (pour standard) ou "--" (dot) */

  /* format de précision par défaut pour l'affichage de XPOS/YPOS dans
     le format F_xy, soient 6 digits par défaut */
  char fmt[17]="%lf %lf\n";
  const int i=Q->i;
  int j=Q->j;

  switch(Q->code){


  case QUERY_INIT:
    //------------------------------------
    // initialise la fonction
    //------------------------------------

    cpt=0;
    last=-1;
    if(CHECK) L0=L1=L2=new_list(); /* initialise la liste */

    switch(FORMAT){

    case F_standard:
      if(HEADER) Header(1);
      sep1="";
      sep2=" ";
      sepe=(DIRECTED)?"->":"-";
      return;

    case F_userdot:
      USERDOT.i=USERDOT.j=-1;
      USERDOT.adj=NULL;
      USERDOT.ptr=NULL;
      // puis pareil que F_dot
    case F_dot:
      if(HEADER) Header(1);
      printf("%sgraph {\n",(DIRECTED)?"di":"");
      if(CAPTION){
	printf("\tgraph [label=\"%s\"];\n",CAPTION);
	free(CAPTION);
	CAPTION=NULL;
      }
      if(DOTSCALE) printf("\tgraph [scale=\"%s\"];\n",DOTSCALE);
      if(VCOLOR&0x10){ /* "list" */
	printf("\tgraph [layout=nop, splines=line];\n");
	printf("\tnode [height=1.0, width=0.4, style=filled, shape=rectangle];\n");
	return;
      }
      if(POS){
	w=PosAspect(Q);
	/* layout=nop: pour dire à dot de ne pas re-calculer les positions */
	printf("\tgraph [layout=nop, splines=line, bb=\"%.2lf,%.2lf,%.2lf,%.2lf\"];\n",
	       w*XMIN-2*VSIZEK*VSIZEXY,w*YMIN-2*VSIZEK*VSIZEXY,
	       w*XMAX+2*VSIZEK*VSIZEXY,w*YMAX+2*VSIZEK*VSIZEXY);
	/* si on ne met pas le "2*" on a des sommets tronqués ... */

	/* affiche éventuellement une sous-grille carrée, avant le
	   graphe pour qu'elle apparaisse en arrière-plan. Le nombre
	   de points vaut XYgrid, Q->n ou sqrt(Q->n) */
	z=Q->n; // 1er sommet disponible après ceux de G
	if(XYgrid<0) XYgrid=(XYtype==XY_PERM)? Q->n : 1+(int)(sqrt(Q->n));
	if(XYgrid>1){ /* il faut au moins deux lignes et colonnes */
	  printf("\n\tsubgraph {\n");
	  printf("\t\tnode [label=\"\", height=0, width=0, shape=point, style=invis];\n");
	  printf("\t\tedge [color=gray, penwidth=0.6];");
	  const double rx=(double)(XMAX-XMIN)/(double)(XYgrid-1); /* pas de la grille en X */
	  const double ry=(double)(YMAX-YMIN)/(double)(XYgrid-1); /* pas de la grille en Y */
	  for(y=0;y<XYgrid;y++){
	      for(x=0;x<XYgrid;x++){
		/* affiche le sommet courant (x,y) ainsi que deux
		   arêtes incidentes vers les voisins (x+1,y) et
		   (x,y+1), s'ils existent */
		printf("\n\t\t%i [pos=\"%lf,%lf\"];",z,
		       w*(XMIN+(double)x*rx),
		       w*(YMIN+(double)y*ry));
		if(x+1<XYgrid) printf("\t%i--%i;",z,z+1); /* arête vers (x+1,y) */
		if(y+1<XYgrid) printf("\t%i--%i;",z,z+XYgrid); /* arête vers (x,y+1) */
		z++; /* prochain sommet */
	      }
	    }
	  printf("\n\t}\n\n");
	}
	if(XYzero){ /* ajoute l'origine (0,0) */
	  printf("%s\tsubgraph {\n",(XYgrid>1)?"":"\n");
	  printf("\t\tzero [label=\"\", pos=\"0,0\", shape=circle, ");
	  printf("height=0.1, width=0.1, color=red, fixedsize=true];\n");
	  printf("u\n;");
	  printf("\n\t}\n\n");
	}
	if(XYborder){ /* ajoute un bord */
	  char *s="rectangle";
	  if((XYtype=XY_DISK)||(XYtype==XY_CIRCLE)||(XYtype==XY_HYPER)) s="circle";
	  printf("%s\tsubgraph {\n",(XYgrid>1)||(XYzero)?"":"\n");
	  printf("\t\tbord [label=\"\", pos=\"0,0\", shape=%s,",s);
	  printf(" height=%lf, width=%lf, color=blue, fixedsize=true];\n",1.0,1.0);
	  printf("u;\n");
	  printf("\n\t}\n\n");
	}
      }
      printf("\tnode [");
      if(LABEL==0 || abs(LABEL)==3) printf("label=\"\", shape=point, "); /* sans label */
      w=POS?VSIZESTD:VSIZEXY; /* taille des sommets */
      printf("height=%.2lf, width=%.2lf];\n",w,w);
      if(strcmp(DOTFILTER,"neato")==0) printf("\tedge [len=%.2lf];\n",LEN);
      sep1=";";
      sep2="; ";
      sepe=(DIRECTED)?"->":"--";
      return;

    case F_html:
      printf("<!DOCTYPE html>\n");
      printf("<html lang=\"fr\">\n");
      printf("<head>\n");
      printf("<meta charset=\"utf-8\" />\n");
      printf("<title>%s</title>\n",MakeCMD(NULL,0,ARGC));
      printf("<script type=\"text/javascript\" src=\""URL_vis_js1"\"></script>\n");
      printf("<script type=\"text/javascript\" src=\""URL_vis_js2"\"></script>\n");
      printf("<style type=\"text/css\">\n");
      printf("html, body { padding: 0; margin: 0; width: 100%%; height: 100%% }\n");
      printf("#G { width: 100%%; height: 100%% }\n");
      printf("</style>\n");
      printf("</head>\n");
      printf("<body>\n");
      printf("<div id=G></div>\n");
      printf("<script>\n");
      printf("var E = new vis.DataSet([\n"); // pour les arêtes
      return;
    }
    return; // termine QUERY_INIT


  case QUERY_END:
    //----------------------------------
    // termine la fonction
    //----------------------------------

    if(CHECK){ /* on crée GF en fonction de la liste L0 */
      free(L1); /* supprime la sentienelle (dernier élément) de L0 */
      if(L0==L1) GF=new_graph(0); /* si premier = dernier alors graphe vide */
      else{
	L2->next=NULL; /* coupe la liste à l'avant dernier élément qui a été supprimer */
	GF=List2Graph(L0,4); /* initialise GF et NF */
      }
      NF=GF->n;
    }
    
    switch(FORMAT){

    case F_standard:
    case F_userdot:
    case F_dot:
      if((VCOLOR&0x10)==0){ /* court-circuite l'affichage du graphe si "-vcolor list" */
	if(cpt>0) printf("%s\n",sep1); /* newline si fini avant la fin de ligne */
	if(FORMAT==F_standard){ /* fin si format standard */
	  if(HEADER) Header(2);
	  return;
	}

	if(POS||LABEL){ // affiche les labels des sommets (voir aussi ComputeName)
	  // NB: si POS=1, alors XPOS/YPOS existent forcément
	  w=PosAspect(Q);
	  printf("\n");
	  for(y=0;y<NF;y++){
	    x=VF[y]; /* le sommet x existe */
	    printf("%i [",V[x]);
	    if(POS) printf("pos=\"%lf,%lf\"",w*XPOS[x],w*YPOS[x]);
	    if(LABEL){
	      VIDE(Q->name);
	      if(abs(LABEL)==1){
		Q->code=QUERY_NAME;
		Q->i=x;
		Q->adj(Q); /* calcule Q->name */
	      }
	      if(POS && abs(LABEL)==3) sprintf(Q->name,"(%g,%g)",XPOS[x],YPOS[x]);
	      printf("%s%slabel=\"%s\"",
		     POS?", ":"",
		     (LABEL<0)? "x": "",
		     NONVIDE(Q->name)? Q->name : "\\N");
	    }
	    printf("];\n");
	  }
	}
	
	if(VSIZE&&(NF>0)){ /* taille en fonction du degré */
	  double alpha,smin;
	  ScanINC(&x,&z,&y); /* x=degmax, z=degmin */
	  smin=POS?VSIZESTD:VSIZEXY;
	  alpha=(x==z)? 0 : smin*(VSIZEK-1)/((double)(x-z));
	  printf("\n");
	  for(y=0;y<NF;y++){
	    x=VF[y]; /* le sommet x existe */
	    w=smin + alpha*(INC[x]-z);
	    printf("%i [height=%.2lf, width=%.2lf];\n",V[x],w,w);
	  }
	}
      } /* fin du if((VCOLOR&0x10)==0) ... */

      if(VCOLOR&&(NF>0)){ /* couleur des sommets */
	color c={0,0,0},*P; /* couleur noire par défaut */
	int *D;

	if(VCOLOR&0x8){ /* option "pal" on initialise la PALETTE */
	  NPAL=(int)strlen(PARAM_PAL);
	  if(NPAL==1) { /* PARAM_PAL="x" alors PARAM_PAL="xx" */
	    PARAM_PAL[1]=PARAM_PAL[0];
	    PARAM_PAL[NPAL=2]='\0';
	  }
	  ALLOC(PALETTE,NPAL); /* PALETTE = tableau de NPAL "color" */
	  for(y=z=0;y<NPAL;y++){ /* z=prochaine entrée libre dans PALETTE */
	    x=(int)(index(COLORCHAR,PARAM_PAL[y])-COLORCHAR); /* x=index de la couleur */
	    x/=sizeof(*COLORCHAR); /* normalement inutile */
	    if((0<=x)&&(x<COLORNB)) PALETTE[z++]=COLORBASE[x]; /* on a trouvé la couleur */
	  }
	  if(z<2){ PALETTE[0]=PALETTE[1]=c; z=2; } /* si pas trouvé au moins 2 couleurs */
	  NPAL=z; /* NPAL=nb de couleurs trouvées */
	}
	
	if(VCOLOR&0x17){ /* fonction de couleur: 1,2,3,4,5 ou "-vcolor list" */
	  if((VCOLOR&0x7)>2){ /* si 3="degm", 4="randg"  ou 5="kcolor" */
	    if((VCOLOR&0x7)==3){ /* si "degm" */
	      int *T=Prune(GF,NULL);
	      D=GreedyColor(GF,T); /* calcule les couleurs selon l'ordre T */
	      y=1+GF->int1; /* y=nb de couleurs */
	      free(T); /* on a plus besoin de T */
	    }
	    if((VCOLOR&0x7)==4){ /* si "randg" */
	      NALLOCZ(int,T,NF,_i);
	      Permute(T,NF); /* T=permutation aléatoire */
	      D=GreedyColor(GF,T); /* calcule les couleurs selon l'ordre T */
	      y=1+GF->int1; /* y=nb de couleurs */
	      free(T); /* on a plus besoin de T */
	    }
	    if((VCOLOR&0x7)==5){ /* si "kcolor" */
              y=MEM(CPARAM,0,int); /* y=nb de couleur */
	      D=kColor(GF,y);
              if(D==NULL){ ALLOCZ(D,NF,0); y=1; } /* pas trouvé -> une seule couleur */
	    }
	  }
	  else{ /* si 1="deg" ou 2="degr" */
	    ScanINC(&x,&z,&y); /* calcule x=degmax, z=degmin */
	    y=x-z+1; /* y=nb a priori de couleurs nécessaires */
	    ALLOCZ(D,NF,INC[VF[_i]]-z);
	    if((VCOLOR&0x7)==2){ /* si "degr" */
	      int *R=SortInt(D,NULL,Q->n,0,&y,SORT_INC_RANK);
	      /* après SortInt(), y=nb de degré différents */
	      free(D);
	      D=R; /* on substitue R à D */
	    }
	  }
	  /* ici D[x]=indice de la couleur du sommet x, et y=nb de couleurs */
	  P=GradColor(PALETTE,NPAL,y); /* calcule une palette P de y couleurs. NB: ici NPAL>1 */
	  printf("\n");
	  if(VCOLOR&0x10){ /* si "-vcolor list" */
	    for(x=0;x<y;x++){
	      c=P[x];
	      for(z=0;z<COLORNB;z++) /* on cherche c dans COLORBASE */
		if((COLORBASE[z].r==c.r)&&(COLORBASE[z].g==c.g)&&(COLORBASE[z].b==c.b)) break;
	      printf("\t%i [pos=\"%i,0\", color=\"#%.02x%.02x%.02x\", label=\"%c\", fontcolor=%s];\n",
		     x,10+28*x,c.r,c.g,c.b,(z<COLORNB)?COLORCHAR[z]:' ',
		     (c.r+c.g+c.b<150)?"white":"black");
	    }
	  }else{ /* si pas "-vcolor list" */
	    for(x=0;x<NF;x++){
	      c=P[D[x]]; /* c=couleur du degré (ou du rang) de x */
	      printf("%i [style=filled, fillcolor=\"#%02x%02x%02x\"];\n",V[VF[x]],c.r,c.g,c.b);
	    }
	  }

	  free(D);
	  free(P);
	}
	
	if(VCOLOR&0x8){
	  free(PALETTE); /* la PALETTE ne sert plus à rien */
	  PALETTE=COLORBASE; /* remet à la valeur par défaut, au cas où */
	  NPAL=COLORNB; /* taille de la palette par défaut */
	}
      }

      printf("}\n"); /* si F_dot on ferme le "}" */
      if(HEADER) Header(2); /* affiche les arêtes */
      return;
      
    case F_no:
      if(HEADER) Header(3);
      return;

    case F_list:
      if(HEADER) Header(3);
      PrintGraphList(GF);
      return;
      
    case F_matrix:
    case F_smatrix:
      if(HEADER) Header(3);
      PrintGraphMatrix(GF);
      return;

    case F_xy:
      if(HEADER) Header(3);
      if((XPOS==NULL)||(YPOS==NULL)) Erreur(8);
      printf("%i\n",NF); /* nombre de sommets final */
      if(ROUND<DBL_DIG){ /* construit le nouveau format */
	int r=imax(ROUND,0); /* si ROUND<0 -> 0 */
	if(ROUND<10){
	  strcpy(fmt,"%.0*lf %.0*lf\n"); /* remplace '*' par r */
	  fmt[3]=fmt[10]=(char)('0'+r);
	}else{
	  strcpy(fmt,"%.0**lf %.0**lf\n"); /* remplace '**' par r */
	  fmt[3]=fmt[11]=(char)('0'+(r/10));
	  fmt[4]=fmt[12]=(char)('0'+(r%10));
	}
      }
      for(y=0;y<NF;y++){
	x=VF[y]; /* le sommet x existe */
	printf(fmt,XPOS[x],YPOS[x]);
      }
      return;
      
    case F_html:
      printf("]);\n"); // fin de E = {...}
      printf("var V = new vis.DataSet([\n");
      for(y=0;y<NF;y++){
	x=VF[y]; /* le sommet x existe */
	VIDE(Q->name);
	if(LABEL){
	  if(abs(LABEL)==1){
	    Q->code=QUERY_NAME;
	    Q->i=x;
	    Q->adj(Q); /* calcule Q->name */
	  }
	  if(POS && abs(LABEL)==3) sprintf(Q->name,"(%g,%g)",XPOS[x],YPOS[x]);
	}
	if(!NONVIDE(Q->name)) sprintf(Q->name,"%i",x);
	printf("{id: %i, %s: \"%s\"},\n",
	       x,
	       (LABEL)?"label":"title",
	       Q->name);
      }
      printf("]);\n");
      printf("var container = document.getElementById('G');\n");
      printf("var data = { nodes: V, edges: E };\n");
      printf("var options = { interaction: { hover: true }");
      if(LABEL<0) printf(", nodes: { shape: 'dot' }");
      printf(" };\n");
      printf("new vis.Network(container, data, options);\n");
      printf("</script>\n");
      printf("</body>\n");
      printf("</html>\n");
      return;

    default: Erreur(5); /* normalement sert à rien */
    }
    return; // termine QUERY_END 

    
  case QUERY_ISOL: j=-1;
  case QUERY_ADJ:
    //-----------------------------------------
    // affichage à la volée: "i-j", "i" ou "-j"
    //-----------------------------------------

    if(CHECK){
      L1=Insert(L2=L1,i,T_NODE); /* on ajoute i */
      if(j>=0) L1=Insert(L2=L1,j,(DIRECTED)?T_ARC:T_EDGE); /* on ajoute ->j ou -j */
      if(VCOLOR&0x10) return; /* ne rien faire d'autre si "-vcolor list".
				 NB: CHECK>0 dans ce cas */
    }

    if((FORMAT==F_standard)||(FORMAT==F_dot)){
      if(j<0) last=-1; /* sommet isolé, donc last sera différent de i */
      if(i!=last){ Q->i=i; printf("%s%s",(cpt==0)?"":sep2,ComputeName(Q)); }
      last=j; /* si i=last, alors affiche -j ou ->j. Si j<0 alors last<0 */
      if(j>=0){ Q->i=j; printf("%s%s",sepe,ComputeName(Q)); }
      if(++cpt==WIDTH){
	cpt=0; last=-1;
	printf("%s\n",sep1);
      }
      return;
    } /* si format matrix, smatrix etc., ne rien faire, c'est fait pqr QUERY_END */
    
    if(FORMAT==F_userdot){ // NB: ici on suppose que WIDTH=1
      if(j<0){ Q->i=i; printf("%s%s\n",ComputeName(Q),sep1); } // sommet isolé
      else{
	Q->code=QUERY_DOT;
	if((USERDOT.i!=i)||(USERDOT.j!=j)) USERDOT.adj(Q); // l'arête n'a pas été calculée
	// ici l'arête i-j vient d'être calculée, il reste à la dessiner
	USERDOT.adj(Q);
      }
      return; 
    }

    if(FORMAT==F_html){ // NB: ici on suppose que WIDTH=1
      printf("{from: %i, to: %i},\n",i,j);
      return;
    }
    
    return; // termine QUERY_ADJ
  }
  return; // termine Out()
}


static inline double CheckProba(const double p){
  /* Vérifie que p est bien une probabilité */
  if((p<0)||(p>1)) Erreur(38);
  return p;
}


/***************************************

           FONCTIONS LIEES
  A L'ANALYSE DE LA LIGNE DE COMMANDE

***************************************/


void Grep(int i){
/*
  Cherche le mot ARGV[i] dans l'aide contenu dans le source du
  programme, puis termine par exit().

  Plusieurs cas peuvent se présenter:

  Cas 0: si un ARGV[j]="-", pour j=0..i, on le saute car cela ne peut
  être ni une option ni un graphe dont on peut trouver l'aide.
  
  Cas 1: ARGV[i] est une option, ou aucun ARGV[j] avec j<i n'est une
  option.  Alors on affiche l'aide allant de "^....ARGV[i]($| )" à
  "^$".

  Cas 2: ARGV[j] est une option mais pas ARGV[i] avec j<i. Dans ce
  cas, on pose mot=ARGV[j]" "ARGV[j+1]" "..." "ARGV[i]. Alors on
  affiche l'aide allant de "[ ]{7}mot($| )" à "^$" ou "^[ ]{7}-fin"
  avec fin défini de sorte que mot ne soit pas un préfixe.

  En fait on bride la recherche de l'option à ARGV[i-1] ou ARGV[i-2]
  seulement, si bien que j=i, i-1 ou i-2.
*/

  int j,k,t;

  // enlève les arguments "-" de la liste des arguments

  for(t=j=0;t<=i;t++)
    if(strcmp(ARGV[j],"-")) j++;
    else ARGV[j]=ARGV[t],i--;
  j=i;

  DEBUG(
	for(t=0;t<=i;t++) printf("%s ",ARGV[t]);
	printf("\n");
	fflush(stdout);
	);
  
  do{ // calcule j=i, i-1 ou i-2
    if((i>0)&&(*ARGV[i-0]=='-')){ j=i-0; break; }
    if((i>1)&&(*ARGV[i-1]=='-')){ j=i-1; break; }
    if((i>2)&&(*ARGV[i-2]=='-')){ j=i-2; break; }
  }while(0);

  // construit mot

  NALLOC(char,mot,CMDMAX); VIDE(mot);
  for(t=j;t<=i;t++){
    strcat(mot,ARGV[t]);
    if(t<i) strcat(mot," ");
  }

  // construit fin (à faire seulement si j<i)
  // si ARGV[j]="-option abc xy"
  // alors fin="-option ([^a]|a[^b]|ab[^c]abc) ([^x]|x[^y])"

  NALLOC(char,fin,CMDMAX); VIDE(fin);
  if(j<i){
    strcat(fin,ARGV[j]);
    strcat(fin," (");
    t=j+1; k=0;
    while(ARGV[t][k]){
      strncat(fin,ARGV[t],k);
      strcat(fin,"[^");
      strncat(fin,ARGV[t]+k,1);
      strcat(fin,"]|");
      k++;
      if(ARGV[t][k]==0){
	if(t==i) break;
	strcat(fin,ARGV[t]);
	strcat(fin,") (");
	t++; k=0;
      }
    }
    fin[strlen(fin)-1]=')';
  }

  // construit la commande s

  NALLOC(char,s,CMDMAX); VIDE(s);
  strcat(s,"sed -n '/*[#] ###/,/### #/p' ");
  strcat(strcat(s,*ARGV),".c|");
  /* rem: sed -E n'est pas standard */

  if(j==i)
    strcat(strcat(strcat(s,"sed -nE '/^[.]{4}"),mot),"($|[ ])/,/^$/p'|");
  else{
    strcat(strcat(strcat(s,"sed -nE '/^[ ]{7}"),mot),"($|[ ])/,/(^$)|(^[ ]{7}");
    strcat(s,strcat(fin,")/p'|tail -r|sed -n '2,$p'|"));
    // les tail -r permettent de supprimer la dernière ligne
    // le awk est pour enlever éventuellement l'avant dernière ligne "...."
    strcat(s,"awk '{if(NR>1)print $0;else if(!match($0,/^[.]{4}/))print $0;}'|");
    strcat(s,"tail -r|");
  }

  strcat(s,"sed -E 's/^[.]{4}/    /g'");
  strcat(s,"| awk '{n++;print $0;}END{if(!n) ");
  strcat(s,"print\"Erreur : argument incorrect.\\nAide sur ");
  strcat(s,mot);
  strcat(s," non trouvée.\\n\";}'");
  printf("\n");
  system(s);

  if(j<i) printf("\n");
  free(s);
  free(mot);
  free(fin);
  exit(EXIT_SUCCESS);  
}


char *GetArgInc(int *i){
/*
  Retourne ARGV[i], s'il existe, puis incrémente i.  Si l'argument
  n'existe pas ou si ARGV[i]="?", on affiche l'aide en ligne sur le
  dernier argument demandé et l'on s'arrête.
*/
  
  if(((*i)==ARGC)||(strcmp(ARGV[*i],"?")==0)) Grep((*i)-1);
  return ARGV[(*i)++];
}


void NextArg(int *i){
/*
  Vérifie si le prochain argument existe et incrémente i. S'il
  n'existe pas ou si c'est "?", une aide sur le dernier argument est
  affichée comme l'aurait fait GetArgInc(). Cette fonction devrait
  être appelée chaque fois que le prochain argument doit être testé
  avec: if EQUAL("..."). Si on ne le fait pas, alors c'est un
  "Segmentation fault".
*/
  (*i)++;
  GetArgInc(i);
  (*i)--;
  return;
}

void CheckHelp(int *i){
/*
  Incrémente i puis vérifie si l'argument est "?". Si tel est le cas,
  une aide est affichée sur ARGV[i] (avant incrémentation).  Cette
  fonction devrait être typiquement appellée lorsqu'on vérifie l'aide
  pour un graphe sans paramètre. Sinon c'est fait par GetArgInc().
*/
  
  (*i)++;
  if(((*i)!=ARGC)&&(strcmp(ARGV[*i],"?")==0)) Grep((*i)-1);
  return;
}


void Help(int i){
/*
  Affiche:
  - l'aide complète si ARGV[i] est "-help" ou "?", ou
  - les paragraphes contenant ARGV[i].
*/
  
  NALLOC(char,s,CMDMAX); VIDE(s);
  strcat(s,"sed -n '/*[#] ###/,/### #/p' "); /* filtre l'aide */
  strcat(strcat(s,*ARGV),".c | "); /* enlève 1ère et dernière ligne */
  strcat(s,"sed -e 's/\\/\\*[#] ###//g' -e 's/### [#]\\*\\///g' ");
  i++;
  if((i==ARGC)||(ARGV[i][0]=='?'))
    strcat(s,"-e 's/^[.][.][.][.][.]*/    /g'|more"); /* aide complète */
  else{
    strcat(s,"|awk 'BEGIN{x=\".....\"}/^[.]{4}./{x=$0} /");
    strcat(strcat(s,ARGV[i]),"/{if(x!=\".....\")print substr(x,5)}'|sort -u");
  }
  system(s);
  free(s);
  exit(EXIT_SUCCESS);
}


void ListGraph(void){
/*
  Affiche les graphes disponibles et puis termine.
*/

  NALLOC(char,s,CMDMAX); VIDE(s);
  strcat(s,"sed -n '/*[#] ###/,/### #/p' ");
  strcat(strcat(s,*ARGV),".c| ");
  strcat(s,"sed -e 's/\\/\\*[#] ###//g' -e 's/### [#]\\*\\///g'| ");
  strcat(s,"grep '^[.][.][.][.][^-.]'| sed 's/^[.][.][.][.]//g'");
  //printf("%s\n",s);
  system(s);
  free(s);
  exit(EXIT_SUCCESS);
}


void Version(void){
/*
  Affiche la version du programme et puis termine.
*/
  NALLOC(char,s,CMDMAX); VIDE(s);
  strcat(s,"sed -n '/*[#] ###/,/### #/p' ");
  strcat(strcat(s,*ARGV),".c| ");
  strcat(s,"sed -n 's/.*[-] v\\([0-9]*[.][0-9]*\\) [-].*/\\1/p'");
  system(s);
  free(s);
  exit(EXIT_SUCCESS);
}


void PipeDot(int j){
/*
  Gère l'option "-format dot<type>".

  Rem: on pourrait utiliser popen() au lieu de réécrire la ligne de
  commande et de lancer system().
*/

  CheckHelp(&j);j--; /* vérifie l'aide, ici ARGV[j]="dot<type>" */
  char *type=strdup(ARGV[j]+3); /* type=<type> */
  ARGV[j][3]='\0'; /* ARGV[j]="dot" plutôt que "dot<type>" */

  /*
    On réécrit la ligne de commande:
    1. en remplaçant "-format dot<type>" par "-format dot"
    2. puis en ajoutant à la fin: "| dot -T<type> -K <filter>"
  */
  
  char *s=MakeCMD(NULL,0,ARGC); // ne pas libérer ce pointeur
  strcat(strcat(strcat(strcat(s,"| dot -T"),type)," -K "),DOTFILTER);
  free(type);
  system(s);
  exit(EXIT_SUCCESS);
}


void Visu(int j){
/*
  Gère l'option "-visu".

  On réécrit la ligne de commande:
  1. en remplaçant "-visu" par "-format userdot" si le FORMAT
     est "userdot" et par "-format dot" sinon
  2. puis en ajoutant à la fin: "| dot -Tpdf -K <filter> -o g.pdf"

  Si le FORMAT est F_no (c'est le cas si l'on a fait "-check maincc"
  ou "loadc" par exemple), alors il y a un problème puisqu'il faut
  qu'un graphe soit généré.
*/
  CheckHelp(&j);j--; /* vérifie l'aide, ici ARGV[j]="-visu" */
  if(FORMAT==F_no) Erreur(24);

  /* on reconstruit dans s la ligne de commande */
  char *s=MakeCMD(NULL,0,j); // ne pas libérer ce pointeur
  strcat(s,"-format ");
  if(FORMAT==F_userdot) strcat(s,"user");
  strcat(s,"dot ");
  MakeCMD(s,j+1,ARGC);
  strcat(strcat(strcat(s,"| dot -Tpdf -K "),DOTFILTER)," -o "xstr(GRAPHPDF));
  system(s);
  exit(EXIT_SUCCESS);
}


void MainCC(int j){
/*
  Gère l'option "-maincc".

  On réécrit la ligne de commande en remplaçant "-maincc" par "-check
  maincc | ./gengraph load - -fast".
*/

  CheckHelp(&j);j--; /* vérifie l'aide, ici ARGV[j]="-maincc" */

  /* on reconstruit dans s la nouvelle ligne de commande */
  char *s=MakeCMD(NULL,0,j); // ne pas libérer ce pointeur
  strcat(s,"-check maincc | ./gengraph load - -fast ");
  MakeCMD(s,j+1,ARGC);
  system(s);
  exit(EXIT_SUCCESS);
}


int gsub(char *s,const char *t,const char *r){
/*
  Remplace, dans la chaîne s, toutes les occurences de t par r et
  renvoie le nombre de remplacements. La chaîne s est modifiée
  directement. Il est donc nécessaire que s ait suffisament de place.
*/

  const int lr=strlen(r); // longueur de t
  const int lt=strlen(t); // longueur de r
  const int d=lt-lr;

  char *p=s+strlen(s)+1; // pointeur sur la fin de s avec son '\0'
  int n=0; // nombre de remplacements effectués

  while((s=strstr(s,t))){
    n++; // une occurrence de plus
    p -= d; // met à jour le pointeur de fin de s
    memmove(s+lr,s+lt,p-s); // déplace la fin
    memcpy(s,r,lr); // copie le remplacement
  }

  return n;
}


/***********************************

           OPTIONS -ALGO

***********************************/


void PrintMorphism(const char *s,const int *P,const int n){
/*
  Affiche le tableau P de n éléments sous la forme:
  i0->j0   i1->j1 ...
  i8->j8   i9->j8 ...
*/
  const int k=8; /* nb de "->" affichés par ligne */
  int i;
  printf("%s",s); /* normalement printf(s) est ok, mais Warning sur certains systèmes */
  for(i=0;i<n;i++)
    printf("%i->%i%s",i,P[i],(((i%k)<(k-1))&&(i<n-1))?"\t":"\n");
  return;
}


/***********************************

           FONCTIONS DE TESTS
              POUR -FILTER

***********************************/


int ftest_minus(graph *G)
/*
  Retourne VRAI ssi G n'est pas dans la famille F, c'est-à-dire G
  n'est isomorphe à aucun graphe de F.
*/
{
  graph *F=MEM(FPARAM,0,graph*);
  if(F==NULL) return (G==NULL);

  int i,*P;
  for(i=0;i<F->f;i++){
    P=Isomorphism(G,F->G[i]);
    free(P);
    if(P!=NULL) return 0;
  }
  
  return 1;
}


int ftest_minus_id(graph *G)
/*
  Retourne VRAI ssi F2 ne contient aucun graphe d'identifiant égale à
  celui de G. La complexité est en O(log|F2|). Il est important que F2
  soit triée par ordre croissant des ID.
*/
{
  graph *F=MEM(FPARAM,0,graph*);
  if(F==NULL) return 0;
  return (bsearch(&G,F->G,F->f,sizeof(graph*),fcmp_graphid)==NULL);
}


int ftest_unique(graph *G)
/*
  Retourne VRAI ssi la sous-famille F allant des indices i+1 à F->f ne contient
  pas G, où i=MEM(FPARAM,0,int).

  Effet de bord: MEM(FPARAM,0,int) est incrémenté.
*/
{
  int i = (MEM(FPARAM,0,int) += 1);
  int *P;

  for(;i<FAMILY->f;i++){
    P=Isomorphism(G,FAMILY->G[i]);
    free(P);
    if(P!=NULL) return 0;
  }

  return 1;
}


int ftest_minor(graph *G)
/*
  Retourne VRAI ssi H est mineur de G.
*/
{
  int *T=Minor(MEM(FPARAM,0,graph*),G);
  if(T==NULL) return 0;
  free(T);
  return 1;
}


int ftest_minor_inv(graph *G)
{
  int *T=Minor(G,MEM(FPARAM,0,graph*));
  if(T==NULL) return 0;
  free(T);
  return 1;
}


int ftest_sub(graph *G)
/*
  Retourne VRAI ssi H est sous-graphe de G avec même nb de sommets.
*/
{
  graph *C=Subgraph(MEM(FPARAM,0,graph*),G);
  if(C==NULL) return 0;
  free_graph(C);
  return 1;
}


int ftest_sub_inv(graph *G)
{
  graph *C=Subgraph(G,MEM(FPARAM,0,graph*));
  if(C==NULL) return 0;
  free_graph(C);
  return 1;
}


int ftest_isub(graph *G)
/*
  Retourne VRAI ssi H est sous-graphe induit de G.
*/
{
  int *C=InducedSubgraph(MEM(FPARAM,0,graph*),G);
  if(C==NULL) return 0;
  free(C);
  return 1;
}


int ftest_isub_inv(graph *G)
{
  int *C=InducedSubgraph(G,MEM(FPARAM,0,graph*));
  if(C==NULL) return 0;
  free(C);
  return 1;
}


int ftest_iso(graph *G)
/*
  Retourne VRAI ssi H est isomorphe à G. Aucun intérêt de faire
  programmer iso-inv.
*/
{
  int *T=Isomorphism(MEM(FPARAM,0,graph*),G);
  if(T==NULL) return 0;
  free(T);
  return 1;
}


int ftest_id(graph *G)
{ return InRange(G->id,FPARAM); }


int ftest_vertex(graph *G)
{ return InRange(G->n,FPARAM); }


int ftest_edge(graph *G)
{ return InRange(NbEdges(G),FPARAM); }


int ftest_degmax(graph *G)
{ return InRange(Degree(G,1),FPARAM); }


int ftest_degmin(graph *G)
{ return InRange(Degree(G,0),FPARAM); }


int ftest_deg(graph *G)
{ int u,b=1,n=G->n;
  for(u=0;(u<n)&&(b);u++) b=InRange(G->d[u],FPARAM);
  return b;
}


int ftest_degenerate(graph *G)
{
  int x;
  int *T=Prune(G,&x);
  free(T);
  return InRange(x,FPARAM);
}


int ftest_gcolor(graph *G)
{
  int *T=Prune(G,NULL);
  int *C=GreedyColor(G,T);
  free(T);
  free(C);
  return InRange(1+G->int1,FPARAM);
}


int ftest_component(graph *G)
{
  param_dfs *p=dfs(G,0,NULL);
  int c=p->nc; /* nb de cc */
  free_param_dfs(p);
  return InRange(c,FPARAM);
}


int ftest_forest(graph *G)
{
  param_dfs *p=dfs(G,0,NULL);
  int c=p->nc; /* nb de cc */
  free_param_dfs(p);
  return InRange(c,FPARAM)&&(NbEdges(G)==G->n-c);
}


int ftest_cutvertex(graph *G)
{
  param_dfs *p=dfs(G,0,NULL);
  int x=p->na;
  free_param_dfs(p);
  return InRange(x,FPARAM);
}


int ftest_biconnected(graph *G)
{
  param_dfs *p=dfs(G,0,NULL);
  int b=(p->nc==1)&&(p->na==0)&&(G->n>2);
  free_param_dfs(p);
  return b;
}


int ftest_ps1xx(graph *G, int version)
{
  path *P=new_path(G,NULL,G->n);
  int v=PS1(G,P,version);
  free_path(P);
  return v;
}


int ftest_ps1(graph *G) { return ftest_ps1xx(G,0); }
int ftest_ps1b(graph *G) { return ftest_ps1xx(G,1); }
int ftest_ps1c(graph *G) { return ftest_ps1xx(G,2); }
int ftest_ps1x(graph *G) { return ftest_ps1xx(G,3); }

int ftest_radius(graph *G)
{
  param_bfs *p;
  const int n=G->n;
  int x=n,u;

  for(u=0;u<n;u++){
    p=bfs(G,u,NULL);
    if(x>=0){
      if(p->n<n) x=-1;
      else x=imin(x,p->radius);
    }
    free_param_bfs(p);
  }
  return InRange(x,FPARAM);
}


int ftest_girth(graph *G)
{
  param_bfs *p;
  int x=1+G->n,u;
  const int n=G->n;
  for(u=0;u<n;u++){
    p=bfs(G,u,NULL);
    if(p->cycle>0) x=imin(x,p->cycle);
    free_param_bfs(p);
  }
  if(x>n) x=-1;
  return InRange(x,FPARAM);
}


int ftest_diameter(graph *G)
{
  param_bfs *p;
  int x=-1,u;
  const int n=G->n;
  for(u=0;u<n;u++){
    p=bfs(G,u,NULL);
    if(p->n==n) x=imax(x,p->radius);
    free_param_bfs(p);
  }
  return InRange(x,FPARAM);
}


int ftest_hyper(graph *G)
{
  param_bfs *p;
  const int n=G->n;
  int u,v,x,y,d1,d2,d3,t,h=0;
  NALLOC(int*,D,n);

  /* calcule la matrice de distance D[][] */
  for(u=0;u<n;u++){
    p=bfs(G,u,NULL); // Dijkstra
    D[u]=p->D;
    if(p->n<n){ /* remplace -1 par +∞ */
      for(v=0;v<n;v++) if(p->D[v]<0) p->D[v]=INT_MAX;
    }
    p->D=NULL;
    free(p); /* efface p, mais pas p->D */
  }

  /* pour tous les quadruplets {u,v,x,y} */
  for(u=0;u<n;u++)
    for(v=u+1;v<n;v++)
      for(x=v+1;x<n;x++)
	for(y=x+1;y<n;y++){
	  d1=D[u][v]+D[x][y];
	  d2=D[u][x]+D[v][y];
	  d3=D[u][y]+D[v][x];
	  if(d1<d2) SWAP(d1,d2,t);
	  if(d1<d3) SWAP(d1,d3,t);
	  if(d2<d3) d2=d3; /* on se fiche de d3 */
	  if(d1-d2>h) h=d1-d2;
	}

  FREE2(D,n); /* efface la matrice de distances */
  if(h==INT_MAX) h=-1; /* cela ne peut jamais arriver */
  return InRange(h,FPARAM);
}


int ftest_tw(graph *G)
{ return InRange(Treewidth(G,1),FPARAM); }


int ftest_tw2(graph *G)
{ return (Treewidth(G,0)<=2); }


int ftest_rename(graph *G)
{
  G->id=SHIFT++;
  return 1;
}


graph *Filter(graph *F,test *f,const int code){
/*
  Etant donnée une famille de graphes et une fonction de test f,
  renvoie une sous-famille de graphes G de F telle que f(G) est
  vraie (si code=0) ou faux (si code=1). Attention! si on libère F,
  cela détruit la sous-famille renvoyée.

  Effet de bord: si PVALUE est vrai, alors dans les graphes filtrés
  G on met dans G->int1 la valeur du paramètre, CVALUE.
*/
  if((F==NULL)||(F->f==0)) return NULL;
  int i,j,n=F->f;

  graph *R=new_graph(0);
  ALLOC(R->G,n); /* a priori R est de même taille que F */

  for(i=j=0;i<n;i++){
    if(f(F->G[i])^code){
      R->G[j++]=F->G[i];
      if(PVALUE) F->G[i]->int1=CVALUE;
    }
  }

  REALLOC(R->G,j);
  R->f=j;
  return R;
}


graph *Graph2Family(graph *G){
/*
  Renvoie une famille composée d'un seul graphe G.
  Effet de bord: met G->id=0.
*/
  graph *F=new_graph(0);
  ALLOC(F->G,1);
  F->G[0]=G;
  F->f=1;
  G->id=0;
  return F;
}


void ApplyFilter(const int code,const int index)
/*
  Applique le filtre FTEST (avec le code=0 ou 1) à la famille de
  graphes FAMILY (voir un graphe seul), et affiche la famille
  résultante. On affiche aussi le nombre de graphes obtenus et la
  ligne de commande (sous forme de commentaire). Si index>=0, alors
  ARGV[index] donne le paramètre.

  Effet de bord: FAMILY est libérée.
*/
{
  graph *R;
  int i;

  if(FAMILY->f==0) FAMILY=Graph2Family(FAMILY); /* transforme en famille si graphe simple */
  R=Filter(FAMILY,FTEST,code); /* calcule la sous-famille */

  printf("// #graphs: %i\n// generated by:",(R==NULL)?0:R->f);
  for(i=0;i<ARGC;i++) printf(" %s",ARGV[i]);
  printf("\n");
  if((index>=0)&&(PVALUE)) /* on affiche la valeur de chaque graphe */
    for(i=0;i<R->f;i++)
      printf("[%i] %s: %i\n",R->G[i]->id,ARGV[index],R->G[i]->int1);
  else PrintGraph(R); /* ou bien on affiche le graphe */

  /* ! aux free() de famille de graphes ! */
  free_graph(FAMILY); /* libère en premier la famille */
  free(R->G); /* libère la sous-famille */
  free(R);
  return;
}


void RS_Start(const char *nom,const char *type,graph *G){
/*
  Partie commune à tous les schémas de routage.
  - nom: est le nom du schéma
  - type: sa catégorie (name-independent ...)
  - G: le graphe sur lequel le schéma doit être appliqué
*/
  printf("\nROUTING SCHEME\n");
  BARRE;
  printf("- name: %s\n",nom);
  printf("- type: %s\n",type);
  printf("- command: %s\n",MakeCMD(NULL,0,ARGC)); // ne pas libérer ce pointeur
  printf("- date: %s\n",DateHeure());
  printf("- seed: %u\n",SEED);
  printf("- time for loading/generating the graph: %s\n",TopChrono(1));
  if(NbEdges(G)<1) Erreur(36); // il faut au moins 1 arête
  param_dfs *X=dfs(G,0,NULL);
  int c=X->nc;
  free_param_dfs(X);
  if(c!=1) Erreur(11); // il doit être connexe
  printf("- checking connectivity: Ok (%s)\n",TopChrono(1));
  int *R=SortGraph(GF,1);
  printf("- checking graph type: simple and undirected (%s)\n",TopChrono(1));
  printf("- #nodes: %i\n",G->n);
  printf("- #edges: %i\n",G->m); // G->m est a jour
  printf("- average degree: %.2lf\n",(double)(G->m<<1)/(double)G->n);
  if(!R[6]) Erreur(31); /* le graphe doit être simple */
  char *s;
  char t[128];
  s="?";
  if(HASH==H_MIX)     s="mix";
  if(HASH==H_PRIME)   s="prime";
  if(HASH==H_SHUFFLE) s="shuffle";
  if(HASH==H_MOD)     s="mod";
  printf("- hash: %s\n",s);
  s="?";
  if(SCENARIO.mode==SC_NONE)   s="none";
  if(SCENARIO.mode==SC_ALL)    s="all";
  if(SCENARIO.mode==SC_NPAIRS) s="n pairs";
  if(SCENARIO.mode==SC_PAIR)   s="pair";
  if(SCENARIO.mode==SC_EDGES)  s="edges";
  if(SCENARIO.mode==SC_ONE)    s="one-to-all";
  if(SCENARIO.mode==SC_UV)     s="u->v";
  if(SCENARIO.mode==SC_UNTIL){ sprintf(t,"until stretch ≥ %g",SCENARIO.stretch);s=t;}
  printf("- scenario: %s\n",s);
  return;
}


/***********************************

               MAIN

***********************************/


int main(const int argc,char* argv[]){

  DEBUG(

	// précaution utile car sinon on peut chercher longtemps une
	// erreur alors qu'elle provient simplement d'effets de bord
	// liée au débugage potentiellement étrange dans certaines
	// situations
	
	fprintf(stderr,"//////////////////////////////////////\n");
	fprintf(stderr,"//!!! WARNING !!! DEBUGGING MODE !!!//\n");
	fprintf(stderr,"//////////////////////////////////////\n");
	
	);
  
  ARGC=argc;
  ARGV=argv;

  if(ARGC==1) Help(0);   /* aide si aucun argument */

  /* initialisations */

  TopChrono(0); /* initialise tous les chronomètres internes */
  VIDE(PARAM_PAL);
  SEED=((getpid()*time(NULL))^clock())&0xFFFF; /* initialise le générateur aléatoire */
  /* On ne prend que 16 bits, plus simple à recopier ... Attention !
     sur certain systèmes clock()=0, et donc pour éviter d'avoir
     toujours SEED=0, il faut éviter de faire '* clock()' ou '&
     clock()' */
  srandom(SEED);
  query *Q=new_query(); /* pour les entrées/sorties d'une fonction d'adjacence */
  Q->adj=ring; /* fonction d'adjacence par défaut */

  int i,j,k,t;

  /******************************************************************

                 ANALYSE DE LA LIGNE DE COMMANDE

    o Il faut éviter de faire des traitements trop coûteux en
      temps/mémoire dans l'analyse de la ligne de commande car on peut
      être amené à la refaire une deuxième fois à cause des alias, les
      options comme -visu, -maincc, -format dot<type> qui causent la
      réécriture puis la ré-analyse de la nouvelle ligne de commande.

    o Il faut éviter d'utiliser random() dans l'analyse de la ligne de
      commande, car si l'option -seed est présente, le comportement
      dépendra de sa position dans la ligne de commande.  Cependant,
      dans certains cas comme "caterpillar" cela est inévitable.

    o Il y a essentiellement quatre cas de figure pour la lecture
      d'arguments ou d'options de la ligne de commande. Il faut
      respecter la structure suivante pour que l'aide en ligne
      fonctionne bien et pour éviter les "Segmentation fault" lors de
      lecture d'argument. Cela arrive dès qu'on essaye de faire
      EQUAL("...") avec une valeur de i en dehors de ARGV[]. En
      général, il faut éviter de faire un "i++" ou un "GetArgInc(&i)"
      suivit d'un EQUAL("...").

      1. option ou graphe avec au moins un argument:

         if EQUAL("kpage"){ i++;
           Q->adj=kpage;
           Q->param[0]=STRTOI(GetArgInc(&i));
           Q->param[1]=STRTOI(GetArgInc(&i));
           ...
           goto fin;
         }

      2. option ou graphe sans argument attendu:

         if EQUAL("tutte"){ CheckHelp(&i);
           Q->adj=tutte;
           goto fin;
         }

      3. option avec une ou plusieurs variantes (de type 1, 2 ou 3):

         if EQUAL("-xy"){ NextArg(&i);
           if EQUAL("load"){ ...; goto fin; }
           if EQUAL("unif"){ ...; goto fin; }
           ...
           Erreur(...);
         }

      4. option avec plusieurs variantes et argument:

         if EQUAL("-xy"){ i++;
           Q->param[0]=STRTOI(GetArgInc(&i));
           i--;NextArg(&i);
           if EQUAL("load"){ ...; goto fin; }
           if EQUAL("unif"){ ...; goto fin; }
           ...
           Erreur(...);
         }

  ******************************************************************/

  i=1; /* on démarre avec le 1er argument */
  while(i<ARGC){

    if EQUAL("-help"){ CheckHelp(&i); i--; }
    if(EQUAL("-help")||EQUAL("?")) Help(i);
    if EQUAL("-list"){ CheckHelp(&i); ListGraph(); }
    if EQUAL("-version"){ CheckHelp(&i); Version(); }
    
    j=i; /* mémorise i pour savoir si on a réussit à lire au moins une
	    option ou un graphe */

    /***********************/
    /* les options -xxx... */
    /***********************/

    if EQUAL("-visu") Visu(i);     /* se termine par system() & exit() */
    if EQUAL("-maincc") MainCC(i); /* se termine par system() & exit() */
    if EQUAL("-not"){ NOT=1-NOT; goto param0; }
    if EQUAL("-permute"){ PERMUTE=1; goto param0; }
    if EQUAL("-undirected"){ DIRECTED=0; LOOP=0; goto param0; }
    if EQUAL("-directed"){ DIRECTED=1; LOOP=1; goto param0; }
    if EQUAL("-noloop"){ LOOP=0; goto param0; }
    if EQUAL("-loop"){ LOOP=1; goto param0; }
    if EQUAL("-header"){ HEADER=1; goto param0; }
    if EQUAL("-vsize"){ VSIZE=1; goto param0; }
    if EQUAL("-fast"){ FAST=1; goto param0; }
    if EQUAL("-seed"){ i++;
	srandom(SEED=STRTOI(GetArgInc(&i)));
	goto fin;
      }
    if EQUAL("-width"){ i++;
	WIDTH=imax(STRTOI(GetArgInc(&i)),0);
	goto fin;
      }
    if EQUAL("-shift"){ i++;
	SHIFT=STRTOI(GetArgInc(&i));
	if(SHIFT<0) Erreur(6);
	goto fin;
      }
    if EQUAL("-pos"){ i++;
	POS=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("-label"){ i++;
	LABEL=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("-dot"){ NextArg(&i);
	if EQUAL("len"){ i++;
	    LEN=STRTOD(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("filter"){ i++;
	    DOTFILTER=GetArgInc(&i); /* pointe sur le nom du filtre */
	    goto fin;
	  }
	if EQUAL("scale"){ i++;
	    DOTSCALE=GetArgInc(&i); /* pointe sur la/les valeurs */
	    goto fin;
	  }
	Erreur(49);
	goto fin;
      }
    if EQUAL("-norm"){ NextArg(&i);
	NORM=NORM_FAIL;
	if EQUAL("L1")    NORM=NORM_L1;
	if EQUAL("L2")    NORM=NORM_L2;
	if EQUAL("Lmax")  NORM=NORM_LMAX;
	if EQUAL("Lmin")  NORM=NORM_LMIN;
	if EQUAL("hyper") NORM=NORM_HYPER;
	if EQUAL("poly"){ i++;
	    NORM=NORM_POLY; 
	    NORM_poly=STRTOI(GetArgInc(&i));
	    if(NORM_poly<3) NORM=NORM_L2;
	    goto fin;
	  }
	if(NORM==NORM_FAIL) Erreur(43);
	i++;
	goto fin;
      }
    if EQUAL("-delv"){ i++;
	DELV=STRTOD(GetArgInc(&i));
	if(DELV>1) DELV=1;
	goto fin;
      }
    if EQUAL("-dele"){ i++;
	DELE=CheckProba(STRTOD(GetArgInc(&i)));
	goto fin;
      }
    if EQUAL("-redirect"){ i++;
	REDIRECT=CheckProba(STRTOD(GetArgInc(&i)));
	goto fin;
      }
    if EQUAL("-star"){ i++;
	STAR=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("-apex"){ i++;
	APEX=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("-caption"){ i++;
	char *c=GetArgInc(&i); // lecture de la légende
	if(CAPTION) free(CAPTION); // si CAPTION déjà allouée
	char *s=strdup(c);
	k=gsub(s,"%SEED","%u");
	if(k>1) Erreur(35);
	if(k==1) asprintf(&CAPTION,s,SEED);
	if(k==0) CAPTION=s;
	goto fin;
      }
    if EQUAL("-variant"){ i++;
	VARIANT=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("-vcolor"){ NextArg(&i);
	/* bits 0-2: codage de la fonction de couleur=1,2,3,4,5
	   bit 3: "pal"
	   bit 4: "list" */
	if EQUAL("deg"){ i++; VCOLOR=(VCOLOR|0x7)^0x7; /* efface les 3 derniers bits */
	    VCOLOR |= 1; goto fin; }
	if EQUAL("degr"){ i++; VCOLOR=(VCOLOR|0x7)^0x7;
	    VCOLOR |= 2; goto fin; }
	if EQUAL("degm"){ i++; VCOLOR=(VCOLOR|0x7)^0x7;
	    VCOLOR |= 3; CHECK=imax(CHECK,CHECK_ON); goto fin; }
	if EQUAL("randg"){ i++; VCOLOR=(VCOLOR|0x7)^0x7;
	    VCOLOR |= 4; CHECK=imax(CHECK,CHECK_ON); goto fin; }
	if EQUAL("kcolor"){ i++; VCOLOR=(VCOLOR|0x7)^0x7;
	    VCOLOR |= 5; CHECK=imax(CHECK,CHECK_ON);
	    if(CPARAM==NULL) ALLOC(CPARAM,PARAMSIZE);
	    MEM(CPARAM,0,int)=STRTOI(GetArgInc(&i)); /* CPARAM[0]=ARGV[i] */
	    goto fin; }
	if EQUAL("pal"){ NextArg(&i); /* teste si arg après pal existe bien */
	    VCOLOR |= 0x8; /* set bit-3 */
	    if(strlen(ARGV[i])>=(int)(sizeof(PARAM_PAL)/sizeof(char))) Erreur(20);
	    strcpy(PARAM_PAL,GetArgInc(&i)); /* PARAM_PAL=ARGV[i] */
	    goto fin; }
	if EQUAL("list"){ i++;
	    VCOLOR |= 0x10;
	    CHECK=imax(CHECK,CHECK_ON);
	    FORMAT = F_dot;
	    goto fin; }
	Erreur(9);
      }
    if EQUAL("-format"){ NextArg(&i);
	FORMAT=-1; /* sentiennelle pour savoir si on a trouvé le FORMAT */
	if EQUAL("standard") FORMAT=F_standard;
	if EQUAL("xy")       FORMAT=F_xy;
	if EQUAL("no")       FORMAT=F_no;
	if EQUAL("userdot"){ FORMAT=F_userdot; WIDTH=1; }
	if EQUAL("html")   { FORMAT=F_html; WIDTH=1; }
	if EQUAL("matrix") { FORMAT=F_matrix; CHECK=imax(CHECK,CHECK_ON);}
	if EQUAL("smatrix"){ FORMAT=F_smatrix;CHECK=imax(CHECK,CHECK_ON);}
	if EQUAL("list")   { FORMAT=F_list;   CHECK=imax(CHECK,CHECK_ON);}
	if EQUAL("vertex"){ i++;
	    FORMAT=F_list;
	    CHECK=imax(CHECK,CHECK_ON);
	    VERTEX0=STRTOI(GetArgInc(&i));
	    i--;
	  }
	if PREFIX("dot"){ /* si "dot" ou "dot<type>" */
	    if EQUAL("dot") FORMAT=F_dot; /* si "dot" seul */
	    else PipeDot(i); /* se termine par system() & exit() */
	  }
	if(FORMAT<0) Erreur(5); /* le format n'a pas été trouvé */
	i++;
	goto fin;
      }
    if EQUAL("-xy"){ NextArg(&i);
	POS=1; /* il faut POS=1 */
	if EQUAL("unif"){ CheckHelp(&i);
	    XYtype=XY_UNIF;
	    goto fin;
	  }
	if EQUAL("load"){ i++;
	    XYtype=XY_FILE;
	    FILEXY=GetArgInc(&i); /* pointe sur le nom du fichier */
	    goto fin;
	  }
	if EQUAL("noise"){ i++;
	    XYnoiser=STRTOD(GetArgInc(&i));
	    XYnoisep=STRTOD(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("box"){ i++;
	    BOXX=STRTOD(GetArgInc(&i));
	    BOXY=STRTOD(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("seed"){ i++;
	    XYseedk=STRTOI(GetArgInc(&i));
	    XYpower=STRTOD(GetArgInc(&i));
	    XYtype=XY_PLAW;
	    goto fin;
	  }
	if EQUAL("hyper"){ i++;
	    XYpower=STRTOD(GetArgInc(&i));
	    XYtype=XY_HYPER;
	    goto fin;
	  }
	if EQUAL("round"){ i++;
	    ROUND=imin(STRTOI(GetArgInc(&i)),DBL_DIG);
	    goto fin;
	  }
	if EQUAL("mesh"){ i++;
	    ROUND=0; /* a priori coordonnées entières dans ce cas */
	    XYtype=XY_MESH;
	    Xmesh=STRTOI(GetArgInc(&i));
	    Ymesh=STRTOI(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("grid"){ i++;
	    XYgrid=STRTOI(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("ratio"){ i++;
	    XYratio=STRTOD(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("vsize"){ i++;
	    XYvsize=STRTOD(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("polygon"){ i++;
	    XYtype=XY_RPOLY;
	    XYpoly=STRTOI(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("surface"){ NextArg(&i);
	    char *s=NULL; // s=signature
	    if EQUAL("square")     s="bb";
	    if EQUAL("plane")      s="bb";
	    if EQUAL("cylinder")   s="hb";
	    if EQUAL("mobius")     s="cb";
	    if EQUAL("torus")      s="hh";
	    if EQUAL("klein")      s="ch";
	    if EQUAL("projective") s="cc";
	    if(s) i++; else s=GetArgInc(&i);
	    int g=0,v; // g=taille de la signature
	    while(s[g]){
	      v=2; // sentinelle
	      if(s[g]=='b') v=0;  // border
	      if(s[g]=='h') v=+1; // handle
	      if(s[g]=='c') v=-1; // crosscap
	      if(v==2) Erreur(45); // caractère non reconnu
	      XYsurface[g++]=v;
	      if(g>SURFACEMAX) Erreur(45); // trop grand
	    }
	    if(g%2) Erreur(45); // problème si g est impaire
	    XYsurfacesize=g;
	    XYtype=XY_RPOLY; // génération dans un polygone
	    XYpoly=2*g; // nombre de cotés du polygone = genre
	    XYratio=1; // polygone régulier
	    goto fin;
	  }
	if EQUAL("zero"){ CheckHelp(&i);
	    XYzero=1;
	    goto fin;
	  }
	if EQUAL("border"){ CheckHelp(&i);
	    XYborder=1;
	    goto fin;
	  }
	if EQUAL("permutation"){ CheckHelp(&i);
	    ROUND=0; /* a priori coordonnées entières dans ce cas */
	    XYtype=XY_PERM;
	    goto fin;
	  }
	if EQUAL("circle"){ CheckHelp(&i);
	    XYtype=XY_CIRCLE;
	    goto fin;
	  }
	if EQUAL("cycle"){ CheckHelp(&i);
	    XYtype=XY_CYCLE;
	    goto fin;
	  }
	if EQUAL("disk") { CheckHelp(&i);
	    XYtype=XY_DISK;
	    goto fin;
	  }
	if EQUAL("convex"){ CheckHelp(&i);
	    XYtype=XY_CONVEX;
	    goto fin;
	  }
	if EQUAL("convex2"){ CheckHelp(&i);
	    XYtype=XY_CONVEX2;
	    goto fin;
	  }
	if EQUAL("unique"){ CheckHelp(&i);
	    XYunique=1;
	    goto fin;
	  }
	Erreur(1); /* l'option après -xy n'a pas été trouvée */
      }
    if EQUAL("-filter"){ NextArg(&i);
	FAMILY=File2Graph(ARGV[i],2); /* lit une famille ou un graphe */
	if(FPARAM==NULL) ALLOC(FPARAM,PARAMSIZE);
	PVALUE=0; /* par défaut, on affiche pas "value" mais les graphes */
	NextArg(&i);
	if EQUAL("not"){ k=1; NextArg(&i); } else k=0;
	 /* vérifie s'il y a bien un autre argument */
	if EQUAL("rename"){
	    FTEST=ftest_rename;
	    i++;
	    SHIFT=STRTOI(GetArgInc(&i));
	    ApplyFilter(0,-1);
	    goto fin;
	  }
	if EQUAL("biconnected"){
	    FTEST=ftest_biconnected;
	  filter0:
	    i++;
	    ApplyFilter(k,-1);
	    goto fin;
	  }
	if EQUAL("id"){
	    FTEST=ftest_id;
	  filter1:
	    i++;
	    ReadRange(GetArgInc(&i),FPARAM);
	    ApplyFilter(k,i-2);
	    goto fin;
	  }
	int c; /* code pour File2Graph() */
	if EQUAL("minor"){
	    FTEST=ftest_minor;
	    c=34; /* détection du shift et charge toujours un graphe */
	  filter2:
	    i++;
	    MEM(FPARAM,0,graph*)=File2Graph(GetArgInc(&i),c);
	    ApplyFilter(k,-1);
	    free_graph(MEM(FPARAM,0,graph*));
	    goto fin;
	  }

	if EQUAL("ps1x"){ i++;
	    FTEST=ftest_ps1x;
	    MEM(CPARAM,0,int)=STRTOI(GetArgInc(&i));
	    for(c=MEM(CPARAM,0,int);c>=1;c--){ /* met les arguments à l'envers */
	      MEM(CPARAM,(2*c-1)*sizeof(int),int)=STRTOI(GetArgInc(&i));
	      MEM(CPARAM,(2*c)*sizeof(int),int)=STRTOI(GetArgInc(&i));
	    }
	    i--;
	    goto filter0;
	  }
	if EQUAL("ps1"){ FTEST=ftest_ps1; goto filter0; }
	if EQUAL("ps1b"){ FTEST=ftest_ps1b; goto filter0; }
	if EQUAL("ps1c"){ FTEST=ftest_ps1c; goto filter0; }
	if EQUAL("tw2"){ FTEST=ftest_tw2; goto filter0; }
	if EQUAL("unique"){ FTEST=ftest_unique; MEM(FPARAM,0,int)=0; goto filter0; }
	if EQUAL("vertex"){ FTEST=ftest_vertex; goto filter1; }
	if(EQUAL("edge")||EQUAL("edges")){ FTEST=ftest_edge; goto filter1; }
	if EQUAL("deg"){ FTEST=ftest_deg; goto filter1; }
	if EQUAL("degenerate"){ FTEST=ftest_degenerate; goto filter1; }
	if EQUAL("degmax"){ FTEST=ftest_degmax; goto filter1; }
	if EQUAL("degmin"){ FTEST=ftest_degmin; goto filter1; }
	if EQUAL("gcolor"){ FTEST=ftest_gcolor; goto filter1; }
	if EQUAL("component"){ FTEST=ftest_component; goto filter1; }
	if EQUAL("cut-vertex"){ FTEST=ftest_cutvertex; goto filter1; }
	if EQUAL("radius"){ FTEST=ftest_radius; goto filter1; }
	if EQUAL("girth"){ FTEST=ftest_girth; goto filter1; }
	if EQUAL("diameter"){ FTEST=ftest_diameter; goto filter1; }
	if EQUAL("hyper"){ FTEST=ftest_hyper; goto filter1; }
	if EQUAL("tw"){ FTEST=ftest_tw; goto filter1; }
	if EQUAL("forest"){ FTEST=ftest_forest; goto filter1; }
	if EQUAL("minor-inv"){ FTEST=ftest_minor_inv; c=34; goto filter2; }
	if EQUAL("sub"){ FTEST=ftest_sub; c=34; goto filter2; }
	if EQUAL("sub-inv"){ FTEST=ftest_sub_inv; c=34; goto filter2; }
	if EQUAL("isub"){ FTEST=ftest_isub; c=34; goto filter2; }
	if EQUAL("isub-inv"){ FTEST=ftest_isub_inv; c=34; goto filter2; }
	if EQUAL("iso"){ FTEST=ftest_iso; c=34; goto filter2; }
	if EQUAL("minus"){ FTEST=ftest_minus; c=2; goto filter2; }
	if EQUAL("minus-id"){ FTEST=ftest_minus_id; c=10; goto filter2; }
	/* alias */
	if EQUAL("connected"){ FTEST=ftest_component;
	    ReadRange("1",FPARAM); goto filter0; }
	if EQUAL("bipartite"){ FTEST=ftest_gcolor;
	    ReadRange("<3",FPARAM); goto filter0; }
	if EQUAL("isforest"){ FTEST=ftest_forest;
	    ReadRange("t",FPARAM); goto filter0; }
	if EQUAL("istree"){ FTEST=ftest_forest;
	    ReadRange("1",FPARAM); goto filter0; }
	if EQUAL("cycle"){ FTEST=ftest_forest;
	    k=1-k; ReadRange("t",FPARAM); goto filter0; }
	if EQUAL("all"){ FTEST=ftest_vertex;
	    ReadRange("t",FPARAM); goto filter0; }
	Erreur(14); /* l'option après -filter n'a pas été trouvée */
      }
    if EQUAL("-check"){ NextArg(&i);
	if(CHECK>CHECK_ON) Erreur(27); /* on ne devrait jamais avoir deux fois -check */
	if(CPARAM==NULL) ALLOC(CPARAM,PARAMSIZE); /* alloue les paramètres */
	if EQUAL("bfs"){
	    CHECK=CHECK_BFS;
	  check0:
	    i++;
	    MEM(CPARAM,0,int)=STRTOI(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("iso"){
	    CHECK=CHECK_ISO;
	  check1:
	    MEM(CPARAM,0,int)=++i;
	    GetArgInc(&i); /* pour vérifier si i existe */
	    goto fin;
	  }
	if(EQUAL("deg")||EQUAL("edge")||EQUAL("edges")){
	  CHECK=CHECK_DEG;
	check_fin:
	  i++;
	  goto fin;
	}
	if EQUAL("paths"){ i++;
	    CHECK=CHECK_PATHS;
	    MEM(CPARAM,0,int)=STRTOI(GetArgInc(&i));
	    MEM(CPARAM,sizeof(int),int)=STRTOI(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("ps1x"){ i++;
	    CHECK=CHECK_PS1x;
	    MEM(CPARAM,0,int)=STRTOI(GetArgInc(&i));
	    for(k=MEM(CPARAM,0,int);k>=1;k--){ /* met les arguments à l'envers */
	      MEM(CPARAM,(2*k-1)*sizeof(int),int)=STRTOI(GetArgInc(&i));
	      MEM(CPARAM,(2*k)*sizeof(int),int)=STRTOI(GetArgInc(&i));
	    }
	    goto fin;
	  }
	if EQUAL("routing"){ NextArg(&i);
	    FORMAT=F_no; /* pour tous les routing schemes */
	    if EQUAL("hash"){ NextArg(&i);
		for(;;){ /* pour faire break (importants à cause des i++) */
		  if EQUAL("mix")    { NextArg(&i); HASH=H_MIX; break; }
		  if EQUAL("prime")  { NextArg(&i); HASH=H_PRIME; break; }
		  if EQUAL("shuffle"){ NextArg(&i); HASH=H_SHUFFLE; break; }
		  if EQUAL("mod")    { NextArg(&i); HASH=H_MOD; break; }
		  Erreur(40); /* option après "hash" non trouvé */
		}
	      }
	    SCENARIO.mode=SC_NONE; /* par défaut aucun scenario */
	    SCENARIO.dist=1; /* par défaut on stocke les distance */
	    if EQUAL("scenario"){ NextArg(&i);
		if EQUAL("nomem"){ NextArg(&i); SCENARIO.dist=0; }
		for(;;){ /* pour faire break (importants à cause des i++) */
		  if EQUAL("none")  { NextArg(&i); SCENARIO.mode=SC_NONE; break; }
		  if EQUAL("all")   { NextArg(&i); SCENARIO.mode=SC_ALL; break; }
		  if EQUAL("npairs"){ NextArg(&i); SCENARIO.mode=SC_NPAIRS; break; }
		  if EQUAL("edges") { NextArg(&i); SCENARIO.mode=SC_EDGES; break; }
		  if EQUAL("one")   { i++;
		      SCENARIO.mode=SC_ONE;
		      SCENARIO.u=STRTOI(GetArgInc(&i));
		      i--;NextArg(&i);
		      break;
		    }
		  if EQUAL("until"){ i++;
		      SCENARIO.mode=SC_UNTIL;
		      SCENARIO.stretch=STRTOD(GetArgInc(&i));
		      i--;NextArg(&i);
		      break;
		    }
		  if EQUAL("pair"){ i++;
		      long p=STRTOL(GetArgInc(&i));
		      SCENARIO.u=abs((int)p);
		      if(p<0L){
			if(-p>(long)INT_MAX) Erreur(46);
			SCENARIO.mode=SC_PAIR;
		      }else{
			SCENARIO.mode=SC_UV;
			SCENARIO.v=STRTOI(GetArgInc(&i));
		      }
		      i--;NextArg(&i);
		      break;
		    }
		  Erreur(41);
		}
	      }
	    if EQUAL("cluster"){ CHECK=CHECK_RS_CLUSTER; goto check0; }
	    if EQUAL("dcr"){ CHECK=CHECK_RS_DCR; goto check0; }
	    if EQUAL("agmnt"){ CHECK=CHECK_RS_AGMNT; goto check0; }
	    if EQUAL("bc"){ CHECK=CHECK_RS_BC; goto check0; }
	    if EQUAL("hdlbr"){ CHECK=CHECK_RS_HDLBR; goto check0; }
	    if EQUAL("tzrplg"){ i++;
		CHECK=CHECK_RS_TZRPLG;
		MEM(CPARAM,0,double)=STRTOD(GetArgInc(&i));
		goto fin;
	      }
	    Erreur(39);
	  }
	if EQUAL("dfs"){ CHECK=CHECK_DFS; goto check0; }
	if EQUAL("bellman"){ CHECK=CHECK_BELLMAN; goto check0; }
	if EQUAL("kcolor"){ CHECK=CHECK_KCOLOR; goto check0; }
	if EQUAL("kcolorsat"){ CHECK=CHECK_KCOLORSAT; FORMAT=F_no; goto check0; }
	if EQUAL("kindepsat"){ CHECK=CHECK_KINDEPSAT; FORMAT=F_no; goto check0; }
	if EQUAL("sub"){ CHECK=CHECK_SUB; goto check1; }
	if EQUAL("isub"){ CHECK=CHECK_ISUB; goto check1; }
	if EQUAL("minor"){ CHECK=CHECK_MINOR; goto check1; }
	if EQUAL("degenerate"){ CHECK=CHECK_DEGENERATE; goto check_fin; }
	if EQUAL("gcolor"){ CHECK=CHECK_GCOLOR; goto check_fin; }
	if EQUAL("ps1"){ CHECK=CHECK_PS1; goto check_fin; }
	if EQUAL("ps1b"){ CHECK=CHECK_PS1b; goto check_fin; }
	if EQUAL("ps1c"){ CHECK=CHECK_PS1c; goto check_fin; }
	if EQUAL("twdeg"){ CHECK=CHECK_TWDEG; goto check_fin; }
	if EQUAL("tw"){ CHECK=CHECK_TW; goto check_fin; }
	if EQUAL("girth"){ CHECK=CHECK_GIRTH; goto check_fin; }
	if EQUAL("info"){ CHECK=CHECK_INFO; FORMAT=F_no; goto check_fin; }
	if(EQUAL("simplify")){ CHECK=CHECK_SIMPLIFY; FORMAT=F_no; goto check_fin; }
	if EQUAL("stretch"){ CHECK=CHECK_STRETCH; goto check_fin; }
	if EQUAL("maincc"){ CHECK=CHECK_MAINCC; FORMAT=F_no; goto check_fin; }
	if(EQUAL("ncc")||EQUAL("connected")){ CHECK=CHECK_NCC; goto check_fin; }
	if EQUAL("diameter"){ CHECK=CHECK_DIAMETER; goto check_fin; }
	if EQUAL("volm"){ CHECK=CHECK_VOLM; goto check_fin; }
	if EQUAL("radius"){ CHECK=CHECK_RADIUS; goto check_fin; }
	Erreur(12); /* l'option après -check n'a pas été trouvée */
      }

    /*******************/
    /* graphes de base */
    /*******************/

    if EQUAL("tutte"){
	Q->adj=tutte;
      param0:
	CheckHelp(&i); /* CheckHelp() au lieu de i++ car aucun paramètre */
	goto fin;
      }
    if EQUAL("prime"){
	Q->adj=prime;
      param1:
	i++;
	Q->param[0]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("arboricity"){
	Q->adj=arboricity;
      param2:
	i++;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("sat"){
	Q->adj=sat;
      param3:
	i++;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("ringarytree"){
	Q->adj=ringarytree;
      param4:
	i++;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=STRTOI(GetArgInc(&i));
	Q->param[3]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("unok"){
	Q->adj=unok; POS=1;
	i++;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=STRTOI(GetArgInc(&i));
	Q->param[3]=STRTOI(GetArgInc(&i)); if(Q->param[3]<0) Q->param[3]=Q->param[1];
	Q->param[4]=STRTOI(GetArgInc(&i)); if(Q->param[4]<0) Q->param[4]=Q->param[2];
	goto fin;
      }
    if EQUAL("ggosset"){
	Q->adj=ggosset;
	i++;
	int d;
	Q->param[0]=STRTOI(GetArgInc(&i)); /* copie p */
	Q->param[1]=STRTOI(GetArgInc(&i)); /* copie k */
	for(k=0;k<(Q->param[1]<<1);k++) Q->param[k+3]=STRTOI(GetArgInc(&i));
	/* il est important de calculer la dimension d pour accélérer l'adjacence */
	for(k=d=0;k<Q->param[1];k++) d += Q->param[(k<<1)+3];
	Q->param[2]=d;
	goto fin;
      }
    if EQUAL("grid"){
	Q->adj=grid;
      param_grid:
	i++;
	int d=STRTOI(GetArgInc(&i));
	if(d<0) Erreur(6);
	if(d>DIMAX) Erreur(4);
	if(d+1>PARMAX){ free(Q->param); ALLOC(Q->param,d+1); }
	for(k=1;k<=d;k++) Q->param[k]=STRTOI(GetArgInc(&i));
	Q->param[0]=d;
	goto fin;
      }
    if EQUAL("hypercube"){
	Q->adj=grid;
      param_hypercube:
	i++;
	int d=STRTOI(GetArgInc(&i));
	if(d<0) Erreur(6);
	if(d>DIMAX) Erreur(4);
	if(d+1>PARMAX){ free(Q->param); ALLOC(Q->param,d+1); }
	for(k=1;k<=d;k++) Q->param[k]=2; // seule différence avec "grid"
	Q->param[0]=d;
	goto fin;
      }
    if EQUAL("ring"){
	Q->adj=ring;
      param_ring:
	i++;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=STRTOI(GetArgInc(&i));
	if(Q->param[1]+2>PARMAX){ free(Q->param); ALLOC(Q->param,Q->param[1]); }
	for(k=0;k<Q->param[1];k++) Q->param[k+2]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("kout"){ i++;
	Q->adj=kout;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=STRTOI(GetArgInc(&i));
	if(Q->param[1]>=Q->param[0]-1) Q->param[1]=Q->param[0]-1;
	if(Q->param[1]<1) Erreur(6);
	goto fin;
      }
    if EQUAL("kneser"){ i++;
	Q->adj=kneser;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=STRTOI(GetArgInc(&i));
	Q->param[0]=imax(Q->param[0],0); /* n>=0 */
	Q->param[1]=imax(Q->param[1],0); /* k>=0 */
	Q->param[1]=imin(Q->param[1],Q->param[0]); /* k<=n */
	goto fin;
      }
    if EQUAL("udg"){ i++;
	Q->adj=udg; POS=1;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->dparam[0]=STRTOD(GetArgInc(&i));
	if(Q->dparam[0]<0) Q->dparam[0]=sqrt(log((double)Q->param[0])/(double)Q->param[0]);
	/* Threshold théorique rc (cf. [EMY07], Theorem 6 avec d=p=2):
	   pour n=10,    rc=0.4798
	   pour n=100,   rc=0.2145
	   pour n=1000,  rc=0.08311
	   pour n=10000, rc=0.03034
	 */
	goto fin;
      }
    if EQUAL("squashed"){ i++;
	Q->adj=squashed;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->dparam[0]=STRTOD(GetArgInc(&i));
	if(Q->dparam[0]<0) Q->dparam[0]=1.0/3; /* proba par défaut */
	Q->dparam[0]=CheckProba(Q->dparam[0]);
	goto fin;
      }
    if EQUAL("rig"){ i++;
	Q->adj=rig;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[0]=imax(Q->param[0],1); /* n>0 */
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[1]=imax(Q->param[1],1); /* k>0 */
	Q->dparam[0]=STRTOD(GetArgInc(&i));
	if(Q->dparam[0]<0){
	  Q->dparam[0]=log((double)Q->param[0])/(double)Q->param[1];
	  if(Q->param[1]>Q->param[0]) Q->dparam[0]=sqrt(Q->dparam[0]/(double)Q->param[0]);
	}
	Q->dparam[0]=CheckProba(Q->dparam[0]); /* proba */
	goto fin;
      }
    if EQUAL("rplg"){ i++;
	Q->adj=rplg;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->dparam[0]=STRTOD(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("thetagone"){ i++;
	Q->adj=thetagone; POS=1;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=STRTOI(GetArgInc(&i));
	Q->dparam[0]=CheckProba(STRTOD(GetArgInc(&i)));
	goto fin;
      }
    if(EQUAL("deltohedron")||EQUAL("trapezohedron")){ i++;
	Q->adj=deltohedron;
	Q->param[0]=(STRTOI(GetArgInc(&i))<<1);
	goto fin;
      }
    if EQUAL("load"){
      param_load:
	i++;
	Q->adj=load;
	Q->sparam=strdup(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("bdrg"){
	Q->adj=bdrg;
      param_bdrg:
	k=++i; while(strcmp(GetArgInc(&i),".")); // NB: i-k=taille séquence avec "."
	if(i-k>PARMAX){ free(Q->param); ALLOC(Q->param,i-k); }
	for(t=k;t<i;t++) Q->param[t-k+1]=STRTOI(ARGV[t]); // lit les valeurs
	Q->param[0]=i-k-1;
	goto fin;
      }
    if EQUAL("rlt"){ i++;
	Q->adj=rlt;
	POS=1; ROUND=0; XYtype=XY_MESH;
	Q->param[0]=Ymesh=STRTOI(GetArgInc(&i)); // nb de colonnes
	Q->param[1]=Xmesh=STRTOI(GetArgInc(&i)); // nb de lignes
	Q->param[2]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("alkane"){ NextArg(&i);
	Q->adj=alkane;
	Q->param[0]=-1;
	if EQUAL("n")      Q->param[0]=ALK_NOR;
	if EQUAL("normal") Q->param[0]=ALK_NOR;
	if EQUAL("cyclo")  Q->param[0]=ALK_CYC;
	if EQUAL("iso")    Q->param[0]=ALK_ISO;
	if EQUAL("neo")    Q->param[0]=ALK_NEO;
	if EQUAL("sec")    Q->param[0]=ALK_SEC;
	if EQUAL("tert")   Q->param[0]=ALK_TER;
	if(Q->param[0]<0) Erreur(6);
	i++;
	Q->param[1]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if SUFFIX("ane"){
	Q->param[0]=ALK_NOR;
	Q->param[1]=-1;
	if SUFFIX("methane") Q->param[1]=1;
	if SUFFIX("ethane")  Q->param[1]=2;
	if SUFFIX("propane") Q->param[1]=3;
	if SUFFIX("butane")  Q->param[1]=4;
	if SUFFIX("pentane") Q->param[1]=5;
	if SUFFIX("hexane")  Q->param[1]=6;
	if SUFFIX("heptane") Q->param[1]=7;
	if SUFFIX("octane")  Q->param[1]=8;
	if SUFFIX("nonane")  Q->param[1]=9;
	if PREFIX("n-")     Q->param[0]=ALK_NOR;
	if PREFIX("cyclo-") Q->param[0]=ALK_CYC;
	if PREFIX("iso-")   Q->param[0]=ALK_ISO;
	if PREFIX("neo-")   Q->param[0]=ALK_NEO;
	if PREFIX("sec-")   Q->param[0]=ALK_SEC;
	if PREFIX("tert-")  Q->param[0]=ALK_TER;
	if SUFFIX("alkane"){ i++; Q->param[1]=STRTOI(GetArgInc(&i)); }
	if(Q->param[1]>=0){ Q->adj=alkane; goto param0; }
      }

    /* graphes de base avec type de paramètres déjà rencontrés */

    if EQUAL("clebsch"){ Q->adj=clebsch; goto param_hypercube; }
    if EQUAL("rpartite"){ Q->adj=rpartite; goto param_grid; }
    if EQUAL("aqua"){ Q->adj=aqua; DIRECTED=1; goto param_grid; }
    if EQUAL("cage"){ Q->adj=cage; goto param_ring; }
    if EQUAL("loadc"){ LOADC=1; FORMAT=F_no; goto param_load; }
    if EQUAL("fdrg"){ Q->adj=fdrg; goto param_bdrg; }
    //
    if EQUAL("icosahedron"){ Q->adj=icosahedron; goto param0; }
    if EQUAL("rdodecahedron"){ Q->adj=rdodecahedron; goto param0; }
    if EQUAL("cuboctahedron"){ Q->adj=cuboctahedron; goto param0; }
    if EQUAL("herschel"){ Q->adj=herschel; goto param0; }
    if EQUAL("goldner-harary"){ Q->adj=goldner_harary; goto param0; }
    if EQUAL("triplex"){ Q->adj=triplex; goto param0; }
    if EQUAL("jaws"){ Q->adj=jaws; goto param0; }
    if EQUAL("starfish"){ Q->adj=starfish; goto param0; }
    if EQUAL("fritsch"){ Q->adj=fritsch; goto param0; }
    if EQUAL("zamfirescu"){ Q->adj=zamfirescu; goto param0; }
    if EQUAL("hatzel"){ Q->adj=hatzel; goto param0; }
    if EQUAL("soifer"){ Q->adj=soifer; goto param0; }
    if EQUAL("poussin"){ Q->adj=poussin; goto param0; }
    if EQUAL("errera"){ Q->adj=errera; goto param0; }
    if EQUAL("kittell"){ Q->adj=kittell; goto param0; }
    if EQUAL("frucht"){ Q->adj=frucht; goto param0; }
    if EQUAL("moser"){ Q->adj=moser; goto param0; }
    if EQUAL("markstrom"){ Q->adj=markstrom; goto param0; }
    if EQUAL("robertson"){ Q->adj=robertson; goto param0; }
    if EQUAL("headwood4"){ Q->adj=headwood4; goto param0; }
    if EQUAL("wiener-araya"){ Q->adj=wiener_araya; goto param0; }
    if EQUAL("hgraph"){ Q->adj=hgraph; goto param0; }
    if(EQUAL("rgraph")||EQUAL("fish")){ Q->adj=rgraph; goto param0; }
    if EQUAL("cricket"){ Q->adj=cricket; goto param0; }
    if EQUAL("moth"){ Q->adj=moth; goto param0; }
    if EQUAL("suzuki"){ Q->adj=suzuki; goto param0; }
    if EQUAL("bull"){ Q->adj=bull; goto param0; }
    if EQUAL("harborth"){ Q->adj=harborth; goto param0; }
    if EQUAL("doily"){ Q->adj=doily; goto param0; }
    if EQUAL("bidiakis"){ Q->adj=bidiakis; goto param0; }
    if EQUAL("schlafli"){ Q->adj=schlafli; goto param0; }
    //
    if EQUAL("gear"){ Q->adj=gear; goto param1; }
    if EQUAL("pstar"){ Q->adj=pstar; Q->param[1]=2; goto param1; }
    if EQUAL("paley"){ Q->adj=paley; goto param1; }
    if EQUAL("comb"){ Q->adj=comb; goto param1; }
    if EQUAL("sunlet"){ Q->adj=sunlet; goto param1; }
    if EQUAL("mycielski"){ Q->adj=mycielski; goto param1; }
    if EQUAL("treep"){ Q->adj=treep; goto param1; }
    if EQUAL("halin"){ Q->adj=halin; goto param1; }
    if EQUAL("windmill"){ Q->adj=windmill; goto param1; }
    if EQUAL("interval"){ Q->adj=interval; goto param1; }
    if EQUAL("circle"){ Q->adj=circle; goto param1; }
    if EQUAL("permutation"){ Q->adj=permutation; goto param1; }
    if EQUAL("pancake"){ Q->adj=pancake; goto param1; }
    if EQUAL("bpancake"){ Q->adj=bpancake; goto param1; }
    if EQUAL("crown"){ Q->adj=crown; goto param1; }
    if EQUAL("shuffle"){ Q->adj=shuffle; goto param1; }
    if EQUAL("flip"){ Q->adj=flip; goto param1; }
    if EQUAL("apollonian"){ Q->adj=apollonian; goto param1; }
    if EQUAL("flower_snark"){ Q->adj=flower_snark; goto param1; }
    if EQUAL("gabriel"){ Q->adj=gabriel; POS=1; goto param1; }
    if EQUAL("sgabriel"){ Q->adj=sgabriel; POS=1; /*FORMAT=F_userdot;*/ goto param1; }
    if EQUAL("rng"){ Q->adj=rng; POS=1; goto param1; }
    if EQUAL("nng"){ Q->adj=nng; POS=1; goto param1; }
    if EQUAL("antiprism"){ Q->adj=antiprism; Q->param[1]=1; goto param1; }
    if EQUAL("butterfly"){ Q->adj=butterfly; goto param1; }
    if EQUAL("matching"){ Q->adj=matching; goto param1; }
    if EQUAL("polygon"){ Q->adj=polygon; goto param1; }
    if EQUAL("cactus"){ Q->adj=cactus; goto param1; }
    if EQUAL("helm"){ Q->adj=helm; goto param1; }
    if EQUAL("haar"){ Q->adj=haar; goto param1; }
    //
    if EQUAL("ktree"){ Q->adj=ktree; Q->param[2]=0; goto param2; }
    if EQUAL("gpetersen"){ Q->adj=gpetersen; goto param2; }
    if EQUAL("debruijn"){ Q->adj=debruijn; goto param2; }
    if EQUAL("kautz"){ Q->adj=kautz; goto param2; }
    if EQUAL("gpstar"){ Q->adj=gpstar; goto param2; }
    if EQUAL("hexagon"){ Q->adj=hexagon; goto param2; }
    if EQUAL("whexagon"){ Q->adj=whexagon; goto param2; }
    if EQUAL("hanoi"){ Q->adj=hanoi; goto param2; }
    if EQUAL("sierpinski"){ Q->adj=sierpinski; goto param2; }
    if EQUAL("banana"){ Q->adj=banana; goto param2; }
    if EQUAL("kpage"){ Q->adj=kpage; goto param2; }
    if EQUAL("line-graph"){ Q->adj=linegraph; goto param2; }
    if EQUAL("linial"){ Q->adj=linial; goto param2; }
    if EQUAL("linialc"){ Q->adj=linialc; goto param2; }
    if EQUAL("expander"){ Q->adj=expander; goto param2; }
    if EQUAL("fan"){ Q->adj=fan; goto param2; }
    if EQUAL("split"){ Q->adj=split; goto param2; }
    if EQUAL("behrend"){ Q->adj=behrend; goto param2; }
    if EQUAL("turan"){ Q->adj=turan; goto param2; }
    //
    if EQUAL("rarytree"){ Q->adj=rarytree; goto param3; }
    if EQUAL("barbell"){ Q->adj=barbell; goto param3; }
    if EQUAL("planar"){ Q->adj=planar; goto param3; }
    if EQUAL("hyperbolic"){ Q->adj=hyperbolic; goto param3; }
    if EQUAL("pat"){ Q->adj=pat; POS=1; goto param3; }
    if EQUAL("uno"){ Q->adj=uno; POS=1; goto param3; }
    if EQUAL("ngon"){ Q->adj=ngon; POS=1; goto param3; }
    //
    if EQUAL("chess"){ Q->adj=chess; goto param4; }

    /********************/
    /* graphes composés */
    /********************/

    if EQUAL("theta"){ i++;
	Q->adj=thetagone; POS=1;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=3;
	Q->param[2]=STRTOI(GetArgInc(&i));
	Q->param[2]=imax(1,Q->param[2]);
	Q->dparam[0]=6.0/Q->param[2];
	goto fin;
      }
    if EQUAL("dtheta"){ i++;
	Q->adj=thetagone; POS=1;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=3;
	Q->param[2]=STRTOI(GetArgInc(&i))/2;
	Q->param[2]=imax(1,Q->param[2]);
	Q->dparam[0]=3.0/Q->param[2];
	goto fin;
      }
    if EQUAL("yao"){ i++;
	Q->adj=thetagone; POS=1;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=0;
	Q->param[2]=STRTOI(GetArgInc(&i));
	Q->param[2]=imax(1,Q->param[2]);
	Q->dparam[0]=2.0/Q->param[2];
	goto fin;
      }
    if EQUAL("path"){ i++;
	Q->adj=grid;
	Q->param[0]=1;
	Q->param[1]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("torus"){ i++;
	Q->adj=grid;
	Q->param[0]=2;
	Q->param[1]=-STRTOI(GetArgInc(&i));
	Q->param[2]=-STRTOI(GetArgInc(&i));
	goto fin;
    }
    if EQUAL("mesh"){ i++;
	Q->adj=grid;
	Q->param[0]=2;
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("mobius"){ i++;
	Q->adj=ring;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=2;
	Q->param[2]=1;
	Q->param[3]=Q->param[0]/2;
	goto fin;
      }
    if EQUAL("ladder"){ i++;
	Q->adj=grid;
	Q->param[0]=Q->param[1]=2;
	Q->param[2]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("johnson"){ i++;
	Q->adj=kneser;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=Q->param[1]-2;
	NOT=1-NOT;
	goto fin;
      }
    if EQUAL("star"){ i++;
	Q->adj=rpartite;
	Q->param[0]=2;
	Q->param[1]=1;
	Q->param[2]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("random"){ i++;
	Q->adj=ring;
	NOT=1-NOT;
	Q->param[0]=STRTOI(GetArgInc(&i));
	DELE=1-STRTOD(GetArgInc(&i));
	Q->param[1]=0;
	goto fin;
      }
    if EQUAL("bipartite"){ i++;
	Q->adj=rpartite;
	Q->param[0]=2;
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("cylinder"){ i++;
	Q->adj=grid;
	Q->param[0]=2;
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=-STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("caterpillar"){ i++;
	Q->adj=grid;
	k=STRTOI(GetArgInc(&i)); /* nb de sommets total */
	STAR=random()%k; /* entre 0...k-1 sommets de deg=1. Active l'opération star() */
	Q->param[0]=1;
	Q->param[1]=k-STAR; /* Q->param[0]=nb de sommets du chemin, qui est >=1 */
	goto fin;
      }
    if EQUAL("centipede"){ i++;
	Q->adj=grid;
	STAR=-1;
	Q->param[0]=1;
	Q->param[1]=STRTOI(GetArgInc(&i)); /* Q->param[0]=nb de sommets du chemin */
	goto fin;
      }
    if EQUAL("sunflower"){ i++;
	Q->adj=cage;
	Q->param[0]=2*STRTOI(GetArgInc(&i));
	Q->param[1]=Q->param[2]=2;
	Q->param[3]=0;
	goto fin;
      }
    if EQUAL("wheel"){ i++;
	Q->adj=ringarytree;
	Q->param[0]=1;
	Q->param[1]=0;
	Q->param[2]=STRTOI(GetArgInc(&i));
	Q->param[3]=2;
	goto fin;
      }
    if EQUAL("tadpole"){ i++;
	Q->adj=barbell;
	Q->param[0]=-STRTOI(GetArgInc(&i));
	Q->param[1]=1;
	Q->param[2]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("pan"){ i++;
	Q->adj=barbell;
	Q->param[0]=-STRTOI(GetArgInc(&i));
	Q->param[1]=1;
	Q->param[2]=1;
	goto fin;
      }
    if EQUAL("web"){ i++;
	Q->adj=ringarytree;
	Q->param[2]=STRTOI(GetArgInc(&i));
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->param[1]=1;
	Q->param[3]=2;
	goto fin;
      }
    if EQUAL("percolation"){ i++;
	Q->adj=udg; POS=1; XYtype=XY_MESH; ROUND=0;
	Xmesh=STRTOI(GetArgInc(&i));
	Ymesh=STRTOI(GetArgInc(&i));
	Q->param[0]=Xmesh*Ymesh; /* nombre de sommets */
	Q->dparam[0]=1; /* rayon */
	DELE=1-CheckProba(STRTOD(GetArgInc(&i))); /* proba existence arête */
	NORM=NORM_L1;
	goto fin;
      }
    if EQUAL("hudg"){ i++; /* paramétrage à revoir */
	Q->adj=udg; POS=1;
	Q->param[0]=STRTOI(GetArgInc(&i));
	Q->dparam[0]=STRTOD(GetArgInc(&i));
	XYtype=XY_HYPER;
	XYpower=Q->dparam[0];
	NORM=NORM_HYPER;
	goto fin;
      }
    if EQUAL("plrg"){ i++;
	Q->adj=bdrg;
	int *S=power_law_seq(STRTOI(GetArgInc(&i)),STRTOD(GetArgInc(&i)),NULL);
	if(S==NULL) Erreur(6);
	Q->param=S; // attention ! Q->param peut être de taille < PARMAX
	goto fin;
      }
    if EQUAL("cubic"){ i++;
	Q->adj=fdrg;
	Q->param[0]=2;
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=3;
	goto fin;
      }
    if EQUAL("regular"){ i++;
	Q->adj=fdrg;
	Q->param[0]=2;
	Q->param[1]=STRTOI(GetArgInc(&i));
	Q->param[2]=STRTOI(GetArgInc(&i));
	goto fin;
      }

    /* graphes composés avec type de paramètres déjà rencontrés */

    if EQUAL("tree"){ Q->adj=arboricity; Q->param[1]=1; goto param1; }
    if EQUAL("cycle"){ Q->adj=ring; Q->param[1]=Q->param[2]=1; goto param1; }
    if EQUAL("rbinary"){ Q->adj=rarytree; Q->param[1]=2; Q->param[2]=0; goto param1; }
    if EQUAL("rbinaryz"){ Q->adj=rarytree; Q->param[1]=2; Q->param[2]=1; goto param1; }
    if (EQUAL("stable")||EQUAL("empty")){ Q->adj=ring; Q->param[1]=0; goto param1; }
    if EQUAL("clique"){ Q->adj=ring; NOT=1-NOT; Q->param[1]=0; goto param1; }
    if EQUAL("outerplanar"){ Q->adj=kpage; Q->param[1]=1; goto param1; }
    if EQUAL("squaregraph"){ Q->adj=planar; Q->param[1]=Q->param[2]=4; goto param1; }
    if EQUAL("prism"){ Q->adj=gpetersen; Q->param[1]=1; goto param1; }
    if EQUAL("d-octahedron"){ Q->adj=matching; NOT=1-NOT; goto param1; }
    if EQUAL("point"){ Q->adj=ring; Q->param[1]=0; POS=1; goto param1; }
    if EQUAL("binary"){ Q->adj=ringarytree;
	Q->param[1]=Q->param[2]=2;Q->param[3]=0; goto param1; }
    if EQUAL("star-polygon"){ Q->adj=ring;
	POS=1; Q->param[1]=Q->param[2]=1; XYtype=XY_DISK; goto param1; }
    if EQUAL("convex-polygon"){ Q->adj=ring;
	POS=1; Q->param[1]=Q->param[2]=1; XYtype=XY_CONVEX; goto param1; }
    if EQUAL("td-delaunay"){ Q->adj=thetagone;
	POS=1; Q->param[1]=Q->param[2]=3; Q->dparam[0]=1; goto param1; }
    //
    if EQUAL("kpath"){ Q->adj=ktree; Q->param[2]=1; goto param2; }
    if EQUAL("kstar"){ Q->adj=ktree; Q->param[2]=2; goto param2; }
    if EQUAL("tw"){ Q->adj=ktree; Q->param[2]=0; DELE=.5; goto param2; }
    if EQUAL("pw"){ Q->adj=ktree; Q->param[2]=1; DELE=.5; goto param2; }
    if EQUAL("knight"){ Q->adj=chess; Q->param[2]=1;Q->param[3]=2; goto param2; }
    if EQUAL("camel"){ Q->adj=chess; Q->param[2]=1;Q->param[3]=3; goto param2; }
    if EQUAL("giraffe"){ Q->adj=chess; Q->param[2]=1;Q->param[3]=4; goto param2; }
    if EQUAL("zebra"){ Q->adj=chess; Q->param[2]=2;Q->param[3]=3; goto param2; }
    if EQUAL("antelope"){ Q->adj=chess; Q->param[2]=2;Q->param[3]=4; goto param2; }
    if EQUAL("lollipop"){ Q->adj=barbell; Q->param[2]=0; goto param2; }
    //
    if EQUAL("arytree"){ Q->adj=ringarytree; Q->param[3]=0; goto param3; }

    /* graphes sans paramètres mais composés, le graphe de base
       contient des paramètres, on doit passer par CheckHelp() */

    // 1 paramètre
    if(EQUAL("cube")||EQUAL("hexahedron")){ Q->adj=crown; Q->param[0]=4; goto param0; }
    if EQUAL("octahedron"){ Q->adj=antiprism; Q->param[0]=3; goto param0; }
    if EQUAL("dodecahedron"){ Q->adj=gpetersen; Q->param[0]=10; Q->param[1]=2; goto param0; }
    if EQUAL("associahedron"){ Q->adj=flip; Q->param[0]=6; goto param0; }
    if EQUAL("tietze"){ Q->adj=flower_snark; Q->param[0]=3; goto param0; }
    if EQUAL("grotzsch"){ Q->adj=mycielski; Q->param[0]=4; goto param0; }
    if EQUAL("egraph"){ Q->adj=comb; Q->param[0]=3; goto param0; }
    // 2 paramètres
    if EQUAL("nauru"){ Q->adj=pstar; Q->param[0]=4;Q->param[1]=2; goto param0; }
    if EQUAL("hajos"){ Q->adj=sierpinski; Q->param[0]=2;Q->param[1]=3; goto param0; }
    if EQUAL("netgraph"){ Q->adj=sierpinski; NOT=1-NOT; Q->param[0]=2;Q->param[1]=3; goto param0; }
    if EQUAL("house"){ Q->adj=grid; NOT=1-NOT; Q->param[0]=1;Q->param[1]=5; goto param0; }
    if EQUAL("tetrahedron"){ Q->adj=ring; NOT=1-NOT; Q->param[0]=4; Q->param[1]=0; goto param0; }
    if EQUAL("claw"){ Q->adj=rpartite; Q->param[0]=2;Q->param[1]=1;Q->param[2]=3; goto param0; }
    if EQUAL("desargues"){ Q->adj=gpetersen; Q->param[0]=10;Q->param[1]=3; goto param0; }
    if EQUAL("durer"){ Q->adj=gpetersen; Q->param[0]=6;Q->param[1]=2; goto param0; }
    if EQUAL("gem"){ Q->adj=fan; Q->param[0]=4;Q->param[1]=1; goto param0; }
    if EQUAL("diamond"){ Q->adj=fan; Q->param[0]=Q->param[1]=2; goto param0; }
    if EQUAL("cross"){ Q->adj=banana; Q->param[0]=1; Q->param[1]=4; goto param0; }
    if(EQUAL("tgraph")||EQUAL("fork")){ Q->adj=banana;Q->param[0]=1;Q->param[1]=3; goto param0; }
    if EQUAL("ygraph"){ Q->adj=banana;Q->param[0]=3;Q->param[1]=1; goto param0; }
    // 3 paramètres
    if EQUAL("petersen"){ Q->adj=kneser; Q->param[0]=5;Q->param[1]=2;Q->param[2]=0; goto param0; }
    if EQUAL("banner"){ Q->adj=barbell; Q->param[0]=-4;Q->param[1]=1;Q->param[2]=1; goto param0; }
    if EQUAL("paw"){ Q->adj=barbell; Q->param[0]=-3;Q->param[1]=1;Q->param[2]=1; goto param0; }
    if EQUAL("theta0"){ Q->adj=barbell; Q->param[0]=Q->param[1]=-5;Q->param[2]=-2; goto param0; }
    if EQUAL("utility"){ Q->adj=rpartite; Q->param[0]=2;Q->param[1]=Q->param[2]=3; goto param0; }
    // 4 paramètres
    if EQUAL("wagner"){ Q->adj=ring;
	Q->param[0]=8; Q->param[1]=2; Q->param[2]=1; Q->param[3]=4; goto param0; }
    if EQUAL("headwood"){ Q->adj=cage;
	Q->param[0]=14; Q->param[1]=2; Q->param[2]=5; Q->param[3]=-5; goto param0; }
    if EQUAL("franklin"){ Q->adj=cage;
	Q->param[0]=12; Q->param[1]=2; Q->param[2]=5; Q->param[3]=-5; goto param0; }
    // 5 paramètres
    if EQUAL("mcgee"){ Q->adj=cage;
	Q->param[0]=24; Q->param[1]=3; Q->param[2]=12; Q->param[3]=7; Q->param[4]=-7;
	goto param0;
      }
    // 6 paramètres
    if EQUAL("dyck"){ Q->adj=cage;
	Q->param[0]=32; Q->param[1]=4;
	Q->param[2]=5; Q->param[3]=0; Q->param[4]=13; Q->param[5]=-13;
	goto param0;
      }
    // 7 paramètres
    if EQUAL("gosset"){ Q->adj=ggosset; // attention! les paramètres sont modifiés (cf. ggosset)
	Q->param[0]=8; Q->param[1]=2;
	Q->param[2]=8; // ajoût du paramètre d=∑_{i=1}^k d_i
	Q->param[3]=2; Q->param[4]=3; Q->param[5]=6; Q->param[6]=-1;
	goto param0;
      }
    // 8 paramètres
    if EQUAL("pappus"){ Q->adj=cage;
	Q->param[0]=18; Q->param[1]=6; Q->param[2]=5; Q->param[3]=7; Q->param[4]=-7;
	Q->param[5]=7; Q->param[6]=-7; Q->param[7]=-5;
	goto param0;
      }
    if EQUAL("tutte-coexter"){ Q->adj=cage;
	Q->param[0]=30; Q->param[1]=6;
	Q->param[2]=-7;Q->param[3]=9;Q->param[4]=13;
	Q->param[5]=-13;Q->param[6]=-9;Q->param[7]=7;
	goto param0;
      }
    if EQUAL("gray"){ Q->adj=cage;
	Q->param[0]=54; Q->param[1]=6; Q->param[2]=7; Q->param[3]=-7;
	Q->param[4]=25;Q->param[5]=-25; Q->param[6]=13; Q->param[7]=-13;
	goto param0;
      }
    // 12 paramètres
    if EQUAL("chvatal"){ Q->adj=cage;
	Q->param[0]=Q->param[1]=12;
	Q->param[2]=Q->param[4]=Q->param[7]=Q->param[10]=Q->param[12]=Q->param[13]=3;
	Q->param[3]=Q->param[5]=Q->param[6]=Q->param[8]=6;
	Q->param[9]=Q->param[11]=-3;
	goto param0;
      }

  fin:
    if(j==i){
      if PREFIX("-") Erreur(2); /* option non trouvée */
      Erreur(10); /* graphe non trouvé */
    }
    
  } /* fin du while(i<ARGC) ... */

  /* options qui ne vont pas ensemble */
  if((STAR)&&(APEX)) Erreur(18); /* on ne peut pas avoir les deux */
  if((LOADC)&&(PERMUTE||NOT)) Erreur(29); /* options incompatibles */
  if((LOADC)&&(CHECK<=CHECK_ON)) Erreur(30); /* manque -check */


  /***********************************

           COEUR DU GENERATEUR

  ***********************************/

  /* si le graphe est apex() ou star() */
  if(APEX){ query *tmp=Q; Q=new_query(); Q->query=tmp; Q->adj=apex; Q->param[0]=APEX; }
  if(STAR){ query *tmp=Q; Q=new_query(); Q->query=tmp; Q->adj=star; Q->param[0]=STAR; }

  /* initialisation du graphe, calcule Q->n avant la suppression
     éventuelle de sommet, lit le graphe Q->G si adj=load, détermine
     (XPOS,YPOS) si graphe géométrique */  

  Q->code=QUERY_INIT;
  if(Q->adj(Q)&&(Q->error)) Erreur(Q->error);
  if(Q->n<0) Q->n=0; /* ne devrait jamais arriver */
 
  if(LOADC){ GF=Q->G; NF=Q->n; Q->adj=NULL; goto check; } /* saute la partie génération d'arêtes */
  if(POS && XPOS==NULL) InitXY(Q); /* il faut déterminer les positions (si pas déjà fait) */
  if(abs(LABEL)==1) PERMUTE=0; /* si on souhaite les labels d'origine, on ne permute rien */

  ALLOC(V,Q->n);       /* V[i]=étiquette du sommet i, -1 si i est supprimé */
  ALLOC(VF,Q->n);      /* VF[j]=indice du j-ème sommet non supprimé */
  NF=InitVertex(Q->n,DELV); /* initialise V[i], VF[i] et renvoie NF=#sommets final */
  ALLOC(INC,Q->n);     /* INC[i]=1 ssi i possède un voisin, 0 si sommet isolé */
  for(j=0;j<NF;j++) INC[VF[j]]=0;
  
  /* constantes pour accélérer les tests de la boucle principale */
  const long seuil_edge=(double)(1-DELE)*(double)RAND_MAX;
  const long seuil_redirect=(double)(REDIRECT)*(double)RAND_MAX;

  /*
    Génère les adjacences i-j en tenant compte des sommets isolés et
    des sommets supprimés. Les sommets isolés sont affichés en dernier
    à cause de l'option -redirect. On a toujours i<j lorsque l'arête
    i-j doit être sortie. Si on a FAST=1, alors on génère le graphe à
    partir de Q->G, s'il existe.
  */

  Q->code=QUERY_INIT;
  Out(Q); /* initialise le format d'affichage */
  
  if(FAST){
    if(Q->G){
      const int redirect=(REDIRECT!=0);
      int d,t;

      /* on ne teste que les arêtes de Q->G */
      for(i=0;i<Q->n;i++)
	if(V[i]>=0) /* si le sommet i existe */
	  for(t=0,d=Q->G->d[i];t<d;t++){
	    j=Q->G->L[i][t]; if((Q->G->sym)&&(j<i)) continue; /* il faut i<j si Q->G symétrique */
	    if((V[j]>=0)&&((!DIRECTED)||(i!=j)||(LOOP))) /* si j existe ... */
	      if(random()<seuil_edge){
		if(redirect){ /* si redirection */
		  j=(random()<seuil_redirect)? random()%Q->n : j;
		  if((V[j]<0)||(j==i)) continue; /* prochain voisin */
		}
		INC[i]++; /* un voisin de plus pour i */
		INC[j]++; /* un voisin de plus pour j */
		Q->code=QUERY_ADJ,Q->i=i,Q->j=j;
		Out(Q); /* sort l'arête i-j avec i<j */
	      }
	  }
    }else Erreur(25); /* -fast mais Q->G n'existe pas */
  }else{
    const int noredirect=(REDIRECT==0);

    /* teste les O(n^2) arcs ou arêtes possibles */
    for(i=0;i<Q->n;i++)       /* pour tous les i */
      if(V[i]>=0){         /* si i existe */
	for(j=(DIRECTED)?0:i+1-LOOP;j<Q->n;j++) /* pour tous les j>i */
	  if((V[j]>=0)&&((!DIRECTED)||(i!=j)||(LOOP))) /* si j existe ... */
	    if(random()<seuil_edge){
	      Q->code=QUERY_ADJ,Q->i=i,Q->j=j;Q->adj(Q);
	      if(Q->a^NOT){
		/* ici l'arête i-j devrait être sortie */
		if(noredirect){ /* si pas de redirection d'arête */
		  INC[i]++; /* un voisin de plus pour i */
		  INC[j]++; /* un voisin de plus pour j */
		  Q->i=i,Q->j=j; // ici Q->code=QUERY_ADJ
		  Out(Q); /* sort l'arête i-j */
		}else{ /* on redirige l'arête i-j vers i-k */
		  k=(random()<seuil_redirect)? random()%Q->n : j;
		  if((V[k]>=0)&&(k!=i)){
		    /* on affiche l'arête que si k existe et si k<>i */
		    /* Attention ! il ne faut toucher ni à i ni à j */
		    INC[i]++; /* un voisin de plus pour i */
		    INC[k]++; /* un voisin de plus pour k */
		    if(k<i){
		      Q->i=k,Q->j=i;
		      Out(Q);
		    }else{
		      Q->i=i,Q->j=k;
		      Out(Q); /* pour avoir i<j */
		    }
		  }
		}
	      }
	    }
      }

  }

  /* affiche les sommets isolés */
  for(Q->i=0;Q->i<Q->n;Q->i++)
    if((V[Q->i]>=0)&&(!INC[Q->i])){ Q->code=QUERY_ISOL; Out(Q); }
  
  /* fin de l'affichage, doit être fait avant adj() QUERY_END) */
  Q->code=QUERY_END; Out(Q); /* NB: calcule GF si CHECK */

  if((CHECK)&&(POS)){ /* mémorise XPOS,YPOS qui vont être supprimés par adj(0,ADJ_END) */
    ALLOCZ(GF->xpos,NF,XPOS[VF[_i]]);
    ALLOCZ(GF->ypos,NF,YPOS[VF[_i]]);
  }

  free(V);
  free(VF);
  free(INC);
  
 /* termine la fonction Q->adj() */
  Q->code=QUERY_END;
  Q->adj(Q);
  Q->adj=NULL;

  /***********************************

           FIN DU GENERATEUR

         Si on a CHECK<>0 alors:
         - le graphe généré est GF
         - son nombre de sommets est NF

  ***********************************/

 check:
  /* NB: dans le cas LOADC on ne fait pas adj(Q) avec
     Q->code=QUERY_END, c'est-à-dire load(Q) car sinon load(Q) va
     libèrer le graphe Q->G. Or on veut conserver GF=Q->G
     précisément. */

  if(CHECK){
    if((GF==NULL)||(GF->n!=NF)) Erreur(32); /* ne devrait jamais arriver */

    switch(CHECK){

    case CHECK_MAINCC:{
      param_dfs *p=dfs(GF,MEM(CPARAM,0,int),NULL);
      int d0,d1,n,c;
      for(i=c=n=d0=d1=0;i<p->nc;i++){ /* détermine la plus grosse cc */
	/* d1-d0=nb de sommets de la cc numéro i */
	if(i+1<p->nc) d1=p->d[p->R[i+1]]; else d1=NF;
	if(d1-d0>n){ n=d1-d0; c=i; } /* n=nb de sommets de cc max */ 
	d0=d1;
      }
      c=p->C[p->R[c]]; /* c=couleur de la composante max */
      NALLOC(int,T,n); /* T=sommet de GF à garder */
      for(i=j=0;i<NF;i++) if(p->C[i]==c) T[j++]=i;
      free_param_dfs(p); /* p ne sert plus à rien */
      graph *C=ExtractSubgraph(GF,T,n,1); /* construit la cc max */
      free(T); /* T ne sert plus à rien */
      CHECK=CHECK_OFF; /* pour ne pas restocker le graphe */
      PrintGraph(C);
      free_graph(C);
    }break;

    case CHECK_BFS:{
      param_bfs *p=bfs(GF,MEM(CPARAM,0,int),NULL);
      int t,c;
      printf("root=%i\n",p->root);
      printf("rad[%i]=%i\n",p->root,p->radius);
      printf("cycle[%i]=%i%s\n",p->root,p->cycle,(p->cycle<0)?" (undefined)":"");
      printf("#vertices traversed=%i\n",p->n);
      printf("distance:\n");
      for(k=0;k<=p->radius;k++){
	printf(" d=%i:",k);
	for(c=t=0;t<NF;t++) if(p->D[t]==k){ c++; printf(" %i",t); }
	printf("  (×%i)\n",c);
      }
      printf("vertices not connected to the source:");
      for(t=c=0;t<NF;t++) if(p->D[t]<0) printf(" %i",t), c++;
      if(c==0) printf(" none\n");
      else printf("  (×%i)\n",c);
      printf("parent:");
      for(k=0;k<NF;k++) printf(" %i",p->P[k]);
      printf("\n");
      free_param_bfs(p);
    }break;

    case CHECK_DFS:
    case CHECK_NCC:{
      if(NF<=0) printf("Empty graph\n");
      int s=(CHECK==CHECK_DFS)? MEM(CPARAM,0,int) : 0;
      if(CHECK==CHECK_DFS) printf("source: %i\n",s);
      param_dfs *p=dfs(GF,s,NULL);
      printf("#component: %i%s\n",p->nc,(p->nc==1)?" (connected)":"");
      printf("#cut-vertex: %i%s\n",p->na,
	     ((p->nc==1)&&(p->na==0)&&(NF>2))?" (biconnected)":"");
      if(CHECK==CHECK_NCC) goto check_ncc;
      printf("root:");
      for(i=0;i<p->nc;i++) printf(" %i",p->R[i]);
      if(p->na) printf("\ncut-vertex:");
      for(i=0;i<NF;i++) if(p->A[i]) printf(" %i",i);
      if(p->nc>1){
	printf("\ncomponent:");
	for(i=0;i<NF;i++) printf(" %i",p->C[i]);
      }
      printf("\nparent:");
      int m,h,c;
      for(i=m=0;i<NF;i++){
	if(p->P[i]<0) printf(" -");
	else printf(" %i",p->P[i]);
	m=imax(m,p->H[i]); // max depth
      }
      printf("\ndepth: %i",m);
      if(m==NF-1) printf(" (hamiltonian path!)");
      PRINTN;
      for(h=0;h<=m;h++){
	printf(" h=%i:",h);
	for(i=c=0;i<NF;i++) if(p->H[i]==h) printf(" %i",i), c++;
	printf("  (×%i)\n",c);
      }
    check_ncc:
      free_param_dfs(p);
    }break;

    case CHECK_BELLMAN:{
      if(!InitWeights(GF,IW_POS)) Erreur(6);
      int u;
      double tk,dk,sk,d,s=0;
      param_bellman *p=Bellman_Ford(GF,MEM(CPARAM,0,int),NULL);
      printf("source=%i\n",p->source);
      for(k=0;k<NF;k++){
	printf("dist[%i]=%lf \tparent=%i",k,p->dist[k],p->parent[k]);
	if(POS){
	  dk=hypot(GF->xpos[p->source]-GF->xpos[k],GF->ypos[p->source]-GF->ypos[k]);
	  sk=(dk==0)? 1 : p->dist[k]/dk;
	  if(sk>s) s=sk,u=k,d=dk,tk=p->dist[k]; // vrai une fois car s=0 au départ
	  printf(" \tstretch=%lf (=%lf/%lf)",sk,dk,p->dist[k]);
	}
	printf("\n");
      }
      if(POS){
	printf("maximum stretch: %lf",s);
	if(d==0) printf("\n"); else printf(" (=%lf/%lf)\n",tk,d);
	printf("maximum stretch pair: %i->%i\n",u,p->source);
	printf("maximum stretch path: %i",u);
	u=p->parent[u];
	while(u>=0){
	  printf("->%i",u);
	  u=p->parent[u];
	}
      }
      printf("\n");
      free_param_bellman(p);
    }break;

    case CHECK_STRETCH:
      if(!InitWeights(GF,IW_GEO)) Erreur(44);
      else{
	int u,v,w,a,b,x,umin,vmin;
	double sv,sw,dw,dv,tv,d,t,s=0;
	double dmin,tmin,smin=-1;
	NALLOC(int,T,NF+1); // T[0..]=chemin du stretch_max (se termine par -1)
	NALLOC(int,C,NF+1); // C[0..]=chemin du stretch_max min (se termine par -1)
	param_bellman *p=new_param_bellman(NF);
	for(v=0;v<NF;v++){ // pour toutes les sources v
	  Bellman_Ford(GF,v,p);
	  sv=0; // sv=stretch max pour la source v
	  for(w=0;w<NF;w++){
	    dw=hypot(GF->xpos[v]-GF->xpos[w],GF->ypos[v]-GF->ypos[w]);
	    sw=(dw==0)? 1 : p->dist[w]/dw; // sw=stretch(w,v)
	    if(sw>sv) sv=sw,u=w,dv=dw,tv=p->dist[w]; // vrai une fois car sv=0 au départ
	  }
	  
	  // ici u->...->v a le stretch_max pour la source v

	  if(sv>s){ // on a trouvé un stretch supérieur -> recopie dans T
	    s=sv,a=u,b=v,d=dv,t=tv;
	    T[w=0]=x=u; while(x>=0) T[++w]=x=p->parent[x]; // termine par -1
	  }
	  if((smin<0)||(sv<smin)){ // stretch_max minimum
	    smin=sv,umin=u,vmin=v,dmin=dv,tmin=tv;
	    C[w=0]=x=u; while(x>=0) C[++w]=x=p->parent[x]; // termine par -1
	  }
	}
	free_param_bellman(p);
	
	printf("maximum stretch: %lf",s);
	if(NF==0){ printf(" (empty graph)\n"); free(T); free(C); break; }
	if(d==0) printf("\n"); else printf(" (=%lf/%lf)\n",tv,d);
	printf("maximum stretch pair: %i->%i\n",a,b);
	printf("maximum stretch path: %i",T[w=0]);
	while(T[++w]>=0) printf("->%i",T[w]);
	printf("\n");
	free(T);
	
	printf("minimum max. stretch: %lf",smin);
	if(dmin==0) printf("\n"); else printf(" (=%lf/%lf)\n",tmin,dmin);
	printf("minimum max. stretch pair: %i->%i\n",umin,vmin);
	printf("minimum max. stretch path: %i",C[w=0]);
	while(C[++w]>=0) printf("->%i",C[w]);
	printf("\n");
	free(C);
      }
      break;

    case CHECK_DEG:
      printf("#edges: %i\n",NbEdges(GF));
      PrintDistribution(GF->d,NF,-1,"degree");
      break;

    case CHECK_DEGENERATE:{
      int *T=Prune(GF,&k);
      printf("Degenerate: %i\n",k);
      for(k=0;k<NF;k++) printf("%i ",T[k]);
      printf("\n");
      free(T);
    }break;

    case CHECK_GCOLOR:{
      int *T=Prune(GF,NULL);
      int *C=GreedyColor(GF,T);
      printf("#colors: %i\n",1+GF->int1);
      PrintMorphism("Coloring (node->color):\n",C,GF->n);
      free(C);
      free(T);
    }break;

    case CHECK_KCOLOR:{
      int k=MEM(CPARAM,0,int);
      int *C=kColor(GF,k);
      if(C==NULL) printf("There is no %i-coloration for this graph.\n",k);
      else{
	printf("#colors: %i\n",1+GF->int1);
	PrintMorphism("Coloring (node->color):\n",C,GF->n);
	free(C);
      }
    }break;

    case CHECK_KCOLORSAT:
      kColorSat(GF,MEM(CPARAM,0,int));
      break;

    case CHECK_KINDEPSAT:
      kIndepSat(GF,MEM(CPARAM,0,int));
      break;

    case CHECK_PS1:  k=0; goto check_ps;
    case CHECK_PS1b: k=1; goto check_ps;
    case CHECK_PS1c: k=2; goto check_ps;
    case CHECK_PS1x: k=3;
    check_ps:;
      path *P=new_path(GF,NULL,NF);
      int v=PS1(GF,P,k);
      printf("#tests: %i\nPS1: %s\n",GF->int1,v?"yes (PS1 for sure)":"no (probably not PS1)");
      free_path(P);
      break;
      
    case CHECK_TWDEG:{
      int *T=Prune(GF,&k);
      printf("treewidth <= %i\n",Treewidth(GF,0));
      printf("treewidth >= %i\n",k);
      free(T);
    }break;

    case CHECK_TW:
      k=Treewidth(GF,1);
      printf("#tests: %i\ntreewidth: %i\n",GF->int1,k);
      break;

    case CHECK_DIAMETER:{
      param_bfs *p;
      int u,x=-1;
      for(u=0;u<NF;u++){
	p=bfs(GF,u,NULL);
	if(p->n<NF) break; /* non connexe */
	x=imax(x,p->radius);
      }
      free_param_bfs(p);
      printf("diameter: ");
      if(x<0) printf("%s","+∞");
      else printf("%i",x);
      printf("\n");
    }break;

    case CHECK_VOLM:{
      // on supprime tous les arcs décroissant
      // on aurait pu ensuite faire un simple DFS
      param_bellman *p=new_param_bellman(NF);
      if(!InitWeights(GF,IW_POS)) Erreur(6);
      HalfGraph(GF,1); // u->v avec u<v, attention! GF est modifié
      //PrintGraphList(GF);
      //PrintGraph(GF);
      NALLOC(int,M,NF); // M[u]=volume monotone de u
      int u,v;
      for(u=0;u<NF;u++){
	Bellman_Ford(GF,u,p);
	// calcule la densité (nombre d'arcs) de la boule couverte
	// par l'arbre de racine u
	M[u]=GF->d[u];
	for(v=0;v<NF;v++) if(p->parent[v]>=0) M[u] += GF->d[v];
      }
      free_param_bellman(p);
      //PRINTT(M,NF);
      PrintDistribution(M,NF,10,"monotonic volume");
    }break;

    case CHECK_RADIUS:{
      param_bfs *p=new_param_bfs();
      int u,x=NF-1;
      p->clean=1;
      for(u=0;u<NF;u++){
	bfs(GF,u,p);
	if(p->n<NF){ x=-1; break; } /* non connexe */
	x=imin(x,p->radius);
      }
      free_param_bfs(p);
      printf("radius: ");
      if(x<0) printf("%s","+∞");
      else printf("%i",x);
      printf("\n");
    }break;

    case CHECK_GIRTH:{
      param_bfs *p=new_param_bfs();
      int u,x=1+NF;
      p->clean=1;
      for(u=0;u<NF;u++){
	bfs(GF,u,p);
	if(p->cycle>0) x=imin(x,p->cycle);
      }
      free_param_bfs(p);
      if(x>NF) x=-1;
      printf("girth: %i%s\n",x,(x<0)?" (undefined)":"");
    }break;

    case CHECK_PATHS:{ /*  sort tous les chemins de x à y */
      path *P=new_path(GF,NULL,NF); /* chemin vide d'au plus NF sommets */
      P->P[0]=MEM(CPARAM,0,int);           /* sommet début */
      P->P[1]=MEM(CPARAM,sizeof(int),int); /* sommet fin */
      if((P->P[0]<0)||(P->P[1]<0)||(P->P[0]>=NF)||(P->P[1]>=NF)) Erreur(37);
      int u,v=NextPath(GF,P,-1); /* initialise le premier chemin */
      while(v){
	for(u=0;u<P->n;u++) printf("%i ",P->P[u]);
	printf("\n");
	v=NextPath(GF,P,0);
      }
      free_path(P);
    }break;

    case CHECK_ISO:{
      char *s=ARGV[MEM(CPARAM,0,int)];
      graph *H=File2Graph(s,34);
      int *P=Isomorphism(GF,H);
      printf("H: %s\n#tests: %i\n",s,H->int1);
      if(P==NULL) printf("Non-isomorphic.\n");
      else{
	PrintMorphism("Isomorphism G->H:\n",P,NF);
	free(P);
      }
      free_graph(H);
    }break;

    case CHECK_SUB:{
      char *s=ARGV[MEM(CPARAM,0,int)];
      graph *H=File2Graph(s,34);
      graph *S=Subgraph(GF,H);
      printf("H: %s\n#tests: %i\n",s,H->int1);
      if(S==NULL) printf("G is not a subgraph of H.\n");
      else{
	printf("Subgraph S of H isomorphic to G:\n");
	PrintGraph(S);
	PrintMorphism("Isomorphism S->G:\n",S->pint1,S->n);
	free_graph(S);
      }
      free_graph(H);
    }break;

    case CHECK_MINOR:{
      char *s=ARGV[MEM(CPARAM,0,int)];
      graph *H=File2Graph(s,34);
      int *C=Minor(H,GF);
      printf("H: %s\n#tests: %i\n",s,H->int1);
      if(C==NULL) printf("H is not a minor of G.\n");
      else{
	int c,u;
	printf("Model of H in G:\n");
	for(c=0;c<H->n;c++){ /* pour chaque couleur c */
	  printf("%i -> {",c);
	  for(u=0;u<NF;u++) /* on affiche les sommets de la couleur c */
	    if(C[u]==c) printf(" %i",u);
	  printf(" }\n");
	}
	free(C);
      }
      free_graph(H);
      }break;

    case CHECK_ISUB:{
      char *s=ARGV[MEM(CPARAM,0,int)];
      graph *H=File2Graph(s,34);
      int *X=InducedSubgraph(H,GF);
      printf("H: %s\n#tests: %i\n",s,GF->int1);
      if(X==NULL) printf("H is not an induced subgraph of G.\n");
      else{
	int u;
	printf("Vertices of the induced subgraph S:");
	for(u=0;u<H->n;u++) printf(" %i",X[u]);
	for(u=0;u<H->n;u++) GF->pint1[u]=X[u];
	PrintMorphism("\nIsomorphism H->S:\n",GF->pint1,H->n);
	free(X);
      }
      free_graph(H);
      }break;

    case CHECK_INFO:{
      printf("- command: %s\n",MakeCMD(NULL,0,ARGC)); // ne pas libérer ce pointeur
      printf("- seed: %u\n",SEED);
      printf("- time to generate or load the graph: %s\n",TopChrono(0));
      int *R=SortGraph(GF,1);
      printf("- time to traverse and sort the graph: %s\n",TopChrono(0));
      if(R){
	printf("- simple and undirected: %s\n",R[6]?"yes":"no");
	printf("- geometric (X,Y positions): %s\n",GF->xpos?"yes":"no");
	printf("- edge-weights: %s\n",GF->W?"yes":"no");
	printf("- #nodes: %s\n",millier(NF));
	printf("- #arcs: %s\n",millier(R[2]));
	printf("- #self-loops: %s\n",millier(R[0]));
	printf("- #multi-arcs: %s\n",millier(R[1]));
	printf("- #asymmetric relations: %s\n",millier(R[3])); 
	printf("- #node IDs < 0: %s\n",millier(R[4])); 
	printf("- #node IDs ≥ n: %s\n",millier(R[5])); 
	printf("- maximum degree: %s\n",millier(R[7]));
	printf("- minimum degree: %s\n",millier(R[8]));
	printf("- #isolated nodes: %s\n",millier(R[9]));
      }else printf("- empty graph\n");
      printf("- memory space: %s bytes\n",millier(SizeOfGraph(GF)));
    }break;

    case CHECK_SIMPLIFY:{
      int u,v,d,i;
      SortGraph(GF,0);
      for(u=0;u<NF;u++){
	d=GF->d[u];
	for(i=0;i<d;i++){
	  v=GF->L[u][i];
	  if(((v>u)||((LOOP)&&(v==u)))&&
	     ((i==0)||(v!=GF->L[u][i-1]))) printf("%i-%i\n",u,v);
	}
      }
      }break;

    case CHECK_RS_CLUSTER:
      RS_Start("cluster",RS_NI_FP,GF);
      
      /* paramètre */
      k=MEM(CPARAM,0,int); /* k=paramètre */
      if(k==-1) k=ceil(sqrt((double)NF)); /* ici k>=1, toujours */
      if(k==-2) k=NF;
      printf("- parameter: %i\n",k); /* ici k>=1, toujours */
      if(k<1) Erreur(6); /* k=0: valeur impossible */
      if(VARIANT) printf("- variant: %i\n",VARIANT);
      BARRE;
      
      /* construit, teste, puis libère les tables */
      rs_cluster_tables *RT=rs_cluster(GF,k); /* construit */
      routing_test(GF,RT,(rt_length)rs_cluster_length,-1,NULL); /* teste */
      free_rs_cluster_tables(RT); /* libère */
      break;

    case CHECK_RS_DCR:{
      RS_Start("dcr",RS_NI_FP,GF);

      /* paramètre */
      k=MEM(CPARAM,0,int); /* k=paramètre */
      if(k==-1) k=ceil(Minimize(func1,&NF,1,NF,0)); /* ici k>=1 */
      if(k==-2) k=NF;
      printf("- parameter: %i\n",k); /* ici k>=1, toujours */
      if(k<1) Erreur(6); /* il faut k>0 */
      if(VARIANT) printf("- variant: %i\n",VARIANT);
      BARRE;

      /* construit, teste, puis libère les tables */
      rs_dcr_tables *RT=rs_dcr(GF,k);
      routing_test(GF,RT,(rt_length)rs_dcr_length,-1,RT->dist);
      free_rs_dcr_tables(RT);
    }break;

    case CHECK_RS_AGMNT:{
      RS_Start("agmnt",RS_NI_FP,GF);

      /* paramètre */
      VARIANT |= 4; // bit-2 a 1 ssi AGMNT
      k=MEM(CPARAM,0,int); /* k=paramètre */
      if(k==-1) k=ceil(Minimize(func1,&NF,1,NF,0)); /* ici k>=1 */
      if(k==-2) k=NF;
      printf("- parameter: %i\n",k); /* ici k>=1, toujours */
      if(k<1) Erreur(6); /* il faut k>0 */
      BARRE;

      /* construit, teste, puis libère les tables */
      rs_dcr_tables *RT=rs_dcr(GF,k);
      routing_test(GF,RT,(rt_length)rs_agmnt_length,-1,RT->dist);
      free_rs_dcr_tables(RT);
    }break;

    case CHECK_RS_TZRPLG:{
      RS_Start("tz rplg",RS_L_FP,GF);

      /* paramètre */
      double t=MEM(CPARAM,0,double); /* t=paramètre (exposant du RPLG) */
      if((VARIANT<2)&&(0<t)&&(t<2)) Erreur(6); /* valeurs impossibles */
      if((VARIANT==2)&&(t<1)) Erreur(6); /* valeurs impossibles */
      printf("- parameter: %g\n",t);
      BARRE;

      /* construit, teste, puis libère les tables */
      rs_tzrplg_tables *RT=rs_tzrplg(GF,t);
      routing_test(GF,RT,(rt_length)rs_tzrplg_length,-1,NULL); /* routage */
      free_rs_tzrplg_tables(RT); /* libère les tables */
    }break;

    case CHECK_RS_BC:{
      RS_Start("bc",RS_L_FP,GF);

      /* paramètre */
      k=MEM(CPARAM,0,int); /* k=paramètre */
      printf("- parameter: %i\n",k); /* ici k>=0, toujours */
      if(k<0) Erreur(6); /* k=0 est possible */
      BARRE;
      
      /* construit, teste, puis libère les tables */
      rs_bc_tables *RT=rs_bc(GF,k); /* construit */
      routing_test(GF,RT,(rt_length)rs_bc_length,-1,RT->dist); /* teste */
      free_rs_bc_tables(RT); /* libère les tables */
    }break;

    case CHECK_RS_HDLBR:{
      RS_Start("hdlbr",RS_NI_FP,GF);
      
      /* paramètre */
      k=MEM(CPARAM,0,int); /* k=paramètre */
      if(k==-1) k=ceil(sqrt((double)NF)); /* ici k>=1, toujours */
      k=imin(k,NF); /* pas plus que le nombre de sommets */
      printf("- parameter: %i\n",k); /* ici k>=1, toujours */
      if(k<=0) Erreur(6); /* k<=0: valeur impossible */
      BARRE;
      
      /* construit, teste, puis libère les tables */
      rs_hdlbr_tables *RT=rs_hdlbr(GF,k);
      routing_test(GF,RT,(rt_length)rs_hdlbr_length,-1,NULL); /* routage */
      free_rs_hdlbr_tables(RT); /* libère les tables */
    }break;

    default: if(CHECK==CHECK_ON) break;
      Erreur(47); // on ne devrait jamais avoir cette erreur
    }// fin du switch(CHECK)
    
    free(CPARAM),CPARAM=NULL; /* supprime les paramètres pour CHECK */
    free(FPARAM),FPARAM=NULL; /* supprime les paramètres pour FILTER */
    if(GF!=Q->G) free_graph(GF),GF=NULL; /* supprime le graphe (évite le double-free) */
  }
  /* fin du "if(CHECK)" */

  free_query(Q),Q=NULL; /* supprime la requête */
  TopChrono(-1); /* libère tous les chronos */
  return 0; /* fin de gengraph */
}


/*# ###
Générateur de graphes - v5.1 - © Cyril Gavoille - Novembre 2017

USAGE

       gengraph [-options] graph [parameters]


DESCRIPTION

       Génère sur la sortie standard un graphe. Par défaut le graphe
       est non orienté et affiché selon une liste d'arêtes (en texte),
       mais d'autres formats sont possibles: liste d'adjacence, format
       dot de GraphViz, xfig ou pdf par exemple. En paramètre figure
       le nom du graphe ainsi que ses paramètres éventuels,
       typiquement le nombre de sommets. La commande appelée seule
       affiche l'aide sur les options du générateur. Si les paramètres
       d'une option ou d'un graphe sont absents ou remplacés par "?",
       une aide spécifique est affichée. Une console supportant l'UTF8
       est préférable.

       Ex: gengraph -help
	   gengraph -list | sort
	   gengraph tree ?
	   gengraph ? arbre
	   gengraph tutte
           gengraph hypercube 8
           gengraph mesh 7 3 -not
	   gengraph mesh 50 50 -dele .5 -maincc -visu
	   gengraph rdodecahedron -visu
           gengraph tree 100 -visu
	   gengraph web 10 3 -visu
	   gengraph gabriel 50 -caption "Grabriel with n=50" -visu
	   gengraph gabriel 2000 -xy seed 1 0.15 -visu
	   gengraph gabriel 700 -xy seed 1 -0.3 -visu
	   gengraph sierpinski 7 3 -visu
	   gengraph udg 400 .1 -visu
	   gengraph udg 400 .1 -xy seed 3 1.5 -visu
	   gengraph udg 400 -1 -vsize -vcolor deg -visu
	   gengraph arytree 6 3 3 -dot filter circo -visu
	   gengraph dyck -dot filter circo -visu
	   gengraph ringarytree 4 2 3 0 -label 1 -visu
	   gengraph arboricity 100 2 -vcolor degr -visu
	   gengraph prime 6 -directed -noloop -visu
	   gengraph aqua 3 3 2 1 -label 1 -dot filter dot -visu
	   echo "0->1->2->0" | gengraph load - -check bfs 0
	   gengraph tutte | gengraph -filter - diameter p
           gengraph rplg 300 3 -maincc -vcolor degr -vcolor pal wz -vsize -visu
           gengraph -xy box 15 15 -xy round 0 -xy grid 16 rng 30 -visu
	   gengraph linial 7 3 -check kcolorsat 3 | ./glucose -model


   LE FORMAT STANDARD

       Le format par défaut (ou standard) est une liste d'arêtes ou de
       chemins écrits en texte simple. Ce format minimaliste est très
       proche de celui du format "dot" de GraphViz.  D'autres formats
       de sortie sont possibles, notamment le format "dot" (voir
       l'option -format). Les sommets sont numérotés consécutivement
       de 0 à n-1 où n est le nombre de sommets présents dans le
       graphe (en fait cela peut être changé avec l'option -shift).
       Une arête entre i et j est représentée par i-j, un arc de i
       vers j par i->j. Les sommets isolés sont simplement représentés
       par le numéro du sommet suivit d'un espace ou d'un retour de
       ligne. Le nombre de sommets du graphe est l'entier le plus
       grand + 1, et s'il y a i-j (ou i->j), alors il existe une arête
       (ou un arc) entre les sommets i et j.

       Pour une représentation plus compacte, les arêtes (ou arcs)
       consécutives d'un chemin du graphe peuvent être regroupées en
       blocs i-j-k-…. Par exemple, les deux arêtes 3-5 et 5-8
       peuvent être regroupées en 3-5-8. Mais ce n'est pas
       obligatoire.  Également, les arêtes (ou arcs) d'une étoile
       peuvent être groupées avec i-(j k …). Par exemple, 3-(5 7 8)
       représente les arêtes 3-5, 3-7 et 3-8. Il n'est pas possible
       cependant de combiner chemins et étoiles, comme 3-(5-7-8) ou
       3-(5-(7 8)). Toutefois 3-5-(7 …) est correct, mais pas 3-(5
       6)-7 ni (3 5)-6. Les sommets isolés et les arêtes (ou les blocs
       d'arêtes) sont séparés par des espaces ou des sauts de ligne.
       Une boucle sur un sommet i est codée par i-i. Les arêtes
       multiples sont codées par la répétition d'une même arête, comme
       par exemple i-j i-j, ou encore i-j-i (même convention pour les
       arcs i->j->i).

       Ex: 0 1-2-3-1            0   1
                                   / \
                                  3---2

       représente un graphe à 4 sommets, composé d'un sommet isolé (0)
       et d'un cycle à trois sommets (1,2,3). Une représentation
       graphique possible est donnée à droite.

       Ex: 4-2-1-0-3-2-5

       représente un graphe à 6 sommets composé d'un cycle de longueur
       4 et de deux sommets de degré 1 attaché à 2. On aurait pu coder
       le même graphe avec l'expression 2-(1 3 4 5) 1-0-3. En voici
       une représentation graphique possible:

       Ex:        1
                 / \
              4-2   0
               / \ /
              5   3

       Plus généralement, une famille de graphes peut être définie en
       précédant chaque graphe par "[n]" où n est un entier unique
       représentant l'identifiant du graphe.

       Ex: [17] 0-1 [22] 0->1->2->0

       représente une famille composée de deux graphes: un chemin à
       deux sommets (d'identifiant 17) ainsi d'un cycle orienté à
       trois sommets (d'identifiant 22).

   COMMENT FONCTIONNE LE GÉNÉRATEUR ?

       Pour chaque graphe une fonction adj(i,j) est définie. Elle
       fournit d'adjacence (0 ou 1) entre les sommets i et j, des
       entiers entre 0 et n-1. Le graphe est affiché en générant
       toutes les paires {i,j} possibles et en appelant adj(i,j) (ou
       tous les couples (i,j) possibles dans le cas orienté). Les
       graphes sont ainsi générés de manière implicite. Les arêtes du
       graphe ne sont pas stockées en mémoire, mais affichées à la
       volée. Ceci permet de générer des graphes de très grande taille
       sans nécessiter O(n²) espace de mémoire centrale. Pour certains
       graphes cependant, comme les arbres, certains graphes
       d'intersections ou les graphes géométriques, une structure de
       données en O(n) peut être utilisée. Pour les formats
       d'affichage liste, matrix, et smatrix une structure de données
       de taille linéaire (en O(n+m) où m est le nombre d'arêtes) est
       utilisée en interne. Ces trois derniers formats sont donc à
       éviter. Pour la génération de très grand graphe, le format
       standard ou dot doit être privilégié.

   COMMANDES EXTERNES

       Le programme fait appel, pour certaines fonctions, aux
       commandes systèmes suivantes qui doivent être installées: sed,
       grep, awk, more, sort, dot.


OPTIONS


....-help [word], ? [word], or [option|graph] ?
....
       Affiche l'aide en ligne qui est contenue dans le fichier source
       du générateur. Pour cela, le code source .c doit être dans le
       même répertoire que l'exécutable. Si "word" est précisé, alors
       les options et noms de graphe contenant "word" sont affichés.
       La variante "[option|graph] ?" affiche une aide détaillée sur
       une option ou un graphe précis.
....
       Ex: gengraph ? arbre
           gengraph ktree ?
	   gengraph ? hedron
	   gengraph ? planaire
....
       La forme ? peut ne pas fonctionner correctement si un fichier
       d'un seul caractère existe dans le répertoire courant (à cause
       de l'interprétation du shell). Il faut alors utiliser '?' au
       lieu de ?.

....-list
....
       Affiche la liste des graphes et leurs paramètres qui sont
       possibles d'être générée. Sont listés d'abord les graphes de
       bases puis les composés. On obtient une aide sur un graphe
       particulier si son nom est suivit de " ?" ou si ses paramètres
       éventuelles sont absents (dans ce cas il doit être le dernier
       mot de la commande).
....
       Ex: gengraph -list | sort
	   gengraph gabriel

....-version
....
       Affiche la version courante du générateur (en fait du programme
       source), un réel > 1. Pour cela, le code source .c doit être
       dans le même répertoire que l'exécutable.

....-directed
....-undirected
....
       L'option -directed produit le graphe orienté en testant les n²
       arcs possibles, l'option -undirected permettant de revenir à la
       situation par défaut, soit en testant les n(n-1)/2 arêtes. Les
       boucles peuvent être testées ou pas grâce aux options -loop et
       -noloop. En format standard ou dot, un arc apparaît comme i->j
       au lieu de i-j ou i--j pour une arête. Tous les graphes ne sont
       pas forcément définis pour fonctionner comme espéré avec
       l'option -directed car certaines fonctions d'adjacence
       supposent en effet que i<j. La plupart des graphes vont
       apparaître comme orientés symétriques.
....
       Ex: gengraph clique 5 -directed
           gengraph cycle 5 -directed -visu

....-noloop
....-loop
....
       L'option -noloop permet de ne pas produire les boucles des
       graphes, alors que -loop les autorise. L'option -noloop est
       celle par défaut pour les graphes non-orientés, et l'option
       -loop celle par défaut pour les graphes orientés. Ces options
       doivent être placées après -(un)directed car modifiée par
       celle-ci.

....-not
....
       Inverse la fonction d'adjacence, et donc affiche le complément
       du graphe. Cette option est prioritaire sur l'option -redirect.

....-dele p
....
       Permet de supprimer chaque arête du graphe générée avec probabilité p.

....-delv p
....
       Similaire à -dele p mais concerne les sommets. Le sommet et ses
       arêtes incidentes sont alors supprimés. Si p est un entier <0,
       alors exactement |p| sommets sont supprimés. Si k sommets sont
       supprimés, alors le nom des sommets restant est dans
       l'intervalle [0,n-k[ où n est le nombre de sommets initial du
       graphe. Les noms des sommets sont donc éventuellement
       renumérotés. Voir aussi les options -permute et -shift. Bien
       sûr la fonction d'adjacence adj(i,j) est appliquée sur les noms
       (i,j) originaux.

....-star n                          [N'EST PLUS EFFECTIF DEPUIS v4.5]
....
       Ajoute n sommets pendant (degré 1) aux sommets du graphe. Si
       n>0, alors n représente le nombre total de sommets ajoutés,
       chacun des n sommets étant connectés aléatoirement uniformément
       aux sommets du graphe original. Si n<0, alors |n| sommets sont
       ajoutés à chacun des sommets du graphe. Cette opération est
       appliquée en priorité et ne peut être appliquée qu'une fois par
       commande. Les options -star et -apex sont incompatibles.

....-apex n                          [N'EST PLUS EFFECTIF DEPUIS v4.5]
....
       Ajoute n sommets universels, donc connectés à tous les sommets
       du graphe. Cette opération est appliquée en priorité. Les
       options -star et -apex sont incompatibles.

....-redirect p
....
       Redirige chaque arête uniformément avec probabilité p. Plus
       précisément, si {i,j} est une arête du graphe original G, alors
       avec probabilité p l'arête affichée est {i,k} au lieu de {i,j}
       où k est un sommet choisi uniformément parmi les sommets du
       graphe G. Si l'arête {i,j} est supprimée par -dele ou si le
       sommet i est supprimé par -delv, la redirection n'a pas lieu.
       Cette option est appliquée après l'option -not. Le graphe G
       tient donc compte de -not avant de rediriger ses arêtes.

....-seed s
....
       Permet d'initialiser le générateur aléatoire avec la graine s,
       permettant de générer plusieurs fois la même suite aléatoire.
       Par défaut, la graine est initialisée par une combinaison du
       numéro de processus de la commande et le temps, donc génère par
       défaut des suites différentes à chaque lancement. Le générateur
       est initialisé lorsque l'option est lue sur la ligne de
       commande. Le comportement du programme peut donc être affecté
       suivant son ordre d'apparition. Cependant le graphe est généré
       après l'analyse de la ligne de commande.

....-width m
....
       Limite à m le nombre d'arêtes et de sommets isolés affichés par
       ligne. Cette option n'a pas de signification particulière en
       dehors des formats standard et dot. Par exemple, -width 1
       affiche une arrête (ou un sommet isolé) par ligne. L'option
       -width 0 affiche tout sur une seule ligne. La valeur par défaut
       est 12.

....-shift s
....
       Permet de renuméroter les sommets à partir de l'entier s
       positif. La valeur par défaut est -shift 0.  L'intérêt de cette
       option est de pouvoir réaliser des unions de graphes simplement
       en renumérotant les sommets et en concaténant les fichiers aux
       formats standard ou list. Cette option n'a pas d'effets pour
       les formats de sortie de type matrice.

....-permute
....
       Permute aléatoirement uniformément le nom des sommets
       lorsqu'ils sont affichés (ou chargés en mémoire dans le cas
       d'options -check). Les numéros restent dans l'intervalle
       initial qui, sauf si l'option -shift a été utilisée, est [0,n[
       où n est le nombre de sommets du graphe réellement généré. Voir
       aussi l'option -label.

....-header
....
       Affiche un préambule donnant certaines informations sur le
       graphe, sous forme de commentaire à la C++ (//). Par défaut
       aucun préambule n'est affiché. Les informations affichées sont:
       - l'heure, la date et la graine du générateur aléatoire
       - la ligne de commande qui a produit la génération du graphe
       - le nombre de sommets, d'arêtes, le degrés maximum et minimum
       Pour les formats standard et dot, le nombre d'arêtes (et les
       degrés min et max) n'est pas déterminé avant l'affichage du
       graphe. Pour cette raison ces nombres ne sont affichés qu'après
       le graphe. Pour n'avoir que les informations sur le graphe,
       utiliser -header avec l'option -format no. Voir aussi -check info.

....-caption title
....
       Permet d'ajouter une légende à un graphe. Cette option n'a
       d'effet qu'avec le format dot et ces variantes. Il est possible
       d'affiche la "seed" avec le format %SEED. On ne peut avoir plus
       d'une occurrence du même format (%SEED) dans cette option. 
....
       Ex: gengraph gabriel 30 -caption ex1 -visu
           gengraph gabriel 30 -caption "Exemple 2" -visu
           gengraph gabriel 30 -caption "graph with seed=%SEED" -visu

....-fast
....
       Génère le graphe sans tester les O(n²) arêtes possibles, mais
       en utilisant directement la liste d'adjacence du graphe
       préalablement générée lors de l'initialisation du graphe, comme
       c'est le cas pour le graphe "load". Le résultat est une
       génération du graphe en temps O(n+m) au lieu de O(n²).
       L'utilisation typique est (voir aussi "loadc"):
....
       Ex: gengraph load file -fast -delv 0.3 -check ncc
....
       Cet exemple permet de calculer le nombre de composantes
       connexes sur un sous-graphe contenu dans un fichier le tout en
       temps linéaire. D'autres graphes peuvent supporter une
       génération rapide si elle est implantée dans l'initialisation
       de la fonction d'adjacence (pour l'instant seul le graphe
       "load" le supporte).  Certaines options, comme -not, n'ont
       alors plus d'effet en présence de -fast. Cependant, -permute,
       -delv, -dele, et d'autres fonctionnent normalement.

....-variant v
....
       Permet de passer un entier v>0 pour contrôler certaines
       fonctionnalités du générateur. Par exemple, "-variant 1 -check
       routing cluster -1" permettra de calculer une variante du
       schéma de routage "cluster".

....-check [parameters]
....
       Stocke en mémoire le graphe généré sous la forme d'une liste
       d'adjacence, et lui applique un algorithme. Le graphe est
       généralement affiché, sauf pour certaine options comme -check
       maincc ou -check routing. Utiliser "-format no" pour ne pas
       afficher le graphe généré. Cette option nécessite un espace
       supplémentaire en O(n+m) pour le stockage du graphe.
....
       -check info
....
          Affiche quelques caractéristiques du graphe, après avoir
          effectué un tri puis un simple parcours de ses listes
          d'adjacences. Indique, par exemple, si le graphe est
          orienté, s'il contient des boucles, des multi-arêtes,
          etc. Le graphe lui-même n'est pas affiché. Indique aussi
          l'occupation mémoire du graphe, le temps de génération ou
          de chargement et le temps de parcours.
....
       -check simplify
....
          Permet de simplifier un graphe en supprimant les boucles et
          multi-arêtes, ceci par un simple parcours du graphe après un
          tri de ses listes d'adjacences. Les arêtes sont affichées au
          format standard une par ligne, comme avec -width 1 -format
          standard. Cette option est un moyen rapide de sortir un
          graphe au "bon" format. Elle est insensible aux options
          d'affichage, en particulier -format. Il est possible de ne
          pas supprimer les boucles en utilisant -loop.
....
          Ex: gengraph loadc G -check simplify > H
	      gengraph load G -fast -permute -check simplify > H
....
	  La différence avec "gengraph load G > H", qui produit aussi
	  une sortie valide, est qu'avec "gengraph loadc G -check
	  simplify > H" le graphe G est lu et affiché en temps
	  quasi-linéaire sans la (re)génération de ses arêtes en temps
	  O(n²).
....
       -check bfs s
....
          Effectue un parcours en largeur d'abord sur le graphe généré
          depuis le sommet s. La distribution des distances depuis s
          est affichée, ainsi que l'arborescence (-1 indique que le
          sommet n'a pas de père). La longueur du plus petit cycle
          passant par s est aussi donnée. Elle vaut -1 s'il n'existe
          pas.
....
       -check bellman s
....
          Calcule les plus courts chemins depuis le sommet s par
          l'algorithme de Bellman-Ford. Si le graphe est géométrique,
          le poids de chaque arête correspond à la distance
          euclidienne entre ses extrémités, sinon il vaut 1 et le
          résultat sera similaire à un bfs. Dans le cas géométrique,
          l'étirement maximum depuis s est calculé, ainsi qu'un chemin
          le réalisant.
....
       -check stretch
....
          Calcule, comme -check bellman s, l'étirement d'un graphe
          géométrique depuis chaque source s. On affiche une source
          atteignant l'étirement maximum, mais aussi une source
          atteignant l'étirement minimum, ainsi qu'un chemin réalisant
          ces étirements. Utilisez l'option "-format no" pour ne pas
          avoir l'affichage de la génération du graphe.
....
       -check volm
....
          Calcule la distribution du volume monotone des sommets. Le
          volume monotone d'un sommet u est le nombre d'arcs du
          sous-graphe des sommets accessibles depuis u dans le graphe
          où seuls les arcs u->v avec u<v ont été gardés. Cette mesure
          dépend de la numérotation des sommets. À utiliser en
          combinaison dans l'option -permute.
....
       -check dfs s
....
          Effectue un parcours en profondeur d'abord de toutes les
          composantes connexes du graphe généré depuis le sommet s. Le
          nombre de composantes ainsi que l'arborescence (-1 indique
          que le sommet n'a pas de père) sont donnés.
....
       -check ncc
       -check connected
....
          Donne le nombre de composantes connexes ainsi que le nombre
          de cut-vertex du graphe. Ces informations sont aussi
          affichées par -check dfs 0.
....
       -check deg
       -check edge
       -check edges
....
          Affiche la distribution des degrés et le nombre d'arêtes du
          graphe.
....
       -check degenerate
....
          Donne la dégénérescence du graphe, ainsi que l'ordre
          d'élimination correspondant des sommets.
....
       -check girth
....
          Donne la maille du graphe dans le cas non orienté. La valeur
          -1 est renvoyée si le graphe est acyclique, et la valeur 0
          dans le cas orienté.
....
       -check diameter
....
          Calcule le diamètre du graphe généré. Affiche +∞ pour un
          graphe non connexe.
....
       -check radius
....
          Calcule le rayon du graphe généré, soit la hauteur du plus
          arbre couvrant. Affiche +∞ pour un graphe non connexe.
....
       -check gcolor
....
          Donne une borne supérieure sur le nombre chromatique du
          graphe en utilisant l'heuristique du degré minimum.
....
       -check kcolor k
....
          Donne une k-coloration du graphe (et la couleur pour chaque
          sommet), si c'est possible. Pour cela une recherche
          exhaustive de toutes les k-colorations est effectuée. Le
          temps est raisonnable si k=3 et n<20.
....
       -check kcolorsat k
....
          Donne une formulation SAT de la k-coloration du graphe. Il
          s'agit de la formulation multivaluée classique, un sommet
          pouvant avoir plusieurs couleurs sans que cela nuise à la
          validité du résultat. Les contraintes sont décrites au
          format Dimacs CNF. On peut alors envoyer le résultat à un
          solveur SAT comme MiniSat ou Glucose. Le graphe n'est pas
          affiché, et donc "-format no" n'est pas nécessaire.
....
          Ex: gengraph linial 6 3 -check kcolorsat 3 | ./glucose -model
....
       -check kindepsat k
....
          Donne une formulation SAT d'un ensemble indépendant de
          taille k du graphe. Les variables i=1 à n indiquent si le
          sommet numéroté i-1 est dans la solution ou pas. Les
          contraintes sont décrites au format Dimacs CNF. On peut
          alors envoyer le résultat à un solveur SAT comme MiniSat ou
          Glucose. Le graphe n'est pas affiché, et donc "-format no"
          n'est pas nécessaire.
....
          Pour le problème clique de taille k, il suffit de chercher
          un ensemble indépendant de taille k pour le complément du
          graphe. Et pour le problème "vertex cover" de taille k,
          c'est un ensemble indépendant de taille n-k sur le
          complémentaire qu'il suffit de chercher.
....
       -check ps1
       -check ps1b
       -check ps1c
       -check ps1x n u_1 v_1 … u_n v_n
....
          Applique le test ps1 ou l'une de ses variantes (voir -filter
          ps1 pour plus de détail sur ce test). Affiche aussi le
          nombre de tests réalisés (nombre de paires de sommets et de
          chemins testés).
....
       -check paths x y
....
          Liste tous les chemins simples entre les sommets x et
          y. N'affiche rien si x et y ne sont pas connectés. L'ordre
          est défini suivant le premier plus court chemins dans
          l'ordre des sommets depuis le sommet x.
....
       -check iso H
....
          Teste si le graphe généré G est isomorphe à H. Si oui,
          l'isomorphisme de G à H est donné. Le nombre de tests
          affichés est le nombre de fois où les graphes sont comparés,
          la comparaison prenant un temps linéaire en la taille des
          graphes. Plus les graphes sont symétriques (comme un cycle
          ou un hypercube), plus le nombre de tests sera important.
....
          Tester l'isomorphisme entre deux cycles de 8 sommets
          étiquetés aléatoirement prends environ 4 mille tests, et
          entre deux cycles de 12 sommets, 30 millions de tests soit
          9" environ. Pour deux arbres à 75 sommets (aléatoires mais
          isomorphes), moins de 20 tests suffisent.
....
       -check sub H
....
          Teste si le graphe généré G est un sous-graphe couvrant de H
          (donc avec le même nombre de sommets). S'ils ont le même
          nombre d'arêtes, le test est équivalent à l'isomorphisme. Le
          nombre de tests est le nombre total de fois où deux graphes
          sont comparés.  On peut tester si H est Hamiltonien en
          prennant pour G un cycle.
....
          Tester un cycle de longueur 12 dans une grille 3x4 prend
          jusqu'à environ 32 millions de tests (parfois bien moins),
          soit au plus 10".
....
       -check minor H
....
          Teste si le graphe G généré contient H comme mineur. Les
	  graphes peuvent être non connexes. S'ils ont le même nombre
	  de sommets le test est équivalent à celui du sous-graphe
	  (voir -check sub). Dans le cas positif, un modèle de H dans
	  G est fourni.
....
          Le principe consiste à contracter des arêtes de G, de toutes
	  les manières possibles, et à tester si H est un sous-graphe
	  du graphe contracté. Le nombre de tests affichés est le
	  nombre de contractions plus le nombre total de tests
	  réalisés par les tests de sous-graphe. Pour H=K₄ il est
	  préférable d'utiliser -check twdeg qui donne < 3 ssi le graphe
	  ne contient pas K₄ comme mineur.
....
       -check twdeg
....
          Donne une borne supérieure et inférieure sur la treewidth du
          graphe. Pour la borne supérieure, on utilise l'heuristique
          du sommet de degré minimum que l'on supprime et dont on
          complète le voisinage par une clique. En cas d'égalité (même
          degré) on sélectionne le sommet dont il faut rajouter le
          moins d'arêtes. La borne inférieure qui est donnée provient
          de la dégénérescence. La treewidth est exacte si 0,1 ou 2
          est retournée. L'algorithme est en O(n²).
....
       -check tw
....
          Calcule la treewidth du graphe en analysant tous les ordres
          d'éliminations. La complexité est donc en n!. Il ne faut
          l'utiliser que si le nombre de sommets est < 12 (Ex:
          gengraph random 12 .5 -check tw donne 5 en environ 750
          millions de tests). Parfois, l'utilisation de -permute peut
          accélérer le traitement, car partir d'un ordre d'élimination
          déjà bon permet d'éliminer rapidement beaucoup d'ordres
          possibles.
....
       -check maincc
....
          Affiche, dans le mode standard seulement, le graphe
          correspondant à la composante connexe ayant le plus grand
          nombre de sommets. Le graphe initial n'est pas affiché. Les
          sommets sont renumérotés si le graphe initial n'était pas
          connexe. Attention ! l'affichage de la composante n'est
          sensible qu'à l'option -width. En particulier il n'est pas
          possible d'afficher la composante dans un autre format
          (-format) ou avec les noms originaux (-label). Cependant,
          avec "-check maincc | ./gengraph load -" on peut afficher le
          graphe dans le format souhaité, ou ajouter -visu. (Voir
          aussi le raccourcis -maincc.) Notez que "-check maincc
          -visu" provoque une erreur, car -visu applique l'option
          "-format dot" incompatible avec -check maincc.
....
       -check routing [hash h] [scenario [nomem] s] scheme [parameters]
....
          Construit les tables de routage pour le graphe selon le
          schéma de routage "scheme", ce schéma pouvant comporter des
          paramètres spécifiques. La sortie consiste en statistiques
          sur les tables (taille, temps de calcul) et le graphe.
          L'option "scenario" permet en plus de tester certains types
          (s) de routage sur le graphe et d'afficher des statistiques
          sur les longueurs de routes générées (dont l'étirement).
          L'option "hash" permet de préciser la fonction de hashage
          (h) appliquée le cas échéant aux sommets souvent utilisé
          dans les schémas de type "name-independent". Le graphe doit
          être connexe et comporter au moins une arête, propriétés qui
          sont toujours testées.
....
          Ex: gengraph -permute rplg 200 2.3 -maincc > G1
	      gengraph loadc G1 -check routing scenario all cluster -1
....
          L'option "nomem" après "scenario" permet d'optimiser la
          mémoire en ne stockant pas les distances calculées pour
          établir l'étirement (par défaut elles le sont). En
          contrepartie le temps de calcul est allongé. Cela peut avoir
          un intérêt si le graphe est très grand (n ≃ 1.000.000
          sommets) et qu'un scenario comme "pair -1000" est testés.
          Les scenarii possibles sont (n=nombre de sommets du graphe):
....
          • scenario none     → aucun routage (scenario par défaut)
          • scenario all      → n(n-1) routage possibles
          • scenario edges    → tous les routages entre voisins
          • scenario npairs   → n paires aléatoires de sommets différents
          • scenario one u    → n-1 routages depuis u (choix aléatoire si -1)
          • scenario pair u v → routage de u à v (choix aléatoire si -1)
          • scenario pair -p  → routage depuis p>0 paires aléatoires
          • scenario until s  → routage jusqu'à l'étirement s ou plus
....
          Les fonctions de hachages h:[0,n[->[0,k[ possibles sont:
	  (shuffle et mod atteignent le nombre de collisions minimum de ⎡ n/k⎤)
....
          • hash prime   → h(x)=((a*x+b)%p)%k où 0<a,b<p sont aléatoires
	                   et p=2^31-1 est premier (hash par défaut)
          • hash mix     → h(x)=mix(a,b,x)%k où a,b sont aléatoires sur 32 bits
	                   et mix() la fonction mélange de Bob Jenkins (2006)
          • hash shuffle → h(x)=π(x)%k où π(x) est une permutation de [0,n[
	                   basée sur deux entiers aléatoires de [0,n[.
          • hash mod     → h(x)=(x+a)%k où a est aléatoire dans [0,k[.
....
       -check routing cluster k
....
          Schéma de routage name-independent "cluster" de paramètre k.
          Un sommet de degré maximum est choisi comme "centre", puis
          k-1 de ces voisins de plus haut degrés sont choisis pour
          former un cluster de taille au plus k. Un arbre BFS est
          enraciné depuis le centre. Chaque sommet possède une boule
          par rayon croissant qui s'arrête avant de toucher un sommet
          du cluster. On route de u vers v d'abord dans la boule de
          u. Sinon on va dans le cluster pour chercher un sommet du
          cluster responsable du hash de v. Une fois atteint on route
          selon l'arbre BFS ou selon un plus court chemin si la
          distance à la destination est ≤ logn/loglogn. Les sommets
          ayant un voisin dans cluster et qui ne sont pas eux-mêmes
          dans le cluster possèdent dans leur table tout leur
          voisinage. Les boules sont optimisées par l'usage d'un
          voisin par défaut.
....
	  Si k=-1, alors k est initialisé à ⎡ √n⎤. Si k=-2 il est
          initialisé à n si bien que le cluster est fixé à tout le
          voisinage du centre. Dans tous les cas, l'étirement est
          toujours ≤ 5. Il est même ≤ 3 si k=1. La taille des tables
          est réduit dans le cas de graphes power-law (comme RPLG). Il
          existe plusieurs variantes (option -variant v) où chacun des
          bits de v mis à 1 est interprété comme suit:
....
          • bit-0: le routage est réalisé sans les boules de
            voisinages ce qui réduit la taille moyenne mais augmente
            l'étirement maximum.
....
          • bit-1: le routage dans le cluster est réalisée dans
            l'étoile couvrant le cluster, sans les autres arêtes du
            cluster.  Cela réduit la taille des tables sans modifier
            l'étirement maximum. Si de plus k=1, le routage est alors
            réalisé via la racine de l'arbre BFS ce qui réduit au
            minimum la taille moyenne des tables (2 en moyenne).
....
          • bit-2: les tables des sommets voisins du cluster sont
	    vides.  L'étirement devient ≤ 7 au lieu de 5 au maximum,
	    mais la taille des tables est réduite.
....
          • bit-3: les tables des sommets voisins du cluster ne
            contiennent que des sommets qui ne sont pas de le cluster.
            L'étirement ≤ 5 est préservé, mais généralement la taille
            maximum des tables est réduite, l'étirement moyen
            augmentant que très légèrement.
....
       -check routing dcr k
....
          Schéma de routage name-independent "dcr" de paramètre k>0
          représentant le nombre de couleurs. C'est une simplification
          du schéma "agmnt". L'étirement est toujours ≤ 5 et le nombre
          d'entrées des tables est en moyenne f(k,n) = 2n/k +
          (k-1)*(H(k)+1) où H(k) ~ ln(k)+0.577… est le k-ième nombre
          harmonique. Le principe du schéma est le suivant.  Chaque
          sommet possède une couleur, un entier aléatoire de [0,k[,
          les sommets landmarks étant ceux de couleur 0. Les boules de
          voisinages des sommets sont définies par volume comme la
          plus petite boule contenant au moins chacune des couleurs
          (ou tous les sommets s'il manque une couleur n'est pas
          représentée), les sommets du dernier niveau étant ordonnés
          par identifiant croissant. Le routage s→t s'effectue selon
          un plus court chemin si t est dans la boule de s ou si t est
          un landmark.  Sinon, on route vers le sommet w de la boule
          de s dont la couleur est égale au hash de t, une valeur
          aussi dans [0,k[. Puis le routage w→t s'effectue dans
          l'arbre BFS enraciné dans le plus proche landmark de s ou de
          t, celui minimisant la distance de w à t.
....
	  Si k=-1, alors k est initialisé à sa valeur optimale
          théorique, celle qui minimise le nombre moyen d'entrées
          f(k,n), valeur calculée numériquement et qui vaut environ k
          ≃ √(n/ln(n))/2, ce qui donne environ 2√(n*ln(n*ln(n)))
          entrées en moyenne.  Si k=-2, le nombre de couleur est
          initialisé à n. Les valeurs de k>n sont possibles, il s'agit
          alors d'un routage de plus court chemins comme pour le cas
          k=n ou k=1. Il existe une variante (-variant 1) lorsque k>1
          qui a pour effet de choisir les landmarks comme les sommets
          de plus haut degré. Plus précisément, les ⎡ n/k⎤ sommets de
          plus haut degré sont coloriés 0 et les autres coloriés
          aléatoirement dans [1,k[. Dans cette variante, la borne sur
          l'étirement est toujours garantie mais plus sur le nombre
          maximum d'entrées.  Cependant, pour certains graphes
          l'étirement moyen est amélioré.
....
       -check routing agmnt k
....
          Schéma de routage name-independent dit "agmnt" du nom de ces
	  auteurs Abraham et al. (2008). C'est la version originale du
	  schéma "dcr" qui diffère par l'algorithme de routage. Comme
	  dans "dcr", le routage s→t s'effectue directement vers t si
	  t est dans la boule de s ou est un landmark. Sinon, vers le
	  sommet w de la boule de s dont la couleur est égale au hash
	  de t. Le routage w→t s'effectue suivant la meilleure des
	  options suivantes: router via un arbre BFS d'un des
	  landmarks; ou bien router via un arbre couvrant la boule
	  d'un sommet s' contenant à la fois w et un sommet x voisin
	  d'un sommet y contenu dans la boule de t (les boules s' et
	  de t, si elles existent, sont dites "contiguës via l'arête
	  x-y"). L'étirement est toujours ≤ 3 et les tables ont le
	  même nombre d'entrées que "dcr", bien que plus complexes. Le
	  temps de calcul des tables est plus important que pour
	  "dcr". Toutes les variantes de "dcr" (-variant, k<0)
	  s'appliquent aussi à "agmnt".
....
       -check routing tzrplg t
....
          Schéma de routage étiqueté inspiré de celui de Thorup &
          Zwick (2001) optimisé pour les Random Power-Law Graphs (voir
          prlg) de paramètre réel t (power-law exponent) et proposé
          par Sommer et al. (2012). L'étirement est toujours ≤ 5. Les
          valeurs de t entre ]0,1.5] sont interdites.  Le schéma
          utilise des sommets landmarks où des arbres BFS sont
          enracinés, ainsi que des boules (de voisinage) définies par
          rayons croissant qui s'arrête avant de toucher un
          landmark. Le routage s'effectue alors en priorité via les
          boules ou alors via le landmark le plus proche de la
          destination (sans raccourcis), information précisée dans
          l'étiquette de la destination. Les landmarks sont les
          sommets de plus haut degré. Par défaut (-variant 0) leur
          nombre vaut:
....
            • si t>1.5, ⎡ n^((t-2)/(2t-3))⎤
            • si t=0,   ⎡ √n⎤
            • si t<0,   |t|
....
	  Si -variant 1 et t>1.5, alors les landmarks sont tous les
	  sommets de degré > n^1/(2t-3). Si -variant 2 et t>0, alors
	  les landmarks sont t sommets choisis aléatoirement.
	  L'étirement est ≤ 3 si un seul landmark est choisi.
....
       -check routing hdlbr k
....
          Schéma de routage name-independent HDLBR selon Tang et
          al. (2013) avec k>0 landmarks qui sont les sommets de plus
          haut degré.  Si k<0, alors k est initialisé à ⎡ √n⎤. Chaque
          sommet qui n'est pas un landmark possède une boule dont le
          rayon est juste inférieure au plus proche landmark. Chaque
          sommet possède sa boule boule inverse (qui peut être vide),
          définie comme l'ensemble des sommets le contenant dans leur
          boule. Chaque landmark a une couleur unique. On route
          directement de u à v si v est dans la boule de u, la boule
          inverse de u, ou si v est un landmark. Sinon, on route vers
          le landmark (selon un plus court chemin) dont sa couleur
          vaut le hash de v. Delà on route suivant un plus court
          chemin vers le plus proche landmark de v, l(v). On utilise
          ensuite le next-hop de l(v) vers v, un sommet nécessairement
          contenant v dans sa boule inverse. A chaque étape du
          routage, si v est dans la boule ou boule inverse du sommet
          courant, on route directement vers celui-ci. La longueur de
          route entre u et v est au plus 2d(u,v)+2r, où r est la
          distance maximum entre deux landmarks. La valeur de r est
          bornée par une constante pour k ≃ √n et pour les Random
          Power-Law Graphs (voir rplg).
....
       -check routing bc k
....
          Schéma de routage étiqueté selon Brady-Cowen (2006).
          L'étirement est additif ≤ 2k, et donc multiplicativement ≤
          2k+1. En particulier, si k=0, il s'agit d'un routage de plus
          court chemin. Il est adapté aux Power-Law Random Graphs
          (voir plrg). Le principe est le suivant: on construit un
          arbre BFS (=T) enraciné dans un sommet (=r) de plus haut
          degré. Le coeur (=C) est la boule de rayon k depuis r. On
          construit une liste (=L) de BFS couvrant G ainsi qu'une
          forêt (=H) de BFS de G comme suit. Au départ, L={T}, et H
          est la forêt T\C. Puis, pour chaque arête {u,v} de G\C\T, on
          vérifie si l'ajout de {u,v} à H crée un cycle ou pas. Si
          c'est non, on met à jour la forêt H en lui ajoutant
          {u,v}. Si c'est oui, on calcule un BFS de G de racine u (ou
          v) qu'on ajoute à L.  Une fois le calcul de H terminé, on
          calcule une forêt BFS couvrante de H qu'on ajoute à L.
          L'algorithme de routage de u à v consiste simplement à
          router dans l'arbre de L qui minimise la distance de u à v.

....-filter family[:range] [not] test [parameters]
....
       Affiche les graphes d'une famille pour lesquels le test est
       vrai (ou faux si "test" est précédé de "not"). Le paramètre
       "family" est un nom de fichier ou "-" pour l'entrée standard.
       La lecture de la famille se fait en temps linéaire.  Il
       contient la famille de graphes (ou un graphe seul) au format
       standard. La première ligne affichée contient le nombre de
       graphe de la famille résultante.  L'affichage de chaque graphe
       est influencé par l'option -width qui doit être placée avant
       -filter. La variante "family:range" permet de sélectionner les
       graphes de la famille dont les identifiants sont spécifiés par
       "range", comme par exemple "family:5-8" qui sélectionne les
       graphes d'identifiant 5,6,7,8. De manière générale, "range" est
       un ensemble de valeurs selon le format "value" décrit ci-après
       (voir aussi -filter F id value). La variante "-:range" est
       possible.
....
       Dans la suite, la présence de "value" dans les paramètres d'une
       option représente un ensemble de valeurs possibles. Par
       exemple, -filter F vertex '>5' filtre les graphes de la famille
       F comportant plus de 5 sommets. De manière générale, "value"
       est une suite d'intervalles d'entiers séparés par des ","
       (interprétées comme "ou"), chaque intervalle étant codé comme
       suit:
....
          • <x       →  valeur inférieure à l'entier x
          • >x       →  valeur supérieure à l'entier x
          • x ou =x  →  valeur égale à l'entier x
          • x-y      →  valeur dans l'ensemble d'entiers {x,…,y}
          • t        →  toujours vrai (intervalle infini)
          • p        →  affiche la valeur plutôt que le graphe
....
       Ex: -filter F vertex '5,7-13,>100'
           -filter F vertex '5-10,p'
           -filter F edge p
           -filter F id 5,7
           -filter F vertex p | head -1
           -filter … p | grep '^\[' | sort -rnk 3 | head
           -filter … p | grep '^\[' | sort -nk 3 | head
....
       Le premier exemple filtre les graphes de la famille F ayant un
       nombre de sommets n vérifiant soit n=5, soit 7 ≤ n ≤ 13, ou
       soit n>100. L'exemple 2 affiche le nombre de sommets des
       graphes ayant entre 5 et 10 sommets.  L'exemple 3 affiche le
       nombre d'arêtes de chaque graphe.  L'exemple 4 affiche les
       graphes d'identifiant 5 et 7 de la famille F. L'exemple 5
       affiche le nombre de graphes de la famille (ce qui correspond à
       la première ligne de commentaires). Les deux derniers exemples
       permettent d'avoir le maximum/minimum de "value".
....
       Si "value" contient le symbole > ou < il est alors préférable
       de mettre des quotes ('>14' par exemple) pour que la commande
       soit correctement interprétée par le shell.
....
       La différence principale avec -check est que le résultat de
       -filter est non verbeux alors que -check, qui ne s'applique pas
       a priori sur des familles de graphes mais sur un graphe seul,
       peut donner des explications sur l'exécution de l'algorithme
       dont la sortie n'est pas forcément un ou une liste de
       graphes. Aussi, avec -check l'algorithme s'applique au graphe
       généré, donc après un temps a priori en O(n²), alors qu'avec
       -filter c'est toujours à partir d'un fichier (ou de l'entrée
       standard), lu en temps linéaire.
....
       -filter F id value
....
          Filtre les graphes de F dont l'identifiant est déterminé par
          value. Cela permet d'extraire un ou plusieurs graphes
          donnés. C'est équivalent à "-filter F:value all".
....
       -filter F rename shift
....
          Affiche tous les graphes de la famille en renumérotant les
          graphes à partir de l'entier "shift".
....
       -filter F vertex value
....
          Filtre les graphes de F ayant un nombre de sommets déterminé
          par value.
....
       -filter F edge value
       -filter F edges value
....
          Filtre les graphes de F d'un nombre d'arêtes déterminé par
          value.
....
       -filter F all (= vertex t)
....
          Affiche tous les graphes de F ce qui permet en particulier
          de savoir combien il y en a en examinant la première ligne
          affichée.
....
       -filter F1 minus F2
....
          Affiche F1\F2, c'est-à-dire tous les graphes de F1 qui ne
          sont pas isomorphes à F2 (si F2 est un graphe) ou à l'un des
          graphes de F2 (dans le cas d'une famille).
....
       -filter F1 minus-id F2
....
          Comme "minus" mais concerne les identifiants: supprime de F1
          les graphes dont l'identifiant existe dans F2 (qui peut être
          un graphe ou une famille de graphes). La complexité est
          environ (|F1|+|F2|)*log|F2|, alors que pour "minus" elle est
          en |F1|*|F2|*T où T est le temps pour décider si deux
          graphes pris dans F1 et F2 sont isomorphes.
....
       -filter F minor[-inv] H
....
          Filtre les graphes de F contenant H comme mineur. La
          variante minor-inv filtre les graphes de F qui sont mineurs
          de H. Si H=K₄, il est préférable d'utiliser -filter tw2.
....
       -filter F sub[-inv] H
....
          Filtre les graphes de F contenant H comme sous-graphe,
          chaque graphe de F devant avoir le même nombre de sommets
          que H. La variante sub-inv filtre les graphes de F qui sont
          un sur-graphe de H.
....
       -filter F isub[-inv] H
....
          Filtre les graphes de F contenant H comme sous-graphe
          induit. La variante isub-inv filtre les graphes de F qui
          sont sous-graphes induits de H.
....
       -filter F iso H
....
          Filtre les graphes de F isomorphes à H.
....
       -filter F degenerate value
....
          Filtre les graphes de F de dégénérescence déterminée par
          value.
....
       -filter F forest value
....
          Filtre les graphes de F qui sont des forêts dont le nombre
          d'arbres est déterminé par value.
....
       -filter F isforest (= forest t)
....
          Filtre les graphes de F qui sont des forêts.
....
       -filter F istree (= forest '=1')
....
          Filtre les graphes de F qui sont des arbres.
....
       -filter F cycle (= not forest t)
....
          Filtre les graphes de F contenant au moins un cycle.
....
       -filter F degmax/degmin value
....
          Filtre les graphes de F de degré maximum (ou minimum)
          déterminé par value.
....
       -filter F deg value
....
          Filtre les graphes de F où tous les sommets ont un degré
          déterminé par value. Ainsi -filter deg 4-7 filtre les
          graphes avec un degré minimum au moins 4 et un degré maximum
          au plus 7.
....
       -filter F gcolor value
....
          Filtre les graphes de F dont le nombre de couleurs obtenu
          selon l'heuristique du degré minimum est déterminé par
          value.
....
       -filter F bipartite (= gcolor <3)
....
          Filtre les graphes de F qui sont bipartis.
....
       -filter F component value
....
          Filtre les graphes de F dont le nombre de composantes
          connexes est déterminé par value.
....
       -filter F connected (= component 1)
....
          Filtre les graphes de F qui sont connexes.
....
       -filter F biconnected
....
          Filtre les graphes de F qui sont 2-connexes. Un graphe G est
          k-connexe s'il n'y a pas d'ensemble avec <k sommets qui
          déconnecte G ou laisse G avec 1 sommet. Un graphe est
          2-connexe s'il est connexe, ne possède pas de sommet
          d'articulation et a plus de 2 sommets. Les cliques de taille
          k+1 sont k-connexes.
....
       -filter F radius value
....
          Filtre les graphes de F dont le rayon est déterminé par
          value. Le rayon est la profondeur du plus petit arbre
          couvrant le graphe. Il vaut -1 si le graphe n'est pas
          connexe.
....
       -filter F girth value
....
          Filtre les graphes de F dont la maille est déterminée par
          value. La maille est la taille du plus petit cycle. Elle
          vaut -1 si le graphe n'a pas de cycle. Elle n'est définie
          que si les graphes sont orientés.
....
       -filter F diameter value
....
          Filtre les graphes de F dont le diamètre est déterminé par
          value. Le diamètre vaut -1 si le graphe n'est pas connexe.
....
       -filter F cut-vertex value
....
          Filtre les graphes de F dont le nombre de sommets
          d'articulations est déterminé par value. Un sommet est un
          point d'articulation si sa suppression augmente le nombre de
          composante connexe. Les sommets de degré 1 ne sont pas des
          point d'articulation. Le graphe est biconnexe ssi value<1 ou
          si le graphe est une clique avec au moins deux sommets. On
          peut tester si un graphe est une clique avec -filter degmin
          ou -filter deg.
....
       -filter F ps1
       -filter F ps1b
       -filter F ps1c
       -filter F ps1x n u_1 v_1 … u_n v_n
....
          Filtre les graphes G de la famille F dont le test ps1 est
          vrai, c'est-à-dire si l'évaluation de la fonction f(G,{})
          décrite ci-après est vraie.
....
	  Soit P un chemin d'un graphe G tel que G\P est connexe. La
          fonction f(G,P) est vraie ssi G est vide (en pratique
          |G|-|P|<3 suffit) ou s'il existe deux sommets x,y de G où y
          n'est pas dans P tels que pour tout chemin Q entre x et y
          dans G "compatible" avec P (c'est-à-dire P et Q
          s'intersectent en exactement un segment) on a les deux
          conditions suivantes: (1) il n'y a pas d'arête entre les
          sommets de P\Q et de G\(Q∪P); et (2) pour toute composante
          connexe C de G\(Q∪P), f(C∪Q,Q) est vraie. Le test est
          optimisé dans un certain nombre de cas, en particulier: les
          arbres (toujours vrai), les cliques (vrai ssi n<5).
....
	  La variante ps1b calcule et affiche de plus un graphe des
          conflits (affichage modifiable par -width), chaque noeud de
          ce graphe correspondant à un argument (C∪Q,Q) évalué à faux
          par f. La valeur (ou code) d'un noeud est 0 (=lourd ou
          faux), 1 (=léger ou vrai) ou - (indéterminée). Suivant
          certaines règles, les valeurs 0 ou 1 sont propagées selon le
          type des arêtes du graphes des conflits.  Résoudre le graphe
          des conflits revient à trouver une affectation des valeurs 0
          ou 1 aux noeuds qui respecte (sans contradiction) toutes les
          règles.
....
	  La fonction f(G,{}) est évaluée à vraie si le graphe des
          conflits n'a pas de solution, c'est-à-dire si une
          contradiction a été découverte ou si pour une paire de
          sommets (x,y) tous ses noeuds sont à 1.
....
	  On affiche le code d'un noeud (0,1,-) ainsi que les sommets
          de sa composante (par ex: [237]).  Les noeuds du graphe des
          conflits sont reliées par des arêtes typées. Les voisins v
          d'un noeud u sont listés avec le type de l'arête, si l'un
          des 4 cas suivants se produit (il n'y a pas d'arête entre u
          et v dans les autres cas):
....
	     v<  (la composante de v est incluse dans celle de u)
	     v>  (la composante de v contient celle de u)
	     v=  (les composantes de u et v sont les mêmes) 
	     v|  (les composantes de u et v sont disjointes) 
....
	  Parmi les règles on trouve par exemple: si deux noeuds du
          graphe des conflits u=(C∪Q,Q) et v=(C'∪Q',Q') sont
          disjoints, c'est-à-dire C n'intersecte pas C', alors seule
          une des expressions f(C∪Q,Q) ou f(C'∪Q',Q') peut être
          fausse, pas les deux. Dit autrement, les composantes de u et
          v ne peuvent pas être "lourdes" (=0) toutes les deux en même
          temps. Et donc, si le code de u est 0, celui de v est
          1. Notons que le code de u et v égale à 1 est compatible
          avec cette règle.
....
	  La variante ps1c est similaire à ps1b sauf que récursivement
          seul le test ps1 est appliqué, et pas ps1b. Le test ps1c est
          plus long que ps1 mais plus rapide que ps1b. La variante
          ps1x est similaire à ps1b sauf que les valeurs v_i sont
          écrites dans le noeuds u_i du graphe des conflits principal
          (pas ceux générés lors des appels récursifs). Plus
          précisément, v_1 (0 ou 1) est écrit dans le noeud u_1, puis
          sa valeur est propagée. Ensuite v_2 est écrit puis propagée,
          etc.
....
          Dans tous les cas, si G n'est pas connexe, le résultat n'est
          pas déterminé.
....
       -filter F tw value
....
          Filtre les graphes de F selon leur treewidth. L'algorithme
          pour le calcul de la treewidth est assez lent. Pour les
          petites valeurs de tw, des alternatives sont possibles (voir
          -check tw et -filter tw2). Pour savoir si un graphe G est de
          treewidth 3 il suffit de savoir si G contient l'un des 4
          mineurs suivants:
....
          !!! echo "[0]"  > F; ./gengraph clique 5 >> F
	      echo "[1]" >> F; ./gengraph wagner >> F
	      echo "[2]" >> F; ./gengraph prism 5 >> F
	      echo "[3]" >> F; ./gengraph hajos >> F ; echo "0-1-2-0" >> F
	      cat G |./gengraph -filter F minor-inv - -format no
....
       -filter F tw2
....
          Affiche les graphes de F de treewidth ≤ 2. L'algorithme est
          en O(n²). Ce test peut être utilisé pour tester (plus
          rapidement qu'avec -filter minor) les graphes sans mineur
          K₄.
....
       -filter F hyper value
....
          Filtre les graphes de F selon leur hyperbolicité. Il s'agit
          de la valeur (entière) maximum, sur tous les quadruplets de
          sommets {u,v,x,y}, de la différence des deux plus grandes
          sommes parmi les sommes de distance : uv+xy, ux+vy et
          uy+vx. La complexité est en O(n⁴).

....-format type
....
       Spécifie le format de sortie. Il est préférable d'utiliser
       cette option en dernier. Les valeurs possibles pour "type"
       sont:
....
       • standard: format standard (liste d'arêtes), c'est le plus compact.
       • list: liste d'adjacence.
       • matrix: matrice d'adjacence.
       • smatrix: matrice supérieure, diagonale comprise.
       • vertex i: liste des voisins du sommet i.
       • dot: format de GraphViz qui est très proche du format standard.
       • dot<type>: dessine le graphe avec GraphViz et converti au format <type>.
       • html: dessin dynamique au format html et vis.js (cf. http://visjs.org).
       • xy: positions X,Y qui ont été utilisées pour le graphe géométrique.
       • no: n'affiche rien, à utiliser en combinaison avec -header ou -check.
....
       Les formats matrix/smatrix/list/vertex nécessitent de stocker
       le graphe en mémoire, donc nécessite un espace en O(n+m), alors
       que le graphe est généré à la volée pour les formats standard,
       dot ou html. Les valeurs les plus utilisées de <type> pour le
       format dot<type> sont: pdf, fig, svg, ps, jpg, gif, png (voir
       man dot ou faire dot -T.).
....
       L'option -format dot<type> est équivalent à "-format dot | dot
       -T<type>". Par conséquent, elle doit donc être utilisée en
       dernier. Le filtre dot utilisé pour dessiner le graphe peut
       être spécifié par l'option -dot filter.  L'affichage des noms
       de sommets est contrôlé par l'option -label.
....
       Remarque: les positions affichées dans le format dot
       ([pos="…"]) diffèrent d'un facteur proportionnel à √n par
       rapport aux positions originales du graphe (qui peuvent être
       affichées par -format xy ou -label -3). Ce facteur permet de
       garder une taille raisonable pour les sommets car sous dot les
       sommets ont une taille fixe minimale.

....-vcolor option [parameters]
....
       Ces options permettent de modifier la couleur des sommets. Ces
       options n'ont d'effets qu'avec le format dot (et ses variantes
       y compris -visu). Par défaut les sommets sont de couleur
       noire. Notez que les attributs par défaut des sommets
       (couleurs, formes, etc.)  peuvent être modifiés directement par
       dot (voir l'option -N de dot). Cependant l'option -vcolor
       permet d'individualiser la couleur d'un sommet, en fonction de
       son degré par exemple. Il peut avoir plusieurs options -vcolor
       pour une même commande.
....
       -vcolor deg[r]
....
          La couleur dépend du degré du sommet (deg) ou du rang du
          degré du sommet (degr). Ainsi, les sommets de plus petit
          degré obtiennent la première couleur de la palette, les
          sommets de plus grand degré la dernière couleur de la
          palette, et les autres sommets une couleur intermédiaire de
          la palette. Donc une seule couleur est utilisée si le graphe
          est régulier.
....
       -vcolor degm
....
          Effectue une coloration propre (deux sommets voisins ont des
          couleurs différentes) suivant l'heuristique du degré
          minimum: récursivement, le sommet de degré minimum obtient
          la plus petite couleur qui n'est pas utilisée par ses
          voisins. Cela donne des colorations avec assez peu de
          couleurs pour les graphes de faible arboricité (planaire,
          tw, pw, kout, expander, …) ou de faible degré. Avec cette
          technique, les graphes bipartis (tree, crown, …) sont
          coloriés avec deux couleurs. Cette option nécessite un
          espace et un temps en O(n+m).
....
       -vcolor randg
....
          Effectue une coloration propre en utilisant un algorithme
          glouton sur un ordre aléatoire des sommets: récursivement,
          le sommet d'indice i obtient la plus petite couleur qui
          n'est pas utilisée par ses voisins d'indice j<i. Cette
          option nécessite un espace et un temps en O(n+m).
....
       -vcolor kcolor k
....
          Effectue une k-coloration propre du graphe, si c'est
          possible. Si cela n'est pas possible, la première couleur
          est appliquée à tous les sommets. L'algorithme (exponentiel)
          est le même que celui utilisé pour -check kcolor.
....
       -vcolor pal grad
....
          Permet de fixer la palette de couleurs utilisée par les
          sommets. Le paramètre "grad" est un mot sur l'alphabet [a-z]
          (sans les guillemets). Les caractères en dehors de cet
          alphabet sont ignorés. Chaque lettre correspond à une
          couleur de base:
....
       !!! a=aquamarine     h=hotpink      o=olive         v=violet
	   b=blue           i=indigo       p=purple        w=white
	   c=cyan           j=orange       q=pink          x=gray
	   d=darkorange     k=khaki        r=red           y=yellow
	   e=chocolate      l=lavender     s=salmon        z=black
	   f=forestgreen    m=magenta      t=teal
	   g=green (lime)   n=navy         u=yellowgreen
....
          La palette est calculée selon une interpolation linéaire
          entre les points définis par le mot "grad". Par exemple, si
          "grad" vaut rb, la palette sera composée d'un dégradé allant
          du rouge (r) au bleu (b). Si "grad" vaut rgbr, le dégradé
          ira du rouge au vert puis au bleu et enfin au rouge. Pour
          avoir une couleur (de base) unique, disons w, sur tous les
          sommets, poser "grad" égale à w. Par exemple, pour avoir
          tous les sommets blancs, on peut faire:
....
          !!! gengraph gabriel 30 -vcolor deg -vcolor pal w -visu
....
          La palette par défaut correspond au mot "grad" suivant:
          redjykugfocatbhsqvmpinzxlw. On peut visualiser la palette
          avec l'option "-vcolor list".
....
       -vcolor list
....
          Produit l'affichage de la palette des couleurs utilisées
          pour un graphe plutôt que le graphe lui-même. Cela permet en
          particulier de savoir combien de couleur ont été utilisées.
          La palette est générée en affichant au format dot un graphe
          particulier où les sommets (représentés par un rectangle)
          sont les couleurs utilisées. Utilisez -visu pour visualiser
          la palette sous forme pdf. Le nom des sommets correspond à
          la lettre de la couleur de base comme spécifié par -vcolor
          pal.
....
          Ex1: gengraph gabriel 50 -vcolor degm -vcolor list
	  (génère la palette utilisée pour ce graphe de Gabriel)
....
          Ex2: gengraph prime 53 -vcolor list
	  (un moyen simple de générer la palette par défaut)
....
          Ex3: gengraph clique 100 -vcolor degm -vcolor pal rb -vcolor list
          (génère un dégradé de 100 couleurs allant du rouge au bleu)

....-vsize
....
       La taille des sommets est proportionnelle à son degré, alors que
       par défaut elle est fixe. Cette option n'a d'effet qu'avec le
       format dot (et ses variantes). Elle est combinable avec
       -vcolor.

....-visu
....
       Crée un fichier "g.pdf" permettant de visualiser le graphe. Il
       s'agit d'un raccourci de l'option "-format dotpdf" qui rajoute
       également la redirection "> g.pdf" en fin de la ligne de
       commande.

....-maincc
....
       Affiche la composante connexe principale du graphe, les sommets
       étant éventuellement renumérotés si le graphe n'est pas
       connexe. C'est un raccourci pour "-check maincc | ./gengraph
       load - -fast". (Voir aussi -check maincc.) Cet affichage est
       réalisé en temps linéaire grâce à l'option -fast. Les options
       placées avant -maincc affectent le graphe initial alors que
       celles placées après affectent la composante principale.  Les
       options ayant un effet pour les formats hors standard (comme
       -vsize ou -visu) ne devraient être placées qu'après cette
       option.

....-dot option [parameters]
....
       Cette option permet de controler la sortie au format dot. Elle
       permet par exemple de modifier le filtre, la longueur des
       arêtes ou l'échelle du dessin.
....
       -dot scale s
....
          Spécifie le facteur d'échelle pour le format dot. Cela
          affecte les coordonnées des sommets et des arêtes, pas des
          étiquettes (sommets ou arêtes). Cela permet d'écarter les
          sommets les uns des autres si nécessaires. Le format pour s
          prend deux formes: x ou x,y pour un facteur d'échelle
          identique ou pas en X et Y.
....
          Ex: gengraph gabriel 10 -label -3 -dot scale 3,2 -visu
....
       -dot len p
....
          Spécifie la longueur des arêtes pour le format dot et le
          filtre "neato". La valeur par défaut est 1, et une valeur
          plus grande (comme 2.5 ou 3) allonge les arêtes et permet
          dans certain cas de mieux visualiser le graphe. C'est
          parfois nécessaire pour éviter l'intersection des sommets
          lorsqu'on utilise -label 1. On peut obtenir le même genre
          d'effet avec -dot scale.
....
       -dot filter f
....
         Spécifie le filtre de GraphViz, c'est-à-dire l'algorithme de
         dessin utilisé par dot. Par défaut, le filtre est
         "neato". Les filtres principaux sont: dot, neato, twopi,
         circo, fdp, sfdp, … Faire "dot -K ." pour afficher les
         filtres disponibles.

....-pos b
....
       Active (b=1) ou désactive (b=0) la génération des positions des
       sommets pour le format dot. Cela sert à indiquer à l'algorithme
       de dessin dot de respecter (b=1) ou pas (b=0) les coordonnées
       des sommets. L'option par défaut est -pos 0, mais cette option
       est activée pour tous les graphes géométriques (udg, gabriel,
       thetagone, …).

....-label b
....
       Active (b≠0) ou désactive (b=0) l'affichage du nom des sommets
       pour les formats dot et standard. Les valeurs possibles de b
       sont b∈{-3,-2,-1,0,1,2,3}. Si b=1, il s'agit du nom original du
       sommet, par exemple un mot binaire pour l'hypercube. Cette
       fonctionnalité n'est pas implémentée pour tous les graphes, le
       nom par défaut étant les entiers de [0,n[ où n est le nombre de
       sommets du graphe généré. L'option -label 1 -visu permet alors
       d'afficher sur le dessin du graphe le nom des sommets. Ils ne
       le sont pas par défaut (b=0).  L'option -label 2 -visu force
       l'affichage des noms sous forme d'entiers de [0,n[ (et non pas
       du nom original). L'option -label 3 permet, dans le cas de
       graphe géométrique (ou si -pos 1), d'afficher les coordonnées
       des points. Si b<0, alors l'effet est similaire à -label |b|
       sauf que le nom du sommet est affiché à coté du sommet et non
       pas sur au center du sommet. L'option -label 1 annule l'option
       -permute, mais -label 2 ne le fait pas. Comme cette option
       influence l'option -format dot<type>, l'option -label devrait
       être placée avant -format.
....
       Ex: gengraph petersen -label 1 -width 1
           gengraph petersen -label 1 -format dot | grep label
           gengraph petersen -label 1 -dot len 2 -visu
           gengraph gabriel 30 -pos 0 -label 1 -visu
	   gengraph gabriel 30 -label -3 -dot scale 4 -xy round 2 -visu

....-norm ℓ [parameters]
....
       Fixe la norme d'un vecteur (x,y) du plan (ou la fonction de
       distance entre deux points du plan) pour l'adjacence de
       certains graphes géométriques (dont udg, gabriel, rng, nng, …).
       Par défaut c'est la norme Euclidienne qui est utilisé. Les
       valeurs possibles pour ℓ sont:
....
         • L1      →  |x|+|y|, distance de Manhattan
         • L2      →  √(x²+y²), norme Euclidienne
         • Lmax    →  max{|x|,|y|}
         • Lmin    →  min{|x|,|y|}
         • poly p  →  distance polygonale de paramètre p
         • hyper   →  distance hyperbolique
....
       Il s'agit de pseudo-norme (ou pseudo-distance) puisque par
       exemple la norme "Lmin" ne vérifie pas l'inégalité
       triangulaire. La norme polygonale est le rayon du cercle
       inscrit dans le polygone régulier convexe à p cotés contenant
       (x,y), le polygone étant centré en (0,0) et orienté de façon à
       avoir son coté le plus à droit vertical. Ainsi "poly 4"
       correspond à la norme "Lmax". Une valeur de p<3 est interprétée
       comme p=+∞ ce qui correspond à la norme Euclidienne. Attention
       ! la norme "poly p" n'est pas toujours symétrique, lorsque p
       est impair par exemple. La norme (ou distance) hyperbolique
       n'est définie que pour des points du disque ouvert unité centré
       en (0,0).

....-xy option [parameters]
....
       Cette option contrôle la façon dont sont générées les
       coordonnées des sommets d'un graphe géométrique. Par défaut les
       positions sont tirées aléatoirement uniformément dans le carré
       [0,1[ × [0,1[, mais cela peut être changé par l'option -xy.
       Notez bien que, même si c'est improbable, deux sommets peuvent
       avoir les mêmes positions (voir l'option -xy unique). Il est
       possible de visualiser les points issus des options -xy (voir
       le graphe "point n").
....
       -xy load file
....
          Charge les positions à partir du fichier "file" ou de
          l'entrée standard si file=-. Cela permet de tester les
          adjacences d'un graphe géométrique à partir de positions
          pré-déterminées. Le format est celui de -format xy.
....
          Ex: gengraph gabriel 10 -xy load file.pos
....
	  Le nombre de sommets du graphe est déterminé par le fichier
          et non par les paramètres du graphe. Cette option n'a
          d'effet que pour les graphes géométriques. La structure du
          fichier texte doit être:
....
          !!!    n
		 x_1 y_1
		 x_2 y_2
		 …
		 x_n y_n
....
	  où n est le nombre de positions. Les positions x_i y_i ne
	  sont pas forcément dans l'intervalle [0,1[. Notez qu'avec
	  l'option -format xy, il est possible d'effectuer la
	  transformation d'un fichier de positions. L'exemple suivant
	  normalise les coordonnées du fichier g.pos dans le carré
	  unité:
....
          Ex: gengraph -xy load g.pos -xy box 1 1 -format xy
....
       -xy box a b
....
          Effectue un redimensionement des positions de sorte quelles
          se situent dans le rectangle [0,a[ × [0,b[. En prenant
          a=b=1, les coordonnées seront renormalisées dans le carré
          [0,1[ × [0,1[. Cette opération est effectuée juste avant la
          génération des arêtes, mais après avoir effectué l'opération
          -xy noise (voir ci-après) et/ou -xy load.
....
       -xy grid p
....
          Ajoute une grille p × p au graphe généré, ce qui est utile
          lorsque les coordonnées des points sont entiers.
          Techniquement, on ajoute au format de sortie dot un
          sous-graphe représentant la grille où les sommets et les
          arêtes sont de couleur grise. Si p<0, alors le paramètre est
          initialisé à 1+⎣ √n⎦ ou bien à n si l'option "-xy
          permutation" est présente, n étant le nombre de sommets du
          graphe. Pour être visible, le nombre de lignes (et de
          colonnes) de la grille générée doit être au moins 2.
....
       -xy zero
....
          Ajoute l'origine (0,0) au dessin qui est représenté par un
          cercle rouge.
....
       -xy vsize f
....
          Facteur de grossissement des sommets pour le format dot. Par
          défaut f=1.
....
       -xy noise r p
....
          Effectue une perturbation aléatoire sur les positions des
	  sommets. Le déplacement de chaque sommet est effectué dans
	  sa boule de rayon r (pour p>0) selon une loi en puissance de
	  paramètre p. Prendre p=0.5 pour une perturbation uniforme
	  dans cette boule, p>0.5 pour une concentration des valeurs
	  vers le centre et p<0.5 pour un écartement du centre. Les
	  valeurs <0 de p donne des écartements au delà du rayon r.
....
          Plus précisément, une direction (angle de 0 à 2π) est
	  choisie aléatoirement uniformément, puis, selon cette
	  direction, un décalage aléatoire est effectué selon une loi
	  en puissance: si x est uniforme dans [0,1[, le décalage sera
	  d(x)=r*x^p.  Après cette opération, il est possible que les
	  points ne soient plus dans le rectangle d'origine, ce qui
	  peut bien sûr être corrigé par -xy box.
....
       -xy seed k p
....
          Génère les points à partir (ou autour) de k>0 graines. Les
	  graines sont choisies uniformément dans le carré [0,1[ ×
	  [0,1[ puis centrées par rapport à leur barycentre. Chaque
	  point est alors tiré aléatoirement autour d'une des graines
	  et à une distance variant selon une loi en puissance (voir
	  -xy noise) de paramètre p et de rayon r ≃ √(ln(k+1)/k). Ce
	  rayon correspond au seuil de connectivité pour un Unit Disk
	  Graph à k sommets dans le carré [0,1[ × [0,1[ (voir udg n
	  r). On peut obtenir une distribution uniforme dans un disque
	  avec -xy seed 1 0.5 (voir aussi -xy disk) sauf que le centre
	  est en (1/2,1/2) au lieu de (0,0) comme avec -xy disk.
....
          Ex: gengraph point 1000 -xy seed 1 1
	      gengraph point 1000 -xy seed 1 0.5
....
       -xy permutation
....
          Génère les points correspondant à une permutation π
          aléatoire uniforme. Le point i aura pour position (i,π(i)).
....
       -xy mesh x y
....
          Génère tous les points de coordonnées entières correspondant
          aux sommets d'une grille de x colonnes et de y lignes.
....
       -xy cycle
....
          Génère les points régulièrement espacés le long d'un cercle
	  de centre (0,0) et de rayon 1. Les points sont ordonnées
	  selon l'angle de leurs coordonnées polaires.
....
          Ex: gengraph cycle 10 -xy cycle -visu
....
       -xy unif
....
          Génère les points aléatoirement uniformément dans le carré
          [0,1[ × [0,1[. C'est la distribution par défaut.
....
       -xy circle
....
          Génère les points aléatoirement uniforme le long d'un cercle
	  de centre (0,0) et de rayon 1. Les points sont ordonnées
	  selon l'angle de leurs coordonnées polaires.
....
       -xy disk
....
          Génère les points aléatoirement uniforme dans le disque
          unité de centre (0,0) triés selon l'angle de leurs
          coordonnées polaires. Cette distribution permet de générer,
          par exemple, un polygone "star-shaped". La distribution est
          similaire à l'option "-xy seed 1 0.5" sauf que les points
          sont ordonnées.
....
          Ex: gengraph cycle 25 -xy disk -visu
....
          Remarque: les points sont générés avant l'application des
          options comme -xy round, -xy noise, ou -xy unique qui
          modifient les coordonnées et qui peuvent donc produire des
          croisements avec le graphe cycle par exemple.
....
       -xy hyper p
....
          Génère les points aléatoires selon une loi exponentielle de
          paramètre p dans le disque unité de centre (0,0). Les poitns
          sont triés selon l'angle de leurs coordonnées polaires.
....
       -xy convex
       -xy convex2
....
          Génère les points aléatoirement en position convexe à
          l'intérieur d'un cercle de rayon 1 et de centre (0,0), ce
          qui peut être modifié par -xy ratio. Ils sont numérotés
          consécutivement selon le parcours de l'enveloppe convexe. On
          les génère comme suit, n étant le nombre points à
          générer. Inductivement, une fois que n-1 points en position
          convexe ont été générés, on choisit un angle 𝛼 aléatoire du
          cercle de rayon 1 et de centre (0,0) supposé à l'intérieur
          du convexe. On détermine ensuite la partie S du segment
          d'angle 𝛼 où chacun des points de S=[a,b[ forment avec les
          n-1 points précédant un convexe. Enfin, on choisit
          aléatoirement un point de S selon la probabilité √|b-a| pour
          obtenir n points en position convexe. Les angles des trois
          premiers points sont choisis parmi trois secteurs non
          adjacents d'angle π/3 si bien que l'origine est toujours à
          l'intérieur de l'ensemble convexe.
....
          La variante -xy convex2 génère également des points
          aléatoires en position convexe, selon la méthode suivante.
          On génère n points aléatoires u_i du carré [0,1[ × [0,1[
          puis on calcule les vecteurs différences v_i ≡ u_{i+1}-u_i
          (mod n). Les vecteurs (dont la somme est nulle) sont ensuite
          triés par angle croissant, puis les points en position
          convexe sont obtenus de proche en proche en ajoutant chacun
          des vecteurs v_i. Cette méthode tend à générer des points
          proches d'un cercle, chaque angle et chaque longueur entre
          deux points consécutifs suivant une loi normale.
....
          Ex: gengraph cycle 25 -xy convex -visu
	      gengraph dtheta 100 6 -xy convex -visu
....
          L'ordre des sommets peut être modifié par certaines options
          (voir la remarque de l'option -xy disk).
....
       -xy polygon p
....
          Génère des points aléatoires uniformément dans un polygone
	  convexe régulier à p≥3 cotés inscrit dans le cercle de
	  centre (0,0) et de rayon 1 de sorte qu'un des cotés du
	  polygone soit vertical. Les sommets ne sont pas
	  spécifiquement ordonnés. Pour une distribution uniforme dans
	  un disque, soit lorsque p=+∞, utiliser -xy disk. L'option
	  pour p=4 est similaire à -xy unif, sauf que pour p=4 la
	  distribution est dans le carré [-c,+c[ × [-c,+c[ où c =
	  cos(π/4) = ½√2 ≃ 0.707… au lieu du carré [0,1[ × [0,1[.
....
       -xy ratio r
....
          Modifie les distributions de points faisant intervenir une
          forme de largeur 1 et de hauteur r, comme: -xy unif, -xy
          circle, -xy cycle, -xy convex, -xy disk, -xy seed. La valeur
          par défaut est r=1. Le réel r>0 est donc le ratio de la
          hauteur par la largeur de la forme. Par exemple, pour la
          distribution par défaut (-xy unif), les points seront
          aléatoires uniformes dans le rectangle [0,1[ × [0,r[. Si la
          forme est un cercle (-xy circle ou -xy disk), alors la forme
          devient une ellipse dont le rayon horizontal est 1 et celui
          vertical r. Dans le cas de -xy seed, les graines sont alors
          générés dans le rectangle [0,1[ × [0,r[.
....
       -xy surface s
....
          Définit la signature s de la surface sur laquelle va être
          construit le graphe géométrique. La surface peut-être
          orientable ou non, avec ou sans bord. Elle est représentée
          par un polygone convexe régulier inscrit dans un cercle de
          rayon 1 et dont les 2|s| cotés sont appariés. Cette option
          se charge également de générer des points aléatoirement
          uniformes sur la surface. La signature est un mot s sur
          l'alphabet {h,c,b} de longueur |s|=2g, où g est le genre de
          la surface, indiquant comment sont appariés les 4g cotés du
          polygone. Chaque coté est apparié avec le coté +2 (le
          suivant du suivant) ou le coté -2 selon l'une des trois
          coutures suivantes:
....
          • h = handle   = couture orientée ou anse
          • c = crosscap = couture non-orientié
          • b = border   = aucune couture
....
          La caractéristique d'Euler de la surface (ou sa courbure)
          vaut 2-|s| = 2-2g. Certaines signatures ont des synonymes.
          Par exemple, "-xy surface torus" est synonyme de "-xy
          surface hh" (voir ci-dessous leurs listes).
....
          Ex: -xy surface bb (ou plane ou square)  → plan réel
              -xy surface hb (ou cylinder) ....... → cylindre
              -xy surface cb (ou mobius) ......... → ruban de Möbius
              -xy surface hh (ou torus) .......... → tore
	      -xy surface ch (ou klein) .......... → bouteille de Klein
              -xy surface cc (ou projective) ..... → plan projectif
	      -xy surface hhhh ................... → double tore
....
	  Cette option active également -xy polygon 4g et -xy ratio 1
	  pour générer des points aléatoires uniformément sur la
	  surface.
....
       -xy round p
....
          Arrondi les coordonnées à 10^-p près. Il faut que p soit un
          entier < DBL_DIG, soit p<15 en général. Donc p=0 arrondi à
          l'entier le plus proche. Cet opérateur est appliqué après
          -xy box. Il sert aussi à préciser le nombre de décimales à
          afficher pour l'option -format xy (par défaut p=6). Par
          exemple, la combinaison -xy box 100 100 -xy round -1 permet
          d'avoir des coordonnées multiples de 10.
....
       -xy unique
....
          Supprime les sommets en double, correspondant aux mêmes
          positions. Cela peut être utile lorsqu'on utilise -xy round
          par exemple. Cette opération est appliquée après toutes les
          autres, notamment après -xy box et -xy round. Ceci est
          réalisé à l'aide d'un tri des points, l'ordre n'est donc pas
          préservé).


   GRAPHES

       Deux types de graphes sont possibles : les graphes de base et
       les graphes composés. Ces derniers sont obtenus en paramétrant
       un graphe de base. Une catégorie importante de graphes sont les
       graphes géométriques (qui peuvent être composés ou de bases).
       L'adjacence est déterminée par les coordonnées associées aux
       sommets. De nombreuses options s'y réfèrent.  Ils activent tous
       par défaut l'option -pos. Les graphes orientés activent quant à
       eux tous l'option -directed.
       


   GRAPHES DE BASE :

....grid k n_1 … n_k
....
       Grille à k dimensions de taille n_1 × … × n_k. Si la taille
       n_i est négative, alors cette dimension est cyclique.  Par
       exemple, "grid 1 -10" donnera un cycle à 10 sommets.

....ring n k c_k … c_k
....
       Anneaux de cordes à n sommets chacun ayant k cordes de longueur
       c_1,…,c_k.

....cage n k c_1 … c_k
....
       Graphe cubique pouvant servir à la construction de graphes
       n-cage, c'est-à-dire aux plus petits graphes cubique à n
       sommets de maille donnée. Ils sont toujours Hamiltoniens. Ils
       peuvent être vus comme des anneaux de cordes irréguliers. Ils
       sont construits à partir d'un cycle de longueur n découpé en
       n/k intervalles de k sommets. Le i-ème sommet de chaque
       intervalle, disons le sommet numéro j du cycle, est adjacent au
       sommet numéro j+c_i du cycle (modulo n). Les valeurs c_i peuvent
       être positives ou négatives. Cette fonction permet aussi de
       construire des graphes avec des sommets de degré 4 comme "cage
       8 2 0 2" (voir aussi le graphe de Chvátal) ou avec des sommets
       de degré 2 comme "cage 4 2 2 0".

....arboricity n k
....
       Graphe d'arboricité k à n sommets aléatoire. Ce graphe est
       composé de l'union de k>0 arbres aléatoires. Il est donc
       toujours connexe. Chacun des arbres est un arbre plan enraciné
       aléatoire uniforme dont les sommets sont permutés
       aléatoirement, sauf le premier arbre dont les sommets sont
       numérotés selon un parcours en profondeur. Ces graphes
       possèdent au plus k(n-1) arêtes, et pour k=1 il s'agit d'un
       arbre.

....rarytree n b z
....
       Arbre b-aire plan aléatoire uniforme à n noeuds internes. Il
       faut b≥2. Il possède bn+1+z sommets, z étant un paramètre
       valant 0 ou 1. La racine est de degré b+z, les autres sommets
       sont de degré b+1 (soit b fils) ou 1 (=feuille). Les sommets
       sont numérotés selon un parcours en profondeur modifié: tous
       les fils du sommet courant sont numérotés avant l'étape de
       récursivité. Si n=1, alors le graphe est une étoile à b+z
       feuilles. Le dessin avec dot (-visu) ne respecte pas le
       plongement de l'arbre.

....ringarytree h k r p
....
       Arbre de hauteur h où chaque noeud interne à exactement k fils,
       le degré de la racine étant de degré r. Lorsque p>0, un chemin
       (si p=1) ou un cycle (si p=2) est ajouté entre les sommets de
       même niveau. Notez que "ringarytree h 1 r 0" génère une étoile
       de degré r où chaque branche est de longueur h. Le nom des
       sommets correspond au chemin depuis la racine.

....kpage n k
....
       Graphe k-pages connexe aléatoire. Un graphe k-page peut être
       représenter en plaçant les sommets le long d'un cercle, en
       dessinant les arêtes comme des segments de droites, et en
       coloriant les arêtes en k>0 couleurs de façon à ce que les
       arêtes de chaque couleur induisent le dessin d'un graphe
       planaire-extérieur. La numérotation des sommets est faite le
       long du cercle. Les graphes 1-page sont les graphes
       planaires-extérieurs, les 2-pages sont les sous-graphes de
       graphes planaires Hamiltoniens. Les graphes planaires de degré
       au plus 4 sont 2-pages, les 3-arbres planaires (ou graphes
       Apolloniens) sont 3-pages, et les cliques avec 2k-1 ou 2k
       sommets des k-pages.
....
       Ces graphes sont construits par le processus aléatoire suivant.
       On génère k graphes planaires-extérieurs aléatoires uniformes
       connexes à n sommets (plan et enraciné) grâce à une bijection
       avec les arbres plans enracinés dont tous les sommets, sauf
       ceux de la dernière branche, sont bicoloriés. On fait ensuite
       l'union de ces k graphes en choisissant aléatoirement la racine
       des arbres, sauf celui du premier planaire-extérieur, ce qui
       correspond à une permutation circulaire des sommets sur la face
       extérieure.

....cactus n
....
       Graphe cactus aléatoire à n sommets. Il s'agit d'arbres de
       cycles, c'est-à-dire de graphes connexes où chaque arête
       appartient à au plus un cycle. Ce sont aussi les graphes
       planaires-extérieurs connexes sans cordes. Ils sont générés à
       partir d'un "outerplanar n" dans lequel les arêtes internes (ou
       cordes) des composantes biconnexes ont été supprimés.

....ktree n k
....
       k-arbre aléatoire à n sommets. Il faut n>k≥0. C'est un graphe
       chordal appelé aussi graphe triangulé (triangulated). Il est
       généré à partir d'un arbre enraciné aléatoire uniforme à n-k
       noeuds de manière similaire à "tree n-k". Cela constitue les
       "sacs" que l'on remplit avec les n sommets comme suit: on met
       k+1 sommets dans le sac racine connecté en clique, puis, selon
       un parcours en profondeur de l'arbre, on met un sommet
       différent pour chacun des autres sacs. Ce sommet est alors
       connectés à exactement k sommets choisis aléatoirement dans le
       sac parent et sont ajoutés à son sac. Lorsque k=1, c'est un
       arbre, et lorsque k=0, c'est un stable.

....kpath n k
....
       k-chemin aléatoire à n sommets. La construction est similaire à
       celle utilisée pour ktree, sauf que l'arbre est un chemin. Ces
       graphes sont des graphes d'intervalles particuliers (voir
       "interval n").

....kstar n k
....
       k-star aléatoire à n sommets. La construction est similaire à
       celle utilisée pour ktree, sauf que l'arbre est une étoile. Ces
       graphes, qui sont des "split graphs", sont composés d'une
       clique à k+1 sommets et de n-k-1 sommets indépendants connectés
       à k sommets aléatoire de la clique. Il est possible d'obtenir
       le graphe "split n k" si à chaque fois les k sommets de la
       clique tirés aléatoires sont toujours les mêmes.

....rig n k p
....
       Graphe d'intersections aléatoire (Uniform Random Intersection
       Graph). Il possède n sommets, chaque sommet u étant représenté
       par un sous-ensemble S(u) aléatoire de {1,…,k} tel que chaque
       élément appartient à S(u) avec probabilité p. Il y a une arête
       entre u et v ssi S(u) et S(v) s'intersectent. La probabilité
       d'avoir une arête entre u et v est donc Pₑ=1-(1-p²)^k, mais les
       arêtes ne sont pas indépendantes (Pr(uv|uw)>Pr(uv)). En
       général, pour ne pas avoir Pₑ qui tend vers 1, on choisit les
       paramètres de façon à ce que kp²<cste.  Lorsque k≥n³, ce modèle
       est équivalent au modèle des graphes aléatoires d'Erdös-Reny
       (voir random n p). Si p<0, alors p est fixée au seuil théorique
       de connectivité, à savoir p=√(ln(n)/(nk)) si k>n et p=ln(n)/k
       sinon.

....apollonian n
....
       Graphe Apollonien aléatoire uniforme à n≥4 sommets. Les graphes
       Apolloniens sont les 3-arbres planaires ou encore les graphes
       planaires maximaux chordaux. Ils sont obtenus en subdivisant
       récursivement un triangle en trois autres. Ils sont
       3-dégénérés, de treewidth 3, et de nombre chromatique 4. La
       distance moyenne est ϴ(logn). Ils sont en bijection avec les
       arbres ternaires à n-3 noeuds internes. Pour n=5, il s'agit
       d'un K₅ moins une arête qu'on peut obtenir aussi avec "split 5
       3".

....polygon n
....
       Triangulation aléatoire uniforme d'un polygone convexe à n≥3
       cotés. Ce sont aussi des graphes planaires-extérieurs maximaux
       aléatoires. Ils sont Hamiltoniens, 2-dégénérés, de treewidth 2,
       et de nombre chromatique 3. Ils sont en bijection avec les
       arbres binaires à n-2 noeuds internes. La numérotation des
       sommets n'est pas cyclique le long du polygone. Ce graphe n'est
       pas un graphe géométrique contrairement à sa variante
       convex-polygon.
....
       Ex: gengraph polygon 20 -dot filter circo -visu

....planar n f d
....
       Graphe planaire aléatoire composé de n faces internes de
       longueur f≥3, les sommets internes étant de degré au moins d et
       ceux de la face externe au moins 2. Ils possèdent entre n+f-1
       et n(f-2)+2 sommets, sont 2-connexes, 2-dégénérés, de maille
       f. Si d>4 alors ils sont d'hyperbolicité O(f). Ils sont
       construits en ajoutant itérativement les faces par le processus
       aléatoire suivant. Au départ, il s'agit d'un cycle de longueur
       f. Pour chaque nouvelle face on ajoute un sommet u que l'on
       connecte à un sommet quelconque du cycle C formant le bord de
       la face extérieure du graphe courant. Puis on ajoute un chemin
       allant de u à un sommet v de C de façon à respecter: 1) la
       contrainte des degrés des sommets qui vont devenir internes; et
       2) la contrainte sur la longueur de la nouvelle face créée. Le
       sommet v est choisit uniformément parmi tous les sommets
       possibles de C respectant les deux contraintes. Si d<0, alors
       on fait comme si d=+∞ (aucun sommet ne pouvant alors être
       interne) et le résultat est un graphe planaire-extérieur
       Hamiltonien, c'est-à-dire 2-connexe. Si f<0, alors chaque face
       créée est de longueur aléatoire uniforme prise dans [3,|f|] au
       lieu d'être de longueur exactement |f|. Si f=d=4, il s'agit
       d'un "squaregraph". Les valeurs d=0,1,2 sont équivalentes.

....hyperbolic p k h
....
       Graphe issu du pavage du plan hyperbolique ou euclidien par des
       polygones réguliers à p≥3 cotés où chaque sommets est de degré
       k≥2. Le graphe est construit par couches successives de
       polygones, le paramètre h≥1 représentant le nombre de couches.
       Lorque h=1, il s'agit d'un seul polygone, un cycle à p sommets.
       Dans tous les cas les graphes sont planaires avec O((pk)^h)
       sommets, ils sont 2-connexes et d'arboricité 2 pour p>3.
       Lorsque p=3, ils sont 3-connexes et d'arboricité 3. Dans être
       les mêmes graphes, il y a des similarités avec les graphes
       "planar n f d". Pour paver le plan hyperbolique, représentable
       sur le disque de Poincaré, il faut 1/p + 1/k < 1/2. Dans ce
       cas, le graphe est d'hyperbolicté O(p). Pour paver le plan
       euclidien il faut prendre p=k=4 (grille carrée) ou p=3, k=6
       (grille triangulaire) ou p=6, k=3 (grille hexagonale). Si k≤3
       et p≤5, alors le graphe existe que pour certaines valeurs de
       h≤3.  Le cas k=3, p=4, h=2 correspond au cube, et k=3, p=5, h=3
       est le graphe dodécahèdre ("dodecahedron"). Si k=2, alors h=1
       et le graphe est un cycle à p sommets.

....rlt p q d
....
       Random Lattice Triangulation. Il s'agit d'un graphe planaire
       aléatoire construit à partir d'une triangulation de la grille p
       × q. Toutes les faces, sauf celle extérieure, sont des
       triangles. Il possède pq sommets et (2p-1)(2q-1)-pq arêtes,
       dont les 2(p+q-2) qui sont sur le bord de la grille et le reste
       est à l'intérieur. Le paramètre d contrôle la longueur des
       arêtes suivant la norme Lmax. Par exemple, si d=1, les arêtes
       seront soit celles de la grille (verticales ou horizontales) ou
       diagonales. Si d<0, alors l'effet est similaire à d=+∞. Si d=0,
       on obtient un stable. Si p=1 ou q=1 (et d≠0), on obtient un
       chemin.
....
       Il est difficile de générer de telles grilles aléatoirement
       uniformément. Il faut pour cela utiliser une technique de flips
       avec une chaîne de Markov dont le mixing time n'est pas
       connue. Il est connue que le milieu de chaque arête e d'une
       triangulation T de grille a pour coordonnées (i+1/2,j) ou
       (i+1/2,j+1/2) où i,j sont des entiers. Et inversement, chaque
       point "milieu" de la grille est coupé par exactement une arête
       de T. Il est aussi connue que si l'on parcoure les points
       milieux de la grille de bas en haut et de gauche à droit, alors
       il n'y qu'au plus deux choix possibles pour l'arête ayant ce
       milieu. Malheureusement, suivre ce parcours et choisir
       aléatoirement l'une ou l'autre de ces arêtes ne donne pas une
       distribution aléatoire uniforme. On propose ici la construction
       d'une triangulation T de manière aléatoire comme suit:
....
       Tant qu'il reste au moins un point milieu faire:
       1. choisir uniformément un point milieu R parmi les points restant
       2. déterminer la liste L des arêtes possibles ayant pour milieu R
          (respectant la planarité et le critère de longueur)
       3. choisir uniformément une arête de L et l'ajouter au graphe

....kneser n k r
....
       Graphe de Kneser généralisé. Le graphe de Kneser K(n,k)
       classique est obtenu avec r=0. Les sommets sont tous les
       sous-ensembles à k éléments de [0,n[ (il faut donc 0≤k≤n).
       Deux sommets sont adjacents ssi leurs ensembles correspondant
       ont au plus r éléments en commun. Le nombre chromatique de
       K(n,k), établit par Lovász, vaut n-2k+2 pour tout n≥2k-1>0. Le
       graphe de Petersen est le graphe K(5,2). Ils ont un lien avec
       les graphes de Johnson J(n,k).

....gpetersen n r
....
       Graphe de Petersen généralisé P(n,r), 0≤r<n/2. Ce graphe
       cubique possède 2n sommets qui sont u_1,…,u_n,v_1,…,v_n.  Les
       arêtes sont, pour tout i: u_i-u_{i+1}, u_i-v_i et v_i-v_{i+r}
       (indice modulo n). Il peut être dessiné tel que toute ses
       arêtes sont de même longueur (unit distance graph). Ce graphe
       est biparti ssi n est pair et r est impair. C'est un graphe de
       Cayley ssi r²=1 (modulo n). P(n,r) est Hamiltonien ssi r≠2 ou
       n≠5 (modulo 6). P(n,r) est isomorphe à P(n,(n-2r+3)/2)).
       P(4,1) est le cube, P(5,2) le graphe de Petersen, P(6,2) le
       graphe de Dürer, P(8,2) le graphe de Möbius-Kantor, P(10,2) le
       dodécaèdre, P(10,3) le graphe de Desargues, P(12,5) le graphe
       de Nauru, P(n,1) un prisme.

....squashed n k p
....
       Squashed Cube aléatoire à n sommets. Il faut 0<k<n et p∈[0,1].
       Les sommets sont des mots aléatoires de k lettres sur {0,1,'*'}
       où p est la probabilité d'obtenir '*'. La probabilité d'obtenir
       0 est la même que celle d'obtenir 1, soit (1-p)/2. Deux sommets
       sont adjacents si la distance de Hamming entre leur mots vaut 1
       avec la convention que la distance à la lettre '*' est nulle.
       Lorsque que p=0, le graphe généré est un sous-graphe
       isométrique de l'hypercube où certains sommets sont dupliqués
       en sommets jumeaux non-adjacents (ce sont les sommets
       correspondant au même mot). En particulier, le graphe est
       biparti et la distance entre deux sommets est données par la
       distance de Hamming entre leur mot, sauf s'ils ont le même
       mot. Si p=-1 alors p est fixée à l'équiprobabilité de chacune
       des lettres, soit p=1/3.

....antiprism n
....
       Graphe composé de deux cycles à n sommets connectés par 2n
       triangles. Le prisme est similaire sauf que pour ce dernier les
       deux cycles sont connectés par des carrés. Il est ainsi
       planaire, 4-régulier, possède 2n sommets et 4n arêtes. C'est
       aussi le dual du trapézoèdre n-gonal.

....rpartite k a_1 … a_k
....
       Graphe k-parti complet K_{a_1,…,a_k}. Ce graphe possède
       n=a_1+…+a_k sommets. Les sommets sont partitionnés en k parts
       comme suit: part 1 = [0,a_1[, part 2 = [a_1,a_1+a_2[, … part k
       = [a_1+…+a_{k-1},n[. Les sommets i et j sont adjacents ssi i et
       j appartiennent à des parts différentes.

....ggosset p k d_1 v_1 … d_1 v_k
....
       Graphe de Gosset généralisé. Les sommets sont tous les vecteurs
       entiers de dimension d = d_1 + … + d_k dont les coordonnées
       comprennent, pour i=1…k, exactement d_i fois la valeur v_i. Il
       existe une arête entre les vecteurs u et v si et seulement le
       produit scalaire entre u et v vaut l'entier p. Des valeurs
       intéressantes sont par exemple: 1 2 2 -1 2 0 ou encore 8 2 2 3
       6 -1 (le graphe de Gosset).

....schlafli
....
       Graphe de Schläfli. Il s'agit du sous-graphe induit par les
       voisins d'un quelconque sommet du graphe de Gosset. Il possède
       27 sommets, 216 arêtes et est 16-régulier. Il est sans-griffe,
       Hamiltonien, de diamètre 2, de maille 3 et de nombre
       chromatique 9.

....crown n
....
       Graphe biparti régulier à 2n sommets où le i-ème sommet de la
       première partie de taille n est voisin au j-ème sommet de la
       seconde partie ssi i≠j. Pour n=3, il s'agit du cycle à 6
       sommets, pour n=4, il s'agit du cube (à 8 sommets).

....split n k
....
       Graphe split (ou fendu) à n sommets et de clique maximum k. Il
       s'agit d'un graphe à n sommets composé d'une clique à k sommets
       et d'un ensemble indépendant de n-k sommets connectés chacun à
       tous ceux de la clique. C'est un graphe triangulé (chordal) et
       un cas particulier de "kstar n k". On peut montrer que presque
       tous les graphes triangulés à n sommets sont des graphes split
       (Bender et al. 1985). Si k≥n-1, alors il s'agit d'une clique,
       et si k=n-2, il s'agit d'une clique moins une arête. Si k=1 il
       s'agit d'une étoile à n-1 branches.

....centipede n
....
       Arbre à 2n sommets et n feuilles en forme de peigne.

....fan p q
....
       Graphe composé d'un chemin à p sommets et de q sommets, chacun
       connectés à tous ceux du chemin. Le graphe classique "fan n"
       correspond à p=n et q=1.

....flip n
....
       Graphe des flips des triangulations d'un polygone convexe à n>2
       sommets. Les sommets, qui sont les triangulations, sont en
       bijection avec des arbres binaires complets à 2n-3 noeuds qui
       sont codés par les mots de Dyck de longueur 2n-4 (que l'on peut
       afficher avec -label 1). Le nombre de sommets est donc C(n-2)
       nombre de Catalan d'ordre n-2. Les adjacences peuvent être vues
       aussi comme des rotations d'arbres. Le diamètre est 2n-10 pour
       n>12. Le nombre chromatique n'est pas connu, on sait pas s'il
       est constant ou pas. Il vaut 3 pour n=5..9, et 4 pour n=10 et
       11.

....interval n
....
       Graphe d'intersection de n intervalles d'entiers aléatoires
       uniformes pris dans [0,2n[. Des graphes d'intervalles peuvent
       aussi être générés par "kpath n k".

....circle n
....
       Graphe d'intersection de n cordes aléatoires d'un cercle. Il
       est réalisé par le graphe d'inclusion de n intervalles
       d'entiers aléatoires uniformes pris dans [0,2n[. Les graphes de
       permutation et les planaires extérieurs sont des exemples de
       circle graphs.

....permutation n
....
       Graphe de permutation sur une permutation aléatoire uniforme
       des entiers de [0,n[.

....prime n
....
       Graphe à n sommets tel que i est adjacent à j ssi i>1 et j
       divisible par i.

....paley n
....
       Graphe de Paley à n sommets. Deux sommets sont adjacents ssi
       leur différence est un carré modulo n. Il faut que n soit la
       puissance d'un nombre premier et que n≡1 (mod 4), mais le
       graphe est aussi défini pour les autres valeurs. Les premières
       valeurs possibles pour n sont: 5, 9, 13, 17, 25, 29, 37, 41,
       49, … Ces graphes sont Hamiltoniens. Si n est simplement
       premier, alors ils sont de plus auto-complémentaires et
       réguliers.  Paley 17 est le plus grand graphe G où ni G ni son
       complémentaire ne contient K₄, d'où Ramsey(4)=18.

....mycielski k
....
       Graphe de Mycielski de paramètre (nombre chromatique) k. C'est
       un graphe sans triangle, k-1 (sommets) connexe, et de nombre
       chromatique k. Le premier graphe de la série est M2 = K2, puis
       on trouve M3=C5, M4 est le graphe de Grötzsch à 11 sommmets.

....windmill n
....
       Graphe composé de n cycles de longueur trois ayant un sommet
       commun.

....barbell n1 n2 p
....
       Graphe des haltères (Barbell Graph) composé de deux cliques de
       n1 et n2 sommets reliées par un chemin de longueur p. Il
       possède n1+n2+p-1 sommets. Si p=0 (p=-1), le graphe est composé
       de deux cliques ayant un sommet (une arête) en commun. Plus
       généralement, si p≤0, le graphe est composé de deux cliques
       s'intersectant sur 1-p sommets.

....chess p q x y
....
       Graphe composé de p x q sommets représentant les cases d'un
       échiquier p x q, deux cases étant connectée s'il existe un
       déplacement d'une case vers l'autres avec un saut de x cases
       selon un axe et y selon un autre. Le "knight graph" classique
       est donc un "chess 8 8 2 3", et "chess n 2 1 0" correspond à
       "ladder n".

....sat n m k
....
       Graphe aléatoire issu de la réduction du problème k-SAT à
       Vertex Cover. Le calcul d'un Vertex Cover de taille minimum
       pour ce graphe est donc difficile pour k>2. Soit F une formule
       de k-SAT avec n>0 variables x_i et m>0 clauses CNF de k>0
       termes.  Le graphe généré par "sat n m k" possède un Vertex
       Cover de taille n+(k-1)m si et seulement si F est satisfiable.
....
       Ce graphe est composé d'une union de n arêtes indépendantes et
       de m cliques à k sommets, plus des arêtes dépendant de F
       connectant certains sommets des cliques aux n arêtes. Les n
       arêtes représentent les n variables, une extrémité pour x_i,
       l'autre pour ¬(x_i). Ces sommets ont des numéros dans [0,2n[,
       x_i correspond au sommet 2i-2 et ¬(x_i) au sommet 2i-1,
       i=1…n. Les sommets des cliques ont des numéros consécutifs ≥ 2n
       et correspondent aux clauses. Le p-ème sommet de la q-ème
       clique (pour p=1…k et q=1…m) est connecté à l'une des
       extrémités de la i-ème arête (pour i=1…n) ssi la p-ème variable
       de la q-ème clause est x_i ou ¬(x_i).
....
       La formule F est construite en choisissant indépendamment et
       uniformément pour chacune des m clauses et chacun des k termes
       une des variables parmi x_1,…,x_n,¬(x_1),…,¬(x_n).  Ainsi
       chaque sommet d'une clique possède exactement un voisin (choisi
       aléatoirement uniforme) parmi les 2n extrémités d'arêtes.

....kout n k
....
       Graphe à n sommets k-dégénéré crée par le processus aléatoire
       suivant: les sommets sont ajoutés dans l'ordre croissant de
       leur numéro, i=0,1,…,n-1. Le sommet i est connecté à d voisins
       qui sont pris aléatoirement uniformément parmi les sommets dont
       le numéro est < i. La valeur d est choisie aléatoirement
       uniformément entre 1 et min{i,k}. Il faut k>0. Le graphe est
       connexe, et pour k=1, il s'agit d'un arbre.

....expander n k
....
       Graphe à n sommets composé de k>0 cycles Hamiltoniens
       aléatoires. Le degré des sommets varie entre 2 et 2k. Il
       possède le cycle 0,1,…,n-1,0 comme cycle Hamiltonien, et a la
       propriété d'expansion à partir de k>3. Plus précisément, avec
       grande probabilité, les valeurs propres de la matrice
       d'adjacence du graphe sont ≤ 2√(2k).

....comb n
....
       Arbre de 2n sommets en forme de peigne, composé d'un chemin à n
       sommets avec un sommet pendant à chacun d'eux.

....sunlet n
....
       Cycle à n sommets avec un sommet pendant à chacun d'eux. Un
       sunlet 3 est parfois appelé netgraph.

....alkane t n
....
       Graphe planaires dont les sommets sont de degré 1 ou 4
       représentant la structure moléculaire d'hydrocarbure alkalin à
       n atomes de carbones.  Le paramètre t (voir ci-dessous les six
       valeurs possibles) contrôle la topologie des liaisons simples
       entre atomes de carbone (C), les atomes d'hydrogènes (H) étant
       des sommets pendants de sorte que chaque atome C soit de degré
       4. Certaines topologies ne sont définies que pour certaines
       valeurs de n. Les alkalins, de formule C_n H_{2n+2} si
       type≠"cyclo", sont des arbres de 3n+2 sommets alors que les
       cycloalkalin, de formule C_n H_{2n} si type="cyclo", ont 3n
       sommets et possède un cycle. Le type "normal" peut être abrégé
       en "n". L'option "-label 1 -visu" permet de distinguer les
       atomes C et H.
....
       !!! t      topologie     n             t     topologie     n
                                                    C─┐
          normal  C─ ⋯ ─C      n≥1           neo      C─ ⋯ ─C    n≥5
                                                    C─┘
                  ┌─C─ ⋯ ─C                           ┌─C─ ⋯ ─C
          cyclo   C       │    n≥3           sec    C─C           n≥6
                  └─C─ ⋯ ─C                           └─C─ ⋯ ─C
                  C─┐                                  ┌─C─ ⋯ ─C
          iso       C─ ⋯ ─C    n≥4           tret   C─C─ ⋯ ─C    n≥7
                  C─┘                                  └─C─ ⋯ ─C
....
       Il est possible d'utiliser les alias suivants:
....
       !!!  n-alkane n ......... (= alkane normal n) 
            cyclo-alkane n ..... (= alkane cyclo n)
	    iso-alkane n ....... (= alkane iso n)
	    neo-alkane n ....... (= alkane neo n)
	    sec-alkane n ....... (= alkane sec n)
	    tret-alkane n ...... (= alkane tret n)
	    methane ............ (= alkane n 1)
	    ethane ............. (= alkane n 2)
	    propane ............ (= alkane n 3)
	    butane ............. (= alkane n 4)
	    pentane ............ (= alkane n 5)
	    hexane ............. (= alkane n 6)
	    heptane ............ (= alkane n 7)
	    octane ............. (= alkane n 8)
	    nonane ............. (= alkane n 9)
....
      Il est possible aussi de combiner les préfixes cyclo-, iso-,
      neo-, sec-, tret- avec les radicaux meth, eth, prop, but, hex,
      hept, oct, non lorsque la condition sur n est satisfaite. Par
      exemple, cyclo-pentane (= alkane cyclo 5) et iso-butane (=
      alkane iso 4).

....icosahedron
....
       Isocahèdre: graphe planaire 5-régulier à 12 sommets. Il possède
       30 arêtes et 20 faces qui sont des triangles. C'est le dual du
       dodécahèdre.

....cuboctahedron
....
       Cuboctaèdre: graphe planaire 4-régulier à 12 sommets. Il
       possède 24 arêtes et 14 faces qui sont des triangles ou des
       carrés. C'est le dual du rhombicdodécahèdre.

....rdodecahedron
....
       Rhombic-dodécaèdre: graphe planaire à 14 sommets avec des
       sommets de degré 3 ou 4. Il possède 21 arêtes et 12 faces qui
       sont des carrés. C'est le dual du cuboctahèdron.

....deltohedron n
....trapezohedron n
....
       Deltoèdre ou trapézoèdre n-gonal: graphe composé de 2n faces en
       forme de cerf-volant (deltoïdes) décalées symétriquement. C'est
       un graphe planaire de 2n+2 sommets et 4n arêtes où toutes les
       faces sont des carrées. C'est aussi le dual de l'antiprisme
       n-gonal. Il s'agit d'un cube si n=3.

....tutte
....
       Graphe de Tutte. C'est un graphe planaire cubique 3-connexe à
       46 sommets qui n'est pas Hamiltonien.

....hgraph
....
       Arbre à six sommets dont quatre feuilles en forme de H.

....rgraph
....fish
....
       Fish Graph, graphe à six sommets en forme de R ou de
       poisson. Il est composé d'un cycle à quatre sommets dont un
       ayant deux sommets pendants.

....cricket
....
       Cricket Graph, graphe à cinq sommets composé d'un triangle où à
       l'un de sommets est attaché deux sommets pendant (de degré 1).

....moth
....
       Moth Graph, graphe à six sommets composé de deux triangles
       partageant une arête et de deux sommets pendant (degré 1)
       attachés à l'un des sommets commun aux deux triangles.

....bull
....
       Bull Graph, graphe à cinq sommets auto-complémentaire en forme
       de A.

....suzuki
....
       Graphe de Suzuki (2010). C'est l'unique graphe 1-planaire à
       n=11 sommets et ayant le nombre optimal d'arêtes, soit 4n-8
       arêtes (ici 36 donc).

....harborth
....
       Graphe de Harborth. C'est un graphe planaire 4-régulier à 52
       sommets qui est distance unitaire aussi appelé graphe allumette
       (voir theta0 et diamond). Il peut ainsi être dessiné sans
       croisement d'arête qui ont toutes la même longueur.

....doily
....
       Graphe Doily (de Payne). C'est un graphe de 15 sommets qui est
       un carré généralisé pouvant être représenté par 15 points et 15
       lignes, avec 3 points par ligne et 3 lignes par point, et sans
       triangle.

....bidiakis
....
       Graphe ou cube de Bidiakis. C'est un graphe planaire cubic à 12
       sommets. On peut le représenter comme un cube où deux faces
       opposées comportent une arête supplémentaire perpendiculaire
       joignant deux bord opposés. Il est Hamiltonien et son nombre
       chromatique est 3.

....herschel
....
       Graphe de Herschel. C'est le plus petit graphe planaire
       3-connexe qui ne soit pas Hamiltonien. Il est biparti, possède
       11 sommets et 18 arêtes.

....goldner-harary
....
       Graphe de Goldner-Haray. C'est le plus petit graphe planaire
       maximal qui ne soit pas Hamiltonien. Il possède 11 sommets et
       donc 27 arêtes (voir aussi Herchel). C'est un 3-arbre planaire
       (voir apollonian).

....fritsch
....
       Graphe de Fritsch. Il est planaire maximal à 9 sommets qui peut
       être vu comme un graphe de Hajós dans un triangle. C'est, avec
       le graphe de Soifer, le plus petit contre-exemple à la
       procédure de coloration de Kempe. C'est le plus petit graphe où
       l'heuristique de degré minimum donne cinq couleurs.

....triplex
....
       Graphe cubic de maille 5 à 12 sommets 1-planaire pouvant être
       dessiné avec seulement deux croisements d'arête. Un des cinq
       graphes (avec le Petersen) a être cycliquement-5-connexe
       (McCuaig 1992).

....jaws
....
       Graphe cubic de maille 5 à 20 sommets qui est un doublecross,
       c'est-à-dire dessinable sur le plan avec deux paires d'arêtes
       se croisant sur la face extérieure. Il est donc
       1-planaire. Tout graphe theta-connecté sans Petersen mais avec
       Jaws comme mineur est un doublecross.

....starfish
....
       Graphe cubic de maille 5 à 20 sommets non-planaire, mais
       peut-être dessiné comme une étoile à cinq branches avec une
       couronne centrale à 15 sommets formant un circulant avec une
       corde de longueur 3. Un graphe theta-connecté (cf. Seymour et
       al. 2015) ssi il ne contient pas de Petersen comme mineur, si
       c'est un graphe apex (planaire plus un sommet), un doublecross
       (voir jaws) ou un starfish.

....soifer
....
       Graphe de Soifer. Il est planaire maximal à 9 sommets. C'est,
       avec le graphe de Fritsch, le plus petit contre-exemple à la
       procédure de coloration de Kempe. C'est le plus petit graphe où
       l'heuristique de degré minimum donne cinq couleurs.

....poussin
....
       Graphe de Poussin. Il est planaire maximal à 15 sommets. C'est
       un contre-exemple à la procédure de coloration de Kempe.

....headwood4
....
       Graphe de Headwood pour la conjecture des 4 couleurs,
       contre-exemple de la preuve de Kempe. Il est planaire maximal
       avec 25 sommets, est de nombre chromatique 4, de diamètre 5, de
       rayon 3 et Hamiltonien.

....errera
....
       Graphe d'Errera. Il est planaire maximal à 17 sommets. C'est un
       contre-exemple à la procédure de coloration de Kempe.

....kittell
....
       Graphe de Kittell. Il est planaire maximal à 23 sommets. C'est
       un contre-exemple à la procédure de coloration de Kempe.

....frucht
....
       Graphe de Frucht. Il est planaire cubique à 12 sommets. Il n'a
       pas de symétrie non-triviale. C'est un graphe de Halin de
       nombre chromatique 3, de diamètre 4 et de rayon 3.

....treep p
....
       Arbre aléatoire à p>2 feuilles sans sommets internes de degré
       deux. Il possède entre p+1 et 2p-2 sommets. Ce graphe est à la
       base de la construction des graphes de Halin.

....halin p
....
       Graphe de Halin aléatoire basé sur un arbre à p>2 feuilles. Il
       possède entre p+1 et 2p-2 sommets. Il est constitué d'un arbre
       sans sommets de degré deux dont les p feuilles sont connectés
       par un cycle (de p arêtes). Ces graphes planaires de degré
       minimum au moins trois sont aussi arête-minimale 3-connexes,
       Hamiltonien (et le reste après la suppression de n'importe quel
       sommet), de treewidth exactement 3 (ils contiennent K₄ comme
       mineur). Ils contiennent toujours au moins trois triangles et
       sont de nombre chromatique 3 ou 4.

....butterfly d
....
       Graphe Butterfly de dimension d. Les sommets sont les paires
       (x,i) où x est un mot binaire de d bits et i un entier de
       [0,d]. Les sommets peuvent être représentés en d+1 niveaux
       chacun de 2^d sommets, les arêtes connectant les niveaux
       consécutifs. Le sommet (x,i) est adjacent à (y,i+1) ssi les
       bits de x sont identiques à ceux de y sauf pour celui de numéro
       i+1 (le bit 1 étant le bit de poids le plus faible). Il possède
       (d+1)*2^d sommets et d*2^(d+1) arêtes, les sommets de niveau 0
       et d étant de degré 2 les autres de degré 4.

....line-graph n k
....
       Line-graphe aléatoire à n sommets et de paramètre k>0 entier
       Plus k est petit, plus le graphe est dense, le nombre d'arêtes
       étant proportionnel à (n/k)². Si k=1, il s'agit d'une clique à
       n sommets. Ces graphes sont obtenus en choisissant, pour chaque
       sommet, deux couleurs de [1,k]. Deux sommets sont alors
       adjacents ssi ils possèdent la même couleur.  Ces graphes sont
       claw-free (sans K_{1,3} induit). Comme les graphes claw-free,
       les line-graphes connexes avec un nombre pair de sommets
       possèdent toujours un couplage parfait. On rappel qu'un graphe
       G est le line-graphe d'un graphe H si les sommets de G
       correspondent aux arêtes de H et où deux sommets de G sont
       adjacents ssi les arêtes correspondantes dans H sont
       incidentes. Pour parle parfois de graphe adjoint.

....shuffle d
....
       Graphe Shuffle-Exchange de dimension d. Les sommets sont les
       mots binaires de d lettres. Le sommet w et w' sont voisins si w
       et w' diffèrent du dernier bit, ou bien si w' peut être obtenu
       par décalage cyclique à droite ou à gauche de w.

....debruijn d b
....
       Graphe de De Bruijn de dimension d≥0 et de base b>0. Il a b^d
       sommets qui sont tous les mots de d lettres sur un alphabet de
       b lettres. Le sommet (x_1,…,x_d) est voisin des sommets
       (x_2,…,x_d,*). Ce graphe est Hamiltonien, de diamètre d et le
       degré de chaque sommet est 2b, 2b-1 ou 2b-2. Pour d=3 et b=2,
       le graphe est planaire.

....kautz d b
....
       Graphe de Kautz de dimension d>0 et de base b>1. Il a
       b*(b-1)^(d-1) sommets qui sont tous les mots de d lettres sur
       un alphabet de b lettres avec la contrainte que deux lettres
       consécutives doivent être différentes. L'adjacence est celle du
       graphe de De Bruijn. C'est donc un sous-graphe induit de De
       Bruijn (debruijn d b). Il est Hamiltonien, de diamètre d et le
       degré de chaque sommet est 2b-2 ou 2b-3. Pour d=b=3 le graphe
       est planaire.

....linial n t
....
       Neighborhood graph des cycles introduit par Linial. C'est le
       graphe de voisinage des vues de taille t d'un cycle orienté
       symétrique à n sommets ayant des identifiants uniques de
       [0,n[. Il faut n≥t>0. Le nombre chromatique de ce graphe est k
       ssi il existe un algorithme distribué qui en temps t-1
       (resp. en temps (t-1)/2 avec t impair) peut colorier en k
       couleurs tout cycle orienté (resp. orienté symétrique) à n
       sommets ayant des identifiants uniques et entiers de [0,n[. Les
       sommets sont les t-uplets d'entiers distincts de [0,n[. Le
       sommet (x_1,…,x_t) est voisin des sommets (x_2,…,x_t,y) où
       y≠x_1 si n>t et y=x_1 si n=t.  C'est un sous-graphe induit de
       linialc n t, et donc du graphe de Kautz (kautz t n) et de De
       Bruijn (debruijn t). Le nombre de sommets est m(m-1)…(m-t+1).
       Certaines propriétés se déduisent du graphe linialc n t.

....linialc m t
....
       Neighborhood graph des cycles colorés.  Il s'agit d'une
       variante du graphe linial n t. La différence est que les
       sommets du cycle n'ont plus forcément des identités uniques,
       mais seulement une m-coloration, m≤n. L'adjacence est
       identique, mais les sommets sont les (2t+1)-uplets
       (x_1,…,x_{2t+1}) d'entiers de [0,m[ tels que x_i≠x_{i+1}. Il
       s'agit donc d'un sous-graphe induit de linialc m t, lui-même
       sous-graphe induit du graphe de Kautz (kautz 2t+1 m) et donc de
       De Bruijn (debruijn 2t+1 m). Le nombre de sommets est
       m(m-1)^{2t} et son degré maximum est 2(m-1). La taille de la
       clique maximum est 3 si m>2 et t>1. Le nombre chromatique de ce
       graphe est 3 pour m=4, 4 pour 5≤m≤24. Pour 25≤m≤70 c'est au
       moins 4 et au plus 5, la valeur exacte n'étant pas connue.

....pancake n
....
       Graphe "pancake" de dimension n. Il a n! sommets qui sont les
       permutations de {1,…,n} et (n-1)-régulier. Une permutation,
       c'est-à-dire un sommet, est voisine de toutes celles obtenues
       en retournant un de ces préfixes. Plus précisément, les sommets
       x=(x_1,…,x_n) et y=(y_1,…,y_n) sont adjacents s'il existe
       un indice k tel que y_i=x_i pour tout i>k et y_i=x_{k-i} sinon.
       Son diamètre, qui est linéaire en n, n'est pas connu
       précisément. Les premières valeurs connues, pour n=1…17,
       sont: 0, 1, 3, 4, 5, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18,
       19. Donc les diamètres 2,6,12 n'existent pas.

....bpancake n
....
       Graphe "burn pancake" de dimension n. Il a n!*2^n sommets qui
       sont les permutations signées de {1,…,n}. Les sommets
       x=(x_1,…,x_n) et y=(y_1,…,y_n) sont adjacents s'il existe
       un indice k tel que y_i=x_i pour tout i>k et y_i=-x_{k-i}
       sinon. Dit autrement la permutation de y doit être obtenue en
       retournant un préfixe de x et en inversant les signes. Par
       exemple, le sommet (+2,-1,-5,+4) est voisin du sommet
       (+5,+1,-2,+4). Comme le graphe pancake, c'est un graphe
       (n-1)-régulier de diamètre linéaire en n.

....gpstar n d
....
       Graphe "permutation star" généralisé de dimension n. Il a n!
       sommets qui sont les permutations de {1,…,n}. Deux sommets sont
       adjacents si leurs permutations diffèrent par d positions. Si
       d<2, il s'agit d'un stable. C'est un graphe régulier.

....pstar n
....
       Graphe "permutation star" de dimension n. Il a n! sommets qui
       sont les permutations de {1,…,n}. Deux sommets sont adjacents
       si une permutation est obtenue en échangeant le premier élément
       avec un autre. Le graphe est (n-1)-régulier. Le graphe est
       biparti et de diamètre ⎣ 3(n-1)/2)⎦. C'est un sous-graphe induit
       d'un "gpstar n 1".

....hexagon p q
....
       Grille hexagonale p x q. C'est un planaire composé de p rangées
       de q hexagones, le tout arrangé comme un nid d'abeille. Ce
       graphe peut aussi être vu comme un mur de p rangées de q
       briques, chaque brique étant représentée par un cycle de
       longueur 6. Il possède (p+1)*(2p+2)-2 sommets et est de degré
       maximum 3. Sont dual est le graphe whexagon.
....
       Ex: gengraph hexagon 20 20 -dele 0.2 -maincc -visu

....whexagon p q
....
       Comme le graphe hexagon p q sauf que chaque hexagone est
       remplacé par une roue de taille 6 (chaque hexagone possède un
       sommet connecté à ses 6 sommets). C'est le dual de l'hexagone.
       Il possède p*q sommets de plus que l'hexagone p q.

....hanoi n b
....
       Graphe de Hanoï généralisé, le graphe classique est obtenu avec
       b=3. Il est planaire avec b^n sommets et est défini de manière
       récursive comme suit. Le niveau n>0 est obtenu en faisant b
       copies du niveau n-1 qui sont connectés comme un cycle par une
       arête, le niveau 0 étant le graphe à un sommet. On obtient le
       graphe de sierpinski n b en contractant ces arêtes là. Il faut
       b≥2 et n≥0. Lorsque n=2, on obtient un sorte de fleur, pour
       n=1 c'est un cycle et pour b=2 il s'agit d'un chemin.

....sierpinski n b
....
       Graphe de Sierpiński généralisé, le graphe classique, le
       triangle Sierpiński qui est planaire, est obtenu avec b=3. Il a
       ((b-2)*b^n+b)/(b-1) sommets et est défini de manière récursive
       comme suit.  Le niveau n est obtenu en faisant b copies du
       niveau n-1 qui sont connectés comme un cycle, le niveau 1 étant
       un cycle de b sommets. Il faut b≥3 et n≥1. A la différence du
       graphe d'Hanoï, les arêtes du cycle sont contractées. Le graphe
       de Hajós est obtenu avec n=2 et b=3. Pour n=1 il s'agit d'un
       cycle.

....banana n k
....
       Arbre à n(k+1)+1 sommets composés de n>0 copies d'étoiles à k>0
       branches connectées par une feuille à un unique sommet.

....moser
....
       Graphe "Moser spindle" découvert par les frères Moser. C'est un
       "unit distance graph" du plan (deux points sont adjacents s'ils
       sont à distance exactement 1) de nombre chromatique 4. Il est
       planaire et possède 7 sommets. On ne connaît pas d'unit
       distance graphe avec un nombre chromatique supérieur.

....markstrom
....
       Graphe de Markström. Il est cubique planaire à 24 sommets. Il
       n'a pas de cycle de longueur 4 et 8.

....robertson
....
       Graphe de Robertson. C'est le plus petit graphe 4-régulier de
       maille 5. Il a 19 sommets, est 3-coloriable et de diamètre 3.

....wiener-araya
....
       Graphe découvert en 2009 par Wiener & Araya. C'est le plus
       petit graphe hypo-Hamiltonien planaire connu, c'est-à-dire
       qu'il n'est pas Hamiltonien mais la suppression de n'importe
       quel sommet le rend Hamiltonien. Il possède 42 sommets, 67
       arêtes, et est de diamètre 7.

....zamfirescu
....
       Graphe de Zamfirescu à 48 sommets découvert en 2007. Il est
       planaire et hypo-Hamiltonien. C'est le second plus petit (voir
       wiener-araya). Il possède 76 arêtes et a un diamètre de 7.

....hatzel
....
       Graphe de Hatzel. Il est planaire, de diamètre 8, possède 57
       sommets et 88 arêtes, et est hypo-Hamiltonien (voir
       wiener-araya). C'était le plus petit planaire hypo-Hamiltonien
       connu avant le graphe de Zamfirescu.

....clebsch n
....
       Graphe de Clebsch d'ordre n. Il est construit à partir d'un
       hypercube de dimension n en ajoutant une arête entre chaque
       paire de sommets opposés, c'est-à-dire à distance n. Le graphe
       classique de Clebsch est réalisé pour n=4 dont le diamètre est
       deux.

....gear n
....
       Graphe planaire à 2n+1 sommets composé d'une roue à n rayons
       ("wheel n") et de n sommets chacun connecté deux sommets
       voisins du bord de la roue. Il est construit à partir du graphe
       "cage 2n 2 2 0" auquel on ajoute un sommet central. Pour n=3,
       c'est le complémentaire de "helm n".

....helm n
....
       Graphe planaire à 2n+1 sommets composé d'une roue à n≥3 rayons
       ("wheel n") et de n sommets pendants connectés au bord de la
       roue. Pour n=3, c'est le complémentaire de "gear n".

....haar n
....
       Graphe de Haar H(n) pour l'entier n>0. C'est un graphe biparti
       régulier possédant 2k sommets où k=1+⎣ log₂(n)⎦ est le nombre
       de bits dans l'écriture binaire de n. Ses sommets sont les
       u_i,v_i pour i=0,…,k-1, et u_i est adjacent à v_{i+j mod k} ssi
       le bit j de n vaut 1 (j=0,…,k-1). Si n est impair, H(n) est
       connexe et de maille 4 ou 6 (il faut n≥1). La valeur maximale
       de n est 2^32-1 = 4.294.967.295 correspondant à un graphe de 64
       sommets. Pour n=37, 69 et 133, on retrouve respectivement les
       graphes de Franklin, de Headwood et de Möbius-Kantor. On a
       aussi H(2^n-1) est le biparti K_{n,n}, H(3*2^{n-2}-1) est le
       "crown n", H(2^{n-1}) est n copies de K_2, H(2^{n-1}+1) est un
       cycle à n sommets, H(2^n+3) est le "mobius 2n+1" et
       H(2^{n-1}+3) est le "prism 2n".

....turan n r
....
       Graphe de Turán à n sommets et r parts. Il s'agit qu'un graphe
       r-parti complet de n sommets avec r-(n mod r) parts de ⎣ n/r⎦
       sommets et (n mod r) parts de ⎡ n/r⎤ sommets. Il faut n≥r>0. On
       peut aussi le définir comme le graphe ayant une arête entre i
       et j ssi |i-j| mod r ≠ 0. C'est la définition utilisée pour
       générer ce graphe. Il est régulier lorsque r divise n. Il
       possède ⎣(r-1)n^2/(2r)⎦ arêtes et est de nombre chromatique
       r. C'est le graphe sans clique de taille r+1 ayant le plus
       grand nombre d'arêtes. Lorsque n=r, il s'agit d'une
       clique. Lorsque n=r+1, il s'agit d'une clique moins une
       arête. Lorsque n=2r, il s'agit du "cocktail party
       graph". Lorsque n=3r, il s'agit du graphe de Moon-Moser, le
       graphe possédant le plus grand nombre de cliques maximales
       (soit 3^(n/3)).Lorsque n=6 et r=3, c'est l'octaèdre.

....flower_snark n
....
       Graphe cubique à 4n sommets construit de la manière suivante:
       1) on part de n étoiles disjointes à 3 feuilles, la i-ème ayant
       pour feuilles les sommets notés u_i,v_i,w_i, i=1…n; 2) pour
       chaque x∈{u,v,w}, x_1-…-x_n induit un chemin; et enfin 3) sont
       adjacents: u_0-u_n, v_0-w_n et w_0-v_n. Pour n>1, ces graphes
       sont non planaires, non Hamiltoniens, 3-coloriables et de
       maille au plus 6. Pour n=1, il s'agit d'un K_{1,3} (claw).

....udg n r
....
       Graphe géométrique aléatoire (random geometric graph) sur n
       points du carré [0,1[ × [0,1[ (distribution par défaut). Deux
       sommets sont adjacents si leurs points sont à distance ≤ r.  Il
       s'agit de la distance selon la norme L2 (par défaut), mais
       cela peut être changée par l'option -norm. Le graphe devient
       connexe avec grande probabilité lorsque r=rc ~ √(ln(n)/n). Si
       r<0, alors le rayon est initialisé à rc. Un UDG (unit disk
       graph) est normalement un graphe d'intersection de disques
       fermés de rayon 1.

....gabriel n
....
       Graphe de Gabriel. Graphe géométrique défini à partir d'un
       ensemble de n points du carré [0,1[ × [0,1[ (distribution par
       défaut). Les points i et j sont adjacents ssi le plus petit
       disque (voir -norm) passant par i et j ne contient aucun autre
       point. Ce graphe est connexe et planaire à condition toutefois
       qu'ils n'existe pas 4 points co-cycliques et par paires
       diamétralement opposées (dans ce cas une clique et des
       croisements d'arêtes apparaissent). C'est un sous-graphe du
       graphe de Delaunay. Son étirement est non borné.

....rng n
....
       Graphe du proche voisinage (Relative Neighborhood Graph).
       Graphe géométrique défini à partir d'un ensemble de n points du
       carré [0,1[ × [0,1[.  Les points i et j sont adjacents ssi il
       n'existe aucun point k tel que max{d(k,i),d(k,j)} < d(i,j) où d
       est la distance (L2 par défaut, voir -norm).  Dit autrement,
       la "lune" définie par i et j doit être vide. Ce graphe est
       planaire et connexe. C'est un sous-graphe du graphe de Gabriel.

....nng n
....
       Graphe du plus proche voisin (Nearest Neighbor Graph). Graphe
       géométrique défini à partir d'un ensemble de n points du carré
       [0,1[ × [0,1[.  Le point i est connecté au plus proche autre
       point (par défaut selon la norme L2, voir -norm). Ce graphe
       est une forêt couvrante du graphe rng de degré au plus 6 (si la
       norme est L2).

....thetagone n p k w
....
       Graphe géométrique défini à partir d'un ensemble de n points du
       carré [0,1[ × [0,1[. En général le graphe est planaire et
       connexe avec des faces internes de longueur au plus p (pour k
       diviseur de p et w=1). On peut interpréter les paramètres comme
       suit: p≥3 est le nombre de cotés d'un polygone régulier, k≥1 le
       nombre d'axes (ou de direction), et w∈[0,1] le cône de
       visibilité. Toute valeur de p<3 est interprétée comme une
       valeur infinie, et le polygone régulier correspondant
       interprété comme un cercle. L'adjacence entre une paire de
       sommets est déterminée en temps O(kn).
....
       Plus formellement, pour tout point u et v, et entier i, on note
       P_i(u,v) le plus petit p-gone (polygone convexe régulier à p
       cotés) passant par u et v dont u est un sommet, et dont le
       vecteur allant de u vers son centre forme un angle de i*2π/k
       avec l'axe des abscisses, intersecté avec un cône de sommet u
       et d'angle w*(p-2)*π/p (w*π si p est infini) et dont la
       bissectrice passe par le centre du p-gone. Alors, u est voisin
       de v s'il un existe au moins un entier i∈[0,k[ tel que
       l'intérieur de P_i(u,v) est vide. La distance entre u et le
       centre du p-gone définit alors une distance (non symétrique) de
       u à v.
....
       Si w=1 (visibilité maximale), P_i est précisément un p-gone. Si
       w=0 (visibilité minimale), P_i se réduit à l'axe d'angle i*2π/k
       pour un entier i. Si w=.5, P_i est un cône formant un angle
       égale à 50% de l'angle défini par deux cotés consécutifs du
       p-gone, ce dernier angle valant (p-2)π/p. Si w=2p/((p-2)k) (ou
       simplement 2/k si p est infini) alors la visibilité correspond
       à un cône d'angle 2π/k, l'angle entre deux axes. Comme il faut
       w≤1, cela implique que k≥2p/(p-2) (k≥2 si p infini). On
       retrouve le Theta_k-Graph pour chaque k≥6 en prenant p=3 et
       w=6/k, le demi-Theta-Graph pour tout k≥3 en prenant p=3 et
       w=3/k, le Yao_k-Graph pour chaque k≥2 en prenant p=0 (infini)
       et w=2/k, et la triangulation de Delaunay si p=0 (infini), k
       très grand et w=1. En fait, ce n'est pas tout-à-fait le graphe
       Yao_k, pour cela il faudrait que u soit le centre du polygone
       (c'est-à-dire du cercle).

....pat p q r
....
       Graphe possédant pqr sommets, issu d'un jeu à un joueur proposé
       par Pat Morin (Barbade, mars 2016). Le jeu se déroule sur une
       grille p×q et comprend r coups. Un coup est un ensemble de
       positions de la grille strictement croissantes (coordonnées en
       x et en y strictement croissantes). De plus, si la position
       (x,y) est jouée alors toutes les positions situées sur la même
       ligne mais avec une abscisse au moins x ou sur la même colonne
       mais avec une ordonnées au moins y sont interdites pour tous
       les coups suivants. Le score est le nombre total de positions
       jouées en r coups. Il s'agit de trouver le score maximum.
       Lorsque r=1, le score maximum vaut min(p,q). Lorsque p=q=n et
       r=2, alors le score maximum vaut ⎣ 4n/3⎦. La question est
       ouverte lorsque r>2, c'est au moins n^1.516 pour r=n où la
       constante vaut log_9(28).
....
       Les sommets du graphes sont les positions dans les r grilles
       p×q et deux sommets sont adjacents les positions sont en
       conflits. Le score du jeu est alors un ensemble indépendant du
       graphe. Si r=1, le graphe est une grille p×q. Ce graphe active
       l'option -pos car un dessin de ce graphe (sous forme de
       grilles) est proposé.
....
       Ex: gengraph pat 4 4 4 -check kindepsat 8 | ./glucose -model

....uno n p q
....
       Line-graphe aléatoire à n sommets issu d'un graphe biparti de
       parts de taille p>0 et q>0. Plus précisément, les sommets sont
       des paires d'entiers aléatoires de [0,p[ × [0,q[, pas
       nécessairement distinctes. Les sommets (x,y) et (x',y') sont
       adjacents ssi x=x' ou y=y'. C'est un sous-graphe induit du
       produit cartésien de deux cliques, K_p × K_q. Ce graphe est
       géométrique. Les sommets représentent aussi des cartes du jeu
       de UNO et les arêtes indiquent si un carte peut être jouée
       consécutivement à une autre.

....unok n p q kp kq
....
       Graphe "uno n p q" particulier où les n points correspondant
       aux sommets sont pris uniformément parmi les ensembles de n
       points distincts de [0,p[ × [0,q[ ayant au plus kp sommets par
       ligne et kq par colonne. Il faut n≤min{p*kp,q*kq} et p, q, kp,
       kq>0. Si kp<0, alors on fait comme si kp=p, de même pour kq=q s
       kq<0.  Contrairement à uno, deux sommets ont toujours des
       coordonnées distinctes. Le graphe résultant est de degré au
       plus kp+kq-2, et est de path-width (et aussi de tree-width) au
       plus celle du produit de clique K_{kp} × K_{kq} soit environ
       kp*kq/2. Le temps de génération des n points est en O(npq)
       contre O(n) pour uno, mais une optimisation (algorithme par
       rejets) fait qu'il est très souvent en seulement O(n+p+q), dans
       les cas peu dense par exemple. Si kp ou kq=1, le graphe est une
       union de cliques, et si kp=kq=2 et n=2p=2q, c'est une union de
       cycles.
....
       Ex: gengraph unok 200 100 100 3 2 -visu

....ngon p c x
....
       Triangulation d'un polygone régulier. Plusieurs types de
       triangulations sont produites suivant la valeur des
       paramètres. Si x&4=0, alors la triangulation a 3p sommets (avec
       p>0) et est composée d'un triangle équilatéral central. Si
       x&4=1, alors la triangulation a 4p sommets et est composée d'un
       carré central avec une diagonale. La triangulation est
       symétrique dans chacun des 3 ou 4 croissants délimités par
       chacune des arêtes du polygone central. Si AB est l'une de ces
       arêtes (A avant B dans le sens direct), alors on note C le
       point de l'arc de cercle de A à B à distance c de A. On doit
       avoir c∈[0,p/2]. Les deux bits de poids faible de x définissent
       comment sont construits les triangulations de l'arc AC et
       CB. Tous les points de AC sont connectés à A si x&1=1 et à C
       sinon. Et, tous les points de BC sont connectés à B si x&2=1 et
       à C sinon. Si x&8=1 alors la triangulation est asymétrique. On
       remplace c par p-c pour un arc AB sur deux.
....
       Si x=-1, alors il s'agit d'une autre triangulation. Elle a p
       sommets et est symétrique par rapport à un axe horizontal
       comprenant trois "fan": un depuis le point 0 vers tous ceux de
       [c,n-c], un depuis c vers tous ceux de [0,c], et enfin un
       depuis n-c vers tous ceux de [n-c,n].
....
       Si x=-2, alors il s'agit de la triangulation récursive à 3p
       sommets. Le paramètre c n'a pas de rôle. Chacun des trois arc
       est coupé en deux récursivement. Si p n'est pas une puissance
       de deux, alors le graphe peut ne pas être une triangulation
       complète, mais le graphe reste cependant planaire.
....
       La triangulation qui expérimentalement minimise le stretch
       maximum (voir -check stretch) est obtenue avec "ngon p 𝛼p 3" où
       𝛼 = 231/512 ≃ 45%. Le stretch maximum est environ 1.455 réalisé
       entre les sommets u=p+30% et v=3p-20%.

....behrend p k
....
       Graphe régulier de pk sommets possédant un très grands nombre
       de cycles de longueur k arête-disjoints où p,k ≥ 2. Si p est
       premier, il en possède exactement p*c! = p^{2-o(1)} où c ~
       log(p)/loglog(p).  Son degré est 2c! si k>2 ou c! si k=2. Le
       graphe est défini que p soit premier ou pas. Il est construit à
       partir de k stables S_0,…,S_{k-1} de chacun p sommets. Chacun
       des cycles de longueur k contient exactement un élément de
       chaque S_j qui est le sommet d'indice i+j*x (mod p) dans S_j
       avec i∈[0,p[ et x∈X où X⊂[0,p[ est un ensemble où k entiers
       quelconques ne sont jamais en progression arithmétique. On
       construit X comme l'ensemble de tous les entiers < p/(k-1)
       s'écrivant sur c chiffres distincts pris dans [0,c[ en base
       ck+1 avec c maximum. Donc |X|=c!. Lorsque p est petit, le
       graphe peut ne pas être connexe.
....
       Par exemple, pour p=421 et k=3 on obtient c=3 et X = { 012,
       021, 102, 120, 201, 210 } (nombres écrits en base ck+1=10). On
       vérifie qu'on a bien p > (k-1)*max{X} = 420. Ce graphe et donc
       12-régulier possède p*k = 1263 sommets et p*c! = 5052 triangles
       arête-disjoints car 421 est premier. La table ci-dessous donne
       en fonction de k et du degré souhaité la plus petite valeur de
       p=p(k) possible. Si p est plus petit que p(k), alors le degré
       sera moindre. Lorsque k=2, le degré est c! au lieu de 2*c!.
....
       !!!          2*2!   2*3!    2*4!       2*5!
              degré    4     12      48        240
    	     --------------------------------------
               p(2)    6    106   2,359     62,811
               p(3)   15    421  13,885    549,921
               p(4)   28  1,054  46,003  2,419,831
....
       Ces graphes sont utilisés en "property testing" pour montrer
       qu'il est difficile de déterminer si un graphe dense possède ou
       pas un cycle de longueur k.

....rplg n t
....
       Random Power-Law Graph. Graphe à n sommets où les degrés des
       sommets suivent une loi de puissance d'exposant t. L'espérance
       du degré du sommet i=0…n-1 est w_i=(n/(i+1))^(1/(t-1)). Il
       s'agit d'une distribution particulière d'un Fixed Degree Random
       Graph où la probabilité d'avoir l'arête i-j est
       min{w_i*w_j/S,1} avec S=∑_k w_k. En général il faut prendre le
       paramètre t comme un réel de ]2,3[, la valeur communément
       observée pour le réseau Internet étant t=2.1.

....bdrg n_1 d_1 … n_k d_k .
....
       Bounded Degree Random Graph. Graphe aléatoire dont la
       distribution des degrés des sommets est fixée par les paires
       (nᵢ,dᵢ) signifiant qu'il y a nᵢ sommets de degré au plus
       dᵢ. Ainsi "bdrg n 3 ." génère un graphe sous-cubic aléatoire à
       n sommets, si n est pair. Les sommets sont dupliqués selon leur
       distribution de degré puis un couplage aléatoire détermine les
       arêtes. Les boucles et les arêtes multiples sont supprimer.  Il
       suit que le degré des sommets ne dépasse pas dᵢ. Ils peuvent
       cependant être inférieurs. Le nombre de sommets est n=∑ᵢ nᵢ et
       le nombre d'arêtes au plus m = ½∑ᵢ (nᵢ*dᵢ). Si cette somme
       n'est pas entière, alors le degré d'un des sommets ayant dᵢ>0
       est diminué d'un. (C'est un sommet avec dᵢ>0 avec le plus grand
       i qui est choisi.)

....fdrg n_1 d_1 … n_k d_k .
....
       Fixed Degree Random Graph. Graphe aléatoire assymptotiquement
       uniforme dont les degrés des sommets sont fixées par les paires
       (nᵢ,dᵢ) signifiant qu'il y a nᵢ sommets de degré dᵢ. La suite
       des degrés doit être graphique, à savoir qu'il existe au moins
       un graphe simple ayant ces degrés (sinon une erreur est
       affichée). Ainsi "fdrg n 3 ." génère un graphe cubic aléatoire
       asymptotiquement uniforme, à condition que n soit pair. La
       complexité est en moyenne O(mΔ+Δ⁴) où m=∑ nᵢdᵢ et Δ=max{dᵢ},
       et pour être asymptotiquement uniforme, il faut Δ=o(m^¼) ou
       Δ=o(√n) pour les graphes réguliers (tous les dᵢ égaux ou k=1).

....matching n
....
       Graphe composé de n arêtes indépendantes, c'est-à-dire de n
       copies de K₂.

....load file[:range]
....loadc file[:range]
....
       Graphe défini à partir du fichier "file" ou de l'entrée
       standard si file vaut "-". Si "file" est une famille de
       graphes, alors il est possible d'utiliser la variante
       "file:range" pour préciser l'identifiant du graphe souhaité
       (sinon c'est le premier graphe de la famille qui sera
       considéré). Le graphe (ou la famille) doit être au format
       standard, les sommets numérotés par des entiers positifs. Les
       caractères situés sur une ligne après "//" sont ignorés, ce qui
       permet de mettre des commentaires.
....
       Le temps et l'espace nécessaire au chargement du graphe sont
       linéaires en la taille du fichier (si "file" est une famille de
       graphes, le fichier est entièrement lu).  Cependant, pour la
       génération à proprement parlée du graphe final, qui peut
       comprendre l'option -not par exemple, toutes les arêtes
       potentielles, soit O(n²), sont passées en revue pour être
       testées. La variante "loadc" (pour "load & check") permet une
       génération plus rapide lorsqu'utilisée avec -check (ou les
       alias utilisant -check, comme -maincc par exemple). Elle permet
       de passer directement de l'étape de chargement du graphe à
       l'étape du test de l'algorithme en sautant la phase de
       génération des arêtes. En contre-partie, le graphe n'est pas
       affiché et les options comme -not, -permute, -delv, -dele,
       etc. n'ont plus d'effet. La variante "loadc file" est environ
       20% plus rapide que "load file -fast".
....
       Pour charger un graphe au format dot on peut utiliser le script
       dot2gen.awk en amont, comme dans l'exemple suivant:
....
       !!! nop file.dot | awk -f dot2gen.awk | ./gengraph load -
....
       Le filtre nop de GraphViz, qui est recommandé mais pas
       nécessaire, permet de standardiser le format dot initial. Il
       transforme par exemple les expressions du type "a--{b;c;}" en
       "a--b;a--c;".
....
       Notez que la suite d'options "load file -fast -format dot<type>"
       permet de convertir "file" au format <type> souhaité. Ce graphe
       active l'option -directed si "file" contient au moins un
       arc. Dans ce cas l'option -undirected n'aura pas d'effet.


   GRAPHES ORIENTÉS :

....aqua n c_1 … c_n
....
       Graphe orienté dont les sommets sont les suites de n entiers
       positifs dont la somme fait c_1 et dont le i-ème élément est au
       plus c_i. Ils représentent les façons de répartir une quantité
       c_1 de liquide dans n récipients de capacité c_1 … c_n. Il y
       a un arc u->v s'ils existent i et j tels que v est le résultat
       du versement du récipient c_i vers le récipient c_j.
....
       Ex: gengraph aqua 3 3 2 1 -label 1 -dot filter dot -visu


   GRAPHES COMPOSÉS :

....mesh p q (= grid 2 p q)
....
       Grille 2D de p x q sommets.

....hypercube d (= grid d 2 … 2)
....
       Hypercube de dimension d.

....path n (= grid 1 n)
....
       Chemin à n sommets.

....cycle n (= ring n 1 1)
....
       Cycle à n sommets.

....torus p q (= grid 2 -p -q)
....
       Tore à p x q sommets.

....stable n (= ring n 0)
....empty n (= ring n 0)
....
       Stable à n sommets.

....clique n (= -not ring n 0)
....
       Graphe complet à n sommets.

....bipartite p q (= rpartite 2 p q)
....
       Graphe biparti complet K_{p,q}.

....utility (= rpartite 2 3 3)
....
       Graphe biparti complet K_{3,3} qui doit son nom au problème de
       la connexion planaire de trois maisons à trois stations (eau,
       gaz, électricité). C'est aussi le graphe de Haar H(7).

....octahedron (= antiprism 3)
....
       Octaèdre: graphe 4-régulier planaire à 6 sommets ayant 8 faces
       triangulaires. Il s'agit de deux pyramides dont la base à 4
       sommets est commune. C'est aussi le graphe de Johnson J(4,2).

....d-octahedron d (= -not matching d)
....
       Octaèdre de dimension d: obtenu à partir d'un octaèdre de
       dimension d-1 auquel on ajoute deux sommets universels,
       l'octaèdre de dimension 1 étant composé d'un stable de deux
       sommets.  L'octaèdre classique est obtenu avec d=3, pour d=2 il
       s'agit d'un carré.

....tetrahedron (= -not ring 4 0)
....
       Tétraèdre: pyramide composée de 4 faces triangulaires. C'est
       aussi une clique à 4 sommets.

....hexahedron (= cube)
....
       Hexaèdre: cube composé de 6 faces carrées.

....associahedron (= flip 6)
....
       Associaèdre (3D): graphe planaire cubique à 14 sommets composé
       de 3 faces carrées et 6 faces pentagonales.

....johnson n k (= -not kneser n k k-2)
....
       Graphe de Johnson J(n,k). Les sommets sont tous les
       sous-ensembles à k éléments de [0,n[ (il faut donc 0≤k≤n). Deux
       sommets sont adjacents ssi leurs ensembles correspondant ont
       k-1 éléments en commun. La distance entre deux sommets est la
       distance de Hamming entre les ensembles correspondant. Ils sont
       réguliers de degré k(n-k), de diamètre min{k,n-k}, de
       sommet-connectivité k(n-k). Ils sont aussi distance
       réguliers. J(n,1) est la clique K_n, J(n,2) est le complément
       du graphe de Kneser K(n,2) et le line-graphe de K_n. En fait,
       tout sous-graphe induit de J(n,2) est un line-graphe. J(4,2)
       est l'octaèdre, J(5,2) le complément du graphe de Petersen.

....claw (= rpartite 2 1 3)
....
       Graphe biparti complet K_{1,3}.

....star n (= rpartite 2 1 n)
....
       Arbre (étoile) à n feuilles et de profondeur 1.

....tree n (= arboricity n 1)
....
       Arbre plan enraciné aléatoire uniforme à n sommets. Les sommets
       sont numérotés selon un parcours en profondeur depuis la racine
       et le long de la face extérieure.

....caterpillar n (= grid 1 n-r -star r, r=random()%n)
....
       Arbre à n sommets dont les sommets internes (de degré > 1)
       induisent un chemin. Il est obtenu à partir d'un chemin de
       longueur n-r (où r est un nombre aléatoire entre 0 et n-1) et
       en appliquant l'option -star r. Si l'option -seed est présente,
       (pour intervenir sur la valeur "r"), il est important qu'elle
       figure avant caterpillar.

....outerplanar n (= kpage n 1)
....
       Graphe planaire-extérieur aléatoire connexe à n sommets (plan
       et enraciné). Ils sont en bijection avec les arbres plans
       enracinés dont tous les sommets, sauf ceux de la dernière
       branche, sont bicoloriés. Les sommets sont numérotés le long de
       la face extérieure. C'est aussi une numérotation selon un
       parcours en profondeur depuis la racine de l'arbre bicolorié.
       Il est aussi possible de générer des graphes
       planaires-extérieurs aléatoires Hamiltoniens, donc 2-connexes,
       avec "planar n f -1" ou "polygon n".

....squaregraph n (= planar n 4 4)
....
       Squaregraph aléatoire à n faces. Ce sont des graphes planaires
       2-connexes dont toutes les faces (sauf l'extérieure) sont des
       carrées. De plus, les sommets des faces internes sont de degré
       au moins 4. Ce sont des sous-graphes de quadrangulations et
       donc des 2-pages.

....random n p (= -not ring n 0 -dele 1-p)
....
       Graphe aléatoire à n sommets et dont la probabilité d'avoir une
       arête entre chaque paire de sommets est p. L'option -dele étant
       déjà présente, il n'est pas conseillé de la réutiliser pour ce
       graphe.

....netgraph (= sierpinski 2 3 -not)
....
       Graphe à 6 sommets composé d'un triangle avec un sommet pendant
       à chacun d'eux. C'est le complémentaire du graphe de Hajós. On
       peut aussi le générer en utilisant "fdrg 3 3 3 1 .".

....sunflower n (= cage 2n 2 2 0)
....
       Tournesol à n pétales. C'est un graphe planaire-extérieur à 2n
       sommets composé d'un cycle de longueur n≥3 où chaque arête
       partage le coté d'un triangle. C'est le graphe "gear n" sans le
       sommet central. Pour n=3, c'est le graphe de Hajós.

....gem (= fan 4 1)
....
       Graphe à 5 sommets composé d'un chemin et d'un sommet universel.

....egraph (= comb 3)
....
       Arbre à 6 sommets et 3 feuilles en forme de E.

....tgraph (= banana 1 3)
....fork (= banana 1 3)
....
       Fork Graph, arbre à 5 sommets dont trois feuilles en forme de
       T.

....ygraph (= banana 3 1)
....
       Arbre à 7 sommets composé d'une étoile à trois branches.

....cross (= banana 1 4)
....
       Cross Graph, arbre à six sommets en forme de croix chrétienne.

....knight p q (= chess p q 1 2)
....
       Graphe des déplacements possible du chevalier dans un échiquier
       p q.

....antelope p q (= chess p q 3 4)
....
       Graphe des déplacements possible d'une antilope dans un
       échiquier p q, une antilope étant une pièce hypothétique se
       déplaçant de 3 cases selon un axe et de 4 selon l'autre.

....camel p q (= chess p q 1 3)
....
       Graphe des déplacements possible d'un chameau dans un échiquier
       p q, un chameau étant une pièce hypothétique se déplaçant de 1
       case selon un axe et 3 de selon l'autre.

....giraffe p q (= chess p q 1 4)
....
       Graphe des déplacements possible d'une giraffe dans un
       échiquier p q, une giraffe étant une pièce hypothétique se
       déplaçant de 1 case selon un axe et de 4 selon l'autre.

....zebra p q (= chess p q 2 3)
....
       Graphe des déplacements possible d'un zèbre dans un échiquier p
       q, un zébre étant une pièce hypothétique se déplaçant de 2
       cases selon un axe et de 3 selon l'autre.

....petersen (= kneser 5 2 0)
....
       Graphe de Kneser particulier. Il est cubique et possède 10
       sommets. Il n'est pas Hamiltonien et c'est le plus petit graphe
       dont le nombre de croisements (crossing number) est 2. C'est le
       complément du line-graphe de K₅.

....tietze (= flower_snark 3)
....
       Graphe de Tietze. Il est cubique avec 12 sommets. Il possède un
       chemin Hamiltonien, mais pas de cycle. Il peut être plongé sur
       un ruban de Möbius, a un diamètre et une maille de 3. Il peut
       être obtenu à partir du graphe de Petersen en appliquant une
       opération Y-Delta.

....mobius-kantor (= gpetersen 8 2)
....
       Graphe de Möbius-Kantor. Graphe cubique à 16 sommets de genre
       1. Il est Hamiltonien, de diamètre 4 et de maille 6. C'est
       aussi le graphe de Haar H(69).

....dodecahedron (= gpetersen 10 2)
....
       Dodécaèdre: graphe planaire cubique à 20 sommets. Il possède 30
       arêtes et 12 faces qui sont des pentagones. C'est le dual de
       l'icosaèdre.

....desargues (= gpetersen 10 3)
....
       Graphe de Desargues. Il est cubique à 20 sommets. Il est
       Hamiltonien, de diamètre 5 et de maille 6.

....durer (= gpetersen 6 2)
....
       Graphe de Dürer. Graphe cubique planaire à 12 sommets de
       diamètre 4 et de maille 3. Il peut être vu comme un cube avec
       deux sommets opposés tronqués (remplacés par un cycle de
       longueur 3).

....prism n (= gpetersen n 1)
....
       Prisme, c'est-à-dire le produit cartésien d'un cycle à n
       sommets et d'un chemin à deux sommets. Pour n=3, c'est un
       graphe de Halin et aussi le complémentaire d'un cycle de
       longueur 6, et pour n=4 il s'agit du cube.

....cylinder p q (= grid p -q)
....
       Produit cartésien d'un chemin à p sommets et d'un cycle à q
       sommets. Cela généralise le prisme (prism n = cylinder n 3). Un
       cube est un "cylinder 2 4".

....nauru (= pstar 4)
....
       Graphe de Nauru. C'est un graphe cubique à 24 sommets. Il
       s'agit d'un graphe "permutation star" de dimension 4. C'est
       aussi un graphe de Petersen généralisé P(12,5).

....headwood (= cage 14 2 5 -5)
....
       Graphe de Headwood. C'est un graphe cubique à 14 sommets, de
       maille 6 et de diamètre 3. C'est le plus petit graphe dont le
       nombre de croisements (crossing number) est 3. C'est aussi le
       graphe de Haar H(69).

....franklin (= cage 12 2 5 -5)
....
       Graphe de Franklin. C'est un graphe cubique à 12 sommets, de
       maille 4 et de diamètre 3. C'est aussi le graphe de Haar H(37).

....dyck (= cage 32 4 5 0 13 -13)
....
       Graphe de Dyck. C'est un graphe cubique 3-connexe biparti à 32
       sommets. C'est le seul graphe cubic à 32 sommets à être
       symétrique, c'est-à-dire qui est à la fois arête et sommet
       transitif. Il est aussi torique, c'est-à-dire de genre 1.

....pappus (= cage 18 6 5 7 -7 7 -7 5)
....
       Graphe de Pappus. C'est un graphe cubique à 18 sommets, de
       maille 6 et de diamètre 4.

....mcgee (= cage 24 3 12 7 -7)
....
       Graphe de McGee. C'est un graphe cubique à 24 sommets, de
       maille 7 et de diamètre 4.

....tutte-coexter (= cage 30 6 -7 9 13 -13 -9 7)
....
       Graphe de Tutte-Coexter appelé aussi 8-cage de Tutte. C'est un
       graphe cubique à 30 sommets, de maille 8 et de diamètre 4.
       C'est un graphe de Levi mais surtout un graphe de Moore,
       c'est-à-dire un graphe d-régulier de diamètre k dont le nombre
       de sommets est 1+d*S(d,k) (d impair) ou 2*S(d,k) (d pair) avec
       S(d,k)=∑_{i=0}^(k-1) (d-1)^i.

....gray (= cage 54 6 7 -7 25 -25 13 -13)
....
       Graphe de Gray. C'est un graphe cubique à 54 sommets qui peut
       être vu comme le graphe d'incidence entre les sommets d'une
       grille 3x3x3 et les 27 lignes droites de la grille. Il est
       Hamiltonien, de diamètre 6, de maille 8, et de genre 7. Il est
       arête-transitif et régulier sans être sommet-transitif.

....grotzsch (= mycielski 4)
....
       Graphe de Grötzsch. C'est le plus petit graphe sans triangle de
       nombre chromatique 4. Il possède 11 sommets et 20 arêtes. Comme
       le graphe de Chvátal, il est non-planaire de diamètre 2, de
       maille 4 et Hamiltonien. C'est le graphe de Mycielskian du
       cycle à 5 sommets.

....hajos (= sierpinski 2 3)
....
       Graphe de Hajós. Il est composé de trois triangles deux à deux
       partageant un sommet distinct. On peut le dessiner comme un
       triangle dans un triangle plus grand. Il est planaire et
       possède 6 sommets. C'est un graphe de Sierpinski ou encore le
       complémentaire d'un "sunlet 3", complémentaire du "netgraph",
       un "sunflower 3" ou encore "cage 6 2 2 0".

....house (= -not grid 1 5)
....
       Graphe planaire à 5 sommets en forme de maison. C'est le
       complémentaire d'un chemin à 5 sommets.

....wagner (= ring 8 2 1 4)
....
       Graphe de Wagner appelé aussi graphe W₈, un cycle à 8 sommets
       où les sommets antipodaux sont adjacents. C'est un graphe
       cubique à 8 sommets qui n'est pas planaire mais sans K₅. C'est
       aussi une échelle de Möbius.

....mobius n (= ring n 2 1 n/2)
....
       Échelle de Möbius, graphe cubic à n sommets obtenu à partir
       d'un cycle à n sommets dont les sommets opposés sont
       adjacents. Lorsque n est pair, il s'agit d'un ruban de Möbius,
       c'est-à-dire d'une échelle dont le premier et dernier barreau
       sont recollés en sens opposé. Pour n≤5, il s'agit d'une clique
       à n sommets. Lorsque n≥5, le graphe n'est plus planaire, et
       pour n=8, il s'agit du graphe de Wagner.

....ladder n (= grid 2 2 n)
....
       Graphe échelle à n barreaux, soit une grille à 2 x n sommets.

....cube (= crown 4)
....
       Hypercube de dimension 3, graphe planaire cubic à 8 sommets où
       toutes les faces sont des rectangles.

....diamond (= fan 2 2)
....
       Clique à quatre sommets moins une arête. C'est un graphe
       allumette, c'est-à-dire planaire et distance unitaire.

....gosset (= ggosset 8 2 2 3 6 -1)
....
       Graphe de Gosset. Il est 27-régulier avec 56 sommets et 756
       arêtes, de diamètre, de rayon et de maille 3. Il est
       27-arête-connexe, 27-sommet-connexe et Hamiltonien. C'est
       localement un graphe de Schläfli, c'est-à-dire que pour tout
       sommet le sous-graphe induit par ses voisins est isomorphe au
       graphe de Schläfli, qui est lui-même localement un graphe de
       Clebsch.

....wheel n (=ringarytree 1 0 n 2)
....
       Roue à n rayons. Graphe planaire à n+1 sommets composé d'un
       cycle à n sommets et d'un sommet universel, donc connecté à
       tous les autres.

....web n r (=ringarytree r 1 n 2)
....
       Graphe planaire à 1+n*r sommets composé d'une étoile à n
       branches de longueur r, les sommets de même niveau étant
       connectés par un cycle. Il généralise "wheel n" (r=1).

....binary h (= ringarytree h 2 2 0)
....
       Arbre binaire complet de profondeur h. Il possède 2^(h+1)-1
       sommets et la racine est de degré deux.

....arytree h k r (= ringarytree h k r 0)
....
       Arbre complet de hauteur h où chaque noeud interne à exactement
       k fils, le degré de la racine étant de degré r.

....rbinary n (= rarytree n 2 0)
....rbinaryz n (= rarytree n 2 1)
....
       Arbre binaire plan aléatoire uniforme à n noeuds internes. Il
       possède 2n-1 sommets (2n pour la variante rbinaryz) numérotés
       selon un parcours en profondeur modifié: tous les fils du
       sommet courant sont numérotés avant l'étape de récursivité. La
       racine est de degré 2 (=rbinary) ou 1 (=rbinaryz). Le dessin
       avec dot (-visu) ne respecte pas le plongement de l'arbre.

....tw n k (= ktree n k -dele .5)
....
       Graphe de largeur arborescente au plus k aléatoire à n
       sommets. Il s'agit d'un k-arbre partiel aléatoire dont la
       probabilité d'avoir une arête est 1/2. L'option -dele étant
       déjà présente, il n'est pas conseillé de la réutiliser pour ce
       graphe.

....pw n k (= kpath n k -dele .5)
....
       Graphe de pathwidth au plus k, aléatoire et avec n sommets.

....tadpole n p (= barbell -n 1 p)
....dragon n p (= barbell -n 1 p)
....
       Graphe à n+p sommets composé d'un cycle à n sommets relié à un
       chemin à p sommets.

....lollipop n p (= barbell n p 0)
....
       Graphe "tapette à mouches" (Lollipop Graph) composé d'une
       clique à n sommets reliée à un chemin de longueur p. Il a n+p
       sommets.

....pan n (= barbell -n 1 1)
....
       Graphe à n+1 sommets composé d'un cycle à n sommets et d'un
       seul sommet pendant.

....banner (= barbell -4 1 1)
....
       Graphe à 5 sommets composé d'un carré et d'un sommet pendant.

....paw (= barbell -3 1 1)
....
       Graphe à 4 sommets composé d'un triangle et d'un sommet
       pendant.

....theta0 (=barbell -5 -5 -2)
....
       Graphe Theta_0. C'est un graphe à 7 sommets série-parallèle
       obtenu à partir d'un cycle de longueur 6 et en connectant deux
       sommets antipodaux par un chemin de longueur 2. C'est un graphe
       allumette, c'est-à-dire planaire et distance unitaire.

....td-delaunay n (= thetagone n 3 3 1)
....
       Triangulation de Delaunay utilisant la distance triangulaire
       (TD=Triangular Distance). Il s'agit d'un graphe planaire défini
       à partir d'un ensemble de n points aléatoires du carré [0,1[ ×
       [0,1[. Ce graphe a un étirement de 2 par rapport à la distance
       euclidienne entre deux sommets du graphe. Ce graphe, introduit
       par Chew en 1986, est le même que le graphe "demi-theta_6", qui
       est un "theta-graph" utilisant 3 des 6 cônes. La dissymétrie
       qui peut apparaître entre le bord droit et gauche du dessin est
       lié au fait que chaque sommet n'a qu'un seul axe dirigé vers la
       droite, alors qu'il y en a deux vers la gauche.

....theta n k (= thetagone n 3 k 6/k)
....
       Theta-graphe à k>0 secteurs réguliers défini à partir d'un
       ensemble de n points du carré [0,1[ × [0,1[. Les sommets u et v
       sont adjacents si le projeté de v sur la bissectrice de son
       secteur est le sommet le plus proche de u. Ce graphe n'est pas
       planaire en général (sauf pour k<3), mais c'est un spanner du
       graphe complet euclidien si k≥6.

....dtheta n k (= thetagone n 3 ⎣ k/2⎦ 6/k)
....
       Demi-Theta-graphe à k≥2 secteurs réguliers défini à partir d'un
       ensemble de n points du carré [0,1[ × [0,1[. La définition est
       similaire au Theta-graphe excepté que seul 1 secteur sur 2 est
       considéré. Il faut k pair. Pour k=2, il s'agit d'un arbre, pour
       k=4, le graphe est de faible tree-width pas toujours connexe.
       Pour k=6, ce graphe coïncide avec le graphe td-delaunay.
....
       Ex: gengraph dtheta 500 6 -visu
           gengraph dtheta 500 4 -pos 0 -visu
           gengraph dtheta 500 2 -pos 0 -visu

....yao n k (= thetagone n 0 k 2/k)
....
       Graphe de Yao à k>0 secteurs réguliers défini à partir d'un
       ensemble de n points du carré [0,1[ × [0,1[. Les sommets u et v
       sont adjacents si v est le sommet le plus proche de u (selon la
       distance euclidienne) de son secteur. Ce graphe n'est pas
       planaire en général, mais c'est un spanner du graphe complet
       euclidien. Le résultat est valide seulement si k≥2. En fait, ce
       n'est pas tout à fait le graphe de Yao (voir thetagone).

....percolation a b p (= udg a*b 1 -norm L1 -xy mesh a b -dele 1-p)
....
       Grille de percolation à coordonnées entières (i,j) de [0,a[ ×
       [0,b[ où p représente la probabilité d'existence de chaque
       arête. La différence avec le graphe "mesh a b -dele 1-p" est
       qu'ici le graphe est géométrique, donc dessiné selon une grille
       si l'option -visu est ajoutée.

....hudg n r (= udg n r -norm hyper -xy hyper r)
....
       Graphe géométrique aléatoire hyperbolique sur n points du
       disque unité. Deux sommets sont adjacents si leurs points sont
       à distance hyperbolique ≤ r. [A FINIR]

....point n (= ring n 0 -pos 1)
....
       Graphe géométrique composé de n points du plan mais sans aucune
       arête (voir "stable"). Ce graphe permet de visualiser la
       distribution des points, par défaut uniforme sur [0,1[ × [0,1[,
       mais qui peut être modifiée avec l'option -xy.
....
       Ex: gengraph point 500 -xy seed 3 2.1 -visu
           gengraph point 1000 -xy seed 3 -0.1 -visu
           gengraph point 1000 -xy disk -visu

....star-polygon n (= ring n 1 1 -xy disk)
....
       Polygone "star-shaped" aléatoire à n cotés contenu dans le
       disque unité (voir -xy ratio).

....convex-polygon n (= ring n 1 1 -xy convex)
....
       Polygone convexe aléatoire à n cotés contenu dans le disque
       unité (voir -xy ratio). Voir aussi le graphe polygon n.

....regular n d (= fdrg n d .)
....
       Graphe d-régulier aléatoire à n sommets asymptotiquement
       uniforme. Il faut que nd soit pair. L'algorithme est de
       complexité O(nd²) et pour être asymptotiquement uniforme, il
       faut d=o(√n). On obtient nécessairement un matching si d=1, un
       stable si d=0, un cycle si d=2 et n<6.

....cubic n (= fdrg n 3 .)
....
       Graphe cubic aléatoire à n sommets asymptotiquement
       uniforme. Il faut que n soit pair.

....plrg n t (= bdrg n_1 d_1 … n_k d_k .)
....
       Power-Law Random Graph. Graphe aléatoire à n sommets dont la
       distribution des degrés suit une loi en puissance d'exposant
       t>0, la probabilité qu'un sommet soit de degré i>0 étant
       proportionnelle à 1/i^t. Plus précisément, la distribution est
       la suivante:
....
         • d_1=1, n_1=⎣ exp(𝛼)⎦ + n-s(𝛼)
         • d_i=i, n_i=⎣ exp(𝛼)/i^t⎦ pour 2≤i≤p(𝛼)
....
       où a est un réel minimisant |n-s(𝛼)| avec p(𝛼)=⎣ exp(𝛼/t)⎦ et
       s(𝛼)=∑_{i=1}^p(𝛼) ⎣ exp(𝛼)/i^t⎦. Ce sont les mêmes graphes que
       ceux générés par Brady-Cowen'06 ou ceux étudiés dans par Lu'01.

.....
HISTORIQUE

       v1.2 octobre 2007:
            - première version

       v1.3 octobre 2008:
            - options: -shift, -width
            - correction d'un bug pour les graphes de permutation
	    - accélération du test d'ajacence pour les arbres, de O(n) à O(1),
              grâce à la représentation implicite
	    - nouveau graphes: outerplanar, sat

       v1.4 novembre 2008:
            - format de sortie: matrix, smatrix, list
            - nouveau graphe: kout
            - correction d'un bug dans l'option -width
	    - correction d'un bug dans la combinaison -format/shift/delv

       v1.5 décembre 2008:
            - correction d'un bug dans tree lorsque n=1

       v1.6 décembre 2009:
            - nouveaux graphes: rpartite, bipartite

       v1.7 janvier 2010:
            - nouveaux graphes: icosa, dodeca, rdodeca, cubocta, geo,
	      wheel, cage, headwood, pappus, mcgee, levi, butterfly,
	      hexagon, whexagone, arytree, binary, ktree, tw, kpath,
	      pw, arboricity, wagner, mobius, tutte-coexter, paley
            - nouveau format de sortie: -format dot
	    - nouvelles options: -header, -h, -redirect, -dotpdf
            - correction d'un bug dans kout, et dans tree lorsque n=0
	    - tree devient un cas particulier d'arboricity.
	    - aide en ligne pour les paramètres des graphes.

       v1.8 juillet 2010:
            - nouveaux graphes: chvatal, grotzsch, debruijn, kautz
	      gpstar, pstar, pancake, nauru, star, udg, gpetersen,
              mobius-kantor, desargues, durer, prism, franklin,
	      gabriel, thetagone, td-delaunay, yao, theta, dtheta
            - suppression du graphe geo (remplacé par udg)
            - nouvelles options: -pos, -norm, -label, -dotfilter
	    - nouvelle famille d'options: -xy file/noise/scale/seed
	    - définition plus compacte dodeca (non explicite)
	    - utilisation du générateur random() plutôt que rand().
	    - correction d'un bug dans "-format standard" qui provoquait une erreur.
	    - correction d'un bug dans kneser pour k=0, n=0 ou k>n/2.
	    - nouveaux formats: -format dot<type>, -format xy
	    - suppression de -dotpdf (qui est maintenant: -format dotpdf)
	    - labeling pour: gpetersen, gpstar, pstar, pancake, interval,
	      permutation

       v1.9 août 2010:
            - renome -h en -list
	    - renome -xy file en -xy load
	    - centrage des positions sur le barycentre des graines (-xy seed)
	    - nouvelles options: -star, -visu, -xy round
	    - les graphes peuvent être stockés en mémoire, sous la forme d'une liste
	      d'adjacence grâce à l'option -check.
	    - généralisation de -delv p avec p<0
	    - nouveaux graphes: caterpillar, hajos, hanoi, sierpinski, sunlet, load
	    - labeling pour: hanoi, sierpinski
	    - aide sur toutes les options (nécessitant au moins un paramètre)
              et non plus seulement pour les graphes
	    - nouvelle famille d'options: -vcolor deg/degr/pal
	    - correction d'un bug pour l'aide dans le cas de commande
	      préfixe (ex: pal & paley)

       v2.0 septembre 2010:
	    - nouvelles options: -vcolor degm/list/randg, -xy unique/permutation,
	      -check bfs, -algo iso/sub
	    - l'option -xy round p admet des valeurs négatives pour p.
	    - les options "load file" et "-xy load file" permettent la
              lecture à partir de l'entrée standard en mettant
              file="-", la lecture de famille de graphes, et supporte les commentaires.
	    - les formats list/matrix/smatrix utilisent un espace
	      linéaire O(n+m) contre O(n²) auparavant.
	    - les sommets sur le bord (graphes géométriques) ne sont plus coupés
	      (bounding-box (bb) plus grandes).
	    - nouveaux graphes: kpage, outerplanar n (=kpage n 1), rng, nng
	      fritsch, soifer, gray, hajos (qui avait été définit mais non
	      implémenté !), crown, moser, tietze, flower_snark, markstrom,
	      clebsch, robertson, kittell, rarytree, rbinary, poussin, errera
	    - les graphes de gabriel (et rng,nng) dépendent maintenant de -norm.
	    - "wheel n" a maintenant n+1 sommets, et non plus n.
	    - aide en ligne améliorée avec "?". Ex: gengraph tutte ? / -visu ?
	    - les options -help et ? permettent la recherche d'un mot clef.
	      Ex: gengraph -help planaire / ? arbre
	    - description plus compacte de tutte (et des graphes à partir d'un tableau)
	    - correction d'un bug pour rpartite (qui ne marchait pas)

       v2.1 octobre 2010:
	    - nouvelles options:
	      -check degenerate/gcolor/edge/dfs/ps1/paths/paths2/iso/sub/minor/isub
	      -filter minor/sub/iso/vertex/edge/degenerate/ps1
	      -filter degmax/degmin/deg/gcolor/component/radius/girth/diameter
	      -filter cut-vertex/biconnected/isub/all/minor-inv/isub-inv/sub-inv
            - suppression de -algo iso/sub: l'option -algo est réservée à la mis
	      au point de -check
	    - extension de -label b à b=2 qui force l'affiche des noms
              sous forme d'entiers même avec -permute.
	    - correction d'un bug pour house (qui ne marchait pas)
	    - nouveau graphe: windmill

       v2.2 novembre 2010:
            - gestion des graphes orientés: lecture d'un fichier de
              graphe (ou d'une famille avec arcs et arêtes)
	    - nouvelles options: -(un)directed, -(no)loop, -check twdeg/tw,
	      -filter tw/id/hyperbol/rename
	    - permet l'affichage de la "value" (p) dans l'option -filter
	    - nouveau graphe: aqua
	    - correction du graphe tutte-coexter et suppression du
              graphe levi (qui en fait était le graphe de tutte-coexter).
	    - généralisation de l'option "load" à load:id family

       v2.3 décembre 2010:
            - nouvelles options: -check ps1bis/edge, -filter ps1bis/tw2
	      -filter minus/minus-id/unique/connected/bipartite/forest
	      -check ps1ter
	    - remplacement de atof()/atoi() par strtod()/strtol() qui
	      sont plus standards.
	    - remplacement de LONG_MAX par RAND_MAX dans les
              expressions faisant intervenir random() qui est de type
              long mais qui est toujours dans [0,2^31[, même si
              sizeof(long)>4. Il y avait un bug pour les architectures
              avec sizeof(long)=8.
	    - nouveau graphe: cylinder
	    - suppression de la variante "load:id" au profit de la
              forme plus générale "file:range" valable pour load, -filter, etc.

       v2.4 janvier 2011:
            - correction d'un bug dans -filter minus-id
	    - correction d'un bug dans rpartite (incorrect à partir de r>5 parts)
	    - correction d'un bug dans whexagon (nb de sommets incorrects)
	    - nouvelles options: -check ps1x/girth, -filter ps1c/ps1x
	    - renomage: ps1bis -> ps1b, ps1ter -> ps1c
	    - nouveau graphe: mycielski
	    - la graphe grotzsch est maintenant défini à partir du graphe
	      mycielski (la définition précédante était fausse)
	    - bug détecté: td-delaunay 500 -check gcolor -format no -seed
              7 | grep '>6' qui donne jusqu'à 7 couleurs; le nb de
              couleurs affichées dans -check gcolor est erroné

       v2.5 mars 2011:
	    - nouveaux graphes: line-graph, claw

       v2.6 juin 2011:
	    - amélioration du test -filter ps1: détection de cliques et d'arbres

       v2.7 octobre 2011:
	    - nouvelle option: -check bellman (pour les géométriques seulement)
	    - ajout des champs xpos,ypos à la structure "graph".
	    - nouveaux graphes: linial, linialc, cube, diamond, theta0,

       v2.8 novembre 2011:
	    - nouveaux graphes: ggosset, gosset, rplg, wiener-araya, headwood4
	    - correction d'un bug pour "-xy seed k n" lorsque k=1.
	    - nouvelles options: -check maincc, -maincc (non documentée)

       v2.9 février 2013:
	    - nouveau graphe: frucht, halin
	    - correction d'un bug pour "-check gcolor" qui ne
	      renvoyait pas le nombre correct de couleurs, et qui de
	      plus n'utilisait pas l'heuristique du degré minimum.
	    - correction d'un bug pour "permutation -label 1"

       v3.0 octobre 2013:
	    - nouveaux graphes: rig, barbell, lollipop
	    - généralisation de l'option -filter forest
	    - nouvelles options: -apex, -filter isforest, -filter istree, -filter cycle
	    - correction d'un bug dans -filter vertex
	    - amélioration de l'aide lors d'erreurs de paramètre comme:
              "-filter F vertex" au lieu de "-filter F vertex n"
	    - amélioration de l'option -header

       v3.1 décembre 2013:
	    - nouveaux graphes: bpancake
	    - légère modification des labels des sommets des graphes pancake, 
	      gpstar et pstar
	    - nouvelles options: -xy grid, -xy vsize
	    - modification de la taille des sommets pour dot permettant de tenir
	      compte de -xy scale.

       v3.2 mai 2014:
            - amélioration du test ps1b (ajoût de règle et réduction
	      du nombre d'indéterminées dans graphes des conflits)

       v3.3 juillet 2014:
            - modification importante du code pour -check ps1
            - modification des graphes linial et linialc
	    - nouvelles options: -check kcolor, -vcolor kcolor, -len, -check kcolorsat

       v3.4 février 2015:
            - documentation et mise au point de l'option -maincc
	    - correction d'un bug lors de la combinaison de "load file" et de "-vcolor pal grad"
	    - correction d'un bug dans la fonction SortInit() qui affectait "-vcolor deg"
	    - correction d'un bug avec l'option -label 1 pour certains graphes (outerplanar ...)
            - création du script dot2gen.awk pour convertir le format dot en format standard
	    - nouvelles options: -fast, -caption
	    - introduction du groupement d'arêtes/arcs i-(j k ...) dans le format standard
	      (il reste un bug si ce format est lu depuis l'entrée standard)

       v3.5 mars 2015:
	    - correction d'un bug pour le dégradé de couleurs avec "-vcolor pal grad"
	    - nouvelles options: -check ncc/connected, -check routing scenario [...] cluster
	    - nouvelle variante de graphe: loadc file
	    - nouveaux graphes: octahedron, turan, hexahedron, tetrahedron, deltohedron,
	      trapezohedron, antiprism, flip, associahedron, shuffle
	    - changement de nom (ajoût du suffixe "hedron") pour: isoca, dodeca, cubocta,
	      rdodeca

       v3.6 juin 2015:
            - nouvelles options: -version, -variant, -check info, -xy mesh,
	      -check routing tzrplg, -check routing dcr
	    - nouveaux graphes: percolation, hgraph, cricket, moth, bull, expander
	    - correction de bug dans "-help ?", dans "-check girth" qui donnait -1 pour
	      un cycle de taille n, dans butterfly (mauvais nombre de paramètres)
	    - vérification que le nombre de sommets n'est pas négatif
	      pour plusieurs graphes dont tree, kout, etc. ce qui
	      pouvait provoquer des erreurs
	    - vérification du format d'entrée pour load/loadc
	    - l'option -check maincc n'affiche plus le graphe initial
	    - description du format des familles de graphes dans l'aide
	    - affiche plus de messages d'erreurs: lecture d'un graphe au mauvais
	      format, mauvaise combinaison de d'options, ...

       v3.7 juillet 2015:
            - nouvelle implémentation des mots de Dyck, qui sont
              maintenant uniformes. En conséquence les graphes
              aléatoires tree, arboricity, rarytree, rbinary, outerplanar,
              kpage … sont générés de manière uniformes.
	    - nouveaux graphes: treep (et simplification de halin),
	      ygraph, ringarytree (et simplification de arytree,
	      binary, wheel), netgraph, egraph, rgraph (=fish),
	      généralisation de barbell, tadpole (=dragon), pan,
	      banner, paw, theta0 (redéfinition), fan, gem, chess,
	      knight, antelope, camel, giraffe, zebra, utility,
	      zamfirescu, hatzel, web, bdrg, fdrg, regular, cubic, plrg
	    - nouvelles options: -check radius, -check diameter
	    - correction d'un bug si combinaison d'options -check ncc -format list
	      (mauvaise gestion de la variable CHECK)
	    - introduction de caractères UTF8 mathématiques dans l'aide

       v3.8 novembre 2015:
            - nouveaux graphes: ladder, matching, d-octahedron, johnson
	    - correction d'un bug pour le graphe kpage (introduit à v3.7)

       v3.9 février 2016:
	    - renomage de l'option "-xy scale" en "-xy box"
	    - correction d'un bug dans rbinary (mauvais alias/définition)
	    - correction d'un bug (introduit à v3.6) dans la fonction
              max(double,double) au lieu de fmax() qui affectait certains
              graphes et options géométriques (rng, -xy grid, percolation, ...)
            - correction d'un bug concernant rarytree lorsque b>2
	    - correction d'un bug dans les statistiques pour -check routing
	    - généralisation du graphe rarytree (augmentation du degré de la racine)
	    - nouveaux graphes: herschel, goldner-harary, rbinaryz,
	      point, empty, apollonian, dyck, cross, tgraph, bidiakis, gear,
	      centipede, harborth, sunflower, planar, squaregraph, polygon
	    - modification de certains graphes (cage, ring, gabriel,
              nng, rng) afin que le test d'adjacence ne dépende plus de N
              mais de PARAM[0] par exemple
	    - introduction du format %SEED pour l'option -caption

       v4.0 mars 2016:
	    - nouveaux graphes: pat, star-polygon, convex-polygon
	    - nouvelles options: -check kindepsat, -xy circle, -xy polar,
	      -xy unif, -xy convex, -xy ratio, -xy zero
	    - aide permettant les préfixes comme dans -check routing cluster

       v4.1 avril 2016:
	    - nouvelles options: -xy convex2, -check routing hdlbr/bc
	    - nouveaux graphes: kstar, split, cactus, squashed
	    - suppression de -check paths2 qui donnait toujours comme -check paths
	    - option -check bellman apparaît dans l'aide

       v4.2 mai 2016:
	    - nouvelles options: -check routing agmnt, -format vertex,
	      -check routing hash mix, -xy unif, -xy polygon
	    - amélioration de l'option -check info
	    - prise en compte de -xy ratio dans -xy unif et -xy seed
	    - affichage de l'écart type dans les stats affiché par -check routing
	    - affichage des erreurs sur stderr plutôt que stdout
	    - nouveau graphe: suzuki

       v4.3 juin 2016:
	    - nouveaux graphes: triplex, jaws, starfish
	    - suppression d'un bug (Abort trap: 6) pour -format dot<type>
	    - modification de l'initialisation du générateur aléatoire qui pouvait
              faire comme -seed 0 sur certains systèmes où clock() est toujours nulle
	    - ajoût de variantes pour -check routing cluster
	    - nouvelles options: -xy cycle, -xy disk, -check stretch, -xy surface
	    - modification de l'option -norm (-norm 2 -> -norm L2, etc.)
	    - renomage de -xy polar en -xy disk

       v4.4 juillet 2016:
	    - nouvelle option: -check simplify
	    - nouveaux graphes: schlafli, doily, fork (=tgraph)
	    - suppression d'un bug pour gosset (mauvais paramétrage)

       v4.5 août 2016:
	    - modification du terminateur de séquence pour bdrg et fdrg: -1 devient '.'
	    - refonte du prototype des fonctions d'adjacence: adj(i,j) -> adj(Q)
	    - désactivation des options -apex, -star (et donc de caterpillar)
	    - finalisation de l'implémentation de fdrg (et donc de regular, cubic)
	    - redéfinition des graphes: netgraph, egraph, sunlet, cross, tgraph, ygraph
	    - nouveaux graphes: comb, alkane (et ses nombreuses variantes), banana
	      rlt, circle
	    - l'option -fast est effective pour fdrg et bdrg
	    - suppression d'un bug pour -xy mesh
	    - correction du graphe "interval" qui renvoyait en fait un "circle"
	    - génération d'une aide html à partir du source: gengraph.html
	    - nouvelle option: -norm poly p

       v4.6 septembre 2016:
	    - nouvelle option: -check routing scenario until s
	    - nouveaux graphes: uno, unok, behrend

       v4.7 janvier 2017:
	    - dans -check dfs/bfs, détection d'une source invalide
	    - dans -check dfs, calcule et affichage la hauteur des sommets
	    - ajoût d'un paramètre pour unok

       v4.8 mars 2017:
	    - calcul dans -check stretch du stretch max. minimum
	    - correction d'un bug dans -xy convex qui pouvait produire des ensembles
	      non convexes et des points en dehors du cercle de rayon 1
	    - nouveaux graphes: ngon

       v4.9 juin 2017:
	    - correction d'un bug pour l'option -xy box et amélioration
	      du dessin de la grille pour -xy grid
	    - généralisation de -check bellman aux graphes non valués
            - renomage des options: -dotfilter -> -dot filter,
	      -len -> -dot len, -filter hyperbol -> -filter hyper
	    - nouvelles options: -norm hyper, -xy hyper, -dot scale,
	      -label b pour b=3 et b<0, -check volm
	    - nouveaux graphes: hudg, hyperbolic

       v5.0 août 2017:
	    - nouveau graphe: helm, haar
	    - nouveau format html

       v5.1 novembre 2017:
            - redéfinition (plus simple) du graphe de Turán
### #*/
