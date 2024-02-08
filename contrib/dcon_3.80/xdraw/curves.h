/******************************************************************************
**  NAME      CURVES.H
**  AUTHOR    Sheryl M. Glasser
**
**  DESCRIPTION
**      structures LOOP and CURVE_SET xdraw.c, xinit.c
**
**  Copyright (c) GlassWare 1994.  All rights reserved.
******************************************************************************/

typedef struct {	/*----- Loop structure */
  char ltype;		/* I,N,X,A,V,H */
  char use_sep;		/* if not '1', use nnode[] instead of count */
  int hcount;		/* # variables as header to loop */
  int ih0;		/* index of first header variable (.in uses) */
  int count;		/* # non-header elements in loop */
  long sep;		/* count to next index at this level */    
  int current;		/* temp: current (fixed) value */
  float *limVals;	/* limVals: pointer to alloc'd 6 floats */
  char *labels;		/* individual family member names from .in */
} LOOP;

typedef struct {	/*----- struct in CURVE_SET for var */
  int index;		/* index from draw.in; if 'i', add 0x1000 */
  int hloop;		/* it's a var in a header; index of loop */
  long off0;		/* CC/z, CG/y: offset to 1st word of variable for contour or y-axis  */
  long dstep;		/* CC/z: step to next y, CG/y: step to next point on curve  */
  long dfaml;		/* CC/z: step to next x, CG/y: step to start of next family member */
  float min, max;	/* (note "CC" is type C with -c, "CG" is type C w/o -c) */
  float avg, rms;
  char *label;
  int nitems;
  XTextItem  ml_text[10];
} CVAR;

typedef struct {
  int ncount;
  long nsep;
} NODE;

struct LABEL {
  int color, boxed;
  float x,y;
  char *text;
  struct LABEL *next;	/* We REALLY want (struct LABEL *), but no compile */
  int nitems;
  XTextItem ml_text[10];
  } ;

#define LBL_WHITE -1
#define LBL_ALIGNL -2
#define LBL_ALIGNC -3
#define LBL_ALIGNR -4

#define MAXLOOP 16

/*-------------------------------------------------------------------------------
|	CURVE_SET
-------------------------------------------------------------------------------*/
typedef struct {	/* Info about what's drawn in each window */
  char gtype;		/* N,G,g,C,S,D,B,K */
  char which;		/* '1' = 1st in the window, '2' for 2nd, 3rd etc. */
  Window window;
  int flags;		/* MARKERS,COLORS,ASPECT,FILL, etc */
  int use_log;
  int lstep;		/* index for next step on curve. G: normally 1.  Use -s */
  int lfaml;		/* index for next family.        G: normally 0.  Use -f */
  int fskip;		/* number of curves to skip, use -j. Eg. -j4 = every 4th curve in family */
  int param[MAXLOOP];	/* >=0 ==>fixed value from '-i', -1=step, -2=family */
  //CVAR ix,iy;		/* index, off0, dstep, dfaml, label, min, max, etc */
  //CVAR iz;		/* for contours only, corresponds to lfaml */
  CVAR x_info, y_info;	/* index, off0, dstep, dfaml, label, min, max, etc */
  CVAR z_info;		/* for contours only, corresponds to lfaml */
  char *title;
  char *subtitle;
  struct LABEL *label;
  int mcount;		/* -1:markers all +; 0:cycle; n>0:n no line*/
  float force;		/* force one specific value for contour */
  float forcemin;	/* G: forcemin, forcemax used for ymin, ymax */
  float forcemax;
  float exclmin;
  float exclmax;
  float xcutoff;	/* N: lower limits in x and y for log */
  float ycutoff;
  float vector_max_len;		// if -1, needs initialization
  float vector_num_points;	// # points scanned to get vector_max_len
  float vector_lscale;
  int use_cutoff;

  int icurve[9];
  int ncurve, hzcurve;		// c:#contour lines; G:1st,2nd family members
  int multiDots;

  int f_block, l_block;		// # of first and last block (from 'b' key; zero means  all)
  int i_block;			// current block being drawn in xcontour
  int iqty, iqty2;		// M: qty being displayed (if vector plot, both qtys)
  int itime;			// time to use for offset in buf[] to current timestep
  int itime_rel, itime_abs;
  int i_single, j_single;
  int vDensity;

  /* Window description */
    int view_number; 
    int axis_width,curves_width; /* width of axis and curves */
  } CURVE_SET;

typedef struct {	/*----- struct for each variable */
  int index;		/* index assigned in draw.in; if 'i', add 0x1000 */
  char *name;
  int nitems;
  XTextItem ml_text[10];
} VARIABLE_ID;

typedef struct {
  char  dataType;	// T if topology start, F = func values for timestep
  int   nx, ny, nt;	// current topology's #coordinates, #timesteps
  float xmin, xmax;
  float ymin, ymax;
  int   iTime_rel;	// time relative to this topology (abs time is i_m2)
  long  iTop;		// offset in xx.bin of this topology
  long  iThis;		// offset in xx.bin of this timestep
} M2_ARRAY;

// array m2_array contains pointers to all topologies (header and timesteps)
// in the xx.bin file.  But array buf only contains one topology 
// (top header + assoc timesteps) or one timestep (top header + timestep).
// n_m2 is total # of timesteps
// i_m2 is pointer in m2_array, i_m2_prev is previous pointer

/*-------------------------------------------------------------------------------
| Structures for Multiblocks
-------------------------------------------------------------------------------*/
typedef struct {	/*----- struct for each XY block */
  int mr;		// counts for this topology
  int mz;
  int nt;
  int mr_prv;		// counts for previous topology
  int mz_prv;
  int nt_prv;
  long offset;         /* in buf */
  float xmin;
  float xmax;
  float ymin;
  float ymax;
} BLOCK_XY;

typedef struct {	/*----- struct for each F block */
  float fmin;
  float fmax;
  long offset;
} BLOCK_F;

/*------------------------------------------------------------------------
|	VIEW
------------------------------------------------------------------------*/
typedef struct			/* In xtools.c and event.c */
  {
    Window window;
    Window work_area;
    GC gc;
    float f[4];			/* x1,y1,x2,y2, as read from cursor */
    float f1[4];
    float old_f[4], old_f1[4];
    int clipped;

    /* Window attributes */
    int window_flags;		/* Visibility,MENU,Status,SIZE,POZITION*/
    int x,y;			/* Position*/
    int width,height;		/* Geometry*/
    float xmin,xmax,ymin,ymax;  /* xmin,xmax,ymin,ymax currently used */
    float old_xmin, old_xmax;	/* save for restore from nan eg */
    float old_ymin, old_ymax;
    byte zoom_limit;
  }
VIEW;

/*-----------------------------------------------------------------
CURVE_SET
* type: N=unused, G=Graph, C=Contour, S=Solid, D=Dots on contour, B=Big
  (some of these are obsolete?
------------------------------------------------------------------*/

#define MARKERS 0x1		/* flag bits */
#define COLORS	0x2		/* Combined graphs in window change color */
#define ASPECT	0x4		/* Graphs in window use data aspect ratio */
#define DOTS	0x8
#define POLAR   0x10
#define FORCED	0x20
#define FILL    0x40
#define ANIMATE 0x80
#define REVERSE 0x100
#define SINGLE  0x200
#define LABELF	0x400
#define SPLINE	0x800
#define E_VIS	0x1000		/* Extrema from visible family members only */
#define E_IGNORE 0x2000		/* Extrema: use / ignore -E option if any */
/*DD*/
#define FILES   0x4000		/* Extrema: use / ignore -E option if any */
/*DD*/

#define NNAME 50
#define LNAME 80
#define MAXWINDOW 16
#define NCHAR (LNAME*(NNAME+MAXWINDOW+1))	/* size of namebuf[] */
#define NCS 90					/* size of curveset[] */

#define LOOPBIT 0x1000
#define INDBITS 0x0FFF

/* Window attributes constant */
#define W_VIS 0x1
#define W_MENU 0x2

#define W_POZ 0x4
#define W_GEO 0x8
/*DD*/
/* Files */
#define LLNAME (LNAME*10)

typedef struct { 
  char name[LLNAME];
  int curve_begin;
  int curve_end;
} BINFILE; 

/*DD*/
#define DATA_LIMITS 1.0E-9	/* was 1.0e-10 */
