/******************************************************************************
**  NAME      BINREAD.C
**  AUTHOR    Sheryl M. Glasser
**
**  DESCRIPTION
**      Read in data from .in and .bin, load structures
**
**  Copyright (c) GlassWare 1993.  All rights reserved.
******************************************************************************/
//------ Outline of function binread ----------------------------------------------------

// Call from init() while reading drawxx.in
// open file;
// allocate buf based on SEEK_END (or realloc);
// read in entire file;
// initialize, e.g. determine (from 1st rec) if need to invert bytes;
// initialize loop, e.g. count=rec=0

// scan thru the buffer {
//   fetch a word *p from buffer, invert bytes if needed
//   if count==0, start a new record: reset count, bookkeeping eg update loop counters
//   transfer back into buf
//   if --count==0, finish current record, e.g. fetch info from header, bookkeeping about loops
//--------------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <fcntl.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#define O_BINARY 0

#ifndef O_RDONLY
#define O_RDONLY 0			/* for penning */
#endif

#include "gendefs.h"
#include "curves.h"
#include "xinit.h"

static void write_smaller_bin(Long *ppp);

static void getlimits_mb(void );
static int  update_outerloop(void);
static void update_nodeloop(int);
static int start_new_rec(unsigned Long count, unsigned Long count0,
		 unsigned Long rec, unsigned Long rec0, int nhead);

static void getlimits_qty(CURVE_SET *cp,int first_block,int last_block,
			  float *min, float *max, int iqt, int itime);

extern int ftype;	/* 0=graph, 1=contour, 2=self-determined, 3,4 */
extern int multi_topology;
extern int M_one_timestep;
static int one_timestep = 0;
static Long current_buf_size;
extern int tellMe;

/* multiblocks variables and structures*/
Long offsetbuf = 0;
BLOCK_XY *xylist = NULL;
BLOCK_F *flist = NULL;
Long offsetf = 0;
int  nqty = -1;
int  ntime = -1;
/* end of multiblocks variables and stroucture*/

extern LOOP loop[];		/* Loop structure of the data */
extern int nloop;
extern int outerloop_added;
static int iloop;
static int counting_outer=0, i_outer=0;
int counting_nrz=0;
extern int param[MAXLOOP];
static int use_hcount=0, stringbytes=0;
static int invert_bytes_flag = -999;
int force_invert_bytes = -1;

int debug_mb=0;			/* 1=output to debugmb.dat,0=no output */
int debug_m=0, debug_m_trc=0;

extern float *buf;
float *end_buf;

static unsigned Long bufsize;
extern char bin_file_name[];

extern NODE *nodelist;
extern int nnode, inode, ivar, ncount_equal;
extern int nnode_r, nnode_t;	// number blocks: "rectangular", "triangular"
#define NBLOCK 20

#define FZ (float)0

//------ Multi-Topology stuff

int debug_M2=0;			// type=M, multi_topology=1
int debug_M2_detail=0;
M2_ARRAY m2_array[1000];
int read_entire_topology = 0;	// 0=read 1 timestep, 1=read entire topology (all timesteps)

int n_m2 = 0;		// number of topologies (items in m2_array[])
int i_m2 = -1;		// current topology
int i_m2_prev = -1;
int t_m2_rel = -1;
int t_m2_abs = -1;
int ntime_m2;
int itime_m2;		// current time relative to current topology

static unsigned Long bufsize_m2;
static unsigned Long size_1stRec = 0;
static unsigned Long size_topRec = 0;

/*=============================================================================
**		READ BIN FILES
**===========================================================================*/
/*-----------------------------------------------------------------------------
|	tell_buffer - from event "!", type 6
-----------------------------------------------------------------------------*/
void tell_buffer()
{
  Long i, i1, i2;
  int n;
  
 TB_AGAIN:
  printf("Move focus to this window.  Enter i1,  n (0 0 to exit): ");
  scanf("%ld %d", &i1, &n);
  if (n > 0) {
    for(i=i1, i2=i1+n; i<i2; i++) {
      //printf("%ld. %ld %.10g\n",		// print Long AND flot versions of value
      //	i, *(Long *)(buf+i), *(float *)(buf+i));
      printf("%ld. %.10g\n", i, *(float *)(buf+i));
    }
    goto TB_AGAIN;
  }
}

/*-----------------------------------------------------------------------------
|	tell_array - from event "!", type 7
-----------------------------------------------------------------------------*/
void tell_array()
{
  FILE *f;
  int i, j, n, itime, iqty, iblock;
  LOOP *pl;
  CURVE_SET *cp;
  M2_ARRAY *mp;
  BLOCK_XY *xyp;
  BLOCK_F *fp;
  extern CURVE_SET *cp2, curveset[];
  extern int nwindow;
  extern VIEW view[];

  f = fopen("arrays.out", "w");

  //------ tell Curveset

  n = cp2 - curveset;
  fprintf(f, "*** Curveset: nwindow=%d, #curveset %d\n", nwindow, n);
  fprintf(f, "  gtype win# flags lstep lfaml fskip       off0  dstep  dfaml\n");
  for(i=0, cp=curveset; i<n; i++, cp++) {
    fprintf(f, "%2d. %c %4d %3x %7d %5d %5d", 
	    i, cp->gtype, cp->view_number, cp->flags, cp->lstep, cp->lfaml, cp->fskip);
    if (cp->gtype != 'G') {
      fprintf(f, " %12ld %7d %5d", 
	      cp->z_info.off0 , cp->z_info.dstep, cp->z_info.dfaml);
      fprintf(f, "  %s on %s vs %s", cp->z_info.label, cp->y_info.label, cp->x_info.label);
    }
    else {
      fprintf(f, " %12ld %7d %5d", 
	      cp->y_info.off0, cp->y_info.dstep, cp->y_info.dfaml);
    }
    fprintf(f, "\n");
  }

  //------ tell Loop

  fprintf(f, "\n*** loop[], nloop=%d\n", nloop);
  fprintf(f, "   ltype use_sep hcount count sep     labels\n");
  for(i=0, pl=loop; i<nloop; i++,pl++) {
    fprintf(f, "%2d. %c    %c       %d      %2d %6ld  %s\n", 
	    i, pl->ltype, pl->use_sep, pl->hcount, pl->count, pl->sep, pl->labels);
  }

  //------ tell M2_Array, xylist, flist

  if (n_m2) {
    fprintf(f, "\n*** m2_array[], nm2=%d\n", n_m2);
    fprintf(f, "type, nx, ny, nt, iTime_rel, iTop, iThis\n");
    for(i=0, mp=m2_array; i<n_m2; i++, mp++) 
      fprintf(f, "%2d. %c %d %d %d %d %ld %ld\n",
	      i, mp->dataType, mp->nx, mp->ny, mp->nt, mp->iTime_rel,
	      mp->iTop, mp->iThis);

    //---- xylist[]: 1 entry per block, plus 0th entry is "all blocks"
    
    fprintf(f, "\n*** xylist[], # xy blocks %d\n", nnode);
    fprintf(f, "mr, mz, nt, offset\n");
    for(i=0, xyp=xylist; i<=nnode; i++, xyp++)
      fprintf(f, "%2d. %d %d %d at %ld\n", i, xyp->mr, xyp->mz, xyp->nt, xyp->offset);

    //---- flist[]: 1 entry per time+qty+block, plus "0-th" block is really for all

    fprintf(f, "\n*** flist[], ntime %d, nqty %d, nblocks %d\n", n_m2, nqty, nnode);
    fprintf(f, "itime, iqty, iblock, offset\n");
    for (itime=0; itime<n_m2; itime++) {
      for (iqty=0; iqty<nqty; iqty++) {
	fprintf(f, "---- time %d, qty %d ------------------------\n", itime, iqty);
	for (iblock=0; iblock<=nnode; iblock++) {
	  fp = flist + iblock + iqty*(nnode+1) + (nnode+1)*nqty*itime;
	  fprintf(f, "%d. %ld\n", iblock, fp->offset);
	}
      }
    }
  }

  fclose(f);
  printf("File arrays.out written\n");
}

/*-----------------------------------------------------------------------------
|	init_type_g
|	n = # "outer" loops (parameterization levels ABOVE node level)
-----------------------------------------------------------------------------*/
void init_type_g(int, int);
void init_type_g( int n , int count)
{
  LOOP *q;
  int i;

  nloop = n + 2;
  nloop += outerloop_added;
  for (q=loop, i=0; i<nloop-2; i++,q++)
  {
    q->ltype = 'I'; q->count = 0; q->use_sep=0;   q->limVals = NULL;
  }
  q->ltype = 'N'; q->count = 0;   q->use_sep='1'; q->limVals=NULL; q++;
  q->ltype = 'V'; q->count=count; q->use_sep='1'; q->limVals=NULL;
  nodelist = (NODE *)malloc( NBLOCK * sizeof(NODE) * n);
  nodelist->ncount = (ftype == 0) ? 1 : 0;
  stringbytes = use_hcount = i_outer = 0;
  counting_outer = 1;
  inode = nloop - 2;
  ivar  = nloop - 1;
  //printf("INIT_TYPE_G, ltype %c, nloop %d, count %d = %d\n", 
  //q->ltype, nloop, count, q->count );
}

#ifdef DEAD_CODE
/*-----------------------------------------------------------------------------
|	initloops_unused -- after read 1st rec
|	* type M2: nloop=4, ltypes I,X,X,V
-----------------------------------------------------------------------------*/

static int initloops_unused(Long *p, Long count)
{
  LOOP *q;
  int nhead, i;

  q = loop;	// loop, nloop etc are all globals
  nloop = 3;
  iloop = 0;
  ivar = -1;
  nnode = 0;

  if (ftype == 0 || ftype == 4)
    {
      init_type_g(ftype == 0 ? 1 : (int)*p, (int)count);
      nhead = (ftype == 0) ? 0 : 1;
    }

  else if (ftype == 1 || ftype == 5 || ftype == 6 )
    {
      if (count==3 && *(p+2)==0 || ftype ==6 )	/* flag for time steps */
        {
	  q->ltype = 'I';			// ltype I
	  q->use_sep = '1';
	  q->count = 0;
	  q->limVals = NULL;
	  q++;
	}

      if (ftype == 6) {
	ntime = 1;					// default values
	nqty =  1;
	nnode_t = 0;
	nnode_r = 1;

	if  (count >= 1) nnode_r = (int)(*p++);		// optionally in 1st rec
	if  (count >= 2) nnode_t = (int)(*p++);
	if  (count >= 3) nqty =    (int)(*p++);
	if  (count == 4) ntime =   (int)(*p++);

	nnode = nnode_r + nnode_t;	   	
	nhead = 2*nnode_r + 3*nnode_t;

	if (xylist==NULL) 
	  xylist = (BLOCK_XY *) malloc((nnode+1)* sizeof (BLOCK_XY));
	if (debug_m || debug_M2) 
	  printf("Initloops, nnode_r=%d, nnode_t=%d, nqty=%d, nnode=%d, xylist alloc for %d\n",
		 nnode_r, nnode_t, nqty,  nnode, (nnode+1));
      }

      for(i=0; i<2; i++,q++) {				// ltype X: X, Y info (R, Z info)
	q->ltype = 'X';
	q->use_sep = '1';
	if (ftype != 6) {
	  q->count = (int)(*p)+1;
	  nhead = 1;
	  if (i==1) p++;
	}
	else {
	  q->count = 0;
	  if (q->limVals==NULL) q->limVals = (float *)malloc(6 * sizeof(float));
	}
      }

      q->ltype = 'V'; 					// ltype V
      q->count = 1; q->use_sep='1'; q->limVals = NULL;

      if(ftype == 5) {
	nhead=2;
	counting_nrz = (q-1)->count * (q-2)->count;
      }

      ivar = q - loop; 
      nloop = ivar + 1;

      if (nloop == 4)
	counting_outer = (q-2)->count * (q-1)->count * q->count;
    }

  else {				// ftype=2,3 (Type I,H): version, nloop, stringbytes*/
    nloop = nhead = (int)*(p+1);
    if (*p==0 && count==3)
      stringbytes = use_hcount = 0;

    else {
      stringbytes = (int)*(p+2); use_hcount = 1;
    }
  }

  return(nhead);
}
#endif		//... of DEAD_CODE for initloops_unused

/*-----------------------------------------------------------------------------
|	initloops
|	* call just after finished reading first record
|	* p0 points to buf[0]
|	* type M2: nloop=4, ltypes I,X,X,V
|	* ditto for type C
-----------------------------------------------------------------------------*/

static int initloops(Long *p, Long count)
{
  LOOP *q, *qr, *qz;
  int nhead;

  q = loop;	// loop, nloop etc are all globals
  nloop = 3;
  iloop = 0;
  ivar = -1;
  nnode = 0;

  if (ftype == 0 || ftype == 4)
    {
      init_type_g(ftype == 0 ? 1 : (int)*p, (int)count);
      nhead = (ftype == 0) ? 0 : 1;
    }

  else if (ftype == 1 || ftype == 5 || ftype == 6 )
    {
      if (count==3 && *(p+2)==0 || ftype == 6 )	/* flag for time steps */
        {
	  q->ltype = 'I';			// enter one ltype = I for eg time
	  q->use_sep = '1';
	  q->count = 0;
	  q->limVals = NULL;
	  q++;
	}

      qz = q; qr = q+1; //q += 2;		// enter two ltype = X for r, z
      qr->ltype = qz->ltype = 'X';
      qr->use_sep  = qz->use_sep  = '1';
      if (ftype == 5 || ftype == 1 )
	{
	  qr->count = (int)(*p++) + 1;		// here's gridsize r, z
	  qz->count = (int)(*p) + 1;
	  nhead=1;
	}

      else     /* ftype=6 - Multi blocks */
	{
	  ntime = 1;
	  nqty =  1;
	  nnode_t = 0;
	  nnode_r = 1;

	  if  (count >= 1) nnode_r = (int)(*p++);
	  if  (count >= 2) nnode_t = (int)(*p++);
	  if  (count >= 3) nqty =    (int)(*p++);
	  if  (count == 4) ntime =   (int)(*p++);

	  nnode = nnode_r + nnode_t;	   	
	  nhead = 2*nnode_r + 3*nnode_t;
	  qr->count = 0;
	  qz->count = 0;
	  if (xylist==NULL) 
	    xylist = (BLOCK_XY *) malloc((nnode+1)* sizeof (BLOCK_XY));
	  if (debug_m || debug_M2) 
	    printf("Initloops, nnode_r=%d, nnode_t=%d, nqty=%d, nnode=%d, xylist alloc for %d\n",
		   nnode_r, nnode_t, nqty,  nnode, (nnode+1));
	 }

      if (qr->limVals==NULL) qr->limVals = (float *)malloc(6 * sizeof(float));
      if (qz->limVals==NULL) qz->limVals = (float *)malloc(6 * sizeof(float));

      q += 2;				// enter one ltype = V
      q->ltype = 'V'; q->count = 1; 
      q->use_sep='1'; q->limVals = NULL;

      if(ftype == 5)
	{
          nhead=2;
          counting_nrz = qr->count * qz->count;
        }

      ivar = q - loop; 
      nloop = ivar + 1;
      if (nloop == 4) 
	counting_outer = qr->count * qz->count * q->count;
    }

  else				// ftype=2,3 (Type I,H): version, nloop, stringbytes*/
    {
      nloop = nhead = (int)*(p+1);
      if (*p==0 && count==3)
	stringbytes = use_hcount = 0;

      else {
	  stringbytes = (int)*(p+2); use_hcount = 1;
	}
    }

  return(nhead);
}

/*-----------------------------------------------------------------------------
|	addloop -- end of header record, enter loop into structure
-----------------------------------------------------------------------------*/
static void addloop(Long *p)
{
  LOOP *q;
  float *f, *f1;
  int i, n;

  if (ftype == 1)			// type C
    {
      f = (float *)p;
      q = loop + nloop - 2;		// loop[2] limits (ltype = X)
      f1 = q->limVals;			/* rmin, rmax, dr_grid for inner loop */
      *f1++ = *f++;
      *f1++ = *f++;
      *f1   = (*(f1-1) - *(f1-2)) / (float)(q->count-1);

      q = loop + nloop - 3;		// loop[1] limits (ltype = X)
      f1 = q->limVals;			/* zmin, zmax dz_grid for outer loop */
      *f1++ = *f++;
      *f1++ = *f++;
      *f1   = (*(f1-1) - *(f1-2)) / (float)(q->count-1);
    }

  if (ftype == 5)			// type N
    {
      for(i=2; i<=3; i++ )
	{
         q=loop+nloop-i;
         f=q->limVals;                   /*rmin,rmax=0,Npoints-inner loop*/
         *f++ = FZ;                    /*zmin,rmin=0,Npoints-outer loop*/
         *f++ = (float)q->count;
         *f= (float)1;                 /* r_step,z_step=1.0 */
        }
     }					/*NN-end-*/

  else if (ftype == 4)			// type L
    {
      for(i=0,q=loop; i<inode; i++,q++)
      {
	q->count = (int)(*p++);
	if (i>0) counting_outer *= q->count;
	else q->count = 0;			/* ALWAYS count outer loop */
      }
    }

  if (ftype != 2) return;
  n=0;
  q = loop + iloop;
  q->ltype = invert_bytes_flag ?			/* 1 char, but must take up 4 bytes */
             *(char *)p : *((char *)p +3);
  p++;

  q->use_sep = '1';
  q->count = (int)(*p++);		/* (type 'N' must include count=0) */
  q->hcount = q->ih0 = 0;

  i_outer = 0;
  if (iloop > 0) q->ih0 = (q-1)->ih0 + (q-1)->hcount;
  else if (iloop==0 && q->count==0) counting_outer=1;
  if (use_hcount) q->hcount = (int)(*p++);

  q->limVals = NULL;
  f1 = (float *)p;

  if (q->ltype=='X')
    {
      f = q->limVals = (float *)malloc(3*sizeof(float));
      n = 2;
      *f++ = *f1++;		/* x1 */
      *f++ = *f1++;		/* x2 */
      *f   = *f1;		/* dx */
      if (*f==FZ && q->count > 1)
	*f = (*(f-1) - *(f-2)) / (q->count-1);
    }
  else if (q->ltype=='A')
    {
      f = q->limVals = (float *)malloc(q->count*sizeof(float));
      n = q->count;
      for(i=0; i<n; i++) *f++ = *f1++;
    }
  else if (q->ltype=='N' && iloop>0) (q-1)->use_sep = 0;
  else if (q->ltype=='N') d_abort("Inappropriate level for Node loop",0,0);
  else if (q->ltype=='V' && ivar<0) ivar = iloop;
  else if (q->ltype=='V') d_abort("Loop contains too many type V",0,0);
  iloop++;
}

/*-----------------------------------------------------------------------------
|	addloop_mb -- mx,my  enter loop into structure for MultiBlocks
-----------------------------------------------------------------------------*/
static void addloop_mb(Long *p,int iblock)
{
  BLOCK_XY *xy;
  int i1,i2;
  int i;
  LOOP *q,*qr,*qz;
  int ivar, nloop, iTop;
  float *f;

  xy = xylist + iblock;			// xylist[iblock] includes ptr into buf, xmin..ymax
  if (multi_topology==2) iTop = *p++;

  i1 = xy->mr = (int)(*p++) + 1;	// # r_blocks
  i2 = xy->mz = (int)(*p++) + 1;	// # z_blocks

  if (multi_topology==2) {
    xy->nt = *p++;			// # timesteps for multi_topology
    xy->mr_prv = *p++;			// previous mr, mz, nt
    xy->mz_prv = *p++;
    xy->nt_prv = *p++;
    if (debug_M2_detail) {
      printf("   addloop_mb, nt=%d, prev=%d, %d, %d\n", xy->nt, xy->mr_prv, xy->mz_prv, xy->nt_prv);
      printf("   addloop_mb, mr=%d, mz=%d\n", xy->mr, xy->mz);
      printf("   addloop_mb, buf size = %ld\n", sizeof(buf));
    }
  }
  //i2 = xy->mz = (iblock<= nnode_r)?(int)(*p++)+1:(int)(*p++);

  xylist->mr += i1;
  xylist->mz += i2;
  xy->offset = offsetbuf;

  if (debug_m) printf("   addloop_mb, block=%d, nnode_r=%d, offset in buf after x,y=%ld",
		       iblock, nnode_r, offsetbuf);
     
  if ( iblock <= nnode_r ) offsetbuf += 2*i1*i2;	// update to skip over coordinate part
  else
    { 
      xy->mz--;
      i2--;
      offsetbuf += 2*i1 + 3*i2;	
    }

  if (debug_m) {
    printf(" --> %ld for %d,%d\n", offsetbuf, xylist->mr, xylist->mz);
    //printf(" --> %ld\n", offsetbuf);
    //printf("   addloop_mb, xylist, xy at %ld %ld; %d, %d\n", 
    //xylist, xy, xylist->mr, xylist->mz);
  }
}

/*-----------------------------------------------------------------------------
|	getminmax - find min and max for x or y or f
-----------------------------------------------------------------------------*/
static void getminmax( float *f, int n, float *min, float *max, char *sText)
{
  int i;
  float fmin,fmax;

  fmin = *f;
  fmax = *f;
  if (debug_M2_detail)printf("GETMINMAX %s, off=%6ld: ", sText, (Long)(f-buf));
  //printf(" 0. f = %.10g\n", *f);

  for (i=1; i<n; i++,f++)
    { 
      //if (*f > fmax) printf(" %d. f = %.10g \n",i,*f);
      if ( *f < fmin ) fmin = *f;
      if ( *f > fmax ) fmax = *f;
     }
   *min = fmin;
   *max = fmax;
   if (debug_M2_detail) printf("%11.4e to %11.4e\n", fmin, fmax);
}

/*-----------------------------------------------------------------------------
|	getlimits_mb
-----------------------------------------------------------------------------*/
static void getlimits_mb( )
{
  float *f;
  float fmin, fmax;
  float xmin_bl, xmax_bl, ymin_bl, ymax_bl;
  
  BLOCK_XY *xy;
  BLOCK_F  *flistp, *ftmp;
  int i,i1,i2,irz;
  int itime,iqt;
  Long off;
  /* debug */
  FILE *fd;
  float *xb,*yb,*p,*p0,*xb0,*yb0;
  int iz,ir;
  char sText[100];

  f=buf;
  if (debug_mb)
    {
      fd=fopen("debugmb.dat","wt");
      printf("**Output to debugmb.dat for MB contour plot\n");
      fprintf(fd,"**getlimits_mb: debug MB Contour Plot \n");
      fprintf(fd,"ntime = %d, nqty = %d, #block == nnode = %d \n\n", ntime, nqty, nnode);
    }

  /*------ Scan all blocks, calculate limits */

  for (i=1; i<=nnode ; i++)
    {
      xy = xylist + i;
      irz=(i <= nnode_r) ? xy->mr * xy->mz : xy->mr;
      f = buf + xy->offset;
      if (debug_M2_detail) printf("GETLIMITS_MB buf=%lx, off=%lx, value = %ld = %.10g\n",
				  buf, xy->offset, *buf, *buf);

      //---- Find xmin, xmax

      getminmax(f, irz, &fmin, &fmax, "X");	/* xmin, xmax for block i */
      xy->xmin = fmin;
      xy->xmax = fmax;
     
      if (debug_mb)
	{
	  fprintf(fd,"**Block %3d", i);
	  if ( i <= nnode_r ) fprintf ( fd," (R-block)\n");
	  else   fprintf ( fd," (T-block)\n");
	  fprintf(fd," Mx = %d, My = %d, Mx * My = %d \n",xy->mr,
		  (i <=nnode_r)?xy->mz:1,irz);
	  fprintf(fd," Min x = %f Max x = %f \n",fmin, fmax);
	}

      if (i == 1) { xmin_bl = fmin; xmax_bl = fmax; }

      if ( fmin < xmin_bl ) xmin_bl = fmin;
      if ( fmax > xmax_bl ) xmax_bl = fmax; 

      //---- Find ymin, ymax
      
      f = addfloat(f,irz);      /*      f+=irz;*/
      getminmax(f, irz, &fmin, &fmax, "Y");
      xy->ymin = fmin;
      xy->ymax = fmax;

      if (debug_mb)
	fprintf(fd," Min y = %f Max y = %f \n",fmin, fmax);

      if (i == 1) { ymin_bl = fmin; ymax_bl = fmax; }

      if ( fmin < ymin_bl ) ymin_bl=fmin;
      if ( fmax > ymax_bl ) ymax_bl=fmax;

      /*------ Tblocks - for each triangle, tell vertices */

      if( (i>nnode_r) && debug_mb)
	{
	  i1=xy->mr;
	  i2=xy->mz;
	  f = addfloat(f,irz);
	  fprintf (fd," mvert = %d, mcells = %d\n", i1,i2);
	  for (iz=1,p0=f,xb0=buf+xy->offset,yb0=xb0+xy->mr;
	       iz<=i2; iz++,f+=1)
	    {
	      //fprintf (fd,"\n triangle =%d,: vertixes: %f %f %f \n",
	      //iz,(*f),(*(f+1)),(*(f+2)) );
	      fprintf (fd,"\n triangle =%d,: vertixes: %d %d %d \n",
		       iz,*(int *)f,*(int *)(f+i2), *(int *)(f+2*i2) );
	    }
	}

      /* f+=irz;*/
    }

  /*------ save XY limits in xylist[0] */

  xylist->xmin = xmin_bl; 
  xylist->xmax = xmax_bl;
  xylist->ymin = ymin_bl;
  xylist->ymax = ymax_bl;

  if (debug_mb)
    {
      fprintf(fd," Min all x = %f Max all x = %f \n",xmin_bl, xmax_bl);
      fprintf(fd," Min all y = %f Max all y = %f \n",ymin_bl, ymax_bl);
    }

  /*------ Find fmin, fmax: Scan time, quantities */

  f = addfloat(buf,offsetf);  /* f= buf + offsetf;*/
  off = 0;

  for (itime=0; itime < ntime ; itime++ )
    { 
        for (i=1; i<=nnode ; i++   )
	  {
           xy=xylist+i;
	   irz=(i <= nnode_r)?xy->mr * xy->mz : xy->mr;
	 /*irz=xy->mr * xy->mz; */

           for (iqt=0; iqt<nqty ; iqt++ )
	     {
	       /*---- get fmin, fmax here, set in flist */

	       sprintf(sText, "F (qty=%d, t=%d)", iqt, itime);
	       getminmax(f, irz, &fmin, &fmax, sText);
	       ftmp = flist + i + iqt*(nnode+1) + (nnode+1)*nqty*itime;
	       ftmp->fmin = fmin;
	       ftmp->fmax = fmax;
	       ftmp->offset = offsetf+off;
	       //printf("GETLIMITS_MB %d, %d, %d: FMAX=%.10g from %lx\n", 
	       //itime, iqt, irz, fmax, f);

	       /*------ Debug print values here */

	       /* Temporary    
		  printf(" Time= %d  Qty =%d \n", itime, iqt );
		  printf(" Block %d \n Min f = %f Max f = %f \n",i,fmin, fmax);
	       */

	       if(debug_mb) {
		 fprintf(fd,"\n**Scanning time=%d, block=%d, qty =%d", itime,i, iqt );

		 if ( i <= nnode_r ) {
		   fprintf ( fd," (R-block)\n");
		   i1 = xy->mr;
		   i2 = xy->mz;
		 }

		 else {
		   fprintf ( fd," (T-block)\n");
		   i1 = xy->mr;
		   i2 =1;
		 }

		 fprintf(fd," Min f = %f Max f = %f \n",fmin, fmax);
		 fprintf(fd," x_off = %d, y_off = %d, f_off =%d \n",
			  xy->offset,xy->offset+i1*i2,offsetf+off);
		 fprintf(fd, "   i    j          x(i,j)          y(i,j)          f(i,j)\n");
		 fprintf(fd, "   -    -          ------          ------          ------\n");

		 for (iz=0,p0=f,xb0=buf+xy->offset,yb0=buf+xy->offset+i1*i2;
		      iz<i2; iz++,p0+=xy->mr,xb0+=xy->mr, yb0+=xy->mr)
		   for (ir=0,p=p0,xb=xb0,yb=yb0; ir<xy->mr; ir++,p+=1,xb+=1,yb+=1)
		     {
		       fprintf (fd,"%4d %4d %15.5g %15.5g %15.5g\n",
				ir,iz,*xb,*yb,*p);
		     }
	       }/* end of debug_mb*/

	       /*------ Update pointer into buf (f = f + irz) */

	       off += irz;
	       f = addfloat(f,irz);
	       
	     } /* end of quantities loop */
	  }/* end of block loop*/
             
    }  /* end of time loop        */

  /*------ get xmin_bl, xmax_bl, fmin, fmax */

  if (debug_mb) fprintf(fd, "\n");

  for (itime=0; itime < ntime ; itime++ )
    for (iqt=0; iqt<nqty ; iqt++ )
      {
        flistp = flist + (nnode+1)*iqt + (nnode+1)*nqty*itime;
        ftmp = flistp + 1;
        for (i=1; i<=nnode ; i++   )
	  {
	    if( i==1)
	      {
		xmin_bl = ftmp->fmin;
		xmax_bl = ftmp->fmax;
	      }
	    else
	      {
		fmin = ftmp->fmin ;
		fmax = ftmp->fmax ;
		if (fmin < xmin_bl) xmin_bl = fmin;
		if (fmax > xmax_bl) xmax_bl = fmax;
	      }
	    ftmp += 1;
	  }			/* ...end of block loop*/

	flistp->fmin = xmin_bl;
	flistp->fmax = xmax_bl;

	if(debug_mb)
	  fprintf(fd," Time = %d, Iqt %d,  Min f = %f Max f = %f \n",
		  itime, iqt, flistp->fmin, flistp->fmax);
      }				/* ...end of quantities loop */
				/* ...end of time loop */
  if (debug_mb) { fprintf(fd, "**done getlimits_mb\n");  fclose(fd); }
}

/*-----------------------------------------------------------------------------
|	getlimits_qty
-----------------------------------------------------------------------------*/
void getlimits_qty(CURVE_SET *cp,int first_block,int last_block,
			  float *min,
			  float *max,
			  int iqt,
			  int itime)
{
  float *f;
  float fmin,fmax;
  
  BLOCK_F  *fl,*ftmp;
  int i;
  if (debug_mb)
    {
      /*fd=fopen("debugmb.dat","wt");*/
      printf("debug getlimits_qty  \n\n");
      printf("Quantaty : %d \nTime: %d \nBlocks: %d - %d \n",iqt,itime,
	     first_block,last_block);
    }
  
  fl = flist+(nnode+1)*iqt+(nnode+1)*nqty*itime;
  if( !first_block || (first_block >last_block) ) /* all blocks*/
    {
      *min=fl->fmin;
      *max=fl->fmax;
      return ;
    };

  ftmp=fl+first_block;
  for (i=first_block; (i<=nnode) || (i<= last_block) ; i++   )
    {
      if( i== first_block)
	{
	  fmin=ftmp->fmin;
	  fmax=ftmp->fmax;
	}
      else
	{
	  if (fmin > ftmp->fmin) fmin= ftmp->fmin ;
	  if (fmax < ftmp->fmax) fmax= ftmp->fmax ;
	}
      ftmp+=1;
    } /* end of block loop*/
  
  /*	     fl->fmin = xmin_bl;
             fl->fmax = xmax_bl;
	     if(debug_mb)
	     fprintf(fd," Time = %d Iqt %d \n Min f = %f Max f = %f \n",itime,
	     iqt,fl->fmin,fl->fmax);
  */
  printf(" Time = %d Iqt %d \n Min f = %f Max f = %f \n",itime,
	  iqt,fmin,fmax);
} 



/*-----------------------------------------------------------------------------
|	update_outerloop -- for ftype == 0 or 4, null record
-----------------------------------------------------------------------------*/
static int update_outerloop()
{
  LOOP *nq;
  int retval;
  i_outer++;
  if (i_outer!=counting_outer) retval=0;
  else { i_outer = 0; retval = 1; }
  if (inode==-1) return(retval);	/* see if exists node loop */
  nq = loop + inode;
  nnode++;

  /*printf("update_outerloop: nnode=%d, inode =%d\n", nnode,inode); */

  if (nnode==1) nq->count = nodelist->ncount;
  if ((nnode % NBLOCK) == 0)
    nodelist = (NODE *)realloc(nodelist, (nnode+NBLOCK)*sizeof(NODE));

  (nodelist + nnode)->ncount = 0;
  return(retval);
}

/*-----------------------------------------------------------------------------
|	update_nodeloop - for ftype=0.  Update # nodes, at end of record
-----------------------------------------------------------------------------*/
static void update_nodeloop( int increase)
{
  //printf (" Update_nodeloop: nnode =%d , ncount =%d , increase = %d\n ",
  //         nnode, (nodelist+nnode)->ncount,increase);
  if (increase) (nodelist + nnode)->ncount += 1;
  else		(nodelist + nnode)->ncount -= 1;
}

/*-----------------------------------------------------------------------------
|	read_m2_array
|	* read data into array m2_array (1 per topology)
|	* this array is loaded before reading in the data file
|	* n_m2 is number of topologies OR number of timesteps
-----------------------------------------------------------------------------*/
int read_m2_array(int fd, unsigned Long count0)
{
  M2_ARRAY *m2p;
  Long buff[20], bval;
  Long *p, iNext, iTop;
  Long skip_coord, timestep_bytes, coord_bytes, total_bytes;
  unsigned Long thisSize;
  int nt, nx, ny, i, nTop, invert_bytes_f;
  int nx_prev, ny_prev, nt_prev;
  int m2_verbose;

  m2_verbose = 0;

  //------ read 1st record: (3 data +2 delim)*4 bytes

  read(fd, (unsigned char *)buff, 20L * sizeof(float));

  invert_bytes_f = ((buff[0] & 0xffff0000) == 0) ? 0 : 1;
  if (invert_bytes_f) {
    for (i=0; i<5; i++) buff[i] = invert_bytes(buff+i, &bval);
  }
  
  count0 = buff[0] / sizeof(int);		// buff[0]=12, count0=3
  size_1stRec = buff[0] + 2L * sizeof(int);	// 3 data words: 12+8 = 20 bytes
  nnode_r = buff[1];
  nqty = buff[3];				// nb variables

  iTop = size_1stRec;				// offset to the next topology block
  ntime_m2 = 0;
  bufsize_m2 = 0;

  size_topRec = (multi_topology==0) ? 4L : 9L ;	// includes rec delimiters: 2+2 or 7+2
  size_topRec *= sizeof(int);			// 16 (M) or 36 (M2)

  if (debug_M2_detail || m2_verbose) 
    printf("read_m2_array, count0=%ld, 1st offset=%ld, nqty=%ld, size_topRec=%ld\n", count0, iTop, nqty, size_topRec);

  //------ scan all topologies: 1 rec in m2_array for each timestep

  p = (Long *)buff;
  for(nTop=0; ; nTop++) {			// scan over topologies

    //------ Seek (offset from beginning), read

    lseek(fd, iTop, SEEK_SET);			// seek next topology record (skip iTop bytes)
    read(fd, p, size_topRec);			// read 16 or 36 bytes (4 or 9 words)
    if (invert_bytes_f) {
      for (i=0; i<9; i++) buff[i] = invert_bytes(buff+i, &bval);	// (only use i=0..4)
    }

    if ((debug_M2_detail || m2_verbose)) {
      printf("Topology %d. Skipped %ld bytes.  Buff = ", nTop, iTop);
      for (i=0; i<9; i++) printf(" %ld", buff[i]);
      printf("\n");
    }

    //------ Type M2, nx = 0 means last topology, break out

    if (multi_topology > 0) {			// type M2 test last topology
      if (*(p+2) == 0) break;			// last topology? (has nx = ny = nt = 0)
      if (*p != 28L) {
	printf("Error reading topologies\n");
	return 1;
      }
    }

    //------ Get nx, ny, nt from buff; coord_bytes, timestep_bytes

    if (multi_topology > 0 ) {
      nx = *(p+2) + 1;				// nx, ny, nt for block (in buff[2] etc)
      ny = *(p+3) + 1;
      nt = *(p+4);				// #timesteps in this topology
      coord_bytes    = (2L + nx*ny*2L) * sizeof(float);
      timestep_bytes = (3L + nqty*(2L + nx * ny)) * sizeof(float);
      skip_coord     = size_topRec + coord_bytes;
    }
    else {
      nx = *(p+1) + 1;				// nx, ny, nt for block
      ny = *(p+2) + 1;
      coord_bytes    = nnode_r * (4L + 2L + nx*ny*2L) * sizeof(float);
      timestep_bytes = nqty * nnode_r * (2L + nx*ny) * sizeof(float);
      skip_coord     = coord_bytes;
      total_bytes = lseek(fd, 0L, SEEK_END);
      nt = (total_bytes - coord_bytes - size_1stRec) / timestep_bytes;
    }

    if ((debug_M2_detail || m2_verbose)) // && iNext==iTop)
      printf("Topology %d. nt=%d, coord_bytes=%d, timestep_bytes=%d\n", nTop, nt, coord_bytes, timestep_bytes);

    thisSize = size_1stRec + skip_coord + timestep_bytes;
    if (read_entire_topology) thisSize += (nt-1) * timestep_bytes;
    if (thisSize > bufsize_m2) bufsize_m2 = thisSize;

    //------ scan all timesteps in the topology

    for(i=0; i<nt; i++) {
      m2p = m2_array + n_m2++;			// m2_array is an array

      m2p->nx = nx;
      m2p->ny = ny;
      m2p->nt = nt;
      m2p->dataType = (i == 0) ? 'T' : 'F';	// T=topology start, F=func values
      m2p->iTime_rel = i;
      //ntime_m2 += m2p->nt;

      //------ .iTop  is offset in .bin of topology (coordinates) for this timestep
      //       .iThis is offset in .bin of data for this timestep

      m2p->iTop = iTop;
      m2p->iThis = iTop + skip_coord + i*timestep_bytes;
      iNext = m2p->iThis + timestep_bytes;

      if (debug_M2_detail || m2_verbose) 
	printf("   %d: nx=%d, ny=%d, nt=%d, iTop=%ld, iThis=%ld, iNext=%ld\n", 
	       n_m2-1, m2p->nx, m2p->ny, m2p->nt, m2p->iTop, m2p->iThis, iNext);
    }

    //------ Slava bug - we have topologies with no timesteps

    if (multi_topology > 0 && nt == 0) {
      iNext = iTop + skip_coord;
    }
    if (multi_topology > 0 && nt > 0 ) {
      nx_prev = nx; ny_prev = ny; nt_prev = nt;
    }

    //------ Get offset for next topology

    iTop = iNext;				// next topology
    if (multi_topology == 0) break;
  }

  if (debug_M2_detail || m2_verbose) 
    printf("Done read_m2_array, ntime=%d, max bufsize = %ld\n\n", n_m2, bufsize_m2);
  return 0;
}

/*-----------------------------------------------------------------------------
|	invert_bytes
-----------------------------------------------------------------------------*/
int invert_bytes(unsigned char *p, unsigned char *q)
{
  int request;
  request = *(Long *)p;
  *(q + 0) = *(p + 3);	/* bytes are written reversed!! */
  *(q + 1) = *(p + 2);
  *(q + 2) = *(p + 1);
  *(q + 3) = *(p + 0);
  request = *(Long *)q;
  return request;
}

/*-----------------------------------------------------------------------------
|	dump_buf
-----------------------------------------------------------------------------*/
void dump_buf(char *label, float *buf0)
{
  float *p;
  Long i;
  FILE *dfile;
  static int inited=0;
  dfile = fopen("dump_buf.out", inited==0? "wt" : "at");
  inited++;
  printf("WRITE TO DUMP_BUF.OUT\n");
  fprintf(dfile, "%d. %s: buf at %lx = %lx\n", inited, label, buf, buf0);
  for (i=0, p=buf; p<end_buf; i++,p++) {
    fprintf(dfile, "%6d. %g\n", i, *p);
  }
  fclose(dfile);
}

/*-----------------------------------------------------------------------------
|	binread. Read from xx.bin data file using handle fd.
|	*sizep = size of buf already allocated, in bytes
|	nh = number of records to be considered as header (G=0, C=2)
|
|	Format of the input data:
|	*  "words" are 4 bytes long: float or Long int
|	*  records are: <nbytes>, ...<nbytes of data>..., <nbytes>
|	*  draw.c had eg. 3 parameters per record, many records per curve,
|	   a curve was terminated by a record with nbytes=0, use size
|	   of file to terminate
|	*  contour.c had 1 parameter, all data in the data record.
|	   other curves followed, but the program did not use.
|	   12, m,n,o, 12, 16, rmin..zmax, 16, m*n*4, data... , m*n*4, ..
|
|	Format of the output data:
|	* header records: <nbytes>, ...<nbytes of data>...
|	* data records: <nparam>,<nbytes>, ...<nbytes*nparam of data>...
|	* if multiple curves, end when nparam=0.
|
|	* contour: 0C, <mr,mz,mpsi>, 10, <rmin..zmax>,
|		   4, <nparam>, data, mparam, data,...
-----------------------------------------------------------------------------*/
float *binread(char *fname, unsigned Long *sizep)
{
  unsigned Long oldsize, newsize, nrec, bytes_scanned, filesize, actual_filesize;

  float temp;
  static float *buf0;

  unsigned Long nchar, ngrid, ncoord, nchar_h, nchar_t;
  unsigned Long rec, rec0, count, count0, request;
  Int i, i2;
  int fd, nhead, oflag, result;
  FILE *file;
  static int nhead0;
  int tellRec;
  int ntry;
  int unfinished;
  Long *p0, *longp;
  unsigned char *p, *q;
  float *f;
  Int *data;

  static int debugg = 500;
  int iError;
  M2_ARRAY *m2p, *m2p_prev;
  int iRead_Top;


  /*------ Open Input .bin File --------------------------*/

  oflag = O_RDONLY;
  if (debug_M2_detail) 
    printf("%d. Read file %s\n", multi_topology, fname);

  fd = open(fname, oflag, 0);				/* eg. bal3.bin */

  if (fd == -1) {
    printf("Abort: Cannot open %s.\n", fname);	// (can have many input files
    exit(0);
  }

  //------ Initialize Multi Topology: load array m2_array[], one for each topology
  //	   get bufsize_m2, size_1stRec; initialize i_m2 (current topology)

  one_timestep = (ftype == 6 && M_one_timestep && 
		  (multi_topology==0));

  if (multi_topology==1 || one_timestep) {
    iError = read_m2_array(fd, 0L);
    if (iError) abort();
    multi_topology = (multi_topology==1) ? 2 : 3;
    i_m2 = 0; i_m2_prev = -1;
    close(fd);
    return;					// ... and return (await another call to binread)
  }

  /*------ determine file size, allocate buffer, and read binary data */

  oldsize = *sizep;
  newsize = lseek(fd, 0L, SEEK_END);	// (warning if not UNIX:  do newsize &= 0xffffL because compiler does cwd?)
  actual_filesize = newsize;
  if (multi_topology>=2) newsize = bufsize_m2;

  filesize = newsize;
  newsize += sizeof(float);		/* one more word, in case need for <nparams> */
  ntry = 0;

  //------ Allocate buffer
  
  if (oldsize == 0) {		// 1st time: allocate
    request = newsize + 8;
    buf0 = (float *) malloc(request);
    if (debug_M2_detail) printf("allocate %ld\n", request);
  }

  //------ Multi Topology, 3rd etc time: already allocated max

  else if (multi_topology>=2) {
    request = current_buf_size;
  }

  //------ Or Realloc buffer for 2nd .bin file

  else {			// 2nd etc time: reallocate
    request = oldsize + newsize + 8;
    if (request>oldsize) {
      buf0 = (float *) realloc(buf, (size_t)request);
      if (debug_M2_detail) printf("3. allocate %ld\n", request);
    }
  }

  current_buf_size = request;
  if (debug_M2) printf("Buf allocated at %lx, size %lx\n", buf0, current_buf_size);

  if (buf0 == 0)		// if allocate/reallocate fails, try again or abort
    {
      printf("Insufficient memory for data.\n");
      exit(1);
    }

  //-----------------------------------------------------------------------------
  // read into buf0 here starting at (buf + 1) 
  // if type=0, could have multiple input files
  // read ENTIRE file, or if multi_topology read next topology or next time step
  // NOTE we do not handle blocks! e.g. there would be 1 "coords" for each block
  //-----------------------------------------------------------------------------

  lseek(fd, 0L, SEEK_SET);				// reset file pointer

  //------ if not multi_topology, read everything

  if (multi_topology < 2) {
     buf = (float *)buf0 + oldsize / sizeof(float);	// append
     p = (unsigned char *)(buf+1);			// load data at (buf+1)
     read(fd, p, newsize);				// read everything at (buf+1)
     //write_smaller_bin((Long *)p);
  }

  //------ if yes multi_topology, read everything except...
  //	    don't read the "timestep header" records (4, 1, 4)

  else {
    buf = buf0;					// overwrite
    m2p = m2_array + i_m2;
    ngrid = m2p->nx * m2p->ny;
    ncoord = ngrid * 2L;		        // #floats in coordinate part of topology
    nchar_t = (2L + ngrid) * nqty * nnode_r;	// # chars in timestep
    nchar_t *= sizeof(float);
    nchar_h = (2L + ncoord)*sizeof(float);	// #chars in header + coordinates in xx.bin
    nchar_h += size_topRec;
    nchar_h *= nnode_r;
    //printf("ngrid %ld, ncoord %ld, nchar_h/4 %ld, nchar_t/4 %ld\n", ngrid, ncoord, nchar_h/4, nchar_t/4);

    iRead_Top = 1;
    if (i_m2_prev >= 0) {
      m2p_prev = m2_array + i_m2_prev;
      if (m2p->iTop == m2p_prev->iTop) iRead_Top = 0;
    }

    //---- Read header part of topology 
    //	   (M2: 5 then 9 + nx*ny*2 + 2), (M: 5 then 4 + nx*ny*2*n_rblock)

    iRead_Top = 1;		// Couldn't make it work to just read new timestep!!

    if (iRead_Top) {
      p = (unsigned char *)(buf+1);		// load data at (buf+1)
      read(fd, p, size_1stRec);			// read 1st record
      lseek(fd, m2p->iTop, SEEK_SET);		// skip to topology header for CURRENT topology
      read(fd, p+size_1stRec, nchar_h);		// read topology header plus coords for topology
      q = p + size_1stRec + nchar_h;		// ready to read timesteps at q
      offsetbuf = 0;
    }						// p, q both unsigned char

    //---- Or initialize to just read the timestep (but I didn't make it work!!)

    else {
      p = q = (unsigned char*)(buf + ncoord);
    }

    //------ Read timestep value(s) here

    lseek(fd, m2p->iThis, SEEK_SET);
    i2 = (read_entire_topology) ? m2p->nt : 1;	// all timesteps or just one?

    for (i=0; i<i2; i++) {		// for each timestep...
      if (multi_topology==2)
	lseek(fd, 12L, SEEK_CUR);	// skip timestep header record (4,1,4)
      read(fd, q, nchar_t);		// read the data for the timestep into buf[]
      q += nchar_t;			// nchar_t includes nqty records
    }
    newsize = q - (unsigned char *)buf0;
    //offsetbuf = 0;
    if (debug_M2_detail) printf("newsize=%ld of %ld\n", newsize, current_buf_size);
  }
  if (debug_M2_detail) printf("data was read at %ld bytes\n", p - (unsigned char *)buf);

  data = (Int *)p;
  if (debug_m) printf("Read %ld bytes, 1st long=%ld\n\n", newsize, *data);

  unfinished = 0;

  //------ Add final null record if missing

  data = (Int *)buf;			/* calc data 2 steps so huge ok */
  data += newsize / sizeof(float) - 2;

  if (*data != 0 || *(data + 1) != 0)	/* check for termination null rec */
    {
      *(data + 2) = 0;			// add missing termination
      *(data + 3) = 0;
      newsize += 8;
      if (ftype == 0) filesize += 8;
    }

  //------- Test if need to invert bytes 
  //        0th word is count, we're taking a chance that it's huge
  //        type G has 1st word float, so we can't use it

  if (invert_bytes_flag == -999) {		
    if (ftype == 0) request = *(Long *)buf;	/* type G use the count (all others are float) */
    else request = *(Long *)(buf+1);
    invert_bytes_flag = (request & 0xffff0000) != 0;
    if (force_invert_bytes != -1)
      invert_bytes_flag = force_invert_bytes;
    //printf("invert_bytes=%d, based on %ld for type %d\n", invert_bytes_flag, request, ftype);
  }

  //----------------------------------------------------------
  // Scan through buffer here
  //----------------------------------------------------------

  rec0 = 0;
  count0 = 0;
  ntry =0;				// number of block for mb
  nhead = 0;
  count = 0;

  if (multi_topology >= 2 && iRead_Top==0) {
    rec = nrec = 3;
    i = m2p->nx * m2p->ny * 2L; //* sizeof(float);	// offset to beginning of new data
    f = (buf + i);					// next data value goes here
    bytes_scanned = (i + 5L + 9L + 2L) * sizeof(float);
    i *= sizeof(float);
  }
  else {
    rec = nrec = 0;
    i = 4;				// offset (bytes) to data as read in from xx.bin
    f = buf;				// pointer (float) to processed data
    bytes_scanned = 0;
  }

  for (; i < newsize; i+=4)
    {
      /*------ Get next value into 'temp' (invert if needed) */

      p = (unsigned char *)buf + i;

      if (invert_bytes_flag) request = invert_bytes(p, (unsigned char *) &temp);
      else temp = *(float *) p;
      bytes_scanned += sizeof(float);
      longp = (Long *)p;		/* debug; current "temp" */

      //if (debug_M2) printf("%ld. %g %ld\n", bytes_scanned, temp, *longp);

      /*------ Start New Record: initialize count, test if finish prev rec, misc debug */

      if (count == 0)
	{
	  count = *(Int *)&temp / 4L;	// get new count for record
	  if (size_1stRec==0) size_1stRec = (count + 2L) * sizeof(float);
	  if (debug_M2_detail)
	    if (rec < 9 || rec==47 || rec==48 || (rec>48 && ((rec-49)%24==0)))
	      printf("rec %5d %5d: count=%6ld,(f-buf)=%8ld,words scanned=%ld, %ld\n", 
		     rec,(rec-49)/24, count, f - buf, bytes_scanned/4L -1, (rec-49)/(24L*43L));

	  nrec++;
	  tellRec=0;
	  if (nrec<=100) tellRec=1;

	  if (debug_m == 1 && tellRec)
	    printf ("binread (Rec Start %d), #data=%ld, off=%ld\n",
		    nrec, count, (Long)(f-buf));					/* ,bytes_scanned, filesize */

	  p0 = (Long *)f;

	  /*---- Previous record might be unfinished */

	  if (bytes_scanned + (count*4L + 4L) > actual_filesize) {
	    if (!unfinished) {
	      //printf("Incomplete record encountered. %ld of %ld bytes scanned, pending %ld \n",filesize, bytes_scanned, (count*4L+4L));
	      unfinished = 1;
	    }

	    if (ftype == 0) {		/* type G: generate a 0,0; binread should then abort */
	      count = 0;
	    }
	    else break;
	  }

	  /*---- update loop counts */

	  if (ftype != 6 && start_new_rec(count,count0,rec,rec0,nhead) == 0) 
	    break;
	  count0 = count;

	  //---- Error Handling (NN)

          if (ftype == 5 && rec >0 && count != counting_nrz)
	    d_abort(" Inappropriate number of variables in .bin file",0,0);

#ifdef HANDLE_ERROR_NOT_YET
          /* Comment out error handling by size of rec  - will be done !!! not yet */
          if (ftype == 6 && rec > 0 && rec > nhead) {
	    if (rec <= nhead + 2*nnode) {
                iblock = ( rec - nhead)/2 ;
                if (count != (xylist[iblock].mr*xylist[iblock].mz))
                   d_abort("Inappropriate number of xy variables in .bin file",0,0);
	    }
            else {
                iblock = ( rec - nhead - 2*nnode  );
                if (count != (xylist[iblock].mr*xylist[iblock].mz))
                   d_abort("Inappropriate number of xy variables in .bin file",0,0);
	    }
	  }
#endif			// ... of HANDLE_ERROR_NOT_YET

	  //------ The NEW record has count=0 (null record terminator)

	  if (count==0 || (count ^ 0xFFFFFFFFL)==0 ) {	// zero-length record (^ is xor)
	    count0 = 0;
	    rec0 = rec;
	    
	    if (count==0) {
	      i+=4; rec++;				// skip end count
	      bytes_scanned +=4;
	    }

	    else {
	      rec++;
	      count=0;
	      while ( (request = *(Long *)(p+4))  == -4) // skip ffff fffc
		{ i+=4; p=(unsigned char *)buf + i;}
	    }
	  }
	  continue;
	}		// .. end of count==0 (new record)

      /*------ Transfer 'temp' into buf[] */

      *f = temp;			// read float data
      f++;
      count--;

      /*------ Value just transferred is last in record */

      if (count == 0)
	{
	  if (rec <= nhead) {			// convenient for breakpoint...
	    count0 = count0;			// ...for header recs only
	  }

	  if (rec==0)					// end of FIRST record: initialize
	    {
	      if (debug_m_trc && tellRec) printf("binread (Rec End),   rec=0\n");
	      if (oldsize==0 || multi_topology>=2) 
		nhead0 = initloops(p0,count0);		// get nnode_r eg
	      else if (ftype == 0) 
		nodelist[nnode].ncount = 1;

	      nhead = nhead0;
	      if (ftype == 5 || ftype == 6) {
		//printf("RESTART BUF POINTER from %lx to %lx\n", f, buf);		
		f = buf;
	      }
            }

	  else if (rec <= nhead)			// end of 2nd etc header records
	    {
	      if (debug_m_trc && tellRec) 
		printf("binread (Rec End),   rec=%ld, nhead=%ld, ntry=%ld, nnode_r=%ld\n",
		       rec, nhead, ntry, nnode_r);

	      //---- End of Multiblock r_block header, odd rec#

              if (ftype == 6  && ntry <= nnode_r && (rec & 1))
		{ 
		  if (debug_m_trc && tellRec) printf("   rblock\n");
		  addloop_mb(p0, ++ntry);		// ntry==iblock
		  if (multi_topology != 2) 	// type M can have 0 or 3 here
		    f = f - 2;			// (don't use "-=" so can search for "f =")
		  else f = f - 7;		// M2 info is in xylist[], don't need topology header
		}

	      //---- End of t_block header, 1st rec of iblock =((rec-nnode*2)+2)/3

	      else if (ftype == 6 && ntry > nnode_r && (rec - nnode_r *2) % 3 == 1)
		{
		  if (debug_m && tellRec) printf("tblock\n");
		  addloop_mb(p0, ++ntry); f-=2;		// ntry == iblock?
		}

	      //---- End of other types of header recs (ftype != 6)
		  
	      else if (ftype != 6)
		{
		  if (oldsize==0 && rec == 1) addloop(p0);

	          if (rec == nhead && ftype != 5) {	// new data will be...
		    f = buf;				// ...at beginning of buf[]

		    if (ivar==-1) d_abort("No type V found in data file",0,0);
		  }
		}
	    }

	  //---- End of records that are NOT header records

	  else if (ftype == 0 || ftype == 4)	// regular data; (generalize type N later)
	    update_nodeloop(1);

	  else if (ftype == 1 && !counting_outer) 
	    break;

	  //---- skip end count

	  i+=4; rec++;
	  bytes_scanned +=4;

	  if (bytes_scanned==filesize) 
	    goto BINREAD_DONE;
	}					// end of "last in record"
    }						// end of "for" loop scanning data

  //------ Done reading in the entire xx.bin file

 BINREAD_DONE:
  end_buf = f;
 
  //------ Multiblocks plot - get limits into xylist[].xmin etc, flist[].fmin etc
  //	   Note: cp's (array curveset) not defined if called from init()

  if( ftype == 6 )
    {
      //printf("%d nodes, %11.4e %11.4e %11.4e\n", nnode, *(buf+77602), *(buf+116403), *(buf+271607));
      ntime = (rec -1 -nhead )/(nnode*nqty);
      loop[0].sep = offsetbuf/2 *nqty*nnode;		// timestep
      if (flist == NULL) 
	flist = (BLOCK_F *) malloc((nnode+1)*sizeof(BLOCK_F)*nqty*ntime);
      offsetf = offsetbuf;
      if (debug_M2_detail) printf("binread, ready to call getlimits_mb flist=%lx, off=%lx\n", flist, offsetf);
      getlimits_mb();
      initialize_vectors(NULL);
      if (debug_m) printf("binread, ntime=%d, sep=%ld, qty=%d, nnode=%d, flist alloc for %d\n",
			  ntime, loop[0].sep, nqty, nnode, (nnode+1)*nqty*ntime);
    }
      
  *sizep = (f - buf0) * 4;
  close(fd);

  if (tellMe) xprintf ("**BINREAD, nhead=%d\n", nhead);
  return((float *)buf0);
}				// ... end binread

/*-----------------------------------------------------------------------------
|	read_topology
|	* If multi_topology, test if next time is a new topology,
|	  call binread() if needed
|	* iTime: absolute time - in 0 .. (n_m2-1)
|	* n_m2 = # of topologies, i_m2 = current topology
|	* m2_array[i].nt: # timesteps in topology i
-----------------------------------------------------------------------------*/
void read_topology(CURVE_SET *cp, int iTime, int iTime_last)
{
  int i, itime1, itime2;
  M2_ARRAY *m2p;

  i_m2_prev = i_m2;
  i_m2 = iTime;
  m2p = m2_array + i_m2;

  cp->itime_abs = iTime;
  cp->itime_rel = cp->itime = m2p->iTime_rel;
  if (debug_M2_detail) 
    printf("read_topology, time abs=%d, rel=%d\n", cp->itime_abs, cp->itime_rel);

  if (read_entire_topology && m2p->iTop == (m2_array+iTime_last)->iTop)
    return;

  bufsize = current_buf_size;
  buf = binread(bin_file_name, &bufsize);
  get_limits(cp);
}

/*-----------------------------------------------------------------------------
|	start_new_rec
-----------------------------------------------------------------------------*/
static int start_new_rec(unsigned Long count, unsigned Long count0,
		unsigned Long rec, unsigned Long rec0, int nhead)
{
  static int need_varcount=0;
  int uses_nullrec;
  int j, retval;

  uses_nullrec = (ftype == 0 || ftype == 4);

  if ((count ^ 0xFFFFFFFFL) == 0)		// ^ is xor
      count = 0;

  if (count==0 && (count0==0 || !uses_nullrec))	/* Double null or single */
    return 0;

  retval = 1;

  //------ ftype 0 or 4: update loop, nnode, nodelist

  if (count==0)					/* single null, ftype 0 or 4 (G,L) */
    {
      if (rec>nhead && rec-rec0 <= 2 && !count0)
						/* if loop header... */
	{					/* ...add to loop[] */
	  for(j=0; j<inode && loop[j].hcount!=0; j++) ;
	  if (j<inode)
	    {
	      loop[j].hcount = (int)count0;
	      if (j==0) loop[j].ih0 = 0;
	      else loop[j].ih0 = loop[j-1].ih0 + loop[j-1].hcount;
	      need_varcount=1;
	    }
	  update_nodeloop(0);
	}
      else					/* else count outer loop */
	{
	  if (update_outerloop()) loop->count++;
	  retval = 2;
	}
    }

  //------ ftype 0 or 4: update loop

  else if (uses_nullrec && need_varcount)	// applies to ftype 0 or 4
    {
      loop[ivar].count = (int)count;
      need_varcount = 0;
    }

  //------ ftype 1, 2, 5: update loop

  else if ((ftype == 1 || ftype == 2 || ftype == 5) && rec>nhead && counting_outer)
    loop->count++;
  
  if (rec==nhead+1 && ftype == 4)		/* L give size of extra loop */
    (loop+nloop-1)->count = (int)count;

  return retval;
}

/*-----------------------------------------------------------------------------
|	loop_structure
-----------------------------------------------------------------------------*/
void loop_structure()
{
  int i, ncount;
  LOOP *lp;
  NODE *np;

  for(i=0,lp=loop+nloop-1; i<nloop; i++,lp--)
    {
      if (i==0) lp->sep = 1;
      else lp->sep = (lp+1)->sep * (lp+1)->count + lp->hcount;
      /*----- Note: any lp->sep above a type Node is invalid! */
      if (lp->ltype=='H')
        {
	  if ((lp+1)->ltype != 'N') lp->sep++;
	  else for(i=0; i<nnode; i++) (nodelist + i)->ncount += 1;
	}
    }

  if (inode != -1)
  {
    for(i=0,lp=loop+inode,np=nodelist,ncount_equal=1; i<nnode; i++,np++)
    {
      if (i==0) ncount = np->ncount;		/* See if equal count */
      else if (np->ncount != ncount) ncount_equal = 0;
      np->nsep = np->ncount * lp->sep;		/* Obtain nsep */
      np->nsep += (lp-1)->hcount;		/* Kludge */
    }
    if (ncount_equal)
    {
      for(i=0,lp=loop; i<inode; i++,lp++) lp->use_sep = '1';
      loop[inode].ltype = 'I';
      inode = -1;
    }
  }
}

/*-----------------------------------------------------------------------------
|       write_smaller_bin
-----------------------------------------------------------------------------*/
void write_smaller_bin(Long *p)
{
  int fd2, value;
  Long temp, nwrite, fbytes, hbytes;
  char fname[100];
  unsigned char c, *q, *qq;
  int i1, i2, nqty, itime, iblock, iqty;

  //value = invert_bytes((unsigned char *)p, (unsigned char *) &temp);
  //temp = *(float *) p;
  //printf("%d\n", temp);

  i1 = 0;
  i2 = 30;
  nqty = i2-i1+1;

  //------ carl2.bin

  strcpy(fname,"/win/xdraw/carl/carl2.bin");
  fd2 = open(fname, O_RDWR | O_CREAT | O_TRUNC, 0777);
  hbytes = (2+3) + 24 * ((2+2) + (2+950));
  hbytes = hbytes * 4;

  q = (unsigned char *)p;
  c = *(q+15);
  *(q+15) = nqty;

  nwrite = write(fd2, q, hbytes);
  printf("%ld bytes written\n", nwrite);
  *(q+15) = c;
  
  //------ Assume time=outer loop, next block, then qtys (are contiguous)
  //	   each block is 19*25 = 475 floats; we write 

  fbytes = (2+475)*4;	// size of one block

  for(itime=0, qq=q; itime < 3; itime++) {
    for (iblock=0; iblock<24; iblock++) {
      for(iqty=0; iqty<43; iqty++) {
	if (iqty < nqty) {
	  nwrite = write(fd2, qq, fbytes);
	  printf("%d %d %d written\n", itime, iblock, iqty);
	}
	qq += fbytes;
      }
    }
  }

  close(fd2);
	 
}

/*-----------------------------------------------------------------------------
|       save_as_g
|	* drawcrash1.in, graphs.bin has 1 block, 11 timesteps, 20 qtys
|	* save x, y, Rho (0), Vy (15) for current t (t=3)
|	* Rho: x_info 0,0,0; y_info 4224,0,0; z_info 67584,64,1
|	* iqty's
-----------------------------------------------------------------------------*/
Long nbyte;
void save_as_g()
{
  int ix, iy;
  float x, y, f1, f2;
  int handle;
  Long off_x, off_y, off_f1, off_f2;
  Long iqty1, iqty2, nx, ny, nxy;
  int temp[3];
  float func[4];

  iqty1 = 0; iqty2 = 15;
  handle = open("bin/crash1/graphs2.bin", O_TRUNC | O_RDWR, S_IRWXU | S_IRWXG | S_IRWXO);
  printf("handle %x\n", handle);
  nx = 64;
  ny = 66;
  nxy = nx * ny;
  temp[0] = (nx * 4) * 4;
  temp[1] = temp[2] = 0;

  for(iy=0; iy<ny; iy++) {
    for(ix=0; ix<nx; ix++) {
      off_x = iy * nx + ix;
      off_y = off_x + nxy;

      off_f1 = 2 * nxy + iqty1 * nxy + iy * nx + ix;
      off_f2 = 2 * nxy + iqty2 * nxy + iy * nx + ix;
      func[0] = buf[off_x];
      func[1] = buf[off_y];
      func[2] = buf[off_f1];
      func[3] = buf[off_f2];

      temp[0] = 4 * 4;
      nbyte = write(handle, temp, 4);
      nbyte = write(handle, func, 16);
      nbyte = write(handle, temp, 4);
    }
    nbyte = write(handle, temp+1, 8);
  }
  close(handle);
  printf("Done\n");
}
