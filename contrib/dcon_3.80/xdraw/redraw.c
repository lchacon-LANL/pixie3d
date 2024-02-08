/******************************************************************************
**  NAME      REDRAW.C
**  AUTHOR    Sheryl M. Glasser
**
**  DESCRIPTION
**      Usage: draw [options] [suffix]
**      options: -d = draw (draw.in), -c = cntour (cont.in)
**
**  Copyright (c) GlassWare 1992.  All rights reserved.
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>

#include "gendefs.h"
#include "curves.h"
#include "xdraw.h"
#include "xinit.h"
#include "spline.h"
#include "xcontour.h"
#include "xedit.h"
#include "xtools.h"
#include "setcolor.h"
#include "menu.h"
#include "ps.h"

extern LOOP loop[];
extern int nloop;
extern CURVE_SET curveset[], *cp;
extern NODE *nodelist;
extern int nnode;			/* TOTAL count of nodelist */
extern int inode;
extern float *buf;
extern Font font_greek,font_latin;
extern XFontStruct *font_struct,*font_struct1;
#define LATIN 1
#define GREEK 2

#define F32 (float)32000.
#define bigf(x,y) x>F32 || x<-F32 || y>F32 || y<-F32

/*=============================================================================
**              SUPPORTING FUNCTIONS - externals from xdraw.c
**===========================================================================*/
extern Display *mydisplay;		/* set via redraw()-->get_expose_info() */
extern Window redraw_win;
extern GC redraw_gc;
extern float clipx1, clipy1, clipx2, clipy2;
extern int monochrome;
extern int wx0, wy0, w, h, wx2, wy2;
extern int nqty_in_win, qtys_in_win[20];
extern FloatSc xscale, xoffset, yscale, yoffset;
extern int debug_scale_evis;

extern int ps_modeon;
extern int verifying;	/* 1 ==> ps file uses actual data */

static char Mark = '.';
static char markerstyle[] = "+x#^";
int DM = 4;

/*-----------------------------------------------------------------------------
|	get_off -- get 1st offset of a variable, k=current family
-----------------------------------------------------------------------------*/
static long get_off(int k, int kl, long off, CVAR *xp)
{
  long j;
  int i, i2;
  int it,ipsi,isurf,nt,npsi,nsurf;
  
  if ((xp->hloop >= 0 && xp->hloop < kl) ||
      k == 0) return(xp->off0);
  j = xp->dfaml;
  if (j<=0)		/* 0 ==> h */
  {
    it = ipsi = isurf = 0;
    nt = npsi = 1;
    nsurf = loop[inode-1].count;
    if (inode==1) isurf = k;
    else if (inode==2)
    {
      npsi = loop[0].count;
      if (kl==1) { ipsi = loop[0].current; isurf = k; }
      else       { ipsi = k; isurf = loop[1].current; }
    }
    else if (inode==3)
    {
      nt = loop[0].count;
      npsi = loop[1].count;
      if (kl==2) { it = loop[0].current; ipsi = loop[1].current; isurf = k; }
      if (kl==1) { it = loop[0].current; ipsi = k; isurf = loop[2].current; }
      if (kl==0) { it = k; ipsi = loop[1].current; isurf = loop[2].current; }
    }
    i2 = it * npsi * nsurf + ipsi * nsurf + isurf;
    for(i=0,j=xp->off0; i<i2; i++) j += (nodelist+i)->nsep;
    return( j );
  }
  else return(off+j);
}

/*-----------------------------------------------------------------------------
|	get_data_pointers
|	* get offset into buf EVEN IF one of ix0,iy0 is index, not data
-----------------------------------------------------------------------------*/
static int get_data_pointers(int k, int first, int reversed,
		       CURVE_SET *cp, LOOP *qs, long *ix0p, long *iy0p)
{
   static long ix0,iy0;
   long off;
   int kl, ik, ns, ncurve0;
   static int ncurve;

   ix0 = *ix0p;
   iy0 = *iy0p;
   kl = cp->lfaml;
   if (inode >= 0)
   {
      if (!first) ncurve0 = ncurve;
      ncurve = get_ncount();		/* #outer loops to step over,depends... */
					/* ...on inode,loop[i].current */
      if (first) off = reversed ? step_over_ncount(0, ncurve) : 0;
      else if (!reversed) off =  step_over_ncount(ncurve0,ncurve);
      else                off = -step_over_ncount(ncurve,ncurve0);
   }
   
   if (cp->x_info.index & LOOPBIT) ix0 = -1L;
   else if (inode==-1)
   {
      if (!reversed) ix0 = get_off(k, kl, ix0, &cp->x_info);
      else for(ik=0; ik<=k; ik++) ix0 = get_off(ik, kl, ix0, &cp->x_info);
   }
   else
   {
      if (first) ix0 = cp->x_info.off0;
      ix0 += off;
   }

   if (cp->y_info.index & LOOPBIT) iy0 = -1L;
   else if (inode==-1)
   {
      if (!reversed) iy0 = get_off(k, kl, iy0, &cp->y_info);
      else for(ik=0; ik<=k; ik++) iy0 = get_off(ik, kl, iy0, &cp->y_info);
   }
   else
   {
      if (first) iy0 = cp->y_info.off0;
      iy0 += off;
   }

   *ix0p = ix0;
   *iy0p = iy0;
   if (inode!=-1 && inode==cp->lstep) ns = (nodelist+ncurve)->ncount;
   else ns = qs->count;
   return ns;
}

/*-----------------------------------------------------------------------------
|	next_data
-----------------------------------------------------------------------------*/
static int next_data(int i, LOOP *qs, long ix0, long iy0, 
               float *fx1u, float *fy1u, float **bxp, float **byp)
{
   float *bx, *by;
   long sep;
   int use_sep, iError;
   double x,y,z;

   iError = 0;
   bx = *bxp;
   by = *byp;
   if (i==0)
   {
     bx = addfloat(buf, ix0);		// bx = buf + ix0
     by = addfloat(buf, iy0);
   }
   use_sep = ( inode==-1 || qs>=(loop+inode) );

   //------ Fetch x

   if (ix0 < 0) *fx1u = get_loopval(qs, i);
   else {
     *fx1u = *bx;
     sep = use_sep ? qs->sep : nodelist[i].nsep;
     bx = addfloat(bx, sep);		// bx += sep
   }

   //------ Fetch y

   if (iy0 < 0) *fy1u = get_loopval(qs, i);
   else {
     *fy1u = *by;
     sep = use_sep ? qs->sep : nodelist[i].nsep;
     by = addfloat(by, sep);		// by += sep
   }

   //------ Convert to Log if appropriate

   if (cp->use_log & 1) { 
     z = (double)(*fx1u);
     //printf("xvalue %lg --> %lg\n", z, x); 
     if (cp->use_cutoff & 1 && z < cp->xcutoff) iError = 1;
     else { x = log10(z); *fx1u = (float)x; }
   }
   if (cp->use_log & 2) { 
     z = (double)(*fy1u);
     //printf("yvalue %lg --> %lg: %d, %f\n", z, y, cp->use_cutoff, cp->ycutoff); 
     if (cp->use_cutoff & 2 && z < cp->ycutoff) iError = 1;
     else { y = log10(z); *fy1u = (float)y; }
   }

   *bxp = bx;
   *byp = by;
   return iError;
}

/*-----------------------------------------------------------------------------
|	shrinkf
|	* fx0..fy1 ok on easy clip test.
|	* If ok, all final values will be positive (i.e. in viewport)
|	* If either x<0, the other must be >=0.  Ditto for y. (easy clip test)	
-----------------------------------------------------------------------------*/
void shrinkf(float,float,float,float,int*,int*);
void shrinkf(float fx1, float fy1, float fx0, float fy0, int *x1, int *y1)
{
  float fdx, fdy, fxc, fyc;
  int new;
  fxc = fx0; fdx = fx1 - fx0;
  fyc = fy0; fdy = fy1 - fy0;
  for(new=0; fxc < FZ || fyc < FZ; new=1 )
    {
      fdx /= (float)2; fxc += fdx;
      fdy /= (float)2; fyc += fdy;
    }
  if (new) { fdx = fx1 - fxc; fdy = fy1 - fyc; }

  for(; bigf(fx1, fy1); )
    {
      fdx /= (float)2; fx1 = fxc + fdx;
      fdy /= (float)2; fy1 = fyc + fdy;
    }
  *x1 = float_to_int(fx1);
  *y1 = float_to_int(fy1);
}

/*-----------------------------------------------------------------------------
|	delay_and_clear -- for animation
-----------------------------------------------------------------------------*/
void delay_and_clear(int, int);
void delay_and_clear(int kk, int ncurve)
{
   clock_t clocktime;
   long now;
   static long last, start;
#define millisec(c) c /= (CLOCKS_PER_SEC / 1000)

   for(last=now=0;;)
   {
      clocktime = clock();
      millisec(clocktime);
      if (kk==0) start = clocktime;
      now = clocktime - start;
      if (kk==0 || (now-last) > ncurve*100) break;
   }
   XClearArea(mydisplay, redraw_win, 0, 0, 0, 0, False);
   draw_frame();
}

/*-----------------------------------------------------------------------------
|	rand_label
-----------------------------------------------------------------------------*/
int rand_label(int ns)
{
  static int now;
  int labeli;
  if (ns==0)
    {
      now = labeli = 2;
    }
  else if (ns<=3)
    {
      labeli = now+1;
      if (labeli==ns) labeli=0; 
      now = labeli;
      }
  else
    {
      now += 3;
      if (now > 10) now -= 10;
      labeli = now * ns / 11;
    }
  return labeli;
}

/*-----------------------------------------------------------------------------
|	clip
|	* for VERY BIG zoomed values, returns clipped, else assumes X handles
-----------------------------------------------------------------------------*/
int clip1(float *, float *, float, float, float, float);
int clip(float *fx0, float *fy0, float *fx1, float *fy1)
{
  int which, ok;
  float x0, y0, dx, dy;
  if (
    (*fx0<clipx1 && *fx1<clipx1) || (*fx0>clipx2 && *fx1>clipx2) ||
    (*fy0<clipy1 && *fy1<clipy1) || (*fy0>clipy2 && *fy1>clipy2)
     ) return 0;
  which = 0;
  if (bigf(*fx0, *fy0)) which |= 1;
  if (bigf(*fx1, *fy1)) which |= 2;
  if (!which) return 1;
  dx = *fx1 - *fx0; x0 = *fx0;
  dy = *fy1 - *fy0; y0 = *fy0;
  if (which | 2      ) ok = clip1(fx1, fy1, x0, y0, dx, dy);
  if (which | 1 && ok) ok = clip1(fx0, fy0, x0, y0, dx, dy);
  return (ok ? 1 : 0);
}

#define newy() dy*(x-x0)/dx+y0
#define newx() dx*(y-y0)/dy+x0

int clip1(float *x1, float *y1, float x0, float y0, float dx, float dy)
{
  float x, y;
  int which;
  x = *x1;
  y = *y1;
  which = 0;
  if      (x < clipx1) { x = clipx1; which = 1; }
  else if (x > clipx2) { x = clipx2; which = 1; }
  else if (y < clipy1) { y = clipy1; which = 2; }
  else if (y > clipy2) { y = clipy2; which = 2; }
  if (which==1)
    {
      y = newy();
      if (y<clipy1)	 y = clipy1;
      else if (y>clipy2) y = clipy2;
      else which = 0;
    }
  else if (which==2)
    {
      x = newx();
      if      (x<clipx1) x = clipx1;
      else if (x>clipx2) x = clipx2;
      else which = 0;
    }
  if	  (which == 1) { x=newx(); if (x<clipx1 || x>clipx2) return 0; }
  else if (which == 2) { y=newy(); if (y<clipy1 || y>clipy2) return 0; }
  *x1 = x;
  *y1 = y;
  return 1;
}

/*-----------------------------------------------------------------------------
|	redraw0.  draw for original option
-----------------------------------------------------------------------------*/
int redraw0(int ic, int iqty_in_win, int nqty_in_win,
	     float xmin, float xmax, float ymin, float ymax)
{
  int x0, y0, x1, y1, lx, ly, jcurve, drawMe;
  int i, j, k, kk, nf, ns, nMarker;
  int labeli, reversed, ilast, moveit, firstmove, wait_for_2nd, flags;
  int nline, align, iError;
  static int spline=0;
  float fx1, fy1, fx1u, fy1u, fx1old, fy1old;
  float fx0, fy0, fx0u, fy0u;
  LOOP *qs, *qf;
  CURVE_SET *cp;
  long ix0, iy0;
  XPoint *points;
  char text[40], *p, *labeltext;
  XTextItem *item;
  int nitems;
  float *bx, *by;
  static int test=0;
  float xtestu, ytestu;
  struct LABEL *plabel;
  int iview;

  /*------ Draw curves */

  iview = getview(redraw_win);
  cp = curveset + ic;
  flags = cp->flags;
  reversed = flags & REVERSE;
  spline = (flags & SPLINE) ? 1 : 0;
  //printf("REDRAW0, x=%g to %g, y= %g to %g\n", xmin, xmax, ymin, ymax);

  for(i=0; i<nloop; i++)
    if ((j=cp->param[i]) >= 0) (loop+i)->current = j;

  if (cp->lfaml >= 0) { qf = loop + cp->lfaml; nf = qf->count; }
  else { qf = NULL; nf = 1; }
  qs = loop + cp->lstep;

  plabel = cp->label;
  rand_label(0);
  if (ps_modeon && flags & LABELF) ps_setlabelfont(0);
  nMarker = strlen(markerstyle);

  //------ Scan through family in curveset cp

  for(kk=0; kk<nf; kk++)
    {
      //printf("**Curve %d\n", kk);
      k = reversed ? nf - 1 - kk : kk;	   
      if (qf) qf->current = k;			/* (NULL if lfaml=-1) */

      //---- ns = number of points in this curve
      //     ix0,iy0 = offset in buf of 1st point of qty+family
      //     eg drawg8.in has 4 qtys, 4 curves in family, 9 points per curve
      //     so for window having x=qty1 (where first qty is qty0)...
      //     ...ns=9 and ix0=1,37,73,109
      ns = get_data_pointers(k, kk==0, reversed, cp, qs, &ix0, &iy0);

      //printf("View %d, Draw family member %d of %d, curveset[%d] (iqty=%d, labelf=%d, ix0,iy0=%ld,%ld)\n", 
      // iview, kk, nf, ic, iqty_in_win, flags & LABELF, ix0, iy0);

      if (flags & SINGLE) {
	for(jcurve=drawMe=0; jcurve<9; jcurve++) {
	  if ( nqty_in_win>1 && cp->icurve[jcurve] == iqty_in_win) drawMe=1;
	  else if (nqty_in_win==1 && cp->icurve[jcurve] == kk) drawMe=1;
	}
	if (!drawMe) continue;
      }

      if (k%cp->fskip !=0) continue;

      if (monochrome)
	setcolor(0, nf, -1);
      else if (nqty_in_win>1 && ((flags & COLORS) || nf==1))
	setcolor(iqty_in_win*nf+k, nf*nqty_in_win, -1);
      else
	setcolor(k, nf, -1);
      
      Mark = (cp->mcount >= 0) ? *(markerstyle + (k%nMarker)) : '+';
      //printf("Mark = '%c' for kk=%d, count=%d\n", Mark, kk, cp->mcount);

      if (flags & ANIMATE) delay_and_clear(kk, cp->ncurve);
      if (flags & FILL) points = (XPoint *)malloc(ns*sizeof(XPoint));

      labeli = rand_label(ns);
      lx = ly = -100;

      fx1u = (xmin+xmax) / 2.;
      fy1u = (ymin+ymax) / 2.; 

      fx1old = fx1u * xscale + xoffset;
      fy1old = fy1u * yscale + yoffset;

      //------ If spline, allocate array (init_spline returns 0 if not fail allocate)

      if (spline==1 && !init_spline(ns,0)) spline=0;

    AGAIN:				/* (Return here 2nd time thru spline)*/
      firstmove = 0;
      wait_for_2nd = 0;
      for(i=0,ilast=-1; i<ns; i++)	/*.......... Step through points */
	{
	  x0 = x1; fx0 = fx1old; fx0u = fx1u;
	  y0 = y1; fy0 = fy1old; fy0u = fy1u;
	  moveit = (!firstmove || i>ilast+1 || (flags & MARKERS));

	  //------ Fetch next data value

	  iError = next_data(i, qs, ix0, iy0, &fx1u, &fy1u, &bx, &by);
	  //printf("%d. ix0=%ld, off=%ld\n", i, ix0, bx-buf);
	  if (iError) {
	    firstmove = 0; wait_for_2nd = 1;
	    continue; 
	  }

	  fx1 = fx1old = xscale * fx1u + xoffset;
	  fy1 = fy1old = yscale * fy1u + yoffset;
	  if (wait_for_2nd) { wait_for_2nd = 0; continue; }

	  if (fx1 != fx1 || fy1 != fy1) {	/* nan detected! */
	    printf("Zoom limits exceeded.\n%c", '\007');
	    return 1;
	  }

	  //---- If spline, 1st time thru build array of internal points

	  if (spline==1)
	    {
	      add_to_spline(i, fx1, fy1);
	      continue;
	    }

	  if (!clip(&fx0, &fy0, &fx1, &fy1)) continue;
	  if (test)
	  {
	     if (fx0u*fx1u < FZ)
		ytestu = fy0u + (fy1u - fy0u) * (-fx0u) / (fx1u - fx0u);
	     else if (fy0u*fy1u < FZ)
		xtestu = fx0u + (fx1u - fx0u) * (-fy0u) / (fx1u - fx0u);
	     test = test;
          }
	     
	  //printf ("Call Float_to_int, %g, %g to  %g, %g, moveit=%d\n", fx0,fy0,fx1,fy1, moveit);

	  x0 = float_to_int(fx0);
	  y0 = float_to_int(fy0);
	  x1 = float_to_int(fx1);	/* (PC, fx1=10**6 eg...) */
	  y1 = float_to_int(fy1);	/*  doesnt to float to int OK) */

	  if (flags & LABELF && i == labeli && cp->label==NULL)
	    {
	      lx = x1; ly = y1;
	    }
	  ilast = i;
	  if (i>0 && x0==x1 && y0==y1) continue;
	  if (ps_modeon && verifying)
	    {
	      ps_flineto(fx1u, fy1u);
	      continue;
	    }

	  if (!ps_modeon && (flags & FILL))	/* for SPLINE, allocate..*/
	    {					/* .. when spline==2 */
	      points[i].x = x1;
	      points[i].y = y1;
	    }
	  if ((flags & MARKERS ||  k < cp->mcount) &&
	      (fx1==fx1old && fy1==fy1old))
	    {
	      drawmarker(x1, y1, fx1, fy1);
	      if ((flags & DOTS) || k < cp->mcount) continue;
	    }
	  if (flags & DOTS)
	    {
	      if (!ps_modeon)
		XDrawPoint(mydisplay, redraw_win, redraw_gc, x1, y1);
	      else
	        {
		  ps_fmoveto(fx1, fy1);
		  ps_flineto(fx1+1,fy1);
		}
	    }
	  else if (i==0) ;

	  //---- draw regular line here (if spline, draws 3 or 5 points per line)

	  else
	    {
	      drawline(fx0, fy0, fx1, fy1, moveit, spline, i);
	      firstmove = 1;
	    }
	}

      if      (spline==1) { spline = 2; fit_spline(); goto AGAIN; }
      else if (spline==2) { spline = 1; zap_spline(); }

      if (ps_modeon && !(flags & FILL)) ps_stroke();

      if ((flags & FILL) && !ps_modeon)
	{
	  XFillPolygon(mydisplay, redraw_win, redraw_gc, points, ns,
		       Complex, CoordModeOrigin);
	  free(points);
	}

      if ((flags & FILL) && ps_modeon) ps_fill();

      if (flags&LABELF && (cp->label==NULL))
	{
	  sprintf(text, "%d", k+iqty_in_win);
	  p = text;
	  nitems =0;
	  draw_label(lx, ly, -1, p,item,nitems);
	}
      else
	{
	  if (flags & LABELF&&(cp->label!=NULL))
	    {
	      plabel = cp->label;
	      if (nqty_in_win>1 && ((flags & COLORS) || nf==1))
		i=iqty_in_win*nf+k;
	      else
		i=k;

	      /*printf ("i=%d\n",i);*/
	      while (plabel !=NULL)
		{
		  plabel=
		    get_label(i, plabel , &lx, &ly, &labeltext, &nline, &align,
			      xscale,xoffset, yscale,yoffset, &nitems,&item);
		  p = labeltext;
		  draw_label(lx, ly, -1, p, item, nitems);
		}
	    }
	}
      if (cp->label != NULL) draw_misc_labels(cp);

    }						//... end of kk loop

  //printf("Done REDRAW0\n");
  return 0;
}				/* redraw0 */

/*-----------------------------------------------------------------------------
|	redraw00.  
|	* called for type G if E_VIS flag is on
|	* scan for xmin ... ymax
|	* doesn't draw, just scans curves and gets limits
-----------------------------------------------------------------------------*/
void redraw00(int ic, int iqty_in_win, int nqty_in_win,
      float *xminp, float *xmaxp, float *yminp, float *ymaxp)
{
  float xmin, xmax, ymin, ymax;
  int i, j, k, nf, ns, reversed, flags, skipMe;
  static int spline=0, any;
  float fx1u, fy1u;
  LOOP *qs, *qf;
  CURVE_SET *cp, *cp0;
  long ix0,iy0;
  float *bx, *by;

  /*------ Initialize for this curve set */

  cp =  curveset + ic;
  flags = cp->flags;
  reversed = flags & REVERSE;
  spline = (flags & SPLINE) ? 1 : 0;

  xmin = *xminp; xmax = *xmaxp;
  ymin = *yminp; ymax = *ymaxp;

  for(i=0; i<nloop; i++)
    if ((j=cp->param[i]) >= 0) (loop+i)->current = j;

  if (cp->lfaml >= 0) { qf = loop + cp->lfaml; nf = qf->count; }
  else { qf = NULL; nf = 1; }
  qs = loop + cp->lstep;
  reversed = 0;

  /*------ Loop through family in this curveset */

  if (iqty_in_win==0) any=0;
  for(k=0; k<nf; k++)
    {
      if (qf) qf->current = k;			/* (NULL if lfaml=-1) */
      ns = get_data_pointers(k, k==0, reversed, cp, qs, &ix0, &iy0);
      //printf("Redraw00, get_data_pointers returns ix0,iy0 = %lx, %lx\n", ix0, iy0);

      if (flags & SINGLE)
	{
	  for(i=0, skipMe=1; i<9; i++) 
	    if (cp->icurve[i] == (k + iqty_in_win)) skipMe = 0;

	  if (debug_scale_evis) printf("redraw00, ic=%d, skipMe=%d\n", ic, skipMe);
	  if (skipMe) continue;
	}
      if (k%cp->fskip !=0) continue;

      if (spline==1 && !init_spline(ns,0)) spline=0;
    AGAIN:				/* (Return here 2nd time thru spline)*/
      for(i=0; i<ns; i++)		/*.......... Step through points */
	{
	  next_data(i, qs, ix0, iy0, &fx1u, &fy1u, &bx, &by);	// get data in bx, by
	  if (any==0 || fx1u < xmin) xmin = fx1u;
	  if (any==0 || fx1u > xmax) xmax = fx1u;
	  if (any==0 || fy1u < ymax) ymax = fy1u;
	  if (any==0 || fy1u > ymin) ymin = fy1u;
	  any = 1;
	  if (spline==1)
	    {
	      add_to_spline(i, fx1u, fy1u);	/* WARNING!! was fx1,fy1 */
	      continue;
	    }

	  else if (spline==2) { /*drawline(fx0, fy0, fx1, fy1, moveit, spline, i);*/ }
	}

      if      (spline==1) { spline = 2; fit_spline(); goto AGAIN; }
      else if (spline==2) { spline = 1; zap_spline(); }
    }

  //if (flags & SINGLE && !skipMe) printf("Done Redraw00, y = %f to %f\n", ymin, ymax);
  *xminp = xmin;
  *xmaxp = xmax;
  *yminp = ymin;
  *ymaxp = ymax;
}

/*-----------------------------------------------------------------------------
|	draw_misc_labels
-----------------------------------------------------------------------------*/
#define HEIGHT 16	/* Empirically determined! */
#define SEP 4
void draw_misc_labels(CURVE_SET *cp)
{
  struct LABEL *next;
  int lx1, ly1, lx2, ly2, nline, align, width;
  char *s, *longtext, *p, *q, delim;
  int i, boxed;
  int nitems;
  XTextItem *item;

  if (!ps_modeon)
    XSetForeground(mydisplay, redraw_gc, white());
  else
    {
    ps_setlabelfont(0);
    ps_color(0.83,0.5,0.0);
  }

  i=0;
  for(next=cp->label; next; )
    {
      boxed = next->boxed;		// only draw ONE boxed label

      next = get_label(LBL_WHITE, next, &lx1, &ly1, &s, &nline, &align,
		       xscale, xoffset, yscale, yoffset,&nitems,&item);
      if (s==NULL) break;
      for(i=0, p=s; i<nline; i++, *q=delim, p=q+2)
	{
	  q = get_newline(p);
	  delim = *q; *q = 0;
	  item->chars = p;
	  item->nchars = strlen(item->chars);
	  //printf("draw_misc_labels %d. %d,%d: %s\n", i, lx1,ly1,p);
	  draw_label(lx1, ly1, i, p,item,nitems);/* (properly corrects ly if boxed) */
	}
      if (align <= LBL_ALIGNL)
	{
	  get_labelbox(&width, &longtext);
	  if (!ps_modeon)
	    {
	      lx2 = lx1 + width;
	      lx1 -= SEP; lx2 += SEP;
	      ly2 = ly1 + (nline-1) * (HEIGHT+SEP) + SEP;
	      ly1 -= (HEIGHT+SEP);
	      XDrawLine(mydisplay, redraw_win, redraw_gc, lx1, ly1, lx2, ly1);
	      XDrawLine(mydisplay, redraw_win, redraw_gc, lx2, ly1, lx2, ly2);
	      XDrawLine(mydisplay, redraw_win, redraw_gc, lx2, ly2, lx1, ly2);
	      XDrawLine(mydisplay, redraw_win, redraw_gc, lx1, ly2, lx1, ly1);
	    }
	  else
	    ps_draw_labelbox(lx1, ly1, nline, longtext);
	}
      //if (boxed) break;
    }
}

/*-----------------------------------------------------------------------------
|	draw_label
-----------------------------------------------------------------------------*/
void draw_label(int lx, int ly, int lineno, char *p, 
		XTextItem *item, int items)
{
  int i;
  int width=0;
  int lang;
  if (p == NULL) return;
  if (!ps_modeon)
    {
      if (lineno >= 0) ly += lineno * (HEIGHT+SEP);	/* XWin 0=UL! */
      if(items)
	{
	  for( i=0; i<items; i++)
	    {
	      XSetFont(mydisplay,redraw_gc,(item[i].delta==LATIN)?font_latin
		       :font_greek);
	      XDrawString(mydisplay, redraw_win, redraw_gc, lx, ly,
			  item[i].chars,item[i].nchars );
	      lx+=XTextWidth((item[i].delta==LATIN)?font_struct:
				font_struct1, item[i].chars,item[i].nchars);
	      
	    }
/*	  XDrawText(mydisplay, redraw_win, redraw_gc, lx, ly, item, items);*/
	}
      else
      XDrawString(mydisplay, redraw_win, redraw_gc, lx, ly, p, strlen(p));
    }
  else 
    {
      if(lineno == -3)
	{
	  if(items)
	    {
	      for( i=0; i<items; i++)
		{
		  ps_draw_label(lx,ly, lineno,
				item[i].chars,item[i].nchars,item[i].delta );
		  lineno=-2;
		  
		}
	      
	    }
	  else
	    ps_draw_label(lx, ly, lineno, p,
			  0 , LATIN);
	  
	}
      else
	if (lineno == -4)
	  {
	    if(items)
	      {

		ps_vertical_string(lx,ly,0,
				   item[i].chars,item[i].nchars,item[i].delta );
		
		for( i=0; i<items; i++)
		  {
		    ps_vertical_string(lx,ly,1,
				       item[i].chars,item[i].nchars,item[i].delta );
		  }
		
		
		ps_vertical_string(lx,ly,-1,
				   item[0].chars,item[0].nchars,item[0].delta );
	      
	      }
	    else
	      ps_draw_label(lx, ly, lineno, p,
			    0 , LATIN);
	  
	  }
      
      
    }
  XSetFont(mydisplay,redraw_gc,font_latin);
}

/*-----------------------------------------------------------------------------
|	drawline
|	* i = 1..ns-1
-----------------------------------------------------------------------------*/
void drawline(float fx0, float fy0, float fx1, float fy1,
	      int moveit, int spline, int i)
{
  int ninterval, j;
  int x0, y0, x1, y1, dx, dy, ok;
  float t, fx, fy, ex, ey;
#define SPLSIZE 8

  if (spline)
    {
      dx = float_to_int(fx1 - fx0); if (dx < 0) dx = -dx;
      dy = float_to_int(fy1 - fy1); if (dy < 0) dy = -dy;
      if (dx < dy) dx = dy;
      ninterval = dx/SPLSIZE + 1;	/* 0..7-->0, 8..15-->1, 16..23-->2 */
    }
  else ninterval = 1;

  x1 = float_to_int(fx0);
  y1 = float_to_int(fy0);
  fx = fx0;
  fy = fy0;
  for (j=1; j<=ninterval; j++)
    {
      x0 = x1;
      y0 = y1;
      ex = fx;
      ey = fy;
    LAST:
      if (j == ninterval)
	{
	  fx = fx1; fy = fy1;
	}
      else
	{
	  /*if (j==1) printf("%d: %d intervals in %g %g...%g %g\n",
			   i,ninterval, fx0,fy0,fx1,fy1);*/
	  t = (float)j / (float)ninterval;
	  ok = eval_spline(i-1, t, &fx, &fy, j==1);
	  if (!ok) { j=ninterval; goto LAST; }
	  /*printf("%g: %g %g\n", t, fx, fy);
	  if (j==ninterval-1) printf("\n");*/
	}
      x1 = float_to_int(fx);
      y1 = float_to_int(fy);
      if (!ps_modeon)
	XDrawLine(mydisplay, redraw_win, redraw_gc, x0, y0, x1, y1);
      else
	{
	  if (moveit) ps_fmoveto(ex, ey);
	  ps_flineto(fx, fy);
	  moveit = 0;
	}
    }
}

/*-----------------------------------------------------------------------------
|	drawmarker -- XWindows uses x1,y1, ps uses fx1,fy1
-----------------------------------------------------------------------------*/
void drawmarker(int x1, int y1, float fx1, float fy1)
{
  //printf("drawmarker %c, mode %d: %d %d\n", Mark, ps_modeon, x1, y1);
  if (!ps_modeon)
    {
      if (Mark == 'x')
	{
	  XDrawLine(mydisplay, redraw_win, redraw_gc, x1-DM, y1-DM, x1+DM, y1+DM);
	  XDrawLine(mydisplay, redraw_win, redraw_gc, x1-DM, y1+DM, x1+DM, y1-DM);
	}
      else if (Mark == '#')
	{
	  XDrawLine(mydisplay, redraw_win, redraw_gc, x1-DM, y1-DM, x1+DM, y1-DM);
	  XDrawLine(mydisplay, redraw_win, redraw_gc, x1+DM, y1-DM, x1+DM, y1+DM);
	  XDrawLine(mydisplay, redraw_win, redraw_gc, x1+DM, y1+DM, x1-DM, y1+DM);
	  XDrawLine(mydisplay, redraw_win, redraw_gc, x1-DM, y1+DM, x1-DM, y1-DM);
	}
      else if (Mark == '^')
	{
	  XDrawLine(mydisplay, redraw_win, redraw_gc, x1-DM, y1+DM, x1+DM, y1+DM);
	  XDrawLine(mydisplay, redraw_win, redraw_gc, x1+DM, y1+DM, x1,    y1-DM);
	  XDrawLine(mydisplay, redraw_win, redraw_gc, x1,    y1-DM, x1-DM, y1+DM);
	}
      else
	{
	  XDrawLine(mydisplay, redraw_win, redraw_gc, x1 - DM, y1, x1 + DM, y1);
	  XDrawLine(mydisplay, redraw_win, redraw_gc, x1, y1 - DM, x1, y1 + DM);
	}
    }
  else
    {
      ps_comment("Marker");
      if (Mark=='x')
	{
	  ps_fmoveto(fx1 - (float)DM, fy1 - (float)DM);
	  ps_flineto(fx1 + (float)DM, fy1 + (float)DM);
	  ps_fmoveto(fx1 - (float)DM, fy1 + (float)DM);
	  ps_flineto(fx1 + (float)DM, fy1 - (float)DM);
	}
      else if (Mark=='#')
	{
	  ps_fmoveto(fx1 - (float)DM, fy1 - (float)DM);
	  ps_flineto(fx1 + (float)DM, fy1 - (float)DM);
	  ps_flineto(fx1 + (float)DM, fy1 + (float)DM);
	  ps_flineto(fx1 - (float)DM, fy1 + (float)DM);
	  ps_flineto(fx1 - (float)DM, fy1 - (float)DM);
	}
      else if (Mark=='^')
	{
	  ps_fmoveto(fx1 - (float)DM, fy1 + (float)DM);
	  ps_flineto(fx1 + (float)DM, fy1 + (float)DM);
	  ps_flineto(fx1            , fy1 - (float)DM);
	  ps_flineto(fx1 - (float)DM, fy1 + (float)DM);
	}
      else
	{
	  ps_fmoveto(fx1 - (float) DM, fy1);
	  ps_flineto(fx1 + (float) DM, fy1);
	  ps_fmoveto(fx1, fy1 - (float) DM);
	  ps_flineto(fx1, fy1 + (float) DM);
	}
    }
}

void setmarkerstyle(char *s)
{
  char *p;
  int i;
  for(i=0,p=markerstyle; i<4; i++)
    {
      if (*s!='+' && *s!='x' && *s!='#' && *s!='^') break;
      *p++ = *s++;
    }
  if (i==0) *p++ = '+';
  *p = 0;
}

/*-----------------------------------------------------------------------------
|	read_flabels
|	* call if a variables line contains backslash (\)
|	* returns pointer to next word after labels
-----------------------------------------------------------------------------*/
char *read_flabels(FILE *frd, char *text, int index)
{
   int ok;
   char *p, *q;
   ok = index & LOOPBIT;
   index &= INDBITS;
   ok = (ok && index < nloop);

   if (ok) loop[index].labels = text;
   for(p=text;;)
   {
      for(; (q=strchr(p,'\\'))!=NULL; p=q+1) *q=0;
      *(q=strchr(p,'\n')) = 0;
      p = ok ? q + 1 : text;
      fgets(p, LNAME, frd);
      for(q=p; *q==' ' || *q=='\t'; q++);
      if (*q != '\\') break;
      strcpy(p, q+1);
   }

   return ok ? p : text;
}

/*-----------------------------------------------------------------------------
|	load_single_label
-----------------------------------------------------------------------------*/
int load_single_label(CURVE_SET *cp, char *string)
{
   int i, j;
   char *p;
   p = loop[cp->lfaml].labels;
   if (!p) return 0;

   *string = '\0';
   for(i=0; *p; i++, p+=strlen(p)+1)
   {
     for(j=0; j<9; j++) {
       if (cp->icurve[j] != -1) {
	 if (*string != '\0') strcat(string, ", ");
	 strcat(string, p);
       }
     }
   }

   return (*string == '\0') ?  0 : 1;
}


