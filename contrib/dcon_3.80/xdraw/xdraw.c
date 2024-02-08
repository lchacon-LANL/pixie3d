/******************************************************************************
**  NAME      XDRAW.C
**  AUTHOR    Sheryl M. Glasser
**
**  DESCRIPTION
**      Usage: draw [options] [suffix]
**      options: -d = draw (draw.in), -c = cntour (cont.in)
**
**  Copyright (c) GlassWare 1992.  All rights reserved.
******************************************************************************/
//------ Outline ---------------------------------------------------------
//main
//  init()
//    readstate() - read xdraw.ini or state.xdraw: various variables eg datapath
//    get_inname() - read command line args; open drawxx.in, read ftype
//    read_asc_header() - drawxx.in next line is "InvertBytes:", "Outer:", "Loop:"
//    read to start of filenames
//    binread() for each file
//	if M or M2, getlimits_mb
//    read variable names from drawxx.in into varid[]
//    read info for each curve from drawxx.in
//      initialize cp->xx, eg cp->title=NULL, cp->label=NULL. vp=&cp->ix
//	arg_suffix() for x - get cp->flags: read 0..9, then suffix '.', '&', 'L'
//	read vp->index etc (eg via parsevar())
//	arg_suffix() for y
//	read vp->index etc
//	dash options, e.g. -x or -y calls axis_label, sets cp->ix.label
//      get_limits()
//        limits(.., xp) where xp= & cp->ix, &cp->iy, &cp->iz
//          i = var_exists(); xp->label=(varid+i)->name
//          loop through data, get xmin, xmax
//  opendisplay()
//  for each cp in curveset[]...
//    parse_title() - e.g. title contains '>'
//    makewindow()
//    expose_event() (#defined as event())
//	redraw()
//
//  process_event (defined as event())
//	redraw()

//------ Summary of Zoom Calculation --------------------------------------
//redraw ()
//  get_box_info
//    get window (bigbox, clipbox)
//    givebox - set global variables zoomdx, fxbox1, fdxbox
//    get xmin..ymax from data
//    modify xmin..ymax for polar, forced
//    get_zoomed_limits
//       get_zoomedclip         -- if v->clipped, return fx1=v->f[0] etc
//       modify xmin..ymax      -- xmax = xmin + fx2*dx; xmin += fx1*dx;
//  getscale		        -- xoffset, xscale: (window/world), use xmax-xmin
//  draw_frame
//  redraw0 or redraw_mb
//
//------ Summary of Label and Ticks ---------------------------------------
//redraw
//  get_box_info()
//  getscale
//  draw_frame
//  redraw0 etc - draw graph
//  axis_details
//    void drawticks()
//      char *tick_info() .... returns "%.0lf" to "%.3lf"
//    draw_label or XDrawString
//    draw text in lower left corner
//    draw_label or XDrawString for y axis
//-------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <fcntl.h>
#include <time.h>
#include <limits.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <X11/keysymdef.h>

#define XDRAW_SOURCE
#include "gendefs.h"
#include "curves.h"
#include "xdraw.h"
#include "xinit.h"
#include "xcontour.h"
#include "xtools.h"
#include "setcolor.h"
#include "menu.h"
#include "ps.h"

extern void read_topology(CURVE_SET *cp, int iTime, int iTime_last);

static void draw_ten(int, int, int, int, float);
static void drawticks(int is_y,
		float xlim[], float xmax, float ylim[], float ymin,
		FloatSc xscale, FloatSc xoffset, FloatSc yscale, FloatSc yoffset);
static char *tick_info(float *, float *, float *, int *, float *, int, int);
//void convert_to_log(float *xmin, float *xmax, int axis, int iview);
static  void getscale(int wx0, int dwx, float xmin, float xmax, int axis,
                      FloatSc *xscale, FloatSc *xoffset);

float *buf;
int xterm;
int default_ncol = 3, figures = 4;
int default_ncurve = 16;
int exitflag = 0;
int redrawflag;

Window dialog_win;
GC dialog_gc;
#define expose_event() event()
#define process_event() event()
int dialogwindow = 0;		/* XTerm can use normal dialog window */

int tellMe=0;
int debug_scale = 0;
int debug_scale_evis=0;

int ps_only =0;
int background =0;
int axis_width = AXIS_WIDTH;
int curves_width = CURVES_WIDTH;
int ps_axis_width = AXIS_WIDTH;
int ps_curves_width = CURVES_WIDTH;
char title[1000];
int titledrawn=0;

LOOP loop[20];
int nloop=0;
int iloop;

NODE *nodelist;
extern BLOCK_XY *xylist;
extern int ntime;
extern int time_stride, fam_stride;
extern VARIABLE_ID varid[];
extern unsigned int xscreen, yscreen;
extern int  nRow, nCol, maxWindow;

int nnode=0;			/* TOTAL count of nodelist */
int nnode_t=0,nnode_r=0;
int inode = -1, ivar = -1;
int ncount_equal = 0;

CURVE_SET curveset[NCS];	/* Describes the plots */
CURVE_SET *cp, *cp2;
int ncset=0;

//static CURVE_SET *auto_cp1, *auto_cp2;
int drawing_all_windows = 0;

/*LOOP loop[20];*/
/*int nloop, iloop;*/
/*int  iloop;*/

int ftype;			/* 0=graph, 1=contour, 2=self-determined */
int multi_topology;		// if Type M2, 1 then 2; else 0 or 0 then 3
int M_one_timestep;		// if Type M and "one_timestep", 1; else 0

extern char bin_file_name[];
extern int n_m2;
extern int ps_modeon;
extern unsigned int xscreen, yscreen;
extern Font font_latin,font_greek;
#define LATIN 1
#define GREEK 2

#define incre(p) p++
#define FP_inc(p)    p++;
#define FP_incn(p,n) p+=n;
#define millisec(c) c /= (CLOCKS_PER_SEC / 1000)

extern int zoomi;
extern VIEW view[];
float data_limits;
unsigned long nRedPixel;
extern Display *mydisplay;
int multiDots;
extern int DM;
static int lsp, wsp;
float floorf(float x);

/*-----------------------------------------------------------------------------
|    main
-----------------------------------------------------------------------------*/
int main(argc, argv)
int argc;
char *argv[];
{
  CURVE_SET *cp;
  int nwin;
  //int ncol, nwinplus;
  char text[200];
  float x,y;
  XColor rCol;

  //------ Read input data, open XTerminal

  init(argc, argv);
  xterm = opendisplay(title);	// sets xscreen, yscreen nRow, nCol, maxWindow
  if (!xterm) exit(0);

  //------ Count nwin

  for(cp=curveset,nwin=0; cp<cp2; cp++)
    if ((cp->which=='1') && (cp->view_number!=-1))
      {
	if ((nwin + dialogwindow) < maxWindow) nwin++;	
	else { cp2 = cp; break; }
      }

  //menu_setwindows(nwin);	// dummy proc (USE_MENU not defined)
				// .. and there was once a menu.c, now gone

  for (cp = curveset; cp < cp2; cp++)
    {
      if ((cp->which=='1') && (cp->view_number !=-1)) {
	parse_title(cp, cp->title, text);
	cp->window = makewindow(text, cp->view_number);
      }
      else 
	cp->window = (cp-1)->window;

      if (cp<cp2-1 && (cp+1)->which=='2') continue;

      /* Color background prepare for xterm */
      if (xterm) {
	XParseColor (mydisplay,
		     DefaultColormap(mydisplay,0), "Red", &rCol);
	XAllocColor (mydisplay, DefaultColormap(mydisplay,0), &rCol);
	nRedPixel = rCol.pixel;
      }

      if (xterm) expose_event();		/* Draw (Expose event) */
      else postscript(cp->title);		/* no xterm: write ps */
    }
  
  if (dialogwindow)				/* UNIX or Windows */
    dialog_win =  makewindow("XDraw Dialogue", WIN_MNGR);
  else 
    init_menus();				/* Motif  menu (dummy) */

  //------ Process events

  if (!ncset) ;			       	/* ncset=#curves=set in xinit.c */
  else while (!exitflag) process_event();	/* process events, Q= exit */

  if (xterm) closedisplay();
}

/*=============================================================================
**		ARRAY INDEXING UTILITIES
**===========================================================================*/
/*-----------------------------------------------------------------------------
|	get_loopval
-----------------------------------------------------------------------------*/
float get_loopval(LOOP *qs, int i)
{
  float *p;
  p = qs->limVals;
  if (qs->ltype=='I') return((float)i);
  else if (qs->ltype=='X') return(*p + i * *(p+2));
  else if (qs->ltype=='A') return(*(p+i));
  else return((float)0);
}

/*=============================================================================
**                    SUPPORTING FUNCTIONS
**===========================================================================*/
Display *mydisplay;		/* set via redraw()-->get_expose_info() */
Window redraw_win;
int myscreen;
GC redraw_gc;
XRectangle clipbox, bigbox;
float clipx1, clipx2, clipy1, clipy2;	/* for xdraw's own clip of BIG numbers*/
unsigned int width, height;
int font_height, monochrome=0;
int wx0, wy0, dwx=-1, dwy=-1, wx2, wy2;
int tickdy=0, tickdx=0;
int frame_sep = 0, border = 0;
int nqty_in_win, qtys_in_win[20], preserve_aspect;
FloatSc xscale, xoffset, yscale, yoffset;
float xmin_data, xmax_data, ymin_data, ymax_data;
static char ps_msg[] = "Generating PostScript file...";

/*-----------------------------------------------------------------------------
|	get_box_info
-----------------------------------------------------------------------------*/
int get_box_info(Window win, int iview, XRectangle *bbx, XRectangle *cbx,
	         float *xminp, float *yminp, float *xmaxp, float *ymaxp)
{
  int icurve, iqty_in_win, i, i0, j, jcurve;
  int wxb, wya, wyb;
  //int border;
  unsigned int depth;
  float xmin, xmax, ymin, ymax, dp;
  float xmin0, xmax0, ymin0, ymax0;
  CURVE_SET *cp, *cq;
  int zoomed;
  static int pass=0;
  XRectangle draw_area;
  double z;

  //------ get size of window (in pixels = screen coordinates)

  if (!ps_modeon) get_properties(win, &width, &height, &depth);
  else { width = 1264; height=984; depth=4; }		/* condor full screen */
  monochrome = (depth==1);
  //printf("get_box_info: width %d, height %d\n", width, height);
  //printf("get_box_info: bits per pixel = %d\n", depth);

  icurve = i0 = xterm ? get_graph(win, 0) : get_graph((Window)0, pass++);
  if (icurve==-1) return(-1);
  cp = curveset + icurve;

  //------ get interior window where graph will be drawn
  //	   'border' gives top or right dist of inner to outer rectangle
  //	   it never changes even if frame is stretched

  if (border == 0)
    border = (int) (.02 * ((width >= height) ? width : height));
  frame_sep = border * 1;

  wx0 = 4 * font_height;		// left of inner frame
  wxb = width - border;			// right of inner frame

  wya = border + 2*font_height;		// bottom of inner frame rel to bottom of window
  wyb = height - border;		// top of inner frame rel to bottom of window
  if (cp->subtitle && !ps_modeon) 
    wyb -= font_height;
  wy0 = height - wyb;			// top of inner frame rel to top of window

  wx0 += frame_sep;
  wxb -= frame_sep;
  wyb -= frame_sep;
  wya += frame_sep;
  wy0 += frame_sep;

  dwx = wxb - wx0;			// width of inner frame
  dwy = wyb - wya;			// height of inner frame
  wx2 = wxb;
  wy2 = wyb;

  //------ bigbox = screen coord starting at 0
  //       clipbox is ditto for inner rectangle; note wy0 is rel to TOP of window

  bigbox.x = bigbox.y = 0;
  bigbox.width = width;
  bigbox.height = height;

  clipbox.x = wx0;	clipx1 = (float)clipbox.x;
  clipbox.y = wy0;	clipy1 = (float)clipbox.y;
  clipbox.width = dwx;	clipx2 = clipx1 + (float)clipbox.width;
  clipbox.height = dwy;	clipy2 = clipy1 + (float)clipbox.height;

  //printf("Get_Box_Info, bigbox = %d, %d, clipbox = %d, %d\n",
  //width, height, dwx, dwy);

  //------ Set global vars zoomdx,zoomdy,fxbox1,fybox1,fdxbox,fdybox

  givebox(wx0, wy0, dwx, dwy, width, height);

  if (cp->gtype == 'B') return(-1);	// it's never this value?

  //------ Get xmin .. ymax for all curves in redraw_win, in world coords)
  //	   (redraw_win is current window - myevent.xexpose.window)
  //	   (need a note for how cq->which is not '1')

  for(cq=cp,i=icurve; cq<cp2 && (i==i0 || cq->which=='2'); i++,cq++)
    {
      if ( ftype != 6 )		// xmin etc for type G etc
	{
	  qtys_in_win[i-i0] = i;
	  if (i==i0 || cq->x_info.min < xmin) xmin = cq->x_info.min;
	  if (i==i0 || cq->x_info.max > xmax) xmax = cq->x_info.max;
	  if (i==i0 || cq->y_info.min < ymax) ymax = cq->y_info.min;
	  if (i==i0 || cq->y_info.max > ymin) ymin = cq->y_info.max;
	  //printf("get_box_info, limits %g,%g to %g,%g\n", xmin,ymin,xmax,ymax);
	}

      else			// xmin etc for type M etc
	{
	  xmin = xylist->xmin;
	  xmax = xylist->xmax;
	  ymin = xylist->ymax;
	  ymax = xylist->ymin;
	  if (debug_scale) 
	    printf("DS 00 get_box_info, initial xmin..ymax = %g, %g to %g, %g\n", xmin, ymin, xmax, ymax);

	  //--------- Set Clip Area for types M, M2

	  if (xmax!=xmin) dp = (float) .05 *(xmax - xmin);
	  else if (xmin)  dp = xmin * FH;
	  else dp = FH;
	  xmin -= dp;
	  xmax += dp;

	  if (ymax!=ymin) dp = (float) .05 *(ymin - ymax);
	  else if (ymin)  dp = ymin * FH;	/* ymin > ymax */
	  else dp = FH;
	  ymin += dp;
	  ymax -= dp;
	}
    }

  nqty_in_win = i - i0;		// see eg ahg2 (bessel.bin), has 3 cp's in one window

  if (cp->flags & POLAR)
    {
      xmax = ymin =  cp->y_info.max;
      xmin = ymax = -xmax;
    }

  if (debug_scale) 
    printf("DS 01 get_box_info, xmin..ymax = %g, %g to %g, %g\n", xmin, ymax, xmax, ymin);

  if (cp->gtype == 'G')					/* set margins */
    {
      if ((cp->forcemin != FZ || cp->forcemax != FZ) &&
	   cp->forcemin != cp->forcemax && !(cp->flags & E_IGNORE))
        {
	  ymax = cp->forcemin;
	  if (cp->forcemax > cp->forcemin)
	    ymin = cp->forcemax;
        }
      if ((cp->flags & E_VIS))	/* && (cp->ncurve != -1 || cp->hzcurve != -1))*/
	{
	  //for (i=0,cq=cp; cq<cp2 && cq->window==cp->window; i++,cq++) 
	  //printf("%d. which=%c, view=%lx\n", i, cq->which, cq->window);

	  if (debug_scale_evis) 
	    printf("E-VIS is set, nqty_in_win=%d, new limits %g,%g to %g,%g\n", 
		   nqty_in_win, xmin, ymin, xmax, ymax);

	  for(i=0; i<9; i++) {
	    if (cp->icurve[i] != -1) {
	      for(iqty_in_win=0; iqty_in_win<nqty_in_win; iqty_in_win++) {
		jcurve = icurve + iqty_in_win;
		//printf("iview=%d, jcurve=%d\n",iview,jcurve);
		redraw00(jcurve, iqty_in_win, nqty_in_win, &xmin, &xmax, &ymin, &ymax);
		//goto E_VIS_DONE;
	      }
	    }
	  }
	  if (debug_scale_evis) 
	    printf("E-VIS is set, nqty_in_win=%d, new limits %g,%g to %g,%g\n", 
		   nqty_in_win, xmin, ymin, xmax, ymax);
	}
    E_VIS_DONE:

      /*------ Log10 Stuff */

      xmin_data = xmin; xmax_data = xmax;
      ymin_data = ymin; ymax_data = ymax;
      //printf("x=%g to %g, y=%g to %g\n", xmin, xmax, ymax, ymin);

      if (cp->use_log & 1) {
	if (cp->use_cutoff & 1 && xmin < cp->xcutoff) xmin = xmin_data = cp->xcutoff;
	if (xmin <= 0.0) {
	  printf("Invalid use of Log X in view %d (X=%g to %g)%c\n", 
		 iview, xmin, xmax, '\007');
	  cp->use_log &= 2;
	}
	else {
	  xmin = (float)log10((double)xmin);
	  xmax = (float)log10((double)xmax);
	}
      }

      if (cp->use_log & 2) {
	if (cp->use_cutoff & 2 && ymax < cp->ycutoff) ymax = ymax_data = cp->ycutoff;
	if (ymin <= 0.0 || ymax<=0.0) {
	  printf("Invalid use of Log Y in view %d (Y=%g to %g)%c\n", iview, ymin, ymax, '\007');
	  cp->use_log &= 1;
	}
	else {
	  ymin = (float)log10((double)ymin);
	  ymax = (float)log10((double)ymax);
	}
      }

      //printf("x=%g to %g, y=%g to %g\n", xmin, xmax, ymax, ymin);

      //------ Still type G - Get dp = margin inside the clipbox

      if (xmax!=xmin) dp = (float) .05 *(xmax - xmin);
      else if (xmin)  dp = xmin * FH;
      else dp = FH;

      //if ((fabs((xmax-xmin))<1.0e-05) && (xmax!=xmin))// obsolete test?
      //dp = (float)(0.05*fabs((xmax-xmin)*xmax));	// we still use for y, below

      if (debug_scale) 
	printf("DS 02 get_box_info, xmin..ymax = %g, %g to %g, %g\n", xmin, ymax, xmax, ymin);

      dp = fabs(dp);
      xmin -= dp;
      xmax += dp;

      if (ymax!=ymin) dp = (float) .05 *(ymin - ymax);
      else if (ymin)  dp = ymin * FH;			/* ymin > ymax */
      else dp = FH;

      if ((fabs((ymin-ymax))<1.0e-05) && (ymin != ymax))
	{dp=0.05*fabs((ymax-ymin)); }

      dp = fabs(dp);
      ymin += dp;
      ymax -= dp;
    }			// done with cp->gtype = 'G'

  /*------ Get new xmin..ymax if zoomed or preserve aspect */

  if (debug_scale) 
    printf("DS 03 get_box_info, xmin..ymax = %g, %g to %g, %g\n", xmin, ymax, xmax, ymin);

  preserve_aspect = cp->flags & ASPECT;
  xmin0 = xmin; xmax0 = xmax;
  ymin0 = ymin; ymax0 = ymax;
  
  zoomed = get_zoomed_limits(cp, &xmin, &ymin, &xmax, &ymax) ;
  //printf("preserve_aspect = %d, zoomed=%d\n", preserve_aspect, zoomed);

  if (debug_scale) 
    printf("DS 04 get_box_info, xmin..ymax = %g, %g to %g, %g, aspect=%d\n", 
	   xmin, ymax, xmax, ymin, preserve_aspect);

  if (preserve_aspect && zoomed != 1) 
    get_aspected_limits(win, dwx, dwy,
	&xmin, &ymin, &xmax, &ymax, xmin0, ymin0, xmax0, ymax0);

  *xminp = xmin; *xmaxp = xmax;
  *yminp = ymin; *ymaxp = ymax;

  //------ if called from show_coord(), bbx and cbx are local XRectangles

  //printf("Get_Box_Info, final xmin..ymax = %g, %g to %g, %g\n", xmin, ymin, xmax, ymax);
  if (bbx != &bigbox)  *bbx = bigbox;
  if (cbx != &clipbox) *cbx = clipbox;
  if (debug_scale) 
    printf("DS 05 get_box_info, xmin..ymax = %g, %g to %g, %g\n", xmin, ymax, xmax, ymin);
  return(icurve);
}

/*-----------------------------------------------------------------------------
|	get_aspected_limits
|	* Only called when zoom flag NOT set
|	* Calculates aspected limits, sets new zoom borders
-----------------------------------------------------------------------------*/
void get_aspected_limits(Window win, int w, int h,
			 float *xmin, float *ymin, float *xmax, float *ymax,
			 float xmin1, float ymin1, float xmax1, float ymax1)
{
  float xc, yc, dx, dy;
  float aspect_d, aspect_w;

  dx = *xmax - *xmin;
  dy = *ymin - *ymax;
  aspect_w = (float)h / (float)w;		/* aspect of box */
  aspect_d = dy / dx;				/* aspect of data */
  if (aspect_d > aspect_w)			/* data is tall */
    {
      xc = (*xmax + *xmin) / 2;
      dx = dy / aspect_w;
      *xmin = xc - dx/2;
      *xmax = xc + dx/2;
    }
  if (aspect_d < aspect_w)			/* data is wide */
    {
      yc = (*ymax + *ymin) / 2;
      dy = dx * aspect_w;
      *ymax = yc - dy/2;
      *ymin = yc + dy/2;
    }
  set_aspected_limits(win, 0, xmin1, xmax1, *xmin, *xmax);
  set_aspected_limits(win, 1, ymin1, ymax1, *ymin, *ymax);
  //printf("get_aspected_limits: x = %f to %f\n", *xmin, *xmax);
}

/*-----------------------------------------------------------------------------
|	get_zoomed_limits
|	* IN:  data's autoscaled unzoomed pointers to xmin..ymax (plus dp)
|	* OUT: adjusted by fractional zoom window fx1..fy2
-----------------------------------------------------------------------------*/
int get_zoomed_limits(CURVE_SET *cp, float *xmin, float *ymin, float *xmax, float *ymax)
{
  float fx1, fy1, fx2, fy2, dx, dy;
  int retval=0, iClipped;
  if (iClipped = get_zoomedclip(cp->window, &fx1, &fy1, &fx2, &fy2))
    {
      dx = *xmax - *xmin;
      dy = *ymin - *ymax;
      *xmax = *xmin + fx2 * dx;      *xmin = *xmin + fx1 * dx;
      *ymax = *ymin - fy1 * dy;      *ymin = *ymin - fy2 * dy;
      //printf("get_zoomed_limits: %d\n", iClipped);
      retval = iClipped;
    }
  return retval;
}

/*-----------------------------------------------------------------------------
|	draw_frame
|	* draw inner rectangle and axes.  Note axes are offset from
|	  inner rectangle by ixoff, iyoff (xoffset, yoffset)
-----------------------------------------------------------------------------*/
void draw_frame()
{
  int ixoff, iyoff;
  float ixmin, ixmax, iymin, iymax;
  int wx1, wy1, w1, h1;
  
  if (!ps_modeon)
    {
      clear_if_needed( redraw_win, redraw_gc, bigbox.width, bigbox.height);

      XSetClipRectangles(mydisplay, redraw_gc, 0, 0, &bigbox, 1, Unsorted);
      wx1 = wx0 - frame_sep; 	// frame is separated from graph...
      wy1 = wy0 - frame_sep; 	// ...this is useful if zoomed
      w1 = dwx + 2*frame_sep;
      h1 = dwy + 2*frame_sep;
      //printf("bigbox:  %2d %2d %d %d\n", bigbox.x, bigbox.y, bigbox.width, bigbox.height);
      //printf("clipbox: %2d %2d %d %d\n", clipbox.x, clipbox.y, clipbox.width, clipbox.height);
      //printf("wx1 etc: %2d %2d %d %d\n", wx1, wy1, w1, h1);
      //printf("\n");

      //XDrawRectangle(mydisplay, redraw_win, redraw_gc, wx0, wy0, w, h);
      XDrawRectangle(mydisplay, redraw_win, redraw_gc, wx1, wy1, w1, h1);
      XSetClipRectangles(mydisplay, redraw_gc, 0, 0, &clipbox, 1, Unsorted);
    }
  else if (ps_modeon)
    {
      //xmessage(wx0 + font_height, wy0 + font_height, ps_msg);
      ps_adjustpage(height, wx0, wy0, dwx, dwy);
      ps_color(FZ, FZ, FZ);
      ps_clip(0, 0, 0, 0);
      ps_rectangle(wx0, wy0, dwx, dwy);
      ps_clip(clipbox.x, clipbox.y, clipbox.width, clipbox.height);

    }

  //------ Type G vertical axis line

  if (cp->gtype != 'G') goto ENDLINE;

  ixmin = (clipx1 - xoffset) / xscale;
  ixmax = (clipx2 - xoffset) / xscale;
  iymin = (clipy2 - yoffset) / yscale;
  iymax = (clipy1 - yoffset) / yscale;

  ixoff = (int) xoffset;
  iyoff = (int) yoffset;
  if ((float)0 >= ixmin && (float)0 <= ixmax)
    {
      if (!ps_modeon)
	//XDrawLine(mydisplay, redraw_win, redraw_gc, 
	//ixoff, wy0, ixoff, wy0 + h);
	XDrawLine(mydisplay, redraw_win, redraw_gc, 
		  ixoff, bigbox.y, ixoff, bigbox.height);
      else {
        ps_linewidth(ps_axis_width);
	ps_line(ixoff, wy0, ixoff, wy0 + dwy);
      }
    }

  //------ Type G horizontal axis line

  if ((float)0 >= iymin && (float)0 <= iymax)
    {
      if (!ps_modeon)
	XDrawLine(mydisplay, redraw_win, redraw_gc, 
		  wx0, iyoff, wx0 + dwx, iyoff);
      else
	ps_line(wx0, iyoff, wx0 + dwx, iyoff);
    }

ENDLINE:
  if (ps_modeon)
    ps_thinline();
}

/*-----------------------------------------------------------------------------
|	redraw. redraws all contents of window on expose events.
-----------------------------------------------------------------------------*/
void redraw()
{
  int iview, iResult, i, j;
  float xmin, xmax, ymin, ymax, xlim2[2], ylim2[2];
  char *p;
  VIEW *v;
  XGCValues rValues;
  int w_Prev, h_Prev;

  /*------ get properties of window, min/max of data in window */
  /*       xmin..ymax are min, max of data plus or minus dp */
  /*       but if log, they are min, max of log10 of data plus or minus dp */

  w_Prev = dwx; h_Prev = dwy;	// need this if window changed and retain aspect
  get_expose_info(&mydisplay, &redraw_win, &redraw_gc, &font_height);

  iview = getview(redraw_win);
  cp = curveset + iview;
  v  = view + iview;
  //printf("** redraw, iview=%d\n", iview);
  //printf("redraw 01, extrema %g, %g\n",cp->iz.min, cp->iz.max);
  //Note extrema for ftype=6 are in BLOCK_F *flist at fmin, fmax

  /*------ get_box_info gets box info AND loads qtys_in_win[], which is indexes for curveset */

  iResult = get_box_info(redraw_win, iview, &bigbox, &clipbox, &xmin, &ymin, &xmax, &ymax);
  if (iResult == -1) return;
  //printf("redraw: w %d --> %d \n", w_Prev, w);

  //convert_to_log(&xmin, &xmax, 1, iview);
  //convert_to_log(&ymin, &ymax, 2, iview);

  //printf("Redraw %d, limits are %.10g,%.10g to %.10g,%.10g\n", 
  //iview, xmin, ymin, xmax, ymax);

  /*------ Data limits */

  //printf("\nredraw, old  xmin..ymax = %.10g to %.10g, %.10g to %.10g\n",
  //v->xmin, v->xmax, v->ymin, v->ymax);
  //printf("redraw, xmin..ymax = %.10g to %.10g, %.10g to %.10g\n",
  //xmin, xmax, ymin, ymax);

  if (iResult == -2) {		/* can occur if DATA_LIMITS test is enabled */
    xmin = v->old_xmin; xmax = v->old_xmax;
    ymin = v->old_ymin; ymax = v->old_ymax;
  }

  v->old_xmin = v->xmin; v->xmin = xmin;
  v->old_ymin = v->ymin; v->ymin = ymin;
  v->old_xmax = v->xmax; v->xmax = xmax;
  v->old_ymax = v->ymax; v->ymax = ymax;

  for(i=0; i<4; i++) {
    v->old_f[i] =  v->f[i];
    v->old_f1[i] = v->f1[i];
  }

  //if (debug_scale) 
  //printf("DS redraw, xmin..ymax = %.10g, %.10g to %.10g, %.10g\n", xmin, xmax, ymin, ymax);
  //printf("w=%d\n", w);		// w is set in get_box_info

  /*------ Get scale factors and offsets */
 REDRAW_AGAIN:
  getscale(wx0, dwx, xmin, xmax, 1, &xscale, &xoffset);
  getscale(wy0, dwy, ymin, ymax, 2, &yscale, &yoffset);

  xoffset += FH;	/* roundoff */
  yoffset += FH;

  /* Axis line width setting */ 
 
  if(!ps_modeon) {
    rValues.line_width = cp->axis_width;
    XChangeGC(mydisplay, redraw_gc, GCLineWidth, &rValues);
  }

  else {		/* yes ps_modeon */
    ps_curves_width = cp->curves_width;
    ps_axis_width = cp->axis_width;
  }

  draw_frame();		/* draw frame, horizontal & vertical lines */

  /* Axis line width setting */ 

  if(!ps_modeon) {
    rValues.line_width = cp->curves_width;
    XChangeGC(mydisplay, redraw_gc, GCLineWidth, &rValues);
    }
  else {
    ps_linewidth(ps_curves_width);
    }

  /*------ Draw Graph: call redraw0, redraw1, or redraw_mb (which later calls redraw1) */

  if (cp->gtype=='G')				//.....Draw Graph
    {
      for(i=0; i<nqty_in_win; i++) {			// nqty_in_win>1 if line in draw.in has several iy's
	//printf("%d. call redraw0(%d)\n",i, qtys_in_win[i]);
	iResult = redraw0(qtys_in_win[i], i, nqty_in_win, xmin, xmax, ymin, ymax);

	/*---- restore if nan encountered - probably no longer occurs */

	if (iResult == 1) {
	  xmin = v->xmin = v->old_xmin;
	  ymin = v->ymin = v->old_ymin;
	  xmax = v->xmax = v->old_xmax;
	  ymax = v->ymax = v->old_ymax;
	  for(i=0; i<4; i++) {
	    v->f[i] = v->old_f[i];
	    v->f1[i] = v->old_f1[i];
	  }
	  XSetForeground(mydisplay, redraw_gc, white());
	  goto REDRAW_AGAIN;
	}
      }
    }

  else if (cp->gtype == 'P')			//.....Dos Test Palette; Unix dummy
    testpalette(wx0, wy0, dwx, dwy);

  else if (ftype == 6) {                        //.....Multiblock Contour
    if (tellMe) xprintf("**redraw Multiblock Contour, i=%d\n", i);
    //printf("Redraw, view has %g, %g to %g, %g\n", v->xmin, v->ymin, v->xmax, v->ymax);
    //printf("Redraw, call redraw_mb using %g, %g to %g, %g\n", xmin, ymin, xmax, ymax);
    redraw_mb(cp,xmin,xmax,ymin,ymax);
  }

  else  {					//.....Contour (and all other ftypes)
    if (tellMe) xprintf("**redraw Contour\n");
    redraw1(cp, xmin, xmax, ymin, ymax);
  }

  //------ Reset clip rectangle for tick marks, axis labels etc

  if (!ps_modeon)
    XSetClipRectangles(mydisplay, redraw_gc, 0, 0, &bigbox, 1, Unsorted);
  else {
    ps_color(FZ, FZ, FZ);
    ps_clip(0, 0, 0, 0);
  }

  //------ draw tick marks on axes
  //       warning! xlim may be changed in axis_details()
  
  xlim2[0] = xmin;
  xlim2[1] = xmax;

  ylim2[0] = ymax;
  ylim2[1] = ymin;

  /*------ Draw labels.
    Each curveset describes one curve; each view might use several.
    Array qtys_in_win[] lists the indices of each curveset for this view.
    It's loaded in get_box_info
  */

  if (ftype != 6) cp = curveset + qtys_in_win[0];
  axis_details(xlim2, xscale, xoffset,	ylim2, yscale, yoffset, cp);

  //if (ps_modeon)
  //   xmessage(wx0 + font_height, wy0 + font_height, ps_msg);
}		// ... end redraw()

#ifdef DEAD_CODE
/*-----------------------------------------------------------------------------
|    convert_to_log
-----------------------------------------------------------------------------*/
void convert_to_log(float *xmin, float *xmax, int axis, int iview)
{
  char c;
  double xxmin, xxmax, x;
  c = (axis == 1) ? 'X' : 'Y';

  if (cp->use_log & axis) {
    xxmin = *xmin;
    xxmax = *xmax;
    x =  (xxmin <= xxmax) ? xxmin : xxmax;
    if (x <= 0.0) {
      printf("Invalid use of Log %c in view %d (Xmin=%lg)\n", c, iview, x);
      cp->use_log &= (3-axis);
    }
    else {
      *xmin = log10(xxmin);
      *xmax = log10(xxmax);
      printf("Log %c: %lg --> %g, %lg --> %g\n",
	     c, xxmin, *xmin, xxmax, *xmax);
    }
  }
}
#endif

/*-----------------------------------------------------------------------------
|    getscale
-----------------------------------------------------------------------------*/
void getscale(int x0, int dx, float xmin, float xmax, int axis,
	      FloatSc *xscale, FloatSc *xoffset)
{
  double xxscale, xxoffset;
  double xxmax, xxmin;
  char c;
  xxmax = xmax;
  xxmin = xmin;
  //---- Log conversion now done in get_box_info
  //if (cp->use_log & axis) { xxmax = log10(xxmax); xxmin = log10(xxmin); }

  xxscale  = (double)dx / (xxmax - xxmin);
  xxoffset = (double)x0 - xxscale * xxmin;
  *xscale  = xxscale;
  *xoffset = xxoffset;
}

/*-----------------------------------------------------------------------------
|    get_graph
-----------------------------------------------------------------------------*/
int get_graph(Window mywindow, int i0)
{
  CURVE_SET *cp;
  int i, n;
  if (mywindow)
    for (i=i0, cp=curveset+i; cp < cp2 && cp->window != mywindow; i++, cp++);

  else
    {
      for (i=0,n=-1,cp=curveset; cp < cp2; i++, cp++)
	{
	  if (cp->which=='1') n++;	// '1' means first in the window
	  if (n==i0) break;
	}
      //i = n;	// NN 99
    }
  return ((cp == cp2) ? -1 : i);
}

#ifdef DEAD_CODE
//***** get N (total) cp for the View **************************************** 

int get_graph1(int view, int i0)
{
  CURVE_SET *cp;
  int i, n;
  if (mywindow)
    for (i=i0, cp=curveset+i; cp < cp2 && cp->window != mywindow; i++, cp++);
  else
    {
      for (i=0,n=-1,cp=curveset; cp < cp2; i++, cp++)
	{
	  if (cp->which=='1') n++;
	  if (n==i0) break;
	}
      i = n;
    }
  return ((cp == cp2) ? -1 : i);
}
#endif

/*-----------------------------------------------------------------------------
|	toggle_markersize
-----------------------------------------------------------------------------*/
void toggle_markersize(CURVE_SET *cp1)
{
  DM++;
  if (DM>8) DM=2;
  redrawflag = 1;
}

/*-----------------------------------------------------------------------------
|	toggle_markers
-----------------------------------------------------------------------------*/
void toggle_markers(CURVE_SET *cp1, int *mark)
{
  CURVE_SET *cp;
  if ((cp1->gtype != 'G') && (cp1->gtype != 'M') && (cp1->gtype != 'V')) return; 

  for(cp=cp1; cp<cp2; cp++)
    {
      if (cp > cp1 && cp->which=='1') break;
      if (cp->flags & MARKERS)
	{ cp->flags &= (~MARKERS); *mark = 0; }
      else
	{ cp->flags |= MARKERS; *mark = 1; }
    }
  redrawflag = 1;
}

/*-----------------------------------------------------------------------------
|	toggle_aspect
-----------------------------------------------------------------------------*/
void toggle_aspect(CURVE_SET *cp1, int *aspect)
{
  CURVE_SET *cp;

  zoomi = getview(cp1->window);
  redrawflag=0;
  for(cp=cp1; cp<cp2; cp++)
    {
      if (cp > cp1 && cp->which == '1') break;

      //----- only toggle unzoomed window
      if (get_zoomedclip(cp->window, NULL,NULL,NULL,NULL) & 1)
	continue;

      //----- aspect off: need set v->clipped to 0
      if (cp->flags & ASPECT) {
	cp->flags &= (~ASPECT); *aspect = 0; 
	newclip(-1, 0, 0);			// -1 means "original" clip
	redrawflag=1; 
      }

      //----- aspect on
      else
	{ cp->flags |= ASPECT; *aspect = 1; redrawflag=1;}
    }
}

/*-----------------------------------------------------------------------------
|	toggle_extrema_use
-----------------------------------------------------------------------------*/
void toggle_extrema_use(CURVE_SET *cp1)
{
  CURVE_SET *cp;
  int single, force, eflag;
  if (cp1->gtype != 'G') return;

  redrawflag=0;
  for(cp=cp1; cp<cp2; cp++)
    {
      if (cp > cp1 && cp->which=='1') break;	// break if new window

      /*----- only toggle if unzoomed ---------*/

      if (get_zoomedclip(cp->window, NULL,NULL,NULL,NULL) & 1)
	continue;
      single = cp->flags & SINGLE;
      force = (cp->forcemin != cp->forcemax);
      if (single || force)
	{
	  redrawflag = 1;
	  eflag = cp->flags & (E_VIS | E_IGNORE);
	  if (single && force)
	    {
	      eflag += E_VIS;
	      if (eflag > E_IGNORE) eflag = 0;
            }
	  else if (single) eflag = E_VIS - eflag;
          else if (force)  eflag = E_IGNORE - eflag;
	  cp->flags = (cp->flags & ~(E_VIS | E_IGNORE)) | eflag;
	  //printf ("toggle_extrema_use, ic = %d, flags = %x\n", (int)(cp-curveset), cp->flags);
	}
    }
}

/*-----------------------------------------------------------------------------
|	toggle_block
|	* f_block = first block + 1; l_block = last block + 1
-----------------------------------------------------------------------------*/
void toggle_block(CURVE_SET *cp, char c)
{
   if ( ftype == 6 )
     {
      redrawflag=1;
      if (c == 'b') {
	if(!cp->f_block) { cp->f_block = cp->l_block = 1; }
	else if (cp->f_block == (nnode)) { cp->f_block = cp->l_block = 1; }
	else { cp->f_block += 1; cp->l_block = cp->f_block;} 
        }

     else if (c=='B') 
       cp->f_block=0;
     }
}

/*-----------------------------------------------------------------------------
|	toggle_fill
-----------------------------------------------------------------------------*/
void toggle_fill(CURVE_SET *cp1)
{
  CURVE_SET *cp;
  for (cp=cp1; cp<cp2; cp++)
    {
      if (cp > cp1 && cp->which=='1') break;
      if (!(cp->flags & DOTS))
	{
	  if (cp->flags & FILL) cp->flags &= ~FILL;
	  else cp->flags |= FILL;
	  redrawflag=1;
	  //printf("TOGGLE %d\n", cp->flags);
	}
    }
}

/*-----------------------------------------------------------------------------
|	toggle_single
|	* c: > < . , / 0 1 ... 9
-----------------------------------------------------------------------------*/
void toggle_single(CURVE_SET *cp1, char c)
{        
  int i, n, ifam, nfam, nt, *which, n_already, ipass;
  int nTime, iTime, iTime_last, initializing;
  CURVE_SET *cp;
  char text[100];

  initializing = 0;
  nTime = ntime;
  if (multi_topology) {		// M, M2 always advance ALL windows
    nTime = n_m2;
    if (c == '/') { c = '.'; initializing = 1; }
    if (c == '>') c = '.';
    if (c == '<') c = ',';
  }

  for(cp=cp1; cp<cp2; cp++)
    {
      iTime = cp->itime;
      if (multi_topology) {
	iTime = cp->itime_abs;
	if (cp->gtype == 'V' && c >= '0' && c <= '9') return;
      }

      if (cp > cp1 && cp->which=='1') break;

      //------ Clear SINGLE flag

      if (c == '0') cp->flags &= ~SINGLE;		/* redraw all curves */

      //------ Draw One, Two .. Nine Curves

      else if (c >= '1' && c <= '9')
	{
	  /*------ Type V or M */

	  if ((cp->gtype == 'V') || (cp->gtype == 'M')  ) {
	    nt = cp->ncurve;
	    redrawflag = 1;

	    if ( !(cp->flags & SINGLE )) {	/* first try  */
	      cp->flags |= SINGLE;
	      cp->i_single = 0;
	      cp->j_single =-1;
	      XClearArea(mydisplay, cp->window, 0, 0, 0, 0, True);
	    }

	    else if (c == '1') {		/* next try */
	      cp->i_single +=1;
	      if (cp->i_single  == nt) cp->i_single = 0;
	    }

	    else if (c == '2') {		/* 2nd coutour (for more, must use cp->icurve[] */
	      cp->j_single +=1;
	      if (cp->j_single  == nt) cp->j_single = 0;
	    }
	  }

	  /*------ Type G or C; note that fam_stride can be negative */

	  else {
	    for(i=0, n_already=0; i<9; i++) {
	      if (!(cp->flags & SINGLE))   cp->icurve[i] = -1;	/* 1st of 1..9: reset all */
	      else if (cp->icurve[i] >= 0) n_already++;		/* count how many already there */
	    }

	    ifam = (int)(c) - (int)('1');			/* 0..8 */
	    nfam = (nqty_in_win>1) ? nqty_in_win : (loop+cp->lfaml)->count;

	    if (cp->icurve[ifam] < 0) {
	      if (n_already >= nfam) return;
	      cp->icurve[ifam] = -fam_stride;
	    }

	    //if (c == '1' && cp->ncurve  < 0) cp->ncurve  = -fam_stride;
	    //if (c == '2' && cp->hzcurve < 0) cp->hzcurve = -fam_stride;

	    //printf("toggle_single %c for ifam=%d, nfam=%d, n_already=%d: %d\n", 
	    //c, ifam, nfam, n_already, cp->icurve[ifam]);
	    cp->flags |= SINGLE;

	    which = &cp->icurve[ifam];
	    if (*which < 0 || *which >= nfam) *which = -fam_stride;
	    ipass = 0;

	  ANOTHER_STEP:
	    *which += fam_stride;
	    //printf("step: %d = %d\n", *which, cp->icurve[ifam]);
	    if (*which >= nfam)	 *which -= nfam;
	    else if (*which < 0) *which += nfam;
	    ipass++;

	    for(i=0; i<9; i++) {
	      if (i != ifam && cp->icurve[i] == *which) {
		if (ipass < 200) goto ANOTHER_STEP;
		else printf("Error - can't find appropriate step\n");
	      }
	    }
	    //printf("which=%d\n", *which);
	  }							/* ... end of type G or C */
	}							/* ... end of 1 .. 9 */

      /*------ Advance Timestep for Focus Window */

      else if ((c == '<' || c == '>') && nloop > 3)
	{
          if(ftype != 6) {
	    for(i=0; i<nloop && (n=cp->param[i])<0 ; i++) ;
	    if (i < nloop) {
	      nt = loop[i].count;
	      if (c == '>') n++; else n--;
	      if (n < 0 ) n = nt - 1;
	      else if (n >= nt) n = 0;
	      cp->param[i] = n;
	      get_limits(cp);		// for advance timestep < or >
	      if (parse_title(cp, cp->title, text))
		XChangeProperty(mydisplay, redraw_win, XA_WM_NAME, XA_WM_NAME,
				8, PropModeReplace, (unsigned char *)text, 1);
	    }
	  }					/* end of !=6 */

	  else {
	    iTime_last = iTime; 
	    if (nTime > 1 && c == '>')
	      iTime = (iTime + time_stride) % nTime;

	    else if (nTime > 1 && c == '<') {
	      nt = iTime - time_stride;
	      iTime = ( nt<0 ) ? (nTime + nt): nt;
	    } 

	    cp->itime = cp->itime_abs = iTime;
	    //printf("toggle_single, cp %d, time %d\n", cp-curveset, cp->itime);
	    //if (multi_topology) read_topology(cp, iTime, iTime_last);
	  }					/* end of ==6 */
	}

      /*------ Advance Timestep for all windows, type M */

      else if (ftype == 6 && ( c == ',' || c == '.') && nTime > 1)
	{
	  iTime_last = iTime; 
	  if (c == '.') {			/* forward */
	    iTime = (iTime + time_stride) % nTime;
	  }

	  else if (c == ',') {			/* backward */
	    iTime = iTime - time_stride;
	    if (iTime < 0) iTime += nTime;
	  }

	  nt = iTime;

	  /*---- Loop through all windows - 'True' triggers expose event */

	  if (!initializing)
	    XSync(mydisplay, True);		/* clear event buffer */
	  drawing_all_windows = 1;

	  for(cp=curveset; cp<cp2; cp++) {
	    cp->itime = cp->itime_abs = nt;
	    if (multi_topology) 		// set itime_rel; call binread if needed
	      read_topology(cp, iTime, iTime_last);
	    if (!initializing)
	      XClearArea(mydisplay, cp->window, 0, 0, 0, 0, True);  // redraw
	  }      
	  redrawflag=0;
	}				// ...end of ftype 6, dot/comma, ntime>1
    }

  /*------ Dot or comma if not type M */

  if(ftype != 6 || (ftype ==6 &&(c != '.'|| c !=',')))
    redrawflag=1;
}

#ifdef DEAD_CODE
/*-----------------------------------------------------------------------------
-----------------------------------------------------------------------------*/
void trigger_next_window()
{
  int nt;
  CURVE_SET *cp;
  cp = auto_cp1; nt = cp->itime;
  cp++; cp->itime = nt;

  auto_cp1++;
  if (auto_cp1 == auto_cp2) drawing_all_windows = 0;
  XClearArea(mydisplay, cp->window, 0, 0, 0, 0, True);
}
#endif

/*-----------------------------------------------------------------------------
|	toggle_flabel
|	* toggle labels for all curves in this window
-----------------------------------------------------------------------------*/
void toggle_flabel(CURVE_SET *cp1)
{
  Window win;
  win = cp1->window;
  for (cp=cp1; cp<cp2 && cp->window==win; cp++) {
    if (cp->flags & LABELF) cp->flags &= ~LABELF;
    else cp->flags |= LABELF;
  }
  redrawflag = 1;
}

/*-----------------------------------------------------------------------------
|	draw_graphtype
-----------------------------------------------------------------------------*/
void draw_graphtype(char *type)
{
  char *p;
  for (cp = curveset; cp < cp2; cp++)
    {
      for (p = type; *p && cp->gtype != *p; p++);
      if (!*p)
	continue;
      XClearArea(mydisplay, cp->window, 0, 0, 0, 0, True);
    }
}

/*-----------------------------------------------------------------------------
|	initialize_vectors
-----------------------------------------------------------------------------*/
void initialize_vectors(CURVE_SET *cp1)
{
  CURVE_SET *cp;
  if (cp1)					// from init() (called AFTER binread)
    initialize_vectors_1(cp1);

  else if (cp2) {				// from binread() (but 1st time cp2 is NULL)
    for(cp=curveset; cp<cp2; cp++) 
      initialize_vectors_1(cp);			// will only act if cp->gtype = 'V'
  }
}

/*=============================================================================
**                LABELS
**===========================================================================*/
/*-----------------------------------------------------------------------------
|	tick_info. computes positions of tick marks
|	* IN:  lim[] = min, max of data.
|	* OUT: new[] = integer start, end of ticks (may be outside)
|	*      majorz = separation, iexp=exponent
-----------------------------------------------------------------------------*/
static char *tick_info(float *lim, float *new, float *dx, 
		   int *iexp, float *offset, int is_y, int use_log)
{
  int i, n, ii, ni, n1, n2;
  int mt=8, m4=4, m8=8, m6=6;
  float ratio, min, max, amin, amax;
  double xmin, xmax, x;
  double t, tl, ti, diff, cutoff, expt, expfac;
  int pass, tellMe, screen_len, window_len;
  static double interval[] = {
    .001, .002, .005, .01,  .02,  .05, .1, .2, .5,
    1.,  2.,  5.,  10., 20., 50., 100., 200., 500., 1000.
    };
  static char *format[] = { "%.3lf", "%.2lf", "%.1lf", "%.0lf" } ;
  float data_lim;
  float tt1, tt2, tt;
  char sText[200];	// for debug
  char axis;

  axis = is_y? 'Y':'X';
  screen_len = is_y? yscreen : xscreen;
  window_len = is_y? bigbox.height : bigbox.width;

#define NI 19		/* WARN! keep as sizeof(interval) */

  tellMe = 0;
  //if (!is_y) tellMe = 1;

  *offset = FZ;
  if (figures < 2) figures = 2;
  if (figures > 5) figures = 5;
  cutoff = 2.5 * pow(10., -(double)figures);
  xmin = lim[0];
  xmax = lim[1];

  //------ find the greater of abs(xmin), abs(xmax); find data limit

  tt1 = (xmin>=0.0) ? xmin : -xmin;
  tt2 = (xmax>=0.0) ? xmax : -xmax;
  tt  = (tt2 < tt1) ? tt1 : tt2;
  tt1 = fabs(xmax - xmin) / 2.0;
  while (tt+tt1 !=tt) tt1 /= 2.0;
  data_lim = tt1*2;

  if (tellMe) printf("************** Ticks %c: Limits %g to %g, data_lim %g\n", axis, xmin, xmax, data_lim);
  
  if(fabs(lim[0]-lim[1]) < data_lim) return (NULL);

  //------ Obtain appropriate limits
  
  pass = 0;
START:  
  pass++;
  if (pass > 10) { 
    printf("Crash in 'ticks' for %c (%.8g to %.8g)\n", axis, lim[0], lim[1]); 
    return(NULL); 
  }

  amin = (float) fabs(xmin);
  amax = (float) fabs(xmax);

  //    if (tellMe && (xmax-xmin)<6) 		// DEBUG!! set breakpoint
  //    tellMe=0;
  //  else if (tellMe)
  //    tellMe = tellMe;

  if (tellMe) {
    printf("pass %d, abs limits %g and %g\n", pass, amin, amax);
    sprintf(sText, "  ");
  }

  if (amin != FZ && amax != FZ)			/* replace small endpoints with 0 */
    {
      ratio = amax / amin;			/* e.g. 1..22, 2..44 etc. --> */
      if (ratio > (float)22) xmin = FZ;		/*      0..22, 0..44 */
      if (ratio < (float)(1./22.)) xmax = FZ;
    }

  //---- find the exponent of the largest limit

  t = (double)((amax > amin) ? amax: amin);	/* compute exponent */
  tl = log10(t);				/* e.g. 80-->1.9 ,.008-->-2.1*/
  expt = floor(tl);				/* e.g. 1.9-->1,  -2.1-->-3*/
  i = (int)expt;				/* exponent! */
  if (i > 0) {
    if (i <= 2) i = 0;				/* reduce exponent */
    //else i -= 2;				/* (can give 6000 = 600 x 10^1) */
    expt = (double)i;
  }

  *iexp = i;
  expfac = pow(10., -expt);			/* e.g. 1-->.1,   -3-->1000*/

  if (tellMe) {
    sprintf(sText, "%s max=%g, log10(max)=%g --> %g, xmin,xmax= %g %g ", 
	    sText, t, tl, expt, xmin, xmax);
    printf("%s\n", sText);
  }

  diff = (float)expfac * (xmax - xmin);		/* see if need offset */
  if (tellMe) printf("   diff=%g, cutoff=%g, expfac=%g\n", diff, cutoff, expfac);

  //------ Calculate offset if needed

  if (diff < cutoff) // && !use_log)
    {
      if (tellMe) printf("Obtain an offset, diff=%g vs %g\n", diff, cutoff);
   NEWOFF:	    
      tl = log10(diff);
      expt = floor(tl);
      expfac = pow(10., -expt);
      t = (double)(xmax+xmin) * expfac / 2.;
      if (tellMe) printf("   tl, expt, expfac, t= %g %g %g %g\n",tl, expt, expfac, t);
      for(;;)
	{
	  if(fabs(t)<data_lim) return NULL;
	  if (t>0. && t<1.) { t *= 10.; expfac *= 10.; }
	  else if (t<0. && t>-1.) { t *= 10.; expfac *= 10.; }
	  else if (t>1000.)  { t /= 10.; expfac /= 10.; }
	  else if (t<-1000.) { t /= 10.; expfac /= 10.; }
	  else break;
	}
      //*offset = (float)(floor(t) / expfac);
      *offset = (float)( t / expfac);
      xmin -= *offset;
      xmax -= *offset;
      goto START;
    }

  //------ Find appropriate interval (1,2,5 etc) for desired significant figures

  min = xmin * (float)expfac;			/* get min > -9.99 and... */
  max = xmax * (float)expfac;			/* ... max <  9.99 */
  if (tellMe) printf("min, max using %lg: %lg --> %g, %lg --> %g\n", expfac, xmin, min, xmax, max);

  i = (figures > 3) ? 0 : 3;
  ratio = (float)window_len / screen_len;	// larger windows can have more tick marks

  //------ Change max# tickmarks from values set in declaration
  //	   m4 - m8 ticks for 1, 10, 100 etc
  //	   m4 - m6 ticks for 0.1, 0.5

  if (ratio > 0.45) {
    mt = 25;
    m4 = 10; m6 = 16; m8 = 21;
  }
  else if (ratio > 0.35) {
    mt = 21;
    m4 = 8; m6 = 10; m8 = 19;
  }
  else if (ratio > 0.25) {
    mt = 15;
    m4 = 6; m6 = 8; m8 = 14;			// 4-->6, 6-->8, 8-->14
  }

  //------ DON'T USE RATIO - same # tickmarks despite axis size

  mt = 8; m4 = 4; m8 = 8; m6 = 6;

  for(ni=ii=0; i<NI; i++)
    {
      modf((double)min/interval[i], &ti); n1 = (int)ti;
      modf((double)max/interval[i], &ti); n2 = (int)ti;
      if (n1==0 && n2==0) break;
      n = n2 - n1;
      //if (!is_y && i>6) 
      //printf("i=%d, n=%d, ni=%d, ii=%d (%.3f)\n", i, n, ni, ii, interval[ii]);

      if (ni==0) ;			// 1st pass: allow anything
      else if (ni>=mt && n<=mt) ;	// ni ticks is too many, n would be ok
      else if (i>=9 && n>=m4 && n<=m8);	// 9 means tick step=1: 4-8 ticks would be ok
      else if (i>=6 && i<9 &&		// 6 means tick step=0.1 ...
	       n>=m4 && n<=m6) ;	// ... 4 to 6 ticks is perfect
      else continue;			/* skip next statements */
      ni = n; ii = i;
    }

  modf((double)min/interval[ii], &ti); n1 = (int)ti;
  modf((double)max/interval[ii], &ti); n2 = (int)ti;
  if (tellMe) printf("modf returns ii=%d, interval=%lg, n1=%d, n2=%d, ratio=%.3f\n",
		     ii, interval[ii], n1, n2, ratio);

  //------ Obtain limits based on value near interval

  new[0] = (float)(n1-1) * (float)(interval[ii] / expfac) ;
  new[1] = (float)(n2+1) * (float)(interval[ii] / expfac) ;
  *dx = (float)(interval[ii] / expfac);
  if (*dx < data_lim) return NULL;

  //------ If number of tick marks not reasonable, make another pass
  //       for eg min/max 7.3 to 9.2, "nearby" ticks are 6.5, 7, 7.5, ..., 8.5, 9, 9.5
  //       the ones that count are 7.5, 8, 8.5, 9 (4 tick marks)

  for(x=new[0]+*offset; x<=new[1]+*offset; x+=*dx)
      if (x >= lim[0] && x <= lim[1]) n++;

  if (tellMe) printf("number ticks=%d (minus 2)\n", n);
  if (n <= 2) goto NEWOFF;
  ii /= 3;
  if (ii > 3) ii = 3;
  return format[ii];
}


/*-----------------------------------------------------------------------------
|   label. places tick marks and numbers on graphs
|   * from redraw
-----------------------------------------------------------------------------*/
void axis_details(float *xlim, FloatSc xscale, FloatSc xoffset,
	   float *ylim, FloatSc yscale, FloatSc yoffset, CURVE_SET *cp)
{
  float xleft, ytop;
  char string[100], *sp, *xname, *yname;
  char text[20];
  int x0, y0, xc, yc, i, j;
  int nitems, iTime;
  XTextItem *item;

  xname = cp->x_info.label;
  yname = cp->y_info.label;
  if (ps_modeon)
    {
      ps_save(1);
      ps_normalline();
      ps_setlabelfont(0);
    }

  /*----------- draw and label tick marks */

  if (!ps_modeon)
    XSetForeground(mydisplay, redraw_gc, white());
  xleft = xlim[0];
  ytop  = ylim[1];

  drawticks(0, xlim, xlim[1], ylim, ylim[0], xscale, xoffset, yscale, yoffset);
  drawticks(1, ylim, ytop,    xlim, xleft,   yscale, yoffset, xscale, xoffset);

  if (ps_modeon) ps_save(0);

  //------------------- draw title on x-axis
  //			Note yscale, yoffset refer to the INNER BOX

  if (ps_modeon) ps_setlabelfont(1);
  sp = xname;
  lsp = strlen(sp);
  xc = (int) (xscale * (xlim[0] + xlim[1]) / 2. + xoffset);
  //y0 = (int) (yscale * ylim[0] + 2.0 * font_height + yoffset);	// font_height=12 always
  y0 = bigbox.height - font_height / 2;
  wsp = textwidth(sp, lsp);
  x0 = xc - wsp / 2;

  if (!ps_modeon)
    {
      if(!cp->x_info.nitems)
	XDrawString(mydisplay, redraw_win, redraw_gc, x0, y0, sp, lsp);

      else
	draw_label(x0, y0, -1, sp, cp->x_info.ml_text, cp->x_info.nitems);

      /*------ Lower left, e.g. "time=.." or "curve# " */

      x0= 10;

      if (ftype == 0 && cp->flags & SINGLE) {
	*text = '\0';
	for(i=0; i<9; i++) {
	  if (cp->icurve[i] >= 0) {
	    if (*text == '\0') sprintf(text, "curve# ");
	    else strcat(text, ", ");
	    sprintf(text+strlen(text), "%d", cp->icurve[i]);
	  }
	}

	lsp = strlen(text);
	if (lsp)
	  XDrawString(mydisplay, redraw_win, redraw_gc, x0, y0, text, lsp);
      }

      if ( ftype == 6 ) {
	iTime = cp->itime;
	//printf("label, cp=%d, time %d\n", cp-curveset, iTime);
	if (multi_topology==2) iTime = cp->itime_abs;
	if (!cp->f_block) 
	  sprintf(text,"%s = %d ", varid->name, iTime);
	else
	  sprintf(text,"%s = %d  Block = %d ", varid->name, iTime, cp->i_block);

	lsp = strlen(text);
	XDrawString(mydisplay, redraw_win, redraw_gc, x0, y0, text, lsp);
      }          
    }

  else					/* yes Postscript */
    if(cp->x_info.nitems) 
      draw_label(x0, y0, -3, sp, cp->x_info.ml_text, cp->x_info.nitems);
    else
      ps_centerstring(xc, y0, 0, sp, lsp);

  /*--------------- draw title on y-axis */

  sp = yname;
  nitems = cp->y_info.nitems;
  item = cp->y_info.ml_text;
  if(!nitems)  lsp = strlen(sp);
  else for (j = 0,lsp=0; j < nitems; j++)
    lsp += item[j].nchars;

  x0 = 10;
  yc = (int) (yscale * (ylim[0] + ylim[1]) / 2.0 + yoffset);
  y0 = yc - lsp * font_height / 2;
  if (y0 < font_height)
    y0 = font_height;
  if (!ps_modeon)
      if(!cp->y_info.nitems)
	  for (i = 0; i < lsp; i++)
	    {
	      wsp = textwidth(sp, 1);		/* width of 1st char */
	      XDrawString(mydisplay, redraw_win, redraw_gc, x0 - wsp / 2, y0, sp++, 1);
	      y0 += font_height;
	    }
      else
	{
	  for (j = 0; j < nitems; j++)
	    {
		  
	      XSetFont(mydisplay,redraw_gc,(item[j].delta==LATIN)?
		       font_latin  :font_greek);
	      sp=item[j].chars;
	      for (i = 0; i < item[j].nchars; i++)
		{
		  wsp = textwidth(sp, 1);		/* width of 1st char */
		  XDrawString(mydisplay, redraw_win, redraw_gc,
				  x0 - wsp / 2, y0, sp++, 1);
		  y0 += font_height;
		}
	    }
	}
  else
    {
    if(cp->y_info.nitems) 
      draw_label(0, yc, -4, sp, cp->y_info.ml_text, cp->y_info.nitems);
    else
      ps_centervertically(0, yc, 1, sp, lsp);
    }

  /*------ draw subtitle */

  if (ps_modeon) return;		/* subtitle for ps: see event(), 'p' */
  if (!parse_subtitle(cp, string)) return;

  sp = string;
  lsp = strlen(sp);
  /* xc = ... already set */
  y0 = font_height;			/* Looks good on condor */
  if (!ps_modeon)
    {
      wsp = textwidth(sp, lsp);
      x0 = xc - wsp / 2;
      XDrawString(mydisplay, redraw_win, redraw_gc, x0, y0, sp, lsp);
    }
  else
    ps_centerstring(xc, y0, 0, sp, lsp);
}

/*-----------------------------------------------------------------------------
|	x_is_simple_power
-----------------------------------------------------------------------------*/
int x_is_simple_power(float x)
{
  double zz, z_int, z_frac, zzfrac;
  zz = (double)x;
  z_int  = floor(zz);
  z_frac = zz - z_int;	// e.g. x=8.5, z_frac=0.5
  zzfrac =  (z_frac <= 0.5) ? z_frac : 1.0 - z_frac;
  //printf("Simple: %.10lg --> %.10lg, %.10lg, zzfrac=%lg\n", zz, z_int, z_frac, zzfrac);
  return (zzfrac < 1.0e-5);
}

/*-----------------------------------------------------------------------------
|	drawticks
-----------------------------------------------------------------------------*/
static void drawticks(int is_y,
	float xlim[], float xmax, float ylim[], float ymin,
       	FloatSc xscale, FloatSc xoffset, FloatSc yscale, FloatSc yoffset)
{
  float xnew[2], offlx, x, x1, dx, dy, xnext;
  float zlim[2], znew[2], dz, offlz, z, basePwr, zpos;
  float xdiff, xfloor, xfrac, xfrac2, xpow;
  int xv, yv, xv1, yv1, xv2, yv2, del, xc, yc, xvmid, yvmid;
  int i, xlast;
  int iexpx, iexpz, new_exp, use_log, axismask;
  char string[40], *format, axisname;
  static float log5 = (float)0;
  double zz, z_frac, z_int, zval, zzp, zbasis, fbetween[10];
  int ibetween, nbetween, ntick, first_tick_drawn, cycle, tell_tick, tell_tick1;

  //------ Log10 Stuff.  Note - I tried lots of different schemes to get
  //       ticks working, code may contain some obsolete variables

  tell_tick  = 0;	// 0=no, 1=x, 2=y, 3=x and y
  tell_tick1 = 0;
  if ((tell_tick & 1) && is_y)  tell_tick = 0;
  if ((tell_tick & 2) && !is_y) tell_tick = 0;

  axismask = is_y ? 2 : 1;
  axisname = is_y ? 'Y' : 'X';
  use_log = cp->use_log & axismask ? 1 : 0;
  basePwr = (float)0;
  zlim[0] = xlim[0]; 
  zlim[1] = xlim[1];
  xdiff = xmax - zlim[0];	// if log, xdiff is # cycles

  if (use_log) {
    xfloor = floorf(zlim[0]);
    xfrac  = zlim[0] - xfloor;
    xfrac2 = (float)pow(10.0, (double)xfrac);
    if (tell_tick1) printf("Drawticks %c, x = %g to %g, diff = %g, floor = %g, frac = %g\n", 
	   axisname, xlim[0], xlim[1], xdiff, xfloor, xfrac2);
    xfrac2 = floorf(xfrac2);
    if (log5==0.0) log5 = (float)log10(5.0);	// 0.69897
    if (floor(xlim[0]) == floor(xlim[1])) {
      if (tell_tick) printf("Log %c Drawticks, one cycle, lim %g, %g\n", 
	    axisname, xlim[0], xlim[1]);
    }
  }

  /*------ Get format for tick labels, e.g. ".0lf" */

  format = tick_info(zlim, xnew, &dx, &iexpx, &offlx, is_y, use_log);
  if (format == NULL) return;

  if (tell_tick) {
    printf("Drawticks %c (use  %s), exp %d, offset %f\n",
	   axisname, format, iexpx, offlx);
    printf("Drawticks %c from %g to %g incr by %g in %g to %g (%g)\n",
	 axisname, xnew[0], xnew[1], dx, zlim[0], xmax, zlim[1]);
    if (use_log)
      printf("Drawticks %c from %lg to %lg in %lg to %lg\n",
	   axisname, pow(10.,(double)xnew[0]), pow(10.,(double)xnew[1]),
	   pow(10.,(double)zlim[0]), pow(10.,(double)xmax));
    }

  //------ Log axis has no more than one cycle, must start "below" the minimum

  if (use_log && floor(xlim[0]) == floor(xlim[1])) {
    iexpx = floor(xlim[0]);
    basePwr = iexpx;
    use_log = 2;
    xnew[0] -=1.0;
    if (tell_tick) printf("Drawticks %c, basePwr = %g\n", axisname, basePwr);
  }

  //------ Initialize for drawing ticks

  dy =(ylim[1] - ylim[0]) * .02;		/* height of ticks in y... */
  del = (int) (yscale * dy + FH);		/* ...and in screen coords */
  /*yv =(int) (yscale * ymin + yoffset);*//* bottom of ticks, screen coord*/
  //yv = is_y ? wx0 : wy0 + dwy;

  if (is_y) { yv = wx0;       if (tickdx==0) tickdx = del; else del = tickdx; }
  else      { yv = wy0 + dwy; if (tickdy==0) tickdy = del; else del = tickdy; }

  /*------ Preliminary scan of tick marks: how many? */

  for (ntick=0, x1 = xnew[0]; x1 <= xnew[1]; x1 += dx) {
    if (fabs((double)(x1/dx)) < .001) x1 = (float)0;
    x = x1 + offlx;
    //printf("Counting, x=%g in %g to %g or %g\n", x, xnew[0], xmax, xnew[1]);
    if (x >= xlim[0] && x <= xmax) {
      if (!use_log || x_is_simple_power(x)) ntick++;
    }
  }

#ifdef DATA_INFO
  xdraw log:  x = 9.33141 to 10.419, floor = 9, frac = 0.331411 = 2.14492
  xdraw log1: x = 6.83285 to 10.5101, floor = 6, frac = 0.832853 = 6.80539
  xdraw log2: x = -5.50129 to -1.08369, floor = -6, frac = 0.498713 = 3.15292
#endif

  first_tick_drawn = 0;
  nbetween = 0;

  if (use_log) {
    if (xdiff > 5.0) nbetween = 0;
    else if (xdiff > 4.5) {
      nbetween = 1; fbetween[0] = 5.;
    }
    else if (xdiff > 3.5) {
      nbetween = 2; fbetween[0] = 2.; fbetween[1] = 5.;
    }
    else if (xdiff > 2.0) {
      nbetween = 3; fbetween[0] = 2.; fbetween[1] = 3.; fbetween[2] = 5.;
    }
    else if (xdiff > 0.5) {
      xnew[1] += 1.0; xmax += 1.0;
      nbetween = 8;
      fbetween[0] = .05;
      for (i=1; i<=5; i++) fbetween[i] = (float)i / (float)10;
      fbetween[6] = .7;
      fbetween[7] = .8;
    }

#ifdef DEAD_CODE
  if (use_log) {
    if (dx > 1.0001) nbetween = 0;
    else if (ntick>=5)      { nbetween = 1; fbetween[0] = 5.; }
    //else if (ntick>=2) { 
    //nbetween = 2; fbetween[0] = 2.; fbetween[1] = 5.; 
    //}
    else if (ntick>=2) { 
      nbetween = 3; fbetween[0] = 2.; fbetween[1] = 3.; fbetween[2] = 5.;
    }
    else if (ntick < 2 && xdiff >= 0.5) {
      //xnew[0] -= 1.0;
      xnew[1] += 1.0;
      xmax += 1.0;
      //nbetween = 9; for (i=1; i<=9; i++) fbetween[i-1] = (float)i / (float) 10;
      nbetween = 8;
      fbetween[0] = .05;
      for (i=1; i<=5; i++) fbetween[i] = (float)i / (float)10;
      fbetween[6] = .7;
      fbetween[7] = .8;
    }
#endif

    xv = (int)fabs((double)xscale);	// sep in pixels, e.g. from 7 to 8 (10**7 to 10**8)

    if (tell_tick) {
      //printf("ntick=%d, nbetween = %d, xdiff = %g, offlx=%g, max = %g\n", 
      //ntick, nbetween, xdiff, offlx, xnew[1]);
      printf("*****Drawticks, sep pixels = %d\n", xv);
    }
  }

  /*------ Scan tick marks for drawing */

  cycle = -1;
  for (x1 = xnew[0]; x1 <= xnew[1]; x1 += dx)
    {
      if (fabs((double)(x1/dx)) < .001) x1 = (float)0;
      x = x1 + offlx;
      ibetween = -1;
      //printf("RESET\n");

    ANOTHER_TICK:
      if (use_log) {
	if (ibetween == -1 && xdiff >= 0.5 && nbetween > 0) {
	  if (!x_is_simple_power(x)) continue;
	  zbasis = (double)x ;	//pow(10.0, (double)x);
	  ibetween = 0;
	}
	//printf("x=%g in %g to %g=%g\n", x, xlim[0], xmax, zlim[1]);
	if (x < xlim[0]) goto NEXT_BETWEEN_TICK;
	if (x > zlim[1]) goto NEXT_BETWEEN_TICK;
      }

      if (x >= xlim[0] && x <= xmax)
        {
	  //------ Position of the tick mark

	  xv = (int)(xscale * x + xoffset);
	  if (use_log && ibetween == 0) {
	    //printf("TICK LOCATION AT %d for x=%g, ibetween=%d\n", xv, x, ibetween);
	  }

	  //------ Label of the tick mark

	  if (use_log) {
	    if (first_tick_drawn == 0) {
	      iexpx = floor(x) + 1;
	      if (xdiff < 0.5) iexpx -= 3;	// e.g. 209 instead of 0.209
	      basePwr = (float)iexpx;
	      first_tick_drawn++;
	    }
	    z = (float)pow(10.0, (double)(x - basePwr));
	    sprintf(string, "%.3g", z);
	    lsp = strlen(string);
	    if (strstr(string, "e+03")==NULL) ;
	    else if (lsp==5 || (lsp==6 && z<0))
	      sprintf(string, "%g", z);
	  }
	  
	  else {
	    sprintf(string, format, (x-offlx) * pow(10., (double) (-iexpx)));
	  }

	  lsp = strlen(string);
	  wsp = textwidth(string, lsp);

	  //------ Coordinates of the tick mark - draw xv1,yv1 to xv2,yv2

	  if (!is_y) {				// x-axis tick
	    xv1 = xv2 = xv;
	    yv1 = yv + frame_sep; 
	    yv2 = yv1 + del;
	    xc = xv - wsp / 2;			// location of label of tick
	    if (xc + wsp > wx2) xc = wx2 - wsp;
	    yc = yv + frame_sep + font_height;
	  }
	  else {				// y-axis tick
	    yv1 = yv2 = xv;
	    xv1 = yv - frame_sep; 
	    xv2 = xv1 + del;
	    xc = xv1 - (wsp + 5);
	    yc = yv1 + (int)(.3 * font_height);
	  }

	  /*------ Draw tick mark */

	  if (!ps_modeon)
	    {
	      XDrawLine(mydisplay, redraw_win, redraw_gc, xv1, yv1, xv2, yv2);
	      XDrawString(mydisplay, redraw_win, redraw_gc, xc,yc, string, lsp);
	      if (tell_tick) printf("Tick %s\n", string);
	      if (use_log && ibetween==7) ibetween++;	/* .7 or .8 but not both */

	      /*------ Log10: next tick "between" */

	    NEXT_BETWEEN_TICK:
	      if (use_log && ibetween >= 0 && ibetween < nbetween) {
		z_frac = log10(fbetween[ibetween++]);
		x = (float)(zbasis + z_frac);
		if (tell_tick) printf("zbasis=%lg, z_frac=%lg,  x=%g, between=%d of %d, min=%g\n", 
		       zbasis, z_frac, x, ibetween, nbetween, xlim[0]);
		goto ANOTHER_TICK;
	      }

	    }

	  else
	    {
	      ps_line(xv1, yv1, xv2, yv2);
	      if (!is_y) ps_centerstring(xv1, yv1, -1, string, lsp);
	      else   ps_centervertically(xv1, yv1, -1, string, lsp);
	    }
        }

      //------ Tick marks between

      if (ibetween>=0) {
	//if (ntick<=3 &&
	//ibetween_log++;
      }
    }

  /*--------- write exponents on x- and y- axis */

  if (iexpx != 0 || offlx != FZ)
    {
      xv = is_y ? 5 : (int) (xscale * xlim[1] + xoffset);
      if (is_y)
        yv = 2.5 * font_height;
      else	      
        yv = (int) (yscale * ylim[0] + yoffset + 2.25*font_height);/* was 2.0*/
      draw_ten(is_y+1, xv, yv, iexpx, offlx);
    }
}

/*-----------------------------------------------------------------------------
|   draw_ten -- draw the power of ten on the axis
|	WARN! for ps, should let IT calculate line width wsp
|	and perhaps want text vertical
-----------------------------------------------------------------------------*/
static void draw_ten(int x_or_y, int xv, int yv, int iexp, float off)
{
  int yvh, n;
  char string[80], *p, *p2;
  
  p = string;
  if (iexp)
    {
      sprintf(p, "x10%d",iexp); p += strlen(p); p2 = p;
    }
  if (off)
    {
      sprintf(p,"+ %g", off);
    }
  n = textwidth(string, strlen(string));
  if (x_or_y == 1) xv -= n;

  if (!ps_modeon)
    {
      XSetFont(mydisplay,redraw_gc,font_latin);
     
      p = string;
      
      yvh = yv - font_height / 2;
      
      if (!strncmp(p,"x10",3))
        {
	  XDrawString(mydisplay, redraw_win, redraw_gc, xv, yv, p, 3);
          n = textwidth(p, 3);    xv += n; p += 3;
	  XDrawString(mydisplay, redraw_win, redraw_gc, xv, yvh, p, p2-p);
	  n = textwidth(p, p2-p); xv += n; p = p2;
        }
      if (off)
	XDrawString(mydisplay, redraw_win, redraw_gc, xv, yv, p,strlen(p));
    }
  else if (iexp || off)
    ps_power(x_or_y, string);
}

float floorf(float x)
{
  float y;
  y=(float)(int)x;
  if(y < 0)y=y-1;
  return(y);
}
