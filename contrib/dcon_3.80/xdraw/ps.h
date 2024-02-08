/*-----------------------------------------------------------------------------
|	PS.H
-----------------------------------------------------------------------------*/
extern	void postscript(char *title);
extern	void postscript1(char *title,char *psname, char *psfont);
extern	void ps_dataname(char *s);
extern  void ps_suffix(char *p);
extern  void ps_init(char *ps_filename,char *ps_Ffontname,char data[]);
extern	void ps_color(float, float, float);
extern	void ps_comment(char *);
extern  void ps_line(int x1,int y1,int x2,int y2);
extern  void ps_stroke(void );
extern  void ps_moveto(int x1,int y1);
extern  void ps_lineto(int x1,int y1);
extern  void ps_fmoveto(float x1,float y1);
extern  void ps_flineto(float x1,float y1);
extern  void ps_rectangle(int x1,int y1,int w,int h);
extern  void ps_clip(int x1,int y1,int w,int h);
extern  void ps_unclip(void );
extern  void ps_thinline(void );
extern  void ps_linewidth(int w);
extern  void ps_normalline(void );
extern  void ps_font(int points );
extern	void ps_setlabelfont(int);
static  void goodstring(char *string,int stringlen);

extern void ps_draw_labelbox(int lx1, int ly1, int nline, char *longtext);
extern  void ps_rshowstring(int xpos,int ypos,char *string,int stringlen);
extern  void ps_showstring(int xpos,int ypos,char *string,int stringlen);
extern  void ps_centerstring(int x,int y,int above,char *s,int n);
extern  void ps_centervertically(int x,int y, int right,char *s,int n);
extern	void ps_power(int x_or_y, char *s2);
extern  void ps_string(char *string,int stringlen);
extern  void ps_rightjustify(int xpos,int ypos,char *string,int stringlen);
extern  void ps_translate(int xpos,int ypos);
extern  void ps_showpage(void );
extern  void ps_adjustpage(int height, int bx1, int by1, int bw, int bh);
extern	void ps_scale(int height, int bx1,int by1,int bdx,int bdy);
extern  void ps_close(void );
extern	void ps_save(int save);
extern	void writetrf(FILE *,FILE *,int,int,int,int,int,int);
extern	void ps_layout(int, int, int);
extern	void ps_country(int);
extern	int get_maxps(void);
extern  void ps_draw_label(int , int , int , char *,int ,int );
extern  void ps_fatline(void);
extern  void ps_fill(void);
extern  void ps_grestore(void);
extern  void ps_vertical_begin(int , int  );
extern void ps_vertical_string( int ,int , int ,
			 char *, int , int );


