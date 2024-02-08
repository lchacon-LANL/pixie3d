//**************************************************************************
//* NAME:	tell.c
//* AUTHOR	Sheryl M. Glasser -- LogiCreativity
//*
//* DESCRIPTION
//* 	Usage e.g.: 'tell big.bin', 'tell big.bin m2', 'tell big.bin m2 w'
//*	Mode: 0=create tell.out, 1=print block delimiters,
//*	      2=create big.bin2 (eg)
//*	see line "type = .." to determine if m or g
//**************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

#ifdef comment
   call addloop_mb for rblocks if ntry<nnode_r and rec is odd
   call addloop_mb for tblocks if ntry>nnode_r and 1st of triangles after rblocks
#endif

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
int write_big(fd, q, nblock, ndata)
	int fd, *q, nblock; long ndata;
{
  if (nblock >= 4 && nblock <=10 ) { }	// skip first timestep
  else {
    if (nblock==2) *(q+3) = *(q+3) - 1;
    write(fd, q-1, (ndata+2L)*sizeof(int));
  }
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
int write_m8(fd, q, nblock, ndata)
     int fd, *q, nblock; long ndata;
{
  long temp[20];
  ndata = *(q-1)/sizeof(long);			// ?? ndata=2 showed up as 17
  if (nblock==1 ) {
    temp[0] = temp[4] = 12L;
    temp[1] = 1L; temp[2] = 0L; temp[3] = 1L;	// rblocks, tblocks, nqty
    write(fd, temp, (3L+2L)*sizeof(int));
  }
  else if (nblock==2) {
    temp[0] = temp[8] = 28L;
    temp[1] = temp[5] = temp[6] = temp[7] = 0L;
    temp[2] = *q; temp[3] = *(q+1); temp[4] = 2L;
    write(fd, temp, (7L+2L)*sizeof(int));
  }
  else if (nblock==4 || nblock==5) {
    temp[0] = temp[2] = 4L;
    temp[1] = 1L;
    write(fd, temp, (1L+2L)*sizeof(int));
    write(fd, q-1, (ndata+2L)*sizeof(int));

    if (nblock==5) {
      temp[0] = temp[8] = 28L;
      temp[1] = temp[2] = temp[3] = temp[4] = 0L;
      temp[5] = 6L; temp[6] = 8L; temp[7] = 2L;
      write(fd, temp, (7L+2L)*sizeof(int));
    }
  }
      
  else {
    write(fd, q-1, (ndata+2L)*sizeof(int));
  }
}


//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
int main(argc, argv)
     int argc;
     char *argv[];
{
  int fd, *q, nseq, nblock, nper_line, mode, nchar, ntime;
  long mr, mz, ir, iz;
  long i,j, ndata;
  unsigned int request, invert_bytes, subtype;
  unsigned long bufSize;
  unsigned char *p, temp[4], b1, b2, b3, b4, type, fname[80];
  float *f;
  FILE *outf;

  //---- Initialize Type (Hard-coded)

  //type = 'm'; subtype = 2; 	// Type M2
  //type = 'm'; subtype = 0;	// Type M
  type = 'g'; subtype = 0;	// Type G

  mode = 0;

  if (mode==0) {
    outf = fopen("tell.out", "w");
    if (outf == NULL) {
      printf("Can't open output file\n"); exit(0);
    }
  }

  nper_line=8;

  *fname = '\0';
  if (argc > 2) {
    for(i=2; i< argc; i++) {
      if (strcmp(argv[i],"w")==0) strcpy(fname,"/win/xdraw/"); 	// prefix /win/xdraw
      if (strcmp(argv[i],"g")==0) { type = 'g'; subtype = 0; }
      if (strcmp(argv[i],"g")==0) { type = 'g'; subtype = 0; }
      if (strcmp(argv[i],"c")==0) { type = 'c'; subtype = 0; }
      if (strcmp(argv[i],"m")==0) { type = 'm'; subtype = 0; }
      if (strcmp(argv[i],"m2")==0){ type = 'm'; subtype = 2; }
      if (strcmp(argv[i],"0")==0) { mode = 0; }
      if (strcmp(argv[i],"1")==0) { mode = 1; }
      if (strcmp(argv[i],"2")==0) { mode = 2; }
    }
  }
  strcat(fname,argv[1]);
  printf("Input file %s\n",fname);

  if (mode==0) {
    fprintf(outf, "TELL.OUT for file %s\n", argv[1]);
    fprintf(outf, "args: ",argc);
    for (i=0; i<argc; i++) {
      fprintf(outf, "%d=%s", i, argv[i]);
      if (i < (argc-1)) fprintf(outf,", ");
    }
    fprintf(outf,"\n");
  }

  fd = open(fname, O_RDONLY, 0);
  if (fd == -1) {
    printf("Can't open file %s: %d\n", fname, errno);
    exit(0);
  }

  bufSize = lseek(fd, 0L, SEEK_END);
  if (mode==0) fprintf(outf, "Buffer size %ld\n\n", bufSize);
  printf("Buffer size %ld\n", bufSize);

  /*------ Allocate buffer, read file into it */

  float *buf0, *buf;
  buf0 = (float*)malloc(bufSize);
  buf = (float *)buf0;
  lseek(fd, 0L, SEEK_SET);
  p = (unsigned char *)(buf);

  read(fd, p, bufSize);
  close(fd);

  if (mode==2) {
    strcat(fname,"2");
    fd = open(fname, O_CREAT | O_RDWR, S_IRWXU|S_IRWXG|S_IRWXO);
    if (fd==-1) {
      fd = open(fname, O_TRUNC | O_RDWR);
      if (fd==-1) {
	printf("Problem opening output file %s\n", fname);
	exit(0);
      }
    }
  }

  /*------ Invert bytes in entire buffer if needed */

  i = (type=='g')? 0 : 1;

  request = *(unsigned int *)(buf+i);
  invert_bytes = (request & 0xffff0000) != 0;
  printf ("invert=%d\n", invert_bytes);
  q = (int *)buf;

  printf("First int in buf is %d = %x hex\n", *q, *q);
  printf("Test int (at buf+%d) is %d = %x hex\n", i,  request, request);
  printf("Size of unsigned int is %d, of int is %d\n", 
	 sizeof(unsigned int), sizeof(int));	// size is 4!!

  printf("invert_bytes? %c\n\n", (invert_bytes==0)?'N':'Y');

  if (invert_bytes) {
    for(i=0; i<bufSize; i+=sizeof(int),p+=sizeof(int)) {
      b1 = *(p+0);
      b2 = *(p+1);
      b3 = *(p+2);
      b4 = *(p+3);
      *(p+3) = b1;
      *(p+2) = b2;
      *(p+1) = b3;
      *(p+0) = b4;
    }
  }

  //------ Initialize variables

  p = (unsigned char *)(buf);
  q = (int *)buf;
  j = 0;
  nseq = 0;
  nblock = 0;
  mr = mz = 0;
  ntime = 0;
  //printf("1st int is %d\n", *q);

  //------ Beginning of next block

NEXT:
  ndata = *q++ / sizeof(int);
  j += sizeof(int);
  printf("Block %d, ndata %d\n", nblock, ndata);
  nblock++;

  if (mode==0) {
    fprintf(outf, "**%d(%d):", ndata, nblock);
    if (type != 'g') fprintf(outf, "\n");
  }
  if (ndata == 0) nseq++;

  //------ Write modified version of xx.bin

  if (mode==2) {
    //write_big(fd, q, nblock, ndata);	// for big.bin
    write_m8(fd, q, nblock, ndata);	// for m8.bin
    q += ndata;
    f = (float *)q;
    j += ndata * sizeof(int);
    goto ENDBLOCK;
  }

  /*------ Output data from buffer */

  f = (float *)q;
  ir = -1; iz = 0;
  for (i=0; i<ndata; i++,f++,j+=sizeof(int))  {

    //---- line breaks for row end in block (C or M)

    if (mz > 0) {
      if ( (type=='m' || type=='c') && nblock>2) {
	ir++;
	if (ir>mr) { ir = 0; iz++; if (mode==0) fprintf(outf, "\n"); }
      }
    }

    //------ M header records

    if (type == 'm' && nblock <= 2) {
      q = (int *)f;
      if (mode==0) fprintf(outf, "%12ld ", *q);

      if (subtype==2) {
	if (nblock==2 && i==1) mr = *q;
	if (nblock==2 && i==2) mz = *q;
      }

      else {
	if (nblock==2 && i==0) mr = *q;
	if (nblock==2 && i==1) mz = *q;
      }
    }

    //------ M topology records?

    else if (type=='m' && (ndata==1 || ndata==7 || ndata==2)) {	// CRUDE implementation of type M2
      if (ir==0) if (mode==0) fprintf(outf, "%4d:", iz);
      q = (int *)f;
      if (mode==0) fprintf(outf, "%12ld ", *q);
    }

    //------ C header records

    else if (type=='c' && nblock <= 2) {
      if (nblock==1) {
	if (i==0) mr = *q;
	if (i==1) mz = *q;
	//if (ir==0) if (mode==0) fprintf(outf, "%4d:", iz);
	q = (int *)f;
	if (mode==0) fprintf(outf, "%12ld ", *q);
      }
      else if (mode==0) 
	fprintf(outf, "%12.5g ", *f);
	
    }

    //---- Data records

    else {
      if (ir==0) if (mode==0) fprintf(outf, "%4d:", iz);
      if (mode==0) fprintf(outf, "%12.5g ", *f);
      //fprintf(outf, "%22.10g ", *f);
    }
  }

 ENDBLOCK:
  q = (int *)f;
  ndata = *q++/sizeof(int);
  if (mode==0) fprintf(outf, "        **%d\n", ndata);
  j += sizeof(int);
  if (j < bufSize) goto NEXT;

DONE:
  if (mode==0) {
    fprintf(outf, "Done.  j=%ld vs %ld, nseq=%d\n", j, bufSize, nseq);
    fprintf(outf, "mr=%d, mz=%d\n", mr, mz);
    fclose(outf); 
    printf("Done.  Output in tell.out\n");
  }
  if (mode==2) 
    close(fd);
}

