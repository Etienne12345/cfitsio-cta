/*  This file, drvrhttp.c contains driver routines for http, ftp and root 
    files. */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.  Users shall not, without prior written   */
/*  permission of the U.S. Government,  establish a claim to statutory     */
/*  copyright.  The Government and others acting on its behalf, shall have */
/*  a royalty-free, non-exclusive, irrevocable,  worldwide license for     */
/*  Government purposes to publish, distribute, translate, copy, exhibit,  */
/*  and perform such material.                                             */


/* Notes on the drivers:

   The ftp driver uses passive mode exclusivly.  If your remote system can't 
   deal with passive mode then it'll fail.  Since Netscape Navigator uses 
   passive mode as well there shouldn't be too many ftp servers which have
   problems.


   The http driver works properly with 301 and 302 redirects.  For many more 
   gory details see http://www.w3c.org/Protocols/rfc2068/rfc2068.  The only
   catch to the 301/302 redirects is that they have to redirect to another 
   http:// url.  If not, things would have to change a lot in cfitsio and this
   was thought to be too difficult.
   
   Redirects look like


   <HTML><HEAD>
   <TITLE>301 Moved Permanently</TITLE>
   </HEAD><BODY>
   <H1>Moved Permanently</H1>
   The document has moved <A HREF="http://heasarc.gsfc.nasa.gov/FTP/software/ftools/release/other/image.fits.gz">here</A>.<P>
   </BODY></HTML>

   This redirect was from apache 1.2.5 but most of the other servers produce 
   something very similiar.  The parser for the redirects finds the first 
   anchor <A> tag in the body and goes there.  If that wasn't what was intended
   by the remote system then hopefully the error stack, which includes notes 
   about the redirect will help the user fix the problem.



   Root protocal doesn't have any real docs, so, the emperical docs are as 
   follows.  

   First, you must use a slightly modified rootd server.  The modifications 
   include implimentation of the stat command which returns the size of the 
   remote file.  Without that it's impossible for cfitsio to work properly
   since fitsfiles don't include any information about the size of the files 
   in the headers.  The rootd server closes the connections on any errors, 
   including reading beyond the end of the file or seeking beyond the end 
   of the file.  The rootd:// driver doesn't reopen a closed connection, if
   the connection is closed you're pretty much done.

   The messages are of the form

   <len><opcode><optional information>

   All binary information is transfered in network format, so use htonl and 
   ntohl to convert back and forth.

   <len> :== 4 byte length, in network format, the len doesn't include the
         length of <len>
   <opcode> :== one of the message opcodes below, 4 bytes, network format
   <optional info> :== depends on opcode

   The response is of the same form with the same opcode sent.  Success is
   indicated by <optional info> being 0.

   Root is a NFSish protocol where each read/write includes the byte
   offset to read or write to.  As a result, seeks will always succeed
   in the driver even if they would cause a fatal error when you try
   to read because you're beyond the end of the file.

   There is file locking on the host such that you need to possibly
   create /usr/tmp/rootdtab on the host system.  There is one file per
   socket connection, though the rootd daemon can support multiple
   files open at once.

   The messages are sent in the following order:

   ROOTD_USER - user name, <optional info> is the user name, trailing
   null is sent though it's not required it seems.  A ROOTD_AUTH
   message is returned with any sort of error meaning that the user
   name is wrong.

   ROOTD_PASS - password, ones complemented, stored in <optional info>. Once
   again the trailing null is sent.  Once again a ROOTD_AUTH message is 
   returned

   ROOTD_OPEN - <optional info> includes filename and one of
     {create|update|read} as the file mode.  ~ seems to be dealt with
     as the username's login directory.  A ROOTD_OPEN message is
     returned.

   Once the file is opened any of the following can be sent:

   ROOTD_STAT - file status and size
   returns a message where <optional info> is the file length in bytes

   ROOTD_FLUSH - flushes the file, not sure this has any real effect
   on the daemon since the daemon uses open/read/write/close rather
   than the buffered fopen/fread/fwrite/fclose.

   ROOTD_GET - on send <optional info> includes a text message of
   offset and length to get.  Return is a status message first with a
   status value, then, the raw bytes for the length that you
   requested.  It's an error to seek or read past the end of the file,
   and, the rootd daemon exits and won't respond anymore.  Ie, don't
   do this.

   ROOTD_PUT - on send <optional info> includes a text message of
   offset and length to put.  Then send the raw bytes you want to
   write.  Then recieve a status message


   When you are finished then you send the message:

   ROOTD_CLOSE - closes the file

   Once the file is closed then the socket is closed.
   
 */

#ifdef HAVE_NET_SERVICES
#include <string.h>
#include "fitsio2.h"

#include <sys/types.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>
#include <setjmp.h>

static jmp_buf env; /* holds the jump buffer for setjmp/longjmp pairs */
static void signal_handler(int sig);


/* Network routine error codes */
#define NET_OK 0
#define NOT_INET_ADDRESS -1000
#define UNKNOWN_INET_HOST -1001
#define CONNECTION_ERROR -1002

/* Network routine constants */
#define NET_DEFAULT 0
#define NET_OOB 1
#define NET_PEEK 2

#define NETTIMEOUT 180 /* in secs */

/* local defines and variables */
#define MAXLEN 1000
#define SHORTLEN 100
static char netoutfile[MAXLEN];


#define ROOTD_USER  2000       /*user id follows */
#define ROOTD_PASS  2001       /*passwd follows */
#define ROOTD_AUTH  2002       /*authorization status (to client) */
#define ROOTD_FSTAT 2003       /*filename follows */
#define ROOTD_OPEN  2004       /*filename follows + mode */
#define ROOTD_PUT   2005       /*offset, number of bytes and buffer */
#define ROOTD_GET   2006       /*offset, number of bytes */
#define ROOTD_FLUSH 2007       /*flush file */
#define ROOTD_CLOSE 2008       /*close file */
#define ROOTD_STAT  2009       /*return rootd statistics */
#define ROOTD_ACK   2010       /*acknowledgement (all OK) */
#define ROOTD_ERR   2011       /*error code and message follow */

typedef struct    /* structure containing disk file structure */ 
{
  int sock;
  long currentpos;
} rootdriver;

static rootdriver handleTable[NIOBUF];  /* allocate diskfile handle tables */

/* static prototypes */

static int NET_TcpConnect(char *hostname, int port);
static int NET_SendRaw(int sock, const void *buf, int length, int opt);
static int NET_RecvRaw(int sock, void *buffer, int length);
static int NET_ParseUrl(const char *url, char *proto, char *host, int *port, 
		 char *fn);
static int CreateSocketAddress(struct sockaddr_in *sockaddrPtr,
			       char *host,int port);
static int ftp_status(FILE *ftp, char *statusstr);
static int http_open_network(char *url, FILE **httpfile, char *contentencoding,
			  int *contentlength);
static int ftp_open_network(char *url, FILE **ftpfile, FILE **command, 
			    int *sock);

static int root_send_buffer(int sock, int op, char *buffer, int buflen);
static int root_recv_buffer(int sock, int *op, char *buffer,int buflen);
static int root_openfile(char *filename, char *rwmode, int *sock);

/*--------------------------------------------------------------------------*/
/* This creates a memory file handle with a copy of the URL in filename. The 
   file is uncompressed if necessary */

int http_open(char *filename, int rwmode, int *handle)
{

  FILE *httpfile;
  char contentencoding[SHORTLEN];
  char newfilename[MAXLEN];
  char errorstr[MAXLEN];
  char recbuf[MAXLEN];
  long len;
  int contentlength;
  int status;
  int closehttpfile = 0;
  int closememfile = 0;

  /* don't do r/w files */
  if (rwmode != 0) {
    ffpmsg("Specify an outfile for r/w access (http_open)");
    goto error;
  }

  /* do the signal handler bits */
  if (setjmp(env) != 0) {
    /* feels like the second time */
    /* this means something bad happened */
    ffpmsg("Timeout (http_open)");
    goto error;
  }

  (void) signal(SIGALRM, signal_handler);
  

  /* Open the network connection */
  alarm(NETTIMEOUT);
  if (http_open_network(filename,&httpfile,contentencoding,
			       &contentlength)) {
    alarm(0);
    /* Try the .gz one */
    strcpy(newfilename,filename);
    strcat(newfilename,".gz");
    alarm(NETTIMEOUT);
    if (http_open_network(newfilename,&httpfile,contentencoding,
			  &contentlength)) {
      /* Now the .Z one */
      alarm(0);
      strcpy(newfilename,filename);
      strcat(newfilename,".Z");
      alarm(NETTIMEOUT);
      if (http_open_network(newfilename,&httpfile,contentencoding,
			    &contentlength)) { 
	alarm(0);
	ffpmsg("Unable to open http file (http_open)");
	goto error;
      }
    }
  }

  closehttpfile++;

  /* Create the memory file */
  if ((status =  mem_create(filename,handle))) {
    ffpmsg("Unable to create memory file (http_open)");
    goto error;
  }

  closememfile++;

  /* Now, what do we do with the file */
  if (!strcmp(contentencoding,"x-gzip") || 
      !strcmp(contentencoding,"x-compress")) {
    /* do the compress dance, which is the same as the gzip dance */
    /* Using the cfitsio routine */

    status = 0;
    /* Ok, this is a tough case, let's be arbritary and say 10*NETTIMEOUT,
       Given the choices for nettimeout above they'll probaby ^C before, but
       it's always worth a shot*/
    
    alarm(NETTIMEOUT*10);
    status = mem_uncompress2mem(filename, httpfile, *handle);
    alarm(0);
    if (status) {
      ffpmsg("Error writing compressed memory file (http_open)");
      goto error;
    }
    
    
  } else {
    /* It's not compressed, bad choice, but we'll copy it anyway */
    if (contentlength % 2880) {
      sprintf(errorstr,"Content-Length not a multiple of 2880 (http_open) %d",
	      contentlength);
      ffpmsg(errorstr);
    }

    /* write a memory file */
    alarm(NETTIMEOUT);
    while(0 != (len = fread(recbuf,1,MAXLEN,httpfile))) {
      alarm(0); /* cancel alarm */
      status = mem_write(*handle,recbuf,len);
      if (status) {
	goto error;
      }
      alarm(NETTIMEOUT); /* rearm the alarm */
    }
  }
  
  fclose(httpfile);

  signal(SIGALRM, SIG_DFL);
  alarm(0);
  return mem_seek(*handle,0);


 error:
  alarm(0); /* clear it */
  if (closehttpfile) {
    fclose(httpfile);
  }
  if (closememfile) {
    mem_close_free(*handle);
  }
  
  signal(SIGALRM, SIG_DFL);
  return (FILE_NOT_OPENED);

}

/*--------------------------------------------------------------------------*/
/* This creates a memory file handle with a copy of the URL in filename.  The
   file must be compressed and is copied to disk first. */

int http_compress_open(char *url, int rwmode, int *handle)
{
  FILE *httpfile;
  FILE *diskfile;
  char contentencoding[SHORTLEN];
  char recbuf[MAXLEN];
  long len;
  int contentlength;
  int status;
  int closehttpfile = 0;
  int closediskfile = 0;
  int closefdiskfile = 0;
  int closememfile = 0;


  /* cfileio made a mistake, should set the netoufile first otherwise 
     we don't know where to write the output file */

  if (!strlen(netoutfile)) 
    {
      ffpmsg
	("Output file not set, shouldn't have happened (http_compress_open)");
      goto error;
    }

  if (rwmode != 0) {
    ffpmsg("Only R/O files for http compressed files(http_compress_open");
    goto error;
  }
  /* do the signal handler bits */
  if (setjmp(env) != 0) {
    /* feels like the second time */
    /* this means something bad happened */
    ffpmsg("Timeout (http_open)");
    goto error;
  }

  signal(SIGALRM, signal_handler);
  
  /* Open the http connectin */
  alarm(NETTIMEOUT);
  if ((status = http_open_network(url,&httpfile,contentencoding,
			       &contentlength))) {
    alarm(0);
    ffpmsg("Unable to open http file (http_compress_open)");
    goto error;
  }

  closehttpfile++;

  /* Better be compressed */
  if (!strcmp(contentencoding,"x-gzip") || 
      !strcmp(contentencoding,"x-compress")) {


    /* Create the new file */
    if ((status =  file_create(netoutfile,handle))) {
      ffpmsg("Unable to create output file (http_compress_open)");
      goto error;
    }
    
    closediskfile++;

    /* write a file */
    alarm(NETTIMEOUT);
    while(0 != (len = fread(recbuf,1,MAXLEN,httpfile))) {
      alarm(0);
      status = file_write(*handle,recbuf,len);
      if (status) {
	ffpmsg("Error writing file (http_compres_open)");
	goto error;
      }
      alarm(NETTIMEOUT);
    }
    file_close(*handle);
    fclose(httpfile);
    closehttpfile--;
    closediskfile--;

    /* File is on disk, let's uncompress it into memory */

    if (NULL == (diskfile = fopen(netoutfile,"r"))) {
      ffpmsg("Unable to reopen disk file (http_compress_open)");
      goto error;
    }
    closefdiskfile++;

    /* Create the memory handle to hold it */
    if ((status =  mem_create(url,handle))) {
      ffpmsg("Unable to create memory file (http_compress_open)");
      goto error;
    }
    closememfile++;

    /* Uncompress it */
    status = 0;
    status = mem_uncompress2mem(url,diskfile,*handle);
    fclose(diskfile);
    closefdiskfile--;
    if (status) {
      ffpmsg("Error writing compressed memory file (http_compress_open)");
      goto error;
    }
      
  } else {
    /* Opps, this should not have happened */
    ffpmsg("Can only compressed files here (http_compress_open)");
    goto error;
  }    
    

  signal(SIGALRM, SIG_DFL);
  alarm(0);
  return mem_seek(*handle,0);

 error:
  alarm(0); /* clear it */
  if (closehttpfile) {
    fclose(httpfile);
  }
  if (closefdiskfile) {
    fclose(diskfile);
  }
  if (closememfile) {
    mem_close_free(*handle);
  }
  if (closediskfile) {
    file_close(*handle);
  } 
  
  signal(SIGALRM, SIG_DFL);
  return (FILE_NOT_OPENED);
}

/*--------------------------------------------------------------------------*/
/* This creates a file handle with a copy of the URL in filename.  The
   file must not be compressed and is copied to disk first. */

int http_file_open(char *url, int rwmode, int *handle)
{
  FILE *httpfile;
  char contentencoding[SHORTLEN];
  char errorstr[MAXLEN];
  char recbuf[MAXLEN];
  long len;
  int contentlength;
  int status;
  int closehttpfile = 0;
  int closefile = 0;


  /* cfileio made a mistake, we need to know where to write the file */
  if (!strlen(netoutfile)) 
    {
      ffpmsg("Output file not set, shouldn't have happened (http_file_open)");
      return (FILE_NOT_OPENED);
    }

  /* do the signal handler bits */
  if (setjmp(env) != 0) {
    /* feels like the second time */
    /* this means something bad happened */
    ffpmsg("Timeout (http_open)");
    goto error;
  }

  signal(SIGALRM, signal_handler);
  
  /* Open the network connection */
  alarm(NETTIMEOUT);
  if ((status = http_open_network(url,&httpfile,contentencoding,
			       &contentlength))) {
    alarm(0);
    ffpmsg("Unable to open http file (http_file_open)");
    goto error;
  }

  closehttpfile++;

  if (!strcmp(contentencoding,"x-gzip") || 
      !strcmp(contentencoding,"x-compress")) {
    /* Opps, this should not have happened */

    ffpmsg("Can't do compressed files here (http_file_open)");
    
    goto error;
  }
    
  /* Create the output file */
  if ((status =  file_create(netoutfile,handle))) {
    ffpmsg("Unable to create output file (http_file_open)");
    goto error;

  }

  /* Give a warning message.  This could just be bad padding at the end
     so don't treat it like an error. */
  closefile++;

  if (contentlength % 2880) {
    sprintf(errorstr,
	    "Content-Length not a multiple of 2880 (http_file_open) %d",
	    contentlength);
    ffpmsg(errorstr);
  }

  /* write a file */
  alarm(NETTIMEOUT);
  while(0 != (len = fread(recbuf,1,MAXLEN,httpfile))) {
    alarm(0);
    status = file_write(*handle,recbuf,len);
    if (status) {
      ffpmsg("Error writing file (http_file_open)");
      goto error;
    }
  }

  
  fclose(httpfile);
  closehttpfile--;

  signal(SIGALRM, SIG_DFL);
  alarm(0);

  file_close(*handle);

  return file_open(netoutfile,rwmode,handle); 


 error:
  alarm(0); /* clear it */
  if (closehttpfile) {
    fclose(httpfile);
  }
  if (closefile) {
    file_close(*handle);
  } 
  
  signal(SIGALRM, SIG_DFL);
  return (FILE_NOT_OPENED);
}

/*--------------------------------------------------------------------------*/
/* This is the guts of the code to get a file via http.  
   url is the input url
   httpfile is set to be the file connected to the socket which you can
     read the file from
   contentencoding is the mime type of the file, returned if the http server
     returns it
   contentlength is the lenght of the file, returned if the http server returns
     it
*/
static int http_open_network(char *url, FILE **httpfile, char *contentencoding,
			  int *contentlength)
{

  int status;
  int sock;
  int tmpint;
  char recbuf[MAXLEN];
  char tmpstr[MAXLEN];
  char errorstr[MAXLEN];
  char proto[SHORTLEN];
  char host[SHORTLEN];
  char fn[MAXLEN];
  char turl[MAXLEN];
  char *scratchstr;
  int port;


  /* Parse the URL apart again */
  strcpy(turl,"http://");
  strcat(turl,url);
  if (NET_ParseUrl(turl,proto,host,&port,fn)) {
    sprintf(errorstr,"URL Parse Error (http_open) %s",url);
    ffpmsg(errorstr);
    return (FILE_NOT_OPENED);
  }

  /* Connect to the remote host */
  sock = NET_TcpConnect(host,port);
  if (sock < 0) {
    ffpmsg("Couldn't connect to host (http_open_network)");
    return (FILE_NOT_OPENED);
  }

  /* Make the socket a stdio file */
  if (NULL == (*httpfile = fdopen(sock,"r"))) {
    ffpmsg ("fdopen failed to convert socket to file (http_open_network)");
    close(sock);
    return (FILE_NOT_OPENED);
  }

  /* Send the GET request to the remote server */
  strcpy(tmpstr,"GET ");
  strcat(tmpstr,fn);
  strcat(tmpstr," http/1.0\n\n");
  status = NET_SendRaw(sock,tmpstr,strlen(tmpstr),NET_DEFAULT);

  /* read the header */
  if (!(fgets(recbuf,MAXLEN,*httpfile))) {
    sprintf (errorstr,"http header short (http_open_network) %s",recbuf);
    ffpmsg(errorstr);
    fclose(*httpfile);
    return (FILE_NOT_OPENED);
  }
  *contentlength = 0;
  contentencoding[0] = '\0';

  /* Our choices are 200, ok, 301, temporary redirect, or 302 perm redirect */
  sscanf(recbuf,"%s %d",tmpstr,&status);
  if (status != 200){
    if (status == 301 || status == 302) {
      /* got a redirect */
      if (status == 301) {
	ffpmsg("Note: Web server replied with a temporary redirect from");
      } else {
	ffpmsg("Note: Web server replied with a redirect from");
      }
      ffpmsg(turl);
      /* now, let's not write the most sophisticated parser here */

      while (fgets(recbuf,MAXLEN,*httpfile)) {
	scratchstr = strstr(recbuf,"<A HREF=\"");
	if (scratchstr != NULL) {
	  /* Ok, we found the beginning of the anchor */
	  scratchstr += 9; /* skip the <A HREF=" bits */
	  scratchstr += 7; /* skip http://, we die if it's really ftp:// */
	  strcpy(turl,strtok(scratchstr,"\""));
	  sprintf(errorstr,"to %s\n",turl);
	  ffpmsg(errorstr);
	  fclose (*httpfile);
	  return 
	    http_open_network(turl,httpfile,contentencoding,contentlength);
	}
      }
      /* if we get here then we couldnt' decide the redirect */
      ffpmsg("but we were unable to find the redirected url in the servers response");
    }
    sprintf(errorstr, 
	    "(http_open_network) Status not 200, was %d\nLine was %s\n",
	    status,recbuf);
    ffpmsg(errorstr);
    fclose(*httpfile);
    return (FILE_NOT_OPENED);
  }

  /* from here the first word holds the keyword we want */
  /* so, read the rest of the header */
  while (fgets(recbuf,MAXLEN,*httpfile)) {
    /* Blank line ends the header */
    if (*recbuf == '\r') break;
    if (strlen(recbuf) > 3) {
      recbuf[strlen(recbuf)-1] = '\0';
      recbuf[strlen(recbuf)-1] = '\0';
    }
    sscanf(recbuf,"%s %d",tmpstr,&tmpint);
    /* Did we get a content-length header ? */
    if (!strcmp(tmpstr,"Content-Length:")) {
      *contentlength = tmpint;
    }
    /* Did we get the content-encoding header ? */
    if (!strcmp(tmpstr,"Content-Encoding:")) {
      if (NULL != (scratchstr = strstr(recbuf,":"))) {
	/* Found the : */
	scratchstr++; /* skip the : */
	scratchstr++; /* skip the extra space */
	strcpy(contentencoding,scratchstr);
      }
    }
  }
  
  /* we're done, so return */
  return 0;
}


/*--------------------------------------------------------------------------*/
/* This creates a memory file handle with a copy of the URL in filename. The 
   file is uncompressed if necessary */

int ftp_open(char *filename, int rwmode, int *handle)
{

  FILE *ftpfile;
  FILE *command;
  int sock;
  char recbuf[MAXLEN];
  long len;
  int status;
  int closememfile;
  int closecommandfile;
  int closeftpfile;

  closememfile = 0;
  closecommandfile = 0;
  closeftpfile = 0;

  /* don't do r/w files */
  if (rwmode != 0) {
    ffpmsg("Specify an outfile for r/w access (ftp_open)");
    return (FILE_NOT_OPENED);
  }

  /* do the signal handler bits */
  if (setjmp(env) != 0) {
    /* feels like the second time */
    /* this means something bad happened */
    ffpmsg("Timeout (http_open)");
    goto error;
  }

  signal(SIGALRM, signal_handler);
  
  /* Open the ftp connetion.  ftpfile is connected to the file port, 
     command is connected to port 21.  sock is the socket on port 21 */

  alarm(NETTIMEOUT);
  if ((status = ftp_open_network(filename,&ftpfile,&command,&sock))) {
    alarm(0);
    ffpmsg("Unable to open ftp file (ftp_open)");
    goto error;
  }

  closeftpfile++;
  closecommandfile++;

  /* create the memory file */
  if ((status = mem_create(filename,handle))) {
    ffpmsg ("Could not create filename to passive port (ftp_open)");
    goto error;
  }
  closememfile++;
  /* This isn't quite right, it'll fail if the file has .gzabc at the end
     for instance */

  /* Decide if the file is compressed */
  if (strstr(filename,".gz") || strstr(filename,".Z")) {
    
    status = 0;
    /* A bit arbritary really, the user will probably hit ^C */
    alarm(NETTIMEOUT*10);
    status = mem_uncompress2mem(filename, ftpfile, *handle);
    alarm(0);
    if (status) {
      ffpmsg("Error writing compressed memory file (ftp_open)");
      goto error;
    }
  } else {
    /* write a memory file */
    alarm(NETTIMEOUT);
    while(0 != (len = fread(recbuf,1,MAXLEN,ftpfile))) {
      alarm(0);
      status = mem_write(*handle,recbuf,len);
      if (status) {
	ffpmsg("Error writing memory file (http_open)");
	goto error;
      }
      alarm(NETTIMEOUT);
    }
  }

  /* close and clean up */
  fclose(ftpfile);
  closeftpfile--;

  NET_SendRaw(sock,"QUIT\n",5,NET_DEFAULT);
  fclose(command);
  closecommandfile--;

  signal(SIGALRM, SIG_DFL);
  alarm(0);

  return mem_seek(*handle,0);

 error:
  alarm(0); /* clear it */
  if (closecommandfile) {
    fclose(command);
  }
  if (closeftpfile) {
    fclose(ftpfile);
  }
  if (closememfile) {
    mem_close_free(*handle);
  }
  
  signal(SIGALRM, SIG_DFL);
  return (FILE_NOT_OPENED);
}



/*--------------------------------------------------------------------------*/
/* This creates a file handle with a copy of the URL in filename. The 
   file must be  uncompressed and is copied to disk first */

int ftp_file_open(char *url, int rwmode, int *handle)
{
  FILE *ftpfile;
  FILE *command;
  char recbuf[MAXLEN];
  long len;
  int sock;
  int status;
  int closeftpfile = 0;
  int closecommandfile = 0;
  int closefile = 0;
  

  /* cfileio made a mistake, need to know where to write the output file */
  if (!strlen(netoutfile)) 
    {
      ffpmsg("Output file not set, shouldn't have happened (ftp_file_open)");
      return (FILE_NOT_OPENED);
    }


  /* do the signal handler bits */
  if (setjmp(env) != 0) {
    /* feels like the second time */
    /* this means something bad happened */
    ffpmsg("Timeout (http_open)");
    goto error;
  }

  signal(SIGALRM, signal_handler);
  
  /* open the network connection to url. ftpfile holds the connection to
     the input file, command holds the connection to port 21, and sock is 
     the socket connected to port 21 */

  alarm(NETTIMEOUT);
  if ((status = ftp_open_network(url,&ftpfile,&command,&sock))) {
    alarm(0);
    ffpmsg("Unable to open http file (ftp_file_open)");
    goto error;
  }
  closeftpfile++;
  closecommandfile++;

  /* Now, what do we do with the file */
  if (strstr(url,".gz") || strstr(url,".Z")) {
    /* Opps, this should not have happened */

    ffpmsg("Can't do compressed files here (ftp_file_open)");
    goto error;
  }
    

  /* Create the output file */
  if ((status =  file_create(netoutfile,handle))) {
    ffpmsg("Unable to create output file (ftp_file_open)");
    goto error;
  }
  closefile++;


  /* write a file */
  alarm(NETTIMEOUT);
  while(0 != (len = fread(recbuf,1,MAXLEN,ftpfile))) {
    alarm(0);
    status = file_write(*handle,recbuf,len);
    if (status) {
      ffpmsg("Error writing file (ftp_file_open)");
      goto error;
    }
    alarm(NETTIMEOUT);
  }

  
  fclose(ftpfile);
  closeftpfile--;
  
  NET_SendRaw(sock,"QUIT\n",5,NET_DEFAULT);
  fclose(command);
  closecommandfile--;

  signal(SIGALRM, SIG_DFL);
  alarm(0);

  file_close(*handle);
  
  return file_open(netoutfile,rwmode,handle);

 error:
  alarm(0); /* clear it */
  if (closeftpfile) {
    fclose(ftpfile);
  }
  if (closecommandfile) {
    fclose(command);
  }
  if (closefile) {
    file_close(*handle);
  } 
  
  signal(SIGALRM, SIG_DFL);
  return (FILE_NOT_OPENED);
}


/*--------------------------------------------------------------------------*/
/* This creates a memory  handle with a copy of the URL in filename. The 
   file must be compressed and is copied to disk first */

int ftp_compress_open(char *url, int rwmode, int *handle)
{
  FILE *ftpfile;
  FILE *command;
  FILE *diskfile;
  char recbuf[MAXLEN];
  long len;
  int status;
  int sock;
  int closeftpfile = 0;
  int closecommandfile = 0;
  int closememfile = 0;
  int closefdiskfile = 0;
  int closediskfile = 0;
  
  
  /* Need to know where to write the output file */
  if (!strlen(netoutfile)) 
    {
      ffpmsg(
	     "Output file not set, shouldn't have happened (ftp_compress_open)");
      return (FILE_NOT_OPENED);
    }
  
  /* do the signal handler bits */
  if (setjmp(env) != 0) {
    /* feels like the second time */
    /* this means something bad happened */
    ffpmsg("Timeout (http_open)");
    goto error;
  }
  
  signal(SIGALRM, signal_handler);
  
  /* Open the network connection to url, ftpfile is connected to the file 
     port, command is connected to port 21.  sock is for writing to port 21 */
  alarm(NETTIMEOUT);
  if ((status = ftp_open_network(url,&ftpfile,&command,&sock))) {
    alarm(0);
    ffpmsg("Unable to open http file (ftp_compress_open)");
    goto error;
  }
  closeftpfile++;
  closecommandfile++;

  /* Now, what do we do with the file */
  if (strstr(url,".gz") || strstr(url,".Z")) {

    /* Create the output file */
    if ((status =  file_create(netoutfile,handle))) {
      ffpmsg("Unable to create output file (ftp_compress_open)");
      goto error;
    }
    closediskfile++;
    
    /* write a file */
    alarm(NETTIMEOUT);
    while(0 != (len = fread(recbuf,1,MAXLEN,ftpfile))) {
      alarm(0);
      status = file_write(*handle,recbuf,len);
      if (status) {
	ffpmsg("Error writing file (ftp_compres_open)");
	goto error;
      }
      alarm(NETTIMEOUT);
    }

    file_close(*handle);
    closediskfile--;
    fclose(ftpfile);
    closeftpfile--;
    /* Close down the ftp connection */
    NET_SendRaw(sock,"QUIT\n",5,NET_DEFAULT);
    fclose(command);
    closecommandfile--;

    /* File is on disk, let's uncompress it into memory */

    if (NULL == (diskfile = fopen(netoutfile,"r"))) {
      ffpmsg("Unable to reopen disk file (ftp_compress_open)");
      return (FILE_NOT_OPENED);
    }
    closefdiskfile++;
  
    if ((status =  mem_create(url,handle))) {
      ffpmsg("Unable to create memory file (ftp_compress_open)");
      goto error;
    }
    closememfile++;

    status = 0;
    status = mem_uncompress2mem(url,diskfile,*handle);
    fclose(diskfile);
    closefdiskfile--;

    if (status) {
      ffpmsg("Error writing compressed memory file (ftp_compress_open)");
      goto error;
    }
      
  } else {
    /* Opps, this should not have happened */
    ffpmsg("Can only compressed files here (ftp_compress_open)");
    goto error;
  }    
    

  signal(SIGALRM, SIG_DFL);
  alarm(0);
  return mem_seek(*handle,0);

 error:
  alarm(0); /* clear it */
  if (closeftpfile) {
    fclose(ftpfile);
  }
  if (closecommandfile) {
    fclose(command);
  }
  if (closefdiskfile) {
    fclose(diskfile);
  }
  if (closememfile) {
    mem_close_free(*handle);
  }
  if (closediskfile) {
    file_close(*handle);
  } 
  
  signal(SIGALRM, SIG_DFL);
  return (FILE_NOT_OPENED);
}


/*--------------------------------------------------------------------------*/
/* Open a ftp connection to filename (really a URL), return ftpfile set to 
   the file connection, and command set to the control connection, with sock
   also set to the control connection */

int ftp_open_network(char *filename, FILE **ftpfile, FILE **command, int *sock)
{
  int status;
  int sock1;
  int tmpint;
  char recbuf[MAXLEN];
  char errorstr[MAXLEN];
  char tmpstr[MAXLEN];
  char proto[SHORTLEN];
  char host[SHORTLEN];
  char *newhost;
  char *username;
  char *password;
  char fn[MAXLEN];
  char *newfn;
  char *passive;
  char *tstr;
  char ip[SHORTLEN];
  char turl[MAXLEN];
  int port;

  /* parse the URL */
  strcpy(turl,"ftp://");
  strcat(turl,filename);
  if (NET_ParseUrl(turl,proto,host,&port,fn)) {
    sprintf(errorstr,"URL Parse Error (ftp_open) %s",filename);
    ffpmsg(errorstr);
    return (FILE_NOT_OPENED);
  }
#ifdef DEBUG
  printf ("proto, %s, host, %s, port %d, fn %s\n",proto,host,port,fn);
#endif
  
  port = 21;
  /* we might have a user name */
  username = "anonymous";
  password = "user@host.com";
  /* is there an @ sign */
  if (NULL != (newhost = strrchr(host,'@'))) {
    *newhost = '\0'; /* make it a null, */
    newhost++; /* Now newhost points to the host name and host points to the 
		  user name, password combo */
    username = host;
    /* is there a : for a password */
    if (NULL != strchr(username,':')) {
      password = strchr(username,':');
      *password = '\0';
      password++;
    }
  } else {
    newhost = host;
  }
  
#ifdef DEBUG
  printf("User %s pass %s\n",username,password); 
#endif
  
  /* Connect to the host on the required port */
  *sock = NET_TcpConnect(newhost,port);
  /* convert it to a stdio file */
  if (NULL == (*command = fdopen(*sock,"r"))) {
    ffpmsg ("fdopen failed to convert socket to stdio file (ftp_open)");
    return (FILE_NOT_OPENED);
    
  }

  /* Wait for the 220 response */
  if (ftp_status(*command,"220 ")) {
    ffpmsg ("error connecting to remote server, no 220 seen (ftp_open)");
    fclose(*command);
    return (FILE_NOT_OPENED);
  }
  
  /* Send the user name and wait for the right response */
  sprintf(tmpstr,"USER %s\n",username);
  status = NET_SendRaw(*sock,tmpstr,strlen(tmpstr),NET_DEFAULT);
  
  if (ftp_status(*command,"331 ")) {
    ffpmsg ("USER error no 331 seen (ftp_open)");
    fclose(*command);
    return (FILE_NOT_OPENED);
    
  }
  
  /* Send the password and wait for the right response */
  sprintf(tmpstr,"PASS %s\n",password);
  status = NET_SendRaw(*sock,tmpstr,strlen(tmpstr),NET_DEFAULT);
  
  if (ftp_status(*command,"230 ")) {
    ffpmsg ("PASS error, no 230 seen (ftp_open)");
    fclose(*command);
    return (FILE_NOT_OPENED);
  }
  
  
  /* now do the cwd command */
  newfn = strrchr(fn,'/');
  if (newfn == NULL) {
    strcpy(tmpstr,"CWD /\n");
    newfn = fn;
  } else {
    *newfn = '\0';
    newfn++;
    if (strlen(fn) == 0) {
      strcpy(tmpstr,"CWD /\n");
    } else {
      /* remove the leading slash */
      if (fn[0] == '/') {
	sprintf(tmpstr,"CWD %s\n",&fn[1]);
      } else {
	sprintf(tmpstr,"CWD %s\n",fn);
      } 
    }
  }
  
#ifdef DEBUG
  printf("CWD command is %s\n",tmpstr);
#endif
  status = NET_SendRaw(*sock,tmpstr,strlen(tmpstr),NET_DEFAULT);
  
  if (ftp_status(*command,"250 ")) {
    ffpmsg ("CWD error, no 250 seen (ftp_open)");
    fclose(*command);
    return (FILE_NOT_OPENED);
  }
  
  if (!strlen(newfn)) {
    ffpmsg("Null file name (ftp_open)");
    fclose(*command);
    return (FILE_NOT_OPENED);
  }
  
  /* Send the retrieve command to see if the file exists*/
  sprintf(tmpstr,"RETR %s\n",newfn);
#ifdef DEBUG
  printf ("Checking to see if %s, exists %d\n",tmpstr,strlen(newfn));
#endif
  status = NET_SendRaw(*sock,tmpstr,strlen(tmpstr),NET_DEFAULT);
  /* We better get a 425 here, we haven't sent a PORT or PASV command
     yet.  If we don't get a 425 then we're hosed and the file
     doesn't exist */
  if (ftp_status(*command,"425 ")) {
    ffpmsg("File doesn't exist on remote server (ftp_open)");
    fclose(*command);
    return (FILE_NOT_OPENED);
  }
  /* we're going to use passive mode here */
  
  status = NET_SendRaw(*sock,"PASV\n",5,NET_DEFAULT);
  if (!(fgets(recbuf,MAXLEN,*command))) {
    ffpmsg ("PASV error (ftp_open)");
    fclose(*command);
    return (FILE_NOT_OPENED);
  }
  
  /*  Passive mode response looks like
      227 Entering Passive Mode (129,194,67,8,210,80) */
  if (recbuf[0] == '2' && recbuf[1] == '2' && recbuf[2] == '7') {
    /* got a good passive mode response, find the opening ( */
    
    if (!(passive = strchr(recbuf,'('))) {
      ffpmsg ("PASV error (ftp_open)");
      fclose(*command);
      return (FILE_NOT_OPENED);
    }
    
    *passive = '\0';
    passive++;
    ip[0] = '\0';
      
    /* Messy parsing of response from PASV *command */
    
    if (!(tstr = strtok(passive,",)"))) {
      ffpmsg ("PASV error (ftp_open)");
      fclose(*command);
      return (FILE_NOT_OPENED);
    }
    strcpy(ip,tstr);
    strcat(ip,".");
    
    if (!(tstr = strtok(NULL,",)"))) {
      ffpmsg ("PASV error (ftp_open)");
      fclose(*command);
      return (FILE_NOT_OPENED);
    }
    strcat(ip,tstr);
    strcat(ip,".");
    
    if (!(tstr = strtok(NULL,",)"))) {
      ffpmsg ("PASV error (ftp_open)");
      fclose(*command);
      return (FILE_NOT_OPENED);
    }
    strcat(ip,tstr);
    strcat(ip,".");
    
    if (!(tstr = strtok(NULL,",)"))) {
      ffpmsg ("PASV error (ftp_open)");
      fclose(*command);
      return (FILE_NOT_OPENED);
    }
    strcat(ip,tstr);
    
    /* Done the ip number, now do the port # */
    if (!(tstr = strtok(NULL,",)"))) {
      ffpmsg ("PASV error (ftp_open)");
      fclose(*command);
      return (FILE_NOT_OPENED);
    }
    sscanf(tstr,"%d",&port);
    port *= 256;
    
    if (!(tstr = strtok(NULL,",)"))) {
      ffpmsg ("PASV error (ftp_open)");
      fclose(*command);
      return (FILE_NOT_OPENED);
    }
    sscanf(tstr,"%d",&tmpint);
    port += tmpint;
    
    /* Always use binary mode */
    sprintf(tmpstr,"TYPE I\n");
    status = NET_SendRaw(*sock,tmpstr,strlen(tmpstr),NET_DEFAULT);
    
    if (ftp_status(*command,"200 ")) {
      ffpmsg ("TYPE I error, 200 not seen (ftp_open)");
      fclose(*command);
      return (FILE_NOT_OPENED);
    }
    
    if (!strlen(newfn)) {
      ffpmsg("Null file name (ftp_open)");
      fclose(*command);
      return (FILE_NOT_OPENED);
    }
    
    /* Send the retrieve command */
    sprintf(tmpstr,"RETR %s\n",newfn);
#ifdef DEBUG
    printf ("Retrieving file %s, %d\n",tmpstr,strlen(newfn));
#endif
    status = NET_SendRaw(*sock,tmpstr,strlen(tmpstr),NET_DEFAULT);

    /* COnnect to the data port */
    sock1 = NET_TcpConnect(ip,port);
    if (NULL == (*ftpfile = fdopen(sock1,"r"))) {
      ffpmsg ("Could not connect to passive port (ftp_open)");
      fclose(*command);
      return (FILE_NOT_OPENED);
    }

    /* now we return */

    return 0;
  }
  
  /* no passive mode */

  NET_SendRaw(*sock,"QUIT\n",5,NET_DEFAULT);
  fclose(*command);
  return (FILE_NOT_OPENED);
}

/*--------------------------------------------------------------------------*/
/* return a socket which results from connection to hostname on port port */
static int NET_TcpConnect(char *hostname, int port)
{
  /* Connect to hostname on port */
 
   struct sockaddr_in sockaddr;
   int sock;
   int stat;
   int val = 1;
 
   CreateSocketAddress(&sockaddr,hostname,port);
   /* Create socket */
   if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
     ffpmsg("Can't create socket");
     return CONNECTION_ERROR;
   }
 
   if ((stat = connect(sock, (struct sockaddr*) &sockaddr, 
		       sizeof(sockaddr))) 
       < 0) {
     close(sock);
     perror("NET_Tcpconnect - Connection error");
     ffpmsg("Can't connect to host, connection error");
     return CONNECTION_ERROR;
   }
   setsockopt(sock, IPPROTO_TCP, TCP_NODELAY, (char *)&val, sizeof(val));
   setsockopt(sock, SOL_SOCKET,  SO_KEEPALIVE, (char *)&val, sizeof(val));

   val = 65536;
   setsockopt(sock, SOL_SOCKET,  SO_SNDBUF,    (char *)&val, sizeof(val));
   setsockopt(sock, SOL_SOCKET,  SO_RCVBUF,    (char *)&val, sizeof(val));
   return sock;
}

/*--------------------------------------------------------------------------*/
/* Write len bytes from buffer to socket sock */
static int NET_SendRaw(int sock, const void *buffer, int length, int opt)
{

  char * buf = (char *) buffer;
 
   int flag;
   int n, nsent = 0;
 
   switch (opt) {
   case NET_DEFAULT:
     flag = 0;
     break;
   case NET_OOB:
     flag = MSG_OOB;
     break;
   case NET_PEEK:            
   default:
     flag = 0;
     break;
   }
 
   if (sock < 0) return -1;
   
   for (n = 0; n < length; n += nsent) {
     if ((nsent = send(sock, buf+n, length-n, flag)) <= 0) {
       return nsent;
     }
#ifdef DEBUG
     printf ("send raw, sent %d bytes\n",nsent);
#endif
   }
#ifdef DEBUG
   printf ("send raw end, sent %d bytes\n",n);
#endif
   return n;
}
 
static int NET_RecvRaw(int sock, void *buffer, int length)
{
  /* Receive exactly length bytes into buffer. Returns number of bytes */
  /* received. Returns -1 in case of error. */


   int nrecv, n;
   char *buf = (char *)buffer;

   if (sock < 0) return -1;
   for (n = 0; n < length; n += nrecv) {
      while ((nrecv = recv(sock, buf+n, length-n, 0)) == -1 && errno == EINTR)
	errno = 0;     /* probably a SIGCLD that was caught */
      if (nrecv < 0)
         return nrecv;
      else if (nrecv == 0)
	break;        /*/ EOF */
   }


   return n;
}
 


/*--------------------------------------------------------------------------*/
/* Yet Another URL Parser 
   url - input url
   proto - input protocol
   host - output host
   port - output port
   fn - output filename
*/

static int NET_ParseUrl(const char *url, char *proto, char *host, int *port, 
		 char *fn)
{
  /* parses urls into their bits */
  /* returns 1 if error, else 0 */

  char *urlcopy, *urlcopyorig;
  char *ptrstr;
  char *thost;
  int isftp = 0;

  /* figure out if there is a http: or  ftp: */

  urlcopyorig = urlcopy = (char *) malloc(strlen(url)+1);
  strcpy(urlcopy,url);

  /* set some defaults */
  *port = 80;
  strcpy(proto,"http:");
  strcpy(host,"localhost");
  strcpy(fn,"/");
  
  ptrstr = strstr(urlcopy,"http:");
  if (ptrstr == NULL) {
    /* Nope, not http: */
    ptrstr = strstr(urlcopy,"root:");
    if (ptrstr == NULL) {
      /* Nope, not root either */
      ptrstr = strstr(urlcopy,"ftp:");
      if (ptrstr != NULL) {
	if (ptrstr == urlcopy) {
	  strcpy(proto,"ftp:");
	  *port = 21;
	  isftp++;
	  urlcopy += 4; /* move past ftp: */
	} else {
	  /* not at the beginning, bad url */
	  free(urlcopyorig);
	  return 1;
	}
      }
    } else {
      if (ptrstr == urlcopy) {
	urlcopy += 5; /* move past root: */
      } else {
	/* not at the beginning, bad url */
	free(urlcopyorig);
	return 1;
      }
    }
  } else {
    if (ptrstr == urlcopy) {
      urlcopy += 5; /* move past http: */
    } else {
      free(urlcopyorig);
      return 1;
    }
  }

  /* got the protocol */
  /* get the hostname */
  if (urlcopy[0] == '/' && urlcopy[1] == '/') {
    /* we have a hostname */
    urlcopy += 2; /* move past the // */
  }
  /* do this only if http */
  if (!strcmp(proto,"http:")) {
    strcpy(host,urlcopy);
    thost = host;
    while (*urlcopy != '/' && *urlcopy != ':' && *urlcopy) {
      thost++;
      urlcopy++;
    }
    /* Now, we should either be at the end of the string, have a /, or have a
       : */
    *thost = '\0';
    if (*urlcopy == ':') {
      /* follows a port number */
      urlcopy++;
      sscanf(urlcopy,"%d",port);
      while (*urlcopy != '/' && *urlcopy) urlcopy++; /* step to the */
    }
  } else {
    /* do this for ftp */
    strcpy(host,urlcopy);
    thost = host;
    while (*urlcopy != '/' && *urlcopy) {
      thost++;
      urlcopy++; 
    }
    *thost = '\0';
    /* Now, we should either be at the end of the string, or have a / */
    
  }
  /* Now the rest is a fn */

  if (*urlcopy) {
    strcpy(fn,urlcopy);
  }
  free(urlcopyorig);
  return 0;
}

/*--------------------------------------------------------------------------*/

/* Small helper functions to set the netoutfile static string */
/* Called by cfileio after parsing the output file off of the input file
   url */

int http_checkfile (char *urltype, char *infile, char *outfile)

{
  char newinfile[MAXLEN];
  FILE *httpfile;
  char contentencoding[MAXLEN];
  int contentlength;
  
  /* default to http://
   */
    
  strcpy(urltype,"http://");

  if (strlen(outfile)) {
    /* there is an output file */
    strcpy(netoutfile,outfile);

    if (!http_open_network(infile,&httpfile,contentencoding,&contentlength)) {
      fclose(httpfile);
      /* It's there, we're happy */
      if (strstr(infile,".gz") || (strstr(infile,".Z"))) {
	/* It's compressed */
	strcpy(urltype,"httpcompress://");
      } else {
	strcpy(urltype,"httpfile://");
      }
      return 0;
    }

    /* Ok, let's try the .gz one */
    strcpy(newinfile,infile);
    strcat(newinfile,".gz");
    if (!http_open_network(newinfile,&httpfile,contentencoding,
			   &contentlength)) {
      fclose(httpfile);
      strcpy(infile,newinfile);
      strcat(outfile,".gz");
      strcat(netoutfile,".gz");
      /* It's there, we're happy, and, it's compressed  */
      strcpy(urltype,"httpcompress://");
      return 0;
    }
    
    /* Ok, let's try the .Z one */
    strcpy(newinfile,infile);
    strcat(newinfile,".Z");
    if (!http_open_network(newinfile,&httpfile,contentencoding,
			   &contentlength)) {
      fclose(httpfile);
      strcpy(infile,newinfile);
      strcat(outfile,".Z");
      strcat(netoutfile,".Z");
      /* It's there, we're happy, and, it's compressed  */
      strcpy(urltype,"httpcompress://");
      return 0;
    }
    

  } 
  return 0;
}
/*--------------------------------------------------------------------------*/
int ftp_checkfile (char *urltype, char *infile, char *outfile)

{
  
  /* default to ftp://
   */
    
  strcpy(urltype,"ftp://");

  if (strlen(outfile)) {
    /* there is an output file */
    strcpy(netoutfile,outfile);
    if (strstr(infile,".gz") || (strstr(infile,".Z"))) {
      /* It's compressed */
      strcpy(urltype,"ftpcompress://");
    } else {
      strcpy(urltype,"ftpfile://");
    }
  } 
  return 0;
}

/*--------------------------------------------------------------------------*/
/* A small helper function to wait for a particular status on the ftp 
   connectino */
static int ftp_status(FILE *ftp, char *statusstr)
{
  /* read through until we find a string beginning with statusstr */
  /* This needs a timeout */

  char recbuf[MAXLEN];
  int len;

  len = strlen(statusstr);
  while (1) {
    if (!(fgets(recbuf,MAXLEN,ftp))) {
      return 1; /* error reading */
    }
#ifdef DEBUG
    printf ("%s",recbuf); 
#endif
    recbuf[len] = '\0'; /* make it short */
    if (!strcmp(recbuf,statusstr)) {
      return 0; /* we're ok */
    }
    if (recbuf[0] > '3') {
      /* oh well, some sort of error */
      return 1; 
    }
  }
}

/*
 *----------------------------------------------------------------------
 *
 * CreateSocketAddress --
 *
 *	This function initializes a sockaddr structure for a host and port.
 *
 * Results:
 *	1 if the host was valid, 0 if the host could not be converted to
 *	an IP address.
 *
 * Side effects:
 *	Fills in the *sockaddrPtr structure.
 *
 *----------------------------------------------------------------------
 */

static int
CreateSocketAddress(
    struct sockaddr_in *sockaddrPtr,	/* Socket address */
    char *host,				/* Host.  NULL implies INADDR_ANY */
    int port)				/* Port number */
{
    struct hostent *hostent;		/* Host database entry */
    struct in_addr addr;		/* For 64/32 bit madness */
    char localhost[MAXLEN];

    strcpy(localhost,host);

    memset((void *) sockaddrPtr, '\0', sizeof(struct sockaddr_in));
    sockaddrPtr->sin_family = AF_INET;
    sockaddrPtr->sin_port = htons((unsigned short) (port & 0xFFFF));
    if (host == NULL) {
	addr.s_addr = INADDR_ANY;
    } else {
        addr.s_addr = inet_addr(localhost);
        if (addr.s_addr == -1) {
            hostent = gethostbyname(localhost);
            if (hostent != NULL) {
                memcpy((void *) &addr,
                        (void *) hostent->h_addr_list[0],
                        (size_t) hostent->h_length);
            } else {
#ifdef	EHOSTUNREACH
                errno = EHOSTUNREACH;
#else
#ifdef ENXIO
                errno = ENXIO;
#endif
#endif
                return 0;	/* error */
            }
        }
    }
        
    /*
     * NOTE: On 64 bit machines the assignment below is rumored to not
     * do the right thing. Please report errors related to this if you
     * observe incorrect behavior on 64 bit machines such as DEC Alphas.
     * Should we modify this code to do an explicit memcpy?
     */

    sockaddrPtr->sin_addr.s_addr = addr.s_addr;
    return 1;	/* Success. */
}

/* Signal handler for timeouts */

static void signal_handler(int sig) {

  switch (sig) {
  case SIGALRM:    /* process for alarm */
    longjmp(env,sig);
    
  default: {
      /* Hmm, shouldn't have happend */
      exit(sig);
    }
  }
}

#endif

/**************************************************************/

/* Root driver */




/*--------------------------------------------------------------------------*/
int root_init(void)
{
    int ii;

    for (ii = 0; ii < NIOBUF; ii++)  /* initialize all empty slots in table */
    {
       handleTable[ii].sock = 0;
       handleTable[ii].currentpos = 0;
    }
    return(0);
}
/*--------------------------------------------------------------------------*/
int root_setoptions(int options)
{
  /* do something with the options argument, to stop compiler warning */
  options = 0;
  return(options);
}
/*--------------------------------------------------------------------------*/
int root_getoptions(int *options)
{
  *options = 0;
  return(0);
}
/*--------------------------------------------------------------------------*/
int root_getversion(int *version)
{
    *version = 10;
    return(0);
}
/*--------------------------------------------------------------------------*/
int root_shutdown(void)
{
  return(0);
}
/*--------------------------------------------------------------------------*/
int root_open(char *url, int rwmode, int *handle)
{
    int ii, status;
    int sock;

    *handle = -1;
    for (ii = 0; ii < NIOBUF; ii++)  /* find empty slot in table */
    {
        if (handleTable[ii].sock == 0)
        {
            *handle = ii;
            break;
        }
    }

    if (*handle == -1)
       return(TOO_MANY_FILES);    /* too many files opened */

    /*open the file */
    if (rwmode) {
      status = root_openfile(url, "update", &sock);
    } else {
      status = root_openfile(url, "read", &sock);
    }
    if (status)
      return(status);
    
    handleTable[ii].sock = sock;
    handleTable[ii].currentpos = 0;
    
    return(0);
}
/*--------------------------------------------------------------------------*/
int root_create(char *filename, int *handle)
{
    int ii, status;
    int sock;

    *handle = -1;
    for (ii = 0; ii < NIOBUF; ii++)  /* find empty slot in table */
    {
        if (handleTable[ii].sock == 0)
        {
            *handle = ii;
            break;
        }
    }

    if (*handle == -1)
       return(TOO_MANY_FILES);    /* too many files opened */

    /*open the file */
    status = root_openfile(filename, "create", &sock);

    if (status) {
      ffpmsg("Unable to create file");
      return(status);
    }
    
    handleTable[ii].sock = sock;
    handleTable[ii].currentpos = 0;
    
    return(0);
}
/*--------------------------------------------------------------------------*/
int root_size(int handle, long *filesize)
/*
  return the size of the file in bytes
*/
{

  int sock;
  int offset;
  int status;
  int op;

  sock = handleTable[handle].sock;

  status = root_send_buffer(sock,ROOTD_STAT,NULL,0);
  status = root_recv_buffer(sock,&op,(char *)&offset, 4);
  *filesize = ntohl(offset);
  
  return(0);
}
/*--------------------------------------------------------------------------*/
int root_close(int handle)
/*
  close the file
*/
{

  int status;
  int sock;

  sock = handleTable[handle].sock;
  status = root_send_buffer(sock,ROOTD_CLOSE,NULL,0);
  close(sock);
  handleTable[handle].sock = 0;
  return(0);
}
/*--------------------------------------------------------------------------*/
int root_flush(int handle)
/*
  flush the file
*/
{
  int status;
  int sock;

  sock = handleTable[handle].sock;
  status = root_send_buffer(sock,ROOTD_FLUSH,NULL,0);
  return(0);
}
/*--------------------------------------------------------------------------*/
int root_seek(int handle, long offset)
/*
  seek to position relative to start of the file
*/
{
  handleTable[handle].currentpos = offset;
  return(0);
}
/*--------------------------------------------------------------------------*/
int root_read(int hdl, void *buffer, long nbytes)
/*
  read bytes from the current position in the file
*/
{
  char msg[SHORTLEN];
  int op;
  int status;
  int astat;

  sprintf(msg,"%ld %ld ",handleTable[hdl].currentpos,nbytes);
  status = root_send_buffer(handleTable[hdl].sock,ROOTD_GET,msg,strlen(msg));
  if (status != strlen(msg)) {
    return (READ_ERROR);
  }
  astat = 0;
  status = root_recv_buffer(handleTable[hdl].sock,&op,(char *) &astat,4);
  if (astat != 0) {
    return (READ_ERROR);
  }
#ifdef DEBUG
  printf("root_read, op %d astat %d\n",op,astat);
#endif
  status = NET_RecvRaw(handleTable[hdl].sock,buffer,nbytes);
  if (status != nbytes) {
    return (READ_ERROR);
  }
  handleTable[hdl].currentpos += nbytes;

  return(0);
}
/*--------------------------------------------------------------------------*/
int root_write(int hdl, void *buffer, long nbytes)
/*
  write bytes at the current position in the file
*/
{

  char msg[SHORTLEN];
  int len;
  int sock;
  int status;
  int astat;
  int op;

  sock = handleTable[hdl].sock;
  sprintf(msg,"%ld %ld ",handleTable[hdl].currentpos,nbytes);

  len = strlen(msg);
  status = root_send_buffer(sock,ROOTD_PUT,msg,len+1);
  if (status != len+1) {
    return (WRITE_ERROR);
  }
  status = NET_SendRaw(sock,buffer,nbytes,NET_DEFAULT);
  if (status != nbytes) {
    return (WRITE_ERROR);
  }
  astat = 0;
  status = root_recv_buffer(handleTable[hdl].sock,&op,(char *) &astat,4);
#ifdef DEBUG
  printf("root_read, op %d astat %d\n",op,astat);
#endif
  if (astat != 0) {
    return (WRITE_ERROR);
  }
  handleTable[hdl].currentpos += nbytes;
  return(0);
}

/*--------------------------------------------------------------------------*/
int root_openfile(char *url, char *rwmode, int *sock)
     /*
       lowest level routine to physically open a root file
     */
{
  
  int status;
  char recbuf[MAXLEN];
  char errorstr[MAXLEN];
  char proto[SHORTLEN];
  char host[SHORTLEN];
  char fn[MAXLEN];
  char turl[MAXLEN];
  int port;
  int op;
  int ii;
  int authstat;
  
  
  /* Parse the URL apart again */
  strcpy(turl,"root://");
  strcat(turl,url);
  if (NET_ParseUrl(turl,proto,host,&port,fn)) {
    sprintf(errorstr,"URL Parse Error (root_open) %s",url);
    ffpmsg(errorstr);
    return (FILE_NOT_OPENED);
  }
  
#ifdef DEBUG
  printf("Connecting to %s on port %d\n",host,port);
#endif
  /* Connect to the remote host */
  *sock = NET_TcpConnect(host,port);
  if (*sock < 0) {
    ffpmsg("Couldn't connect to host (http_open_network)");
    return (FILE_NOT_OPENED);
  }
  
  /* get the username */
  if (NULL != getenv("ROOTUSERNAME")) {
    strcpy(recbuf,getenv("ROOTUSERNAME"));
  } else {
    printf("Username: ");
    fgets(recbuf,MAXLEN,stdin);
    recbuf[strlen(recbuf)-1] = '\0';
  }
  
  status = root_send_buffer(*sock, ROOTD_USER, recbuf,strlen(recbuf));
  if (status < 0) {
    ffpmsg("error talking to remote system on username ");
    return (FILE_NOT_OPENED);
  }
  
  status = root_recv_buffer(*sock,&op,(char *)&authstat,4);
  if (!status) {
    ffpmsg("error talking to remote system on username");
    return (FILE_NOT_OPENED);
  }
  
#ifdef DEBUG
  printf("op is %d and authstat is %d\n",op,authstat);
#endif
  
  if (op != ROOTD_AUTH) {
    ffpmsg("ERROR on ROOTD_USER");
    ffpmsg(recbuf);
    return (FILE_NOT_OPENED);
  }
  

  /* now the password */
  if (NULL != getenv("ROOTPASSWORD")) {
    strcpy(recbuf,getenv("ROOTPASSWORD"));
  } else {
    printf("Password: ");
    fgets(recbuf,MAXLEN,stdin);
    recbuf[strlen(recbuf)-1] = '\0';
  }
  /* ones complement the password */
  for (ii=0;ii<strlen(recbuf);ii++) {
    recbuf[ii] = ~recbuf[ii];
  }
  
  status = root_send_buffer(*sock, ROOTD_PASS, recbuf, strlen(recbuf));
  if (status < 0) {
    ffpmsg("error talking to remote system sending password");
    return (FILE_NOT_OPENED);
  }
  
  status = root_recv_buffer(*sock,&op,(char *)&authstat,4);
  if (status < 0) {
    ffpmsg("error talking to remote system acking password");
    return (FILE_NOT_OPENED);
  }
  
#ifdef DEBUG
  printf("op is %d and authstat is %d\n",op,authstat);
#endif
  if (op != ROOTD_AUTH) {
    ffpmsg("ERROR on ROOTD_PASS");
    ffpmsg(recbuf);
    return (FILE_NOT_OPENED);
  }
  
  /* now the file open request */
  strcpy(recbuf,fn);
  strcat(recbuf," ");
  strcat(recbuf,rwmode);

  status = root_send_buffer(*sock, ROOTD_OPEN, recbuf, strlen(recbuf));
  if (status < 0) {
    ffpmsg("error talking to remote system on open ");
    return (FILE_NOT_OPENED);
  }

  status = root_recv_buffer(*sock,&op,(char *)&authstat,4);
  if (status < 0) {
    ffpmsg("error talking to remote system on open");
    return (FILE_NOT_OPENED);
  }

#ifdef DEBUG
  printf("op is %d and recbuf is %d\n",op,authstat);
#endif
  
  if ((op != ROOTD_OPEN) && (authstat != 0)) {
    ffpmsg("ERROR on ROOTD_OPEN");
    ffpmsg(recbuf);
    return (FILE_NOT_OPENED);
  }

  return 0;

}

static int root_send_buffer(int sock, int op, char *buffer, int buflen)
{
  /* send a buffer, the form is
     <len>
     <op>
     <buffer>

     <len> includes the 4 bytes for the op, the length bytes (4) are implicit


     if buffer is null don't send it, not everything needs something sent */

  int len;
  int status;

  int hdr[2];

  len = 4;

  if (buffer != NULL) {
    len += buflen;
  }
  
  hdr[0] = htonl(len);

#ifdef DEBUG
  printf("len sent is %x\n",hdr[0]);
#endif

  hdr[1] = htonl(op);
#ifdef DEBUG
  printf("op sent is %x\n",hdr[1]);
#endif
  

#ifdef DEBUG
  printf("Sending op %d and length of %d\n",op,len);
#endif

  status = NET_SendRaw(sock,hdr,sizeof(hdr),NET_DEFAULT);
  if (status < 0) {
    return status;
  }
  if (buffer != NULL) {
    status = NET_SendRaw(sock,buffer,buflen,NET_DEFAULT);
  }
  return status;
}
  
static int root_recv_buffer(int sock, int *op, char *buffer, int buflen)
{
  /* recv a buffer, the form is
     <len>
     <op>
     <buffer>

  */

  int recv = 0;
  int len;
  int status;
  char recbuf[MAXLEN];

  status = NET_RecvRaw(sock,&len,4);
#ifdef DEBUG
  printf("Recv: status from rec is %d\n",status);
#endif
  if (status < 0) {
    return status;
  }
  recv += status;

  len = ntohl(len);
#ifdef DEBUG
  printf ("Recv: length is %d\n",len);
#endif

  /* ok, have the length, recive the operation */
  len -= 4;
  status = NET_RecvRaw(sock,op,4);
  if (status < 0) {
    return status;
  }

  recv += status;

  *op = ntohl(*op);
#ifdef DEBUG
  printf ("Recv: Operation is %d\n",*op);
#endif
  
  if (len > MAXLEN) {
    len = MAXLEN;
  }

  if (len > 0) { /* Get the rest of the message */
    status = NET_RecvRaw(sock,recbuf,len);
    if (len > buflen) {
      len = buflen;
    }
    memcpy(buffer,recbuf,len);
    if (status < 0) {
      return status;
    }
  } 

  recv += status;
  return recv;
}
