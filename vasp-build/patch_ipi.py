#!/usr/bin/env python
import re
import os
from argparse import ArgumentParser
from pathlib import Path
from copy import copy
import shutil

fpath = Path(__file__)
fdir = fpath.parent

def replace(s):
    return s

def insert(s):
    return r"\1" + s + r"\3"

patch_fsockets_f90 = """!F90 ISO_C_BINGING wrapper for socket communication.

!Copyright (C) 2013, Michele Ceriotti

!Permission is hereby granted, free of charge, to any person obtaining
!a copy of this software and associated documentation files (the
!"Software"), to deal in the Software without restriction, including
!without limitation the rights to use, copy, modify, merge, publish,
!distribute, sublicense, and/or sell copies of the Software, and to
!permit persons to whom the Software is furnished to do so, subject to
!the following conditions:

!The above copyright notice and this permission notice shall be included
!in all copies or substantial portions of the Software.

!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
!IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
!CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
!TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
!SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


!Contains both the functions that transmit data to the socket and read the data
!back out again once finished, and the function which opens the socket initially.

!Functions:
!   open_socket: Opens a socket with the required host server, socket type and
!      port number.
!   write_buffer: Writes a string to the socket.
!   read_buffer: Reads data from the socket.

   MODULE F90SOCKETS
   USE ISO_C_BINDING
   
   IMPLICIT NONE

  INTERFACE writebuffer
      MODULE PROCEDURE writebuffer_s, &
                       writebuffer_d, writebuffer_dv, &
                       writebuffer_i
                       
  END INTERFACE 

  INTERFACE readbuffer
      MODULE PROCEDURE readbuffer_s, &
                       readbuffer_dv, readbuffer_d, &
                       readbuffer_i
                       
  END INTERFACE 

  INTERFACE
    SUBROUTINE open_csocket(psockfd, inet, port, host) BIND(C, name="open_socket")
      USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: psockfd, inet, port
    CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: host

    END SUBROUTINE open_csocket

    
    SUBROUTINE writebuffer_csocket(psockfd, pdata, plen) BIND(C, name="writebuffer")
      USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: psockfd
    TYPE(C_PTR), VALUE                       :: pdata
    INTEGER(KIND=C_INT)                      :: plen

    END SUBROUTINE writebuffer_csocket       

    SUBROUTINE readbuffer_csocket(psockfd, pdata, plen) BIND(C, name="readbuffer")
      USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: psockfd
    TYPE(C_PTR), VALUE                       :: pdata
    INTEGER(KIND=C_INT)                      :: plen

    END SUBROUTINE readbuffer_csocket   
  END INTERFACE

   CONTAINS
   
   SUBROUTINE open_socket(psockfd, inet, port, host)      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: inet, port
      INTEGER, INTENT(OUT) :: psockfd
      CHARACTER(LEN=1024), INTENT(IN) :: host
      CHARACTER(LEN=1,KIND=C_CHAR) :: chost(1024)

      CALL fstr2cstr(host, chost)
      CALL open_csocket(psockfd, inet, port, host)
   END SUBROUTINE

   SUBROUTINE fstr2cstr(fstr, cstr, plen)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: fstr
      CHARACTER(LEN=1,KIND=C_CHAR), INTENT(OUT) :: cstr(:)
      INTEGER, INTENT(IN), OPTIONAL :: plen
      
      INTEGER i,n
      IF (PRESENT(plen)) THEN
         n = plen
         DO i=1,n
            cstr(i) = fstr(i:i)
         ENDDO
      ELSE
         n = LEN_TRIM(fstr)
         DO i=1,n
            cstr(i) = fstr(i:i)
         ENDDO
         cstr(n+1) = C_NULL_CHAR
      END IF
   END SUBROUTINE

  SUBROUTINE writebuffer_d (psockfd, fdata)
      USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    REAL(KIND=8), INTENT(IN)                :: fdata

    REAL(KIND=C_DOUBLE), TARGET              :: cdata

      cdata = fdata
      CALL writebuffer_csocket(psockfd, c_loc(cdata), 8)
  END SUBROUTINE

  SUBROUTINE writebuffer_i (psockfd, fdata)
      USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd, fdata

    INTEGER(KIND=C_INT), TARGET              :: cdata

      cdata = fdata
      CALL writebuffer_csocket(psockfd, c_loc(cdata), 4)
  END SUBROUTINE

  SUBROUTINE writebuffer_s (psockfd, fstring, plen)
      USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    CHARACTER(LEN=*), INTENT(IN)             :: fstring
    INTEGER, INTENT(IN)                      :: plen

    INTEGER                                  :: i
    CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cstring(plen)

      DO i = 1,plen
         cstring(i) = fstring(i:i)
      ENDDO
      CALL writebuffer_csocket(psockfd, c_loc(cstring(1)), plen)
  END SUBROUTINE

  SUBROUTINE writebuffer_dv(psockfd, fdata, plen)
      USE ISO_C_BINDING  
    INTEGER, INTENT(IN)                      :: psockfd, plen
    REAL(KIND=8), INTENT(IN), TARGET        :: fdata(plen)

      CALL writebuffer_csocket(psockfd, c_loc(fdata(1)), 8*plen)
  END SUBROUTINE

  SUBROUTINE readbuffer_d (psockfd, fdata)
      USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    REAL(KIND=8), INTENT(OUT)               :: fdata

    REAL(KIND=C_DOUBLE), TARGET              :: cdata

      CALL readbuffer_csocket(psockfd, c_loc(cdata), 8)
      fdata=cdata
  END SUBROUTINE

  SUBROUTINE readbuffer_i (psockfd, fdata)
      USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    INTEGER, INTENT(OUT)                     :: fdata

    INTEGER(KIND=C_INT), TARGET              :: cdata

      CALL readbuffer_csocket(psockfd, c_loc(cdata), 4)
      fdata = cdata
  END SUBROUTINE

  SUBROUTINE readbuffer_s (psockfd, fstring, plen)
      USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    CHARACTER(LEN=*), INTENT(OUT)            :: fstring
    INTEGER, INTENT(IN)                      :: plen

    INTEGER                                  :: i
    CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cstring(plen)

      CALL readbuffer_csocket(psockfd, c_loc(cstring(1)), plen)
      fstring=""   
      DO i = 1,plen
         fstring(i:i) = cstring(i)
      ENDDO
  END SUBROUTINE

  SUBROUTINE readbuffer_dv(psockfd, fdata, plen)
      USE ISO_C_BINDING  
    INTEGER, INTENT(IN)                      :: psockfd, plen
    REAL(KIND=8), INTENT(OUT), TARGET       :: fdata(plen)

      CALL readbuffer_csocket(psockfd, c_loc(fdata(1)), 8*plen)
  END SUBROUTINE
  END MODULE

"""

patch_socket_c = """/* A minimal wrapper for socket communication.

Copyright (C) 2013, Joshua More and Michele Ceriotti

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


Contains both the functions that transmit data to the socket and read the data
back out again once finished, and the function which opens the socket initially.
Can be linked to a FORTRAN code that does not support sockets natively.

Functions:
   error: Prints an error message and then exits.
   open_socket_: Opens a socket with the required host server, socket type and
      port number.
   write_buffer_: Writes a string to the socket.
   read_buffer_: Reads data from the socket.
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <netdb.h>

void open_socket(int *psockfd, int* inet, int* port, const char* host)
/* Opens a socket.

Note that fortran passes an extra argument for the string length, but this is
ignored here for C compatibility.

Args:
   psockfd: The id of the socket that will be created.
   inet: An integer that determines whether the socket will be an inet or unix
      domain socket. Gives unix if 0, inet otherwise.
   port: The port number for the socket to be created. Low numbers are often
      reserved for important channels, so use of numbers of 4 or more digits is
      recommended.
   host: The name of the host server.
*/

{
   int sockfd, ai_err;

   if (*inet>0)
   {  // creates an internet socket
      
      // fetches information on the host      
      struct addrinfo hints, *res;  
      char service[256];
   
      memset(&hints, 0, sizeof(hints));
      hints.ai_socktype = SOCK_STREAM;
      hints.ai_family = AF_INET;
      hints.ai_flags = AI_PASSIVE;

      sprintf(service,"%d",*port); // convert the port number to a string
      ai_err = getaddrinfo(host, service, &hints, &res); 
      if (ai_err!=0) { perror("Error fetching host data. Wrong host name?"); exit(-1); }

      // creates socket
      sockfd = socket(res->ai_family, res->ai_socktype, res->ai_protocol);
      if (sockfd < 0) { perror("Error opening socket"); exit(-1); }
    
      // makes connection
      if (connect(sockfd, res->ai_addr, res->ai_addrlen) < 0) 
      { perror("Error opening INET socket: wrong port or server unreachable"); exit(-1); }
      freeaddrinfo(res);
   }
   else
   {  
      struct sockaddr_un serv_addr;

      // fills up details of the socket addres
      memset(&serv_addr, 0, sizeof(serv_addr));
      serv_addr.sun_family = AF_UNIX;
      strcpy(serv_addr.sun_path, "/tmp/ipi_");
      strcpy(serv_addr.sun_path+9, host);
      // creates a unix socket
  
      // creates the socket
      sockfd = socket(AF_UNIX, SOCK_STREAM, 0);

      // connects
      if (connect(sockfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) 
      { perror("Error opening UNIX socket: path unavailable, or already existing"); exit(-1); }
   }


   *psockfd=sockfd;
}

void writebuffer(int *psockfd, const char *data, int* plen)
/* Writes to a socket.

Args:
   psockfd: The id of the socket that will be written to.
   data: The data to be written to the socket.
   plen: The length of the data in bytes.
*/

{
   int n;
   int sockfd=*psockfd;
   int len=*plen;

   n = write(sockfd,data,len);
   if (n < 0) { perror("Error writing to socket: server has quit or connection broke"); exit(-1); }
}


void readbuffer(int *psockfd, char *data, int* plen)
/* Reads from a socket.

Args:
   psockfd: The id of the socket that will be read from.
   data: The storage array for data read from the socket.
   plen: The length of the data in bytes.
*/

{
   int n, nr;
   int sockfd=*psockfd;
   int len=*plen;

   n = nr = read(sockfd,data,len);

   while (nr>0 && n<len )
   {  nr=read(sockfd,&data[n],len-n); n+=nr; }

   if (n == 0) { perror("Error reading from socket: server has quit or connection broke"); exit(-1); }
}



"""


patch_makefile = [
    # 0: Remove accidental socket.o <- socket.f90 dependency
    """## sockets.o should be generated by sockets.c, not corresponding f90 src
F90SRC_all=F90SRC
F90SRC=$(filter-out sockets.f90, $(F90SRC_all))
""",
    # 1: Update rsync
    "rsync -u $(SRCDIR)/*.F $(SRCDIR)/*.inc $(SRCDIR)/*.f90 $(SRCDIR)/*.c .\n",
    # 2: Actual sockets.o <- sockets.c
    """\nsockets.o: sockets.c
\t$(CC_LIB) -c sockets.c
""",
]

patch_reader_F = [
    # 0: Add ihost etc to the incar params
    "     &       ,IHOST,PORT,INET &\n",
    # 1: ipi type configs
    """! -ipi- port, socket type, host
      INTEGER   PORT, INET
      CHARACTER(1024) IHOST
! -ipi- end of addition
""",
    # 2: Add NSW setting for ibrion 23
    "      IF (IBRION==23) NSW=100000 ! default large steps for ipi driver\n",
    # 3: Add ISIF setting for ibrion 23
    """      ! ipi driver mode default calculate stress and allow cell to change
      IF (IBRION==23) ISIF=3
""",
    # No need to update IWAVEPR=1 since it's included in IBRION>0 condition
    # 4: Additional setting for IBRION == 23 (VASP5)
    """! read ipi host & port 
      IF (IBRION==23) THEN
        SZNAM='localhost'
        CALL RDATAB(LOPEN,INCAR,IU5,'IHOST','=','#',';','S', &
        &            IDUM,RDUM,CDUM,LDUM,SZNAM,N,256,IERR)
        IF ((IERR/=0).AND.(IERR/=3)) THEN
            IF (IU0>=0) &
            WRITE(IU0,*)'Error reading item ''IHOST'' from file INCAR.'
            GOTO 150
        ENDIF
        CALL STRIP(SZNAM,N,'A')
        IHOST=SZNAM
        CALL XML_INCAR('IHOST','S',IDUM,RDUM,CDUM,LDUM,SZNAM,N)
        PORT=23333
        CALL RDATAB(LOPEN,INCAR,IU5,'PORT','=','#',';','I', &
        &            PORT,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
        IF (((IERR/=0).AND.(IERR/=3)).OR. &
        &                    ((IERR==0).AND.(N<1))) THEN
            IF (IU0>=0) &
            WRITE(IU0,*)'Error reading item ''PORT'' from file INCAR.'
            GOTO 150
        ENDIF
        CALL XML_INCAR('PORT','I',PORT,RDUM,CDUM,LDUM,CHARAC,N)
        INET=1
        CALL RDATAB(LOPEN,INCAR,IU5,'INET','=','#',';','I', &
        &            INET,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
        IF (((IERR/=0).AND.(IERR/=3)).OR. &
        &                    ((IERR==0).AND.(N<1))) THEN
            IF (IU0>=0) &
            WRITE(IU0,*)'Error reading item ''INET'' from file INCAR.'
            GOTO 150
        ENDIF
        CALL XML_INCAR('INET','I',INET,RDUM,CDUM,LDUM,CHARAC,N)
      ENDIF
""",
    # 5: same as 4 but for VASP 6 (use PROCESS_INCAR), much shorter
    """     
! read ipi host & port 
      IF (IBRION==23) THEN
        SZNAM='localhost'
        CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'IHOST', SZNAM, 1024, IERR, WRITEXMLINCAR)
        CALL STRIP(SZNAM,N,'A')
        IHOST=SZNAM

        PORT=23330
        CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'PORT', PORT, IERR, WRITEXMLINCAR)

        INET=1
        CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'INET', INET, IERR, WRITEXMLINCAR)
      ENDIF
""",
    # 6: Disable symmetry for IBRION 23
    "      IF (IBRION==23) ISYM=0 ! default no symm for ipi driver mode\n",
]



def patch_txt(old_file, pattern, patch_content, replace_func):
    """Backup the old file and write the patch"""
    old_file = Path(old_file)
    # Replace all content of file if not exists or pattern is None
    if (not old_file.is_file()) or (pattern is None):
        txt_new = patch_content
    else:
        txt = open(old_file, "r").read()
        repl = replace_func(patch_content)
        matches = re.findall(pattern, txt, flags=(re.MULTILINE | re.IGNORECASE))
        assert len(matches) == 1
        # print(matches[0])
        txt_new = re.sub(pattern, repl, txt, count=1, flags=(re.MULTILINE | re.IGNORECASE))
        backup_file = old_file.with_suffix(old_file.suffix + ".bak")
        shutil.copy(old_file, backup_file)
    with open(old_file, "w") as fd:
        fd.write(txt_new)
    return


patches = []
patches.append(
    {
        "name": "fsockets.f90",
        "desc": "Create new fsockets.f90",
        "pattern": None,
        "patch_content": patch_fsockets_f90,
        "replace_func": None,
    }
)
patches.append(
    {
        "name": "makefile",
        "desc": "Remove unwanted sockets.o <- sockets.f90 dependency",
        "pattern": r"(^F90SRC\+=main\.f90$\n)([\s\S]*?)(^\s*?$\n)",
        "patch_content": patch_makefile[0],
        "replace_func": insert,
    }
)
patches.append(
    {
        "name": "makefile",
        "desc": "Update rsync",
        "pattern": r"rsync\s+-u\s+\$\(SRCDIR\)\/\*\.F\s+\$\(SRCDIR\)\/\*\.inc.*?\.\s*$\n",
        "patch_content": patch_makefile[1],
        "replace_func": replace,
    }
)
patches.append(
    {
        "name": "makefile",
        "desc": "Actual sockets.o <- sockets.c",
        "pattern": r"(^main\.o \:[\s\S]*?$\n^[\s\S]*?$\n)([\s\S]*?)(^\n)",
        "patch_content": patch_makefile[2],
        "replace_func": insert,
    }
)
patches.append(
    {
        "name": "reader.F",
        "desc": "0: Add ihost etc to the incar params",
        "pattern": r"(SUBROUTINE READER[\s\S]*?#endif\s*?$\n)([\s\S]*?)(^\s*?&\s*?\))",
        "patch_content": patch_reader_F[0],
        "replace_func": insert,
    }
)
# VASP 5 / 6 specific
patches.append(
    {
        "name": "reader.F",
        "desc": "1: ipi type configs",
        "pattern": r"(SUBROUTINE READER[\s\S]*?#endif\s*?$\n)([\s\S]*?)(^\s*?&\s*?\))",
        "patch_content": patch_reader_F[1],
        "replace_func": insert,
    }
)
patches.append(
    {
        "name": "reader.F",
        "desc": "2: Add NSW setting for ibrion 23",
        "pattern": r"(IF \(IBRION.*\) NSW=1\s*$\n)([\s\S]*?)(^\s*CALL PROCESS_INCAR\(.*NSW.*$)",
        "patch_content": patch_reader_F[2],
        "replace_func": insert,
    }
)
patches.append(
    {
        "name": "reader.F",
        "desc": "3: Add ISIF setting for ibrion 23",
        "pattern": r"(IF \(IBRION.*\) ISIF=0\s*$\n)([\s\S]*?)(^\s*CALL PROCESS_INCAR\(.*ISIF)",
        "patch_content": patch_reader_F[3],
        "replace_func": insert,
    }
)
patches.append(
    {
        "name": "reader.F",
        "desc": "5: same as 4 but for VASP 6 (use PROCESS_INCAR), much shorter",
        "pattern": r"(IWAVPR=MOD[\s\S]*#endif\s*\n)([\s\S]*?)(\n^! switch on symmetry)",
        "patch_content": patch_reader_F[5],
        "replace_func": insert,
    }
)
patches.append(
    {
        "name": "reader.F",
        "desc": "6: Disable symmetry for IBRION 23",
        "pattern": r"(^! switch on symmetry[\s\S]*?ISYM=2\s*?$\n)([\s\S]*?)(^\s*CALL PROCESS_INCAR\(.*ISYM)",
        "patch_content": patch_reader_F[6],
        "replace_func": insert,
    }
)




def main():
    parser = ArgumentParser(
        description="Patch VASP source code files for better interactive mode integration"
    )
    parser.add_argument(
        "src", type=str, help="Path to VASP source code, e.g. vasp.X.Y.Z/src"
    )
    args = parser.parse_args()
    src = Path(os.path.expanduser(args.src)).resolve()
    for patch in patches:
        old_file = src / patch["name"]
        patch_content = patch["patch_content"]
        pattern = patch["pattern"]
        func = patch["replace_func"]
        print(patch["desc"])
        print(f"Making patch for {old_file.as_posix()}")
        patch_txt(old_file, pattern, patch_content, func)
    print("Success")
    return


if __name__ == "__main__":
    main()
