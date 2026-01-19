/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

int main(int argc, char **argv) {
    int key = 0, port = 0, sfd = 0;
    FILE *escript = NULL;
    struct sockaddr_in sa;

    if (argc < 4) {
        fprintf(stderr, "Missing minimum arguments: %s port key cmd [args]\n",
                argv[0]);
        return 1;
    }
    key = atoi(argv[2]);
    port = atoi(argv[1]);
    
    
    sa.sin_family = AF_INET;
    sa.sin_port = htons(port);
    sa.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
    memset(sa.sin_zero, '\0', sizeof(sa.sin_zero));
    
    sfd = socket(PF_INET, SOCK_STREAM, 0);
    if (sfd < 0) {
        perror("overlord socket creation failed");
        return 1;
    }

    if (connect(sfd, (struct sockaddr*)&sa, sizeof(sa)) < 0) {
        perror("overlord connect() call failed");
        return 1;
    }
    
    escript = fdopen(sfd, "w");
    if (escript == NULL) {
        perror("overlord failed to open file descriptor for writes");
        return 1;
    }
    if (fwrite(&key, sizeof(int), 1, escript) != 1) {
        fprintf(stderr, "overlord failed to initialise communication with escript\n");
        return 1;
    }
        
    fflush(escript);
    execvp(argv[3], argv+3);
    perror("overlord exec failed");
    return 1;
}

