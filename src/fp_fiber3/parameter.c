/*---------------------------------------------------------------------------*/

/*  parameter.c  */

/*  Copyright (C) 2005 Takashi Uneyama

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions
    are met:

      1. Redistributions of source code must retain the above copyright
         notice, this list of conditions and the following disclaimer.
      2. Redistributions in binary form must reproduce the above copyright
         notice, this list of conditions and the following disclaimer in the
         documentation and/or other materials provided with the distribution.
      3. The name of the author may not be used to endorse or promote
         products derived from this software without specific prior written
         permission.

    THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
    WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */

/*---------------------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include "parameter.h"

#define BUFFER_SIZE 0x100 //256

/*---------------------------------------------------------------------------*/
static void search_and_read_parameter(FILE *fp,char *name,input_parameter *input_parameters);

/*---------------------------------------------------------------------------*/
void read_input_parameters(FILE *fp,input_parameter *input_parameters)
{
    char name[BUFFER_SIZE];

    while(fscanf(fp,"%s",name) != EOF)
    {
        search_and_read_parameter(fp,name,input_parameters);
        // printf("complete search parameters in %s \n",name);
    }
}

/*---------------------------------------------------------------------------*/
void search_and_read_parameter(FILE *fp,char *name,input_parameter *input_parameters)
{
    while(input_parameters->name != NULL)
    {
        if(strcmp(input_parameters->name,name) == 0)
        {
            fscanf(fp,input_parameters->specifier,input_parameters->pointer);
            return;
        }
        input_parameters++;
    }
    fprintf(stderr,"warning: unknown input parameter \"%s\" has found and ignored.\n",name);
}

/*---------------------------------------------------------------------------*/



