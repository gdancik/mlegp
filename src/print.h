
#ifndef __print__
#define __print__

#include <stdarg.h>

#ifdef __useR__
#include <R_ext/Print.h>
#include <R.h>
#else 
#include <stdio.h>
#endif

void printout(const char *str, ...);
void printerr(const char *str, ...);

#endif


