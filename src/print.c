
#include "print.h"

void printout(const char *str, ...) {
	va_list argp;
	va_start(argp, str);
	#ifdef __useR__
	Rvprintf(str, argp);
	R_FlushConsole();
	#else
	vfprintf(stdout, str, argp);
	fflush(stdout);
	#endif
 	va_end(argp);
}

void printerr(const char *str, ...) {
	va_list argp;
	va_start(argp, str);
	#ifdef __useR__
	REvprintf(str, argp);
	R_FlushConsole();
	#else
	//vfprintf(stderr, str, argp);
	//fflush(stderr);
	#endif
 	va_end(argp);
}

