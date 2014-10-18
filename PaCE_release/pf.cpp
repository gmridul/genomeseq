#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>

void main_constructor(void) 
    __attribute__((no_instrument_function, constructor));

void main_destructor(void)
    __attribute__((no_instrument_function, destructor));

void __cyg_profile_func_enter(void*, void*)
    __attribute__((no_instrument_function));

void __cyg_profile_func_exit(void*, void*)
    __attribute__ ((no_instrument_function));

static FILE *pf_fp__;

static short int B_ALN__ = 1 + sizeof(void*);
static unsigned int B_SIZE__ = 0;
static char* B_BUF__ = 0;
static unsigned int B_POS__ = 0;

void main_constructor(void) {
    char fn[] = "trace.0000000000";
    char* env = 0;

    pid_t pid = getpid();
    sprintf(fn + 6, "%d", pid);

    pf_fp__ = fopen(fn, "w");

    if (pf_fp__ == NULL) {
	printf("error: could not create %s\n", fn);
	exit(-1);
    }
    
    env = getenv("GLOW_BUFFER");
    B_SIZE__ = 65536;
    B_SIZE__ = env ? (atoi(env) <= B_SIZE__ ? B_SIZE__ : atoi(env)) : B_SIZE__;
    B_BUF__ = (char*)malloc(B_SIZE__ + B_ALN__);
}

void main_destructor(void) {
    if (B_POS__) fwrite(B_BUF__, B_POS__, 1, pf_fp__);
    free(B_BUF__);
    fclose(pf_fp__);
}

void __cyg_profile_func_enter(void* fun, void* caller) {
    B_BUF__[B_POS__] = '>';
    memcpy(B_BUF__ + B_POS__ + 1, &fun, sizeof(void*));
    B_POS__ += B_ALN__;
    if (B_SIZE__ <= B_POS__) {
	fwrite(B_BUF__, B_POS__, 1, pf_fp__);
	B_POS__ = 0;
    }
}

void __cyg_profile_func_exit(void* fun, void* caller) {
    B_BUF__[B_POS__] = '<';
    memcpy(B_BUF__ + B_POS__ + 1, &fun, sizeof(void*));
    B_POS__ += B_ALN__;
    if (B_SIZE__ <= B_POS__) {
	fwrite(B_BUF__, B_POS__, 1, pf_fp__);
	B_POS__ = 0;
    }
}
