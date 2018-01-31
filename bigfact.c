#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "bint.h"

typedef struct thread_PartialFactInfo {
    uint       threads;
    uint       operand;
    bint_exp_t offset;
} thread_PartialFactInfo;

void* thread_partialFact(void* _info) {
    thread_PartialFactInfo* info = (thread_PartialFactInfo*) _info;
    uint       threads = info->threads;
    uint       operand = info->operand;
    bint_exp_t offset  = info->offset;
    free(_info);

    bint_t ret = bint_fromWord(1);
    for (uint i = 1 + offset; i <= operand; i += threads) {
        bint_mulWord(ret, i, 0);
    }

    return (void*) ret;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "Requires at least one argument.\n");
        fprintf(stderr, "Format is: bigfact [operand] [thread count]\n");
        return 1;
    }
    uint operand = atoi(argv[1]);
    uint threads = 4;
    if (argc >= 3) threads = atoi(argv[2]);

    pthread_t* threadPool = malloc(sizeof(pthread_t) * threads);
    for (int i = 0; i < threads; ++i) {
        thread_PartialFactInfo* info = malloc(sizeof(thread_PartialFactInfo));
        info->threads = threads;
        info->operand = operand;
        info->offset  = i;

        pthread_create(&threadPool[i], NULL, thread_partialFact, (void*) info);
    }

    bint_t ret;
    pthread_join(threadPool[0], (void**) &ret);

    for (int i = 1; i < threads; ++i) {
        bint_t partial;
        pthread_join(threadPool[i], (void**) &partial);

        bint_t prod = bint_mulThreaded(ret, partial, threads);
        bint_destroy(ret);
        bint_destroy(partial);

        ret = prod;
    }

    bint_printThreaded(ret, threads);

    bint_destroy(ret);
    free(threadPool);

    return 0;
}
