#include <stdlib.h>
#include <stdio.h>
#include <math.h>   // M_LN10, M_LN2
#include <string.h> // memcpy, memset
#include <pthread.h>

typedef unsigned long long ull;
typedef unsigned int uint;

const uint WORD_LENGTH = 64;
const uint HALF_LENGTH = 32;
const ull  HALF_SIZE   = 1L << 32;

const uint DECIMAL_LENGTH = 19;
const ull  DECIMAL_SIZE   = 1e19;

typedef struct BigInt {
    ull* values;
    uint maxExp;
} BigInt;

BigInt* bint_fromWord(ull value);
BigInt* bint_clone(BigInt* num);
void bint_destroy(BigInt* num);
void bint_print(BigInt* num);

BigInt* bint_addWord(BigInt* lhs, ull rhsVal, uint rhsExp); /* inplace */
BigInt* bint_mulWord(BigInt* lhs, ull rhsVal, uint rhsExp); /* inplace */
BigInt* bint_divWord(BigInt* lhs, ull rhsVal, ull* rem);    /* inplace */

BigInt* bint_add(BigInt* lhs, BigInt* rhs);            /* inplace */
BigInt* bint_sub(BigInt* lhs, BigInt* rhs, int* neg);  /* inplace */
BigInt* bint_mul(BigInt* lhs, BigInt* rhs, uint threads);
BigInt* bint_mulClassical(BigInt* lhs, BigInt* rhs);
BigInt* bint_mulKaratsuba(BigInt* lhs, BigInt* rhs);

BigInt* bint_leftWordShift(BigInt* lhs, uint rhs);  /* inplace */
BigInt* bint_rightWordShift(BigInt* lhs, uint rhs); /* inplace */

BigInt* hi(BigInt* num, uint cutoff) {
    if (cutoff >= num->maxExp) {
        return bint_fromWord(0);
    }

    BigInt* ret = malloc(sizeof(BigInt));
    uint testSize = num->maxExp - cutoff - 1;
    ret->maxExp = (num->maxExp < testSize) ? num->maxExp : testSize;
    ret->values = malloc(sizeof(ull) * (ret->maxExp + 1));
    memcpy(ret->values, num->values + cutoff + 1, sizeof(ull) * (ret->maxExp + 1));

    return ret;
}

BigInt* lo(BigInt* num, uint cutoff) {
    if (cutoff >= num->maxExp) {
        return bint_clone(num);
    }
    BigInt* ret = malloc(sizeof(BigInt));
    ret->maxExp = (num->maxExp < cutoff) ? num->maxExp : cutoff;
    ret->values = malloc(sizeof(ull) * (ret->maxExp + 1));
    memcpy(ret->values, num->values, sizeof(ull) * (ret->maxExp + 1));

    return ret;
}

ull bigMul(ull lhs, ull rhs, ull* _over) {
    ull ret, over;
    asm ("movq %2, %%rax;"
         "mulq %3;"
         "movq %%rax, %0;"
         "movq %%rdx, %1;"
         : "=r" (ret), "=r" (over)
         : "r"  (lhs), "r"  (rhs)
         : "%rdx", "%rax"
    );

    if (_over) *_over = over;
    return ret;
}

ull bigDiv(ull lhsHi, ull lhsLo, ull rhs, ull* _rem) {
    ull ret, rem;
    asm ("movq %2, %%rdx;"
         "movq %3, %%rax;"
         "divq %4;"
         "movq %%rax, %0;"
         "movq %%rdx, %1;"
         : "=r" (ret),   "=r" (rem)
         : "r"  (lhsHi), "r"  (lhsLo), "r" (rhs)
         : "%rdx", "%rax"
    );
    if (_rem) *_rem = rem;
    return ret;
}

BigInt* bint_fromWord(ull value) {
    BigInt* ret = malloc(sizeof(BigInt));
    ret->values = malloc(sizeof(ull));
    ret->maxExp = 0;

    ret->values[0] = value;

    return ret;
}

BigInt* bint_clone(BigInt* num) {
    BigInt* ret = malloc(sizeof(BigInt));
    ret->values = malloc(sizeof(ull) * (num->maxExp + 1));
    ret->maxExp = num->maxExp;

    memcpy(ret->values, num->values, sizeof(ull) * (num->maxExp + 1));

    return ret;
}

void bint_destroy(BigInt* num) {
    free(num->values);
    free(num);
}

void bint_print(BigInt* num) {
    uint decimalMaxExp = (uint) ((double) num->maxExp * WORD_LENGTH * M_LN2 / M_LN10 / DECIMAL_LENGTH) + 1;

    ull* decimalValues = malloc(sizeof(ull) * (decimalMaxExp + 1));

    BigInt* runningDividend = bint_clone(num);
    for (uint i = 0; i <= decimalMaxExp; ++i) {
        ull rem;
        bint_divWord(runningDividend, DECIMAL_SIZE, &rem);
        decimalValues[i] = rem;

        if (runningDividend->values[0] == 0 && runningDividend->maxExp == 0){
            decimalMaxExp = i;
            break;
        }
    }
    bint_destroy(runningDividend);

    printf("%llu", decimalValues[decimalMaxExp]);
    for (uint i = decimalMaxExp; i != 0;) {
        --i;
        printf("%019llu", decimalValues[i]);
    }
    printf("\n");

    free(decimalValues);
}

BigInt* bint_addWord(BigInt* lhs, ull rhsVal, uint rhsExp) {
    uint sumMaxExp = lhs->maxExp + 1;
    if (sumMaxExp < rhsExp) sumMaxExp = rhsExp;
    lhs->values = realloc(lhs->values, sizeof(ull) * (sumMaxExp + 1));

    memset(lhs->values + lhs->maxExp + 1, 0, (sumMaxExp - lhs->maxExp) * sizeof(ull));

    lhs->values[rhsExp] += rhsVal;
    if (lhs->values[rhsExp] < rhsVal) {
        // Overflow
        for (uint i = rhsExp+1; i <= sumMaxExp; ++i) {
            lhs->values[i] += 1;
            if (lhs->values[i] != 0) break;
        }
    }

    lhs->maxExp = sumMaxExp;
    while (lhs->maxExp != 0 && lhs->values[lhs->maxExp] == 0) lhs->maxExp -= 1;

    return lhs;
}

BigInt* bint_mulWord(BigInt* lhs, ull rhsVal, uint rhsExp) {
    uint prodMaxExp = lhs->maxExp + rhsExp + 1;
    lhs->values = realloc(lhs->values, sizeof(ull) * (prodMaxExp + 1));

    memset(lhs->values + lhs->maxExp + 1, 0, (prodMaxExp - lhs->maxExp) * sizeof(ull));

    lhs->maxExp = prodMaxExp;

    for (uint i = prodMaxExp - rhsExp; i != 0;) {
        --i;
        ull lhsVal = lhs->values[i];
        ull over;
        lhs->values[i+rhsExp] = bigMul(lhsVal, rhsVal, &over);
        if (over != 0) {
            // Overflow
            bint_addWord(lhs, over, i+rhsExp+1);
        }
    }
    for (uint i = 0; i < rhsExp; ++i) {
        lhs->values[i] = 0;
    }

    while (lhs->maxExp != 0 && lhs->values[lhs->maxExp] == 0) lhs->maxExp -= 1;

    return lhs;
}

BigInt* bint_divWord(BigInt* lhs, ull rhsVal, ull* _rem) {
    ull rem = 0;
    for (uint i = lhs->maxExp + 1; i != 0;) {
        --i;

        ull divHi = rem;
        ull divLo = lhs->values[i];

        lhs->values[i] = bigDiv(divHi, divLo, rhsVal, &rem);
    }

    while (lhs->maxExp != 0 && lhs->values[lhs->maxExp] == 0) lhs->maxExp -= 1;

    if (_rem) *_rem = rem;
    return lhs;
}

BigInt* bint_add(BigInt* lhs, BigInt* rhs) {
    uint sumMaxExp = (lhs->maxExp > rhs->maxExp) ? lhs->maxExp : rhs->maxExp;
    ++sumMaxExp;

    lhs->values = realloc(lhs->values, sizeof(ull) * (sumMaxExp + 1));
    memset(lhs->values + lhs->maxExp + 1, 0, (sumMaxExp - lhs->maxExp) * sizeof(ull));

    rhs->values = realloc(rhs->values, sizeof(ull) * (sumMaxExp + 1));
    memset(rhs->values + rhs->maxExp + 1, 0, (sumMaxExp - rhs->maxExp) * sizeof(ull));

    ull loops = sumMaxExp + 1;
    asm ("movq $0, %%rax;"
         "movq %2, %%rcx;"
         "clc;"
         "bint_add_loop%=:"
             "movq (%1, %%rax, 8), %%rbx;"
             "adcq %%rbx, (%0, %%rax, 8);"
             "incq %%rax;"
         "loopq bint_add_loop%=;"
         :
         : "r" (lhs->values), "r" (rhs->values), "r" (loops)
         : "%rax", "%rbx", "%rcx"
        );

    lhs->maxExp = sumMaxExp;
    if (lhs->values[sumMaxExp] == 0) lhs->maxExp -= 1;

    return lhs;
}

BigInt* bint_sub(BigInt* lhs, BigInt* rhs, int* neg) {
    if (lhs->maxExp > rhs->maxExp) {
        rhs->values = realloc(rhs->values, sizeof(ull) * (lhs->maxExp + 1));
        memset(rhs->values + rhs->maxExp + 1, 0, (lhs->maxExp - rhs->maxExp) * sizeof(ull));
    }

    ull carrySet;
    ull loops = lhs->maxExp + 1;
    asm ("movq $0, %%rax;"
         "movq %3, %%rcx;"
         "clc;"
         "bint_sub_loop%=:"
             "movq (%2, %%rax, 8), %%rbx;"
             "sbbq %%rbx, (%1, %%rax, 8);"
             "incq %%rax;"
         "loopq bint_sub_loop%=;"
         "movq $0, %%rax;"
         "adcq $0, %%rax;"
         "movq %%rax, %0;"
         : "=r" (carrySet)
         : "r"  (lhs->values), "r" (rhs->values), "r" (loops)
         : "%rax", "%rbx", "%rcx"
        );

    if (carrySet) {
        for (uint i = 0; i <= lhs->maxExp; ++i) {
            lhs->values[i] = ~(lhs->values[i]);
        }
        bint_addWord(lhs, 1, 0);
    }

    if (neg) *neg = carrySet;
    return lhs;
}

typedef struct thread_PartialMulInfo {
    uint threads;
    uint offset;
    BigInt* lhs;
    BigInt* rhs;
} thread_PartialMulInfo;

void* thread_partialMul(void* _info) {
    thread_PartialMulInfo* info = (thread_PartialMulInfo*) _info;
    uint threads = info->threads;
    uint offset  = info->offset;
    BigInt* lhs  = info->lhs;
    BigInt* rhs  = info->rhs;
    free(_info);

    uint maxExp = lhs->maxExp + rhs->maxExp + 1;

    BigInt* ret = malloc(sizeof(BigInt));
    ret->values = malloc(sizeof(ull) * (maxExp + 1));
    ret->maxExp = 0;
    ret->values[0] = 0;

    for (uint i = offset; i <= rhs->maxExp; i += threads) {
        BigInt* step = bint_clone(lhs);
        bint_mulWord(step, rhs->values[i], i);
        bint_add(ret, step);
        bint_destroy(step);
    }

    return (void*) ret;
}

BigInt* bint_mul(BigInt* lhs, BigInt* rhs, uint threads) {
    pthread_t* threadPool = malloc(sizeof(pthread_t) * threads);
    for (int i = 0; i < threads; ++i) {
        thread_PartialMulInfo* info = malloc(sizeof(thread_PartialMulInfo));
        info->threads = threads;
        info->offset  = i;
        info->lhs     = lhs;
        info->rhs     = rhs;

        pthread_create(&threadPool[i], NULL, thread_partialMul, (void*) info);
    }

    BigInt* ret;
    pthread_join(threadPool[0], (void**) &ret);

    for (int i = 1; i < threads; ++i) {
        BigInt* partial;
        pthread_join(threadPool[i], (void**) &partial);

        bint_add(ret, partial);
        bint_destroy(partial);
    }

    free(threadPool);
    return ret;
}

BigInt* bint_mulClassical(BigInt* lhs, BigInt* rhs) {
    BigInt* ret = bint_fromWord(0);
    ret->values = realloc(ret->values, sizeof(ull) * (lhs->maxExp + lhs->maxExp + 1));

    for (uint i = 0; i <= lhs->maxExp; ++i) {
        for (uint j = 0; j <= rhs->maxExp; ++j) {
            ull prod, over;
            prod = bigMul(lhs->values[i], rhs->values[j], &over);
            bint_addWord(ret, prod, i + j);
            if (over != 0) bint_addWord(ret, over, i + j + 1);
        }
    }

    return ret;
}

BigInt* bint_mulKaratsuba(BigInt* lhs, BigInt* rhs) {
    const uint CUTOFF = 10;
    if (lhs->maxExp < CUTOFF || rhs->maxExp < CUTOFF) {
        return bint_mulClassical(lhs, rhs);
    }
    BigInt* ret;

    uint lhsHalfExp = (lhs->maxExp + 1) / 2 - 1;
    uint rhsHalfExp = (rhs->maxExp + 1) / 2 - 1;
    uint halfExp = (lhsHalfExp > rhsHalfExp) ? lhsHalfExp : rhsHalfExp;

    BigInt* lhsLo = lo(lhs, halfExp);
    BigInt* lhsHi = hi(lhs, halfExp);
    BigInt* rhsLo = lo(rhs, halfExp);
    BigInt* rhsHi = hi(rhs, halfExp);

    BigInt *p0, *p1, *p2;

    p0 = bint_mulKaratsuba(lhsLo, rhsLo);
    p2 = bint_mulKaratsuba(lhsHi, rhsHi);

    ret = bint_clone(p2);
    bint_leftWordShift(ret, 2 * (halfExp + 1));
    bint_add(ret, p0);

    bint_add(lhsHi, lhsLo);
    bint_add(rhsHi, rhsLo);

    p1 = bint_mulKaratsuba(lhsHi, rhsHi);
    bint_sub(p1, p0, NULL);
    bint_sub(p1, p2, NULL);

    bint_leftWordShift(p1, halfExp + 1);
    bint_add(ret, p1);

    bint_destroy(lhsLo);
    bint_destroy(lhsHi);
    bint_destroy(rhsLo);
    bint_destroy(rhsHi);
    bint_destroy(p0);
    bint_destroy(p1);
    bint_destroy(p2);

    return ret;
}

BigInt* bint_leftWordShift(BigInt* lhs, uint rhs) {
    lhs->maxExp += rhs;
    lhs->values = realloc(lhs->values, sizeof(ull) * (lhs->maxExp + 1));

    for (uint i = lhs->maxExp + 1; i != rhs;) {
        --i;
        lhs->values[i] = lhs->values[i - rhs];
    }

    for (uint i = rhs; i != 0;) {
        --i;
        lhs->values[i] = 0;
    }

    return lhs;
}

typedef struct thread_PartialFactInfo {
    uint threads;
    uint offset;
    uint operand;
} thread_PartialFactInfo;

void* thread_partialFact(void* _info) {
    thread_PartialFactInfo* info = (thread_PartialFactInfo*) _info;
    uint threads = info->threads;
    uint offset  = info->offset;
    uint operand = info->operand;
    free(_info);

    BigInt* ret = bint_fromWord(1);
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
        info->offset  = i;
        info->operand = operand;

        pthread_create(&threadPool[i], NULL, thread_partialFact, (void*) info);
    }

    BigInt* ret;
    pthread_join(threadPool[0], (void**) &ret);

    for (int i = 1; i < threads; ++i) {
        BigInt* partial;
        pthread_join(threadPool[i], (void**) &partial);

        BigInt* prod = bint_mulKaratsuba(ret, partial);
        bint_destroy(ret);
        bint_destroy(partial);

        ret = prod;
    }

    bint_print(ret);

    bint_destroy(ret);
    free(threadPool);

    return 0;
}
