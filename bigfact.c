#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>

typedef unsigned long long ull;
typedef unsigned int uint;

const uint LENGTH      = 64;
const uint HALF_LENGTH = 32;

ull  HALF_SIZE  = 1L << 32;

typedef struct BigInt {
    ull* values;
    uint maxExp;
} BigInt;

BigInt* ull2BigInt(ull value);
BigInt* cloneBigInt(BigInt* num);
void destroyBigInt(BigInt* num);
void printBigInt(BigInt* num);

BigInt* bint_addWord(BigInt* lhs, ull rhsVal, uint rhsExp);  /* inplace */
BigInt* bint_mulWord(BigInt* lhs, ull rhsVal, uint rhsExp);  /* inplace */
BigInt* bint_divHalfWord(BigInt* lhs, unsigned long rhsVal, ull* rem);

BigInt* bint_add(BigInt* lhs, BigInt* rhs); /* inplace */
BigInt* bint_mul(BigInt* lhs, BigInt* rhs, uint threads);

ull hi(ull op) {
    return op >> HALF_LENGTH;
}

ull lo(ull op) {
    return ((1L << HALF_LENGTH) - 1 & op);
}

BigInt* ull2BigInt(ull value) {
    BigInt* ret = malloc(sizeof(BigInt));
    ret->values = malloc(sizeof(ull));
    ret->maxExp = 0;

    ret->values[0] = value;

    return ret;
}

BigInt* cloneBigInt(BigInt* num) {
    BigInt* ret = malloc(sizeof(BigInt));
    ret->values = malloc(sizeof(ull) * (num->maxExp + 1));
    ret->maxExp = num->maxExp;

    for (uint i = 0; i<= num->maxExp; ++i) {
        ret->values[i] = num->values[i];
    }

    return ret;
}

void destroyBigInt(BigInt* num) {
    free(num->values);
    free(num);
}

void printBigInt(BigInt* num) {
    const uint DECIMAL_LENGTH = 9;
    const ull  DECIMAL_SIZE   = 1e9;
    uint decimalMaxExp = (uint) ((double) num->maxExp * LENGTH * M_LN2 / M_LN10 / DECIMAL_LENGTH) + 2;

    unsigned long* decimalValues = malloc(sizeof(unsigned long) * (decimalMaxExp + 1));

    BigInt* runningDividend = cloneBigInt(num);
    for (uint i = 0; i <= decimalMaxExp; ++i) {
        ull rem;
        BigInt* quot = bint_divHalfWord(runningDividend, DECIMAL_SIZE, &rem);
        decimalValues[i] = rem;
        destroyBigInt(runningDividend);
        runningDividend = quot;

        if (runningDividend->values[0] == 0 && runningDividend->maxExp == 0){
            decimalMaxExp = i;
            break;
        }
    }
    destroyBigInt(runningDividend);

    printf("%llu", decimalValues[decimalMaxExp]);
    for (uint i = decimalMaxExp; i != 0;) {
        --i;
        printf("%09llu", decimalValues[i]);
    }
    printf("\n");

    free(decimalValues);
}

BigInt* bint_addWord(BigInt* lhs, ull rhsVal, uint rhsExp) {
    uint sumMaxExp = lhs->maxExp + 1;
    if (sumMaxExp < rhsExp) sumMaxExp = rhsExp;
    lhs->values = realloc(lhs->values, sizeof(ull) * (sumMaxExp + 1));

    for (uint i = lhs->maxExp + 1; i <= sumMaxExp; ++i) lhs->values[i] = 0;

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

ull mulOverflow(ull a, ull b) {
    ull overflow;
    ull x;

    x = lo(a) * lo(b);
    x = lo(a) * hi(b) + hi(x);
    overflow = hi(x);

    x = hi(a) * lo(b) + lo(x);
    overflow += hi(a) * hi(b) + hi(x);

    return overflow;
}

BigInt* bint_mulWord(BigInt* lhs, ull rhsVal, uint rhsExp) {
    uint prodMaxExp = lhs->maxExp + rhsExp + 1;
    lhs->values = realloc(lhs->values, sizeof(ull) * (prodMaxExp + 1));

    for (uint i = lhs->maxExp + 1; i <= prodMaxExp; ++i) {
        lhs->values[i] = 0;
    }
    lhs->maxExp = prodMaxExp;

    for (uint i = prodMaxExp - rhsExp; i != 0;) {
        --i;
        ull lhsVal = lhs->values[i];
        ull prod   = lhsVal * rhsVal;
        if (rhsVal != 0 && prod/rhsVal != lhsVal) {
            // Overflow
            bint_addWord(lhs, mulOverflow(lhsVal, rhsVal), i+rhsExp+1);
        }
        lhs->values[i+rhsExp] = prod;
    }
    for (uint i = 0; i < rhsExp; ++i) {
        lhs->values[i] = 0;
    }

    while (lhs->maxExp != 0 && lhs->values[lhs->maxExp] == 0) lhs->maxExp -= 1;

    return lhs;
}

BigInt* bint_divHalfWord(BigInt* lhs, unsigned long rhsVal, ull* rem) {
    BigInt* quot = malloc(sizeof(BigInt));
    quot->values = malloc(sizeof(ull) * (lhs->maxExp + 1));
    quot->maxExp = lhs->maxExp;

    quot->values[quot->maxExp] = 0;

    ull remainder = 0;
    for (uint i = lhs->maxExp + 1; i != 0;) {
        --i;

        ull dividend, quotient;

        dividend  = (remainder << HALF_LENGTH) + hi(lhs->values[i]);
        quotient  = dividend / rhsVal;
        remainder = dividend % rhsVal;

        quotient <<= HALF_LENGTH;

        dividend  = (remainder << HALF_LENGTH) + lo(lhs->values[i]);
        quotient += dividend / rhsVal;
        remainder = dividend % rhsVal;

        quot->values[i] = quotient;
    }

    while (quot->maxExp != 0 && quot->values[quot->maxExp] == 0) quot->maxExp -= 1;

    if (rem) *rem = remainder;
    return quot;
}

BigInt* bint_add(BigInt* lhs, BigInt* rhs) {
    for (uint i = 0; i <= rhs->maxExp; ++i) {
        bint_addWord(lhs, rhs->values[i], i);
    }
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
        BigInt* step = cloneBigInt(lhs);
        bint_mulWord(step, rhs->values[i], i);
        bint_add(ret, step);
        destroyBigInt(step);
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
        destroyBigInt(partial);
    }

    free(threadPool);
    return ret;
}

typedef struct thread_PartialFactInfo {
    uint threads;
    uint offset;
    uint value;
} thread_PartialFactInfo;

void* thread_partialFact(void* _info) {
    thread_PartialFactInfo* info = (thread_PartialFactInfo*) _info;
    uint threads = info->threads;
    uint offset  = info->offset;
    uint value   = info->value;
    free(_info);

    BigInt* ret = ull2BigInt(1);
    for (uint i = 1 + offset; i <= value; i += threads) {
        bint_mulWord(ret, i, 0);
    }

    return (void*) ret;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "Requires at least one argument\n");
        return 1;
    }
    uint value = atoi(argv[1]);
    uint threads = 4;
    if (argc >= 3) threads = atoi(argv[2]);

    pthread_t* threadPool = malloc(sizeof(pthread_t) * threads);
    for (int i = 0; i < threads; ++i) {
        thread_PartialFactInfo* info = malloc(sizeof(thread_PartialFactInfo));
        info->threads = threads;
        info->offset  = i;
        info->value   = value;

        pthread_create(&threadPool[i], NULL, thread_partialFact, (void*) info);
    }

    BigInt* ret;
    pthread_join(threadPool[0], (void**) &ret);

    for (int i = 1; i < threads; ++i) {
        BigInt* partial;
        pthread_join(threadPool[i], (void**) &partial);

        BigInt* prod = bint_mul(ret, partial, threads);
        destroyBigInt(ret);
        destroyBigInt(partial);

        ret = prod;
    }

    printBigInt(ret);

    destroyBigInt(ret);
    free(threadPool);

    return 0;
}
