#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>
#include <math.h>   // M_LN10, M_LN2
#include <string.h> // memcpy, memset
#include <pthread.h>

#define bool int
#define false 0
#define true !(false)

typedef unsigned long long ull;
typedef unsigned int uint;

const uint WORD_LENGTH = 64;
const ull  WORD_MAX    = ULLONG_MAX;
const uint HALF_LENGTH = 32;
const ull  HALF_SIZE   = 1L << 32;

const uint DECIMAL_LENGTH = 19;
const ull  DECIMAL_SIZE   = 1e19;

typedef struct BigInt {
    ull* values;
    uint length;
} BigInt;

BigInt* bint_fromWord(ull value);
BigInt* bint_clone(BigInt* num);
void bint_destroy(BigInt* num);
void bint_print(BigInt* num);
void bint_printThreaded(BigInt* num, uint threads);

BigInt* bint_addWord(BigInt* lhs, ull rhsVal, uint rhsExp);            /* inplace */
BigInt* bint_subWord(BigInt* lhs, ull rhsVal, uint rhsExp, int* neg);  /* inplace */
BigInt* bint_mulWord(BigInt* lhs, ull rhsVal, uint rhsExp);            /* inplace */
BigInt* bint_divWord(BigInt* lhs, ull rhsVal, ull* rem);               /* inplace */

BigInt* bint_add(BigInt* lhs, BigInt* rhs);                            /* inplace */
BigInt* bint_sub(BigInt* lhs, BigInt* rhs, int* neg);                  /* inplace */
BigInt* bint_mul(BigInt* lhs, BigInt* rhs, uint threads);
BigInt* bint_div(BigInt* lhs, BigInt* rhs, BigInt** rem);

BigInt* bint_mulClassical(BigInt* lhs, BigInt* rhs);
BigInt* bint_mulKaratsuba(BigInt* lhs, BigInt* rhs);

BigInt* bint_leftWordShift(BigInt* lhs, uint rhs);                     /* inplace */
BigInt* bint_rightWordShift(BigInt* lhs, uint rhs);                    /* inplace */

BigInt* bint_radMulWord(BigInt* lhs, ull rhsVal, ull rhsExp, ull rad); /* inplace */
BigInt* bint_radAddWord(BigInt* lhs, ull rhsVal, ull rhsExp, ull rad); /* inplace */

BigInt* bint_radAdd(BigInt* lhs, BigInt* rhs, ull rad);                /* inplace */

BigInt* hi(BigInt* num, uint cut) {
    if (cut >= num->length) {
        return bint_fromWord(0);
    }

    BigInt* ret = malloc(sizeof(BigInt));
    uint testSize = num->length - cut;
    ret->length = (num->length < testSize) ? num->length : testSize;
    ret->values = malloc(sizeof(ull) * ret->length);
    memcpy(ret->values, num->values + cut, sizeof(ull) * ret->length);

    return ret;
}

BigInt* lo(BigInt* num, uint cut) {
    if (cut >= num->length) {
        return bint_clone(num);
    }
    BigInt* ret = malloc(sizeof(BigInt));
    ret->length = (num->length < cut) ? num->length : cut;
    ret->values = malloc(sizeof(ull) * ret->length);
    memcpy(ret->values, num->values, sizeof(ull) * ret->length);

    return ret;
}

ull bigMul(ull lhs, ull rhs, ull* _over) {
    ull ret, over;
    asm ("mov %2, %%rax;"
         "mul %3;"
         "mov %%rax, %0;"
         "mov %%rdx, %1;"
         : "=r" (ret), "=r" (over)
         : "r"  (lhs), "r"  (rhs)
         : "%rdx", "%rax"
    );

    if (_over) *_over = over;
    return ret;
}

ull bigDiv(ull lhsHi, ull lhsLo, ull rhs, ull* _rem, int* overflow) {
    if (rhs <= lhsHi) {
        if (overflow) *overflow = true;
        if (_rem) _rem = 0;
        return ULLONG_MAX;
    }
    ull ret, rem;
    asm ("mov %2, %%rdx;"
         "mov %3, %%rax;"
         "div %4;"
         "mov %%rax, %0;"
         "mov %%rdx, %1;"
         : "=r" (ret),   "=r" (rem)
         : "r"  (lhsHi), "r"  (lhsLo), "r" (rhs)
         : "%rdx", "%rax"
    );
    if (overflow) *overflow = false;
    if (_rem) *_rem = rem;
    return ret;
}

BigInt* bint_fromWord(ull value) {
    BigInt* ret = malloc(sizeof(BigInt));
    ret->values = malloc(sizeof(ull));
    ret->length = 1;

    ret->values[0] = value;

    return ret;
}

BigInt* bint_clone(BigInt* num) {
    BigInt* ret = malloc(sizeof(BigInt));
    ret->length = num->length;
    ret->values = malloc(sizeof(ull) * ret->length);
    memcpy(ret->values, num->values, sizeof(ull) * ret->length);

    return ret;
}

void bint_destroy(BigInt* num) {
    free(num->values);
    free(num);
}

typedef struct _thread_print_info {
    BigInt* value;
    ull*    decimalValues;
} _thread_print_info;

void* _thread_print(void* _info) {
    _thread_print_info info = * (_thread_print_info*) _info;
    BigInt* value         = info.value;
    ull*    decimalValues = info.decimalValues;
    free(_info);

    uint i = 0;
    while (true) {
        ull rem;
        bint_divWord(value, DECIMAL_SIZE, &rem);
        decimalValues[i] = rem;

        if (value->values[0] == 0 && value->length == 1) break;
        ++i;
    }

    bint_destroy(value);
    return NULL;
}

void bint_printThreaded(BigInt* num, uint threads) {
    uint decimalLength = (uint) ((double) num->length * WORD_LENGTH * M_LN2 / M_LN10 / DECIMAL_LENGTH) + 1;
    ull* decimalValues = malloc(sizeof(ull) * decimalLength);

    uint decimalLengthPerThread = decimalLength / 4;

    BigInt maxValuePerThread;
    maxValuePerThread.length = 1;
    maxValuePerThread.values = malloc(sizeof(ull) * maxValuePerThread.length);
    maxValuePerThread.values[0] = 1;
    for (uint i = 0; i < decimalLengthPerThread; ++i) {
        bint_mulWord(&maxValuePerThread, 1e19, 0);
    }
    while (maxValuePerThread.values[maxValuePerThread.length - 1] == 0) maxValuePerThread.length -= 1;

    BigInt* runningDividend = bint_clone(num);
    pthread_t* threadPool = malloc(sizeof(pthread_t) * threads);
    for (uint i = 0; i < threads; ++i) {
        _thread_print_info* info = malloc(sizeof(_thread_print_info));
        info->decimalValues = decimalValues + decimalLengthPerThread * i;

        BigInt* oldDividend = runningDividend;
        runningDividend = bint_div(runningDividend, &maxValuePerThread, &info->value);
        bint_destroy(oldDividend);
        pthread_create(&threadPool[i], NULL, _thread_print, (void*) info);
    }
    free(maxValuePerThread.values);

    for (uint i = 0; i < threads; ++i) pthread_join(threadPool[i], NULL);
    free(threadPool);

    while (decimalLength > 1 && decimalValues[decimalLength - 1] == 0) decimalLength -= 1;

    printf("%llu", decimalValues[decimalLength - 1]);;
    for (uint i = decimalLength - 1; i != 0;) {
        --i;
        printf("%019llu", decimalValues[i]);
    }
    printf("\n");

    free(decimalValues);
}

void bint_print(BigInt* num) {
    uint decimalLength = (uint) ((double) num->length * WORD_LENGTH * M_LN2 / M_LN10 / DECIMAL_LENGTH) + 1;
    ull* decimalValues = malloc(sizeof(ull) * decimalLength);

    BigInt* runningDividend = bint_clone(num);
    for (uint i = 0; i < decimalLength; ++i) {
        ull rem;
        bint_divWord(runningDividend, DECIMAL_SIZE, &rem);
        decimalValues[i] = rem;

        if (runningDividend->values[0] == 0 && runningDividend->length == 1){
            decimalLength = i+1;
            break;
        }
    }
    bint_destroy(runningDividend);

    printf("%llu", decimalValues[decimalLength - 1]);;
    for (uint i = decimalLength - 1; i != 0;) {
        --i;
        printf("%019llu", decimalValues[i]);
    }
    printf("\n");

    free(decimalValues);
}

BigInt* bint_addWord(BigInt* lhs, ull rhsVal, uint rhsExp) {
    if (lhs->length < rhsExp + 1) {
        lhs->values = realloc(lhs->values, sizeof(ull) * (rhsExp + 1));
        memset(lhs->values + lhs->length, 0, sizeof(ull) * (rhsExp + 1 - lhs->length));
        lhs->length = rhsExp + 1;
        lhs->values[rhsExp] = rhsVal;

        return lhs;
    }

    lhs->values[rhsExp] += rhsVal;
    bool over = false;
    if (lhs->values[rhsExp] < rhsVal) {
        over = true;
        for (uint i = rhsExp + 1; i < lhs->length; ++i) {
            lhs->values[i] += 1;
            if (lhs->values[i] != 0) {
                over = false;
                break;
            }
        }
    }

    if (over) {
        lhs->length += 1;
        lhs->values = realloc(lhs->values, sizeof(ull) * lhs->length);
        lhs->values[lhs->length - 1] = 1;
    }

    return lhs;
}

BigInt* bint_subWord(BigInt* lhs, ull rhsVal, uint rhsExp, int* _neg) {
    uint diffLength = lhs->length;
    if (lhs->length < rhsExp + 1) {
        diffLength = rhsExp + 1;
        lhs->values = realloc(lhs->values, sizeof(ull) * diffLength);
        memset(lhs->values + lhs->length, 0, (diffLength - lhs->length) * sizeof(ull));
    }

    lhs->values[rhsExp] -= rhsVal;
    bool neg = false;
    if (lhs->values[rhsExp] > rhsVal) {
        neg = true;
        for (uint i = rhsExp + 1; i < diffLength; ++i) {
            lhs->values[i] -= 1;
            if (lhs->values[i] != WORD_MAX) {
                neg = false;
                break;
            }
        }
    }

    if (neg) {
        lhs->length = diffLength;
        for (uint i = 0; i < diffLength; ++i) {
            lhs->values[i] = ~lhs->values[i];
        }
        bint_addWord(lhs, 1, 0);
    }

    while (lhs->length > 1 && lhs->values[lhs->length - 1] == 0) lhs->length -= 1;

    if (_neg) *_neg = neg;
    return lhs;
}

BigInt* bint_mulWord(BigInt* lhs, ull rhsVal, uint rhsExp) {
    BigInt carry;
    carry.length = lhs->length + 1;
    carry.values = malloc(sizeof(ull) * carry.length);
    carry.values[0] = 0;

    for (uint i = 0; i < lhs->length; ++i) {
        ull over;
        lhs->values[i] = bigMul(lhs->values[i], rhsVal, &over);
        carry.values[i+1] = over;
    }

    bint_add(lhs, &carry);
    while (lhs->length > 1 && lhs->values[lhs->length - 1] == 0) lhs->length -= 1;
    bint_leftWordShift(lhs, rhsExp);

    free(carry.values);
    return lhs;
}

BigInt* bint_divWord(BigInt* lhs, ull rhsVal, ull* _rem) {
    ull rem = 0;
    for (uint i = lhs->length; i != 0;) {
        --i;

        ull divHi = rem;
        ull divLo = lhs->values[i];

        lhs->values[i] = bigDiv(divHi, divLo, rhsVal, &rem, NULL);
    }

    while (lhs->length > 1 && lhs->values[lhs->length - 1] == 0) lhs->length -= 1;

    if (_rem) *_rem = rem;
    return lhs;
}

BigInt* bint_add(BigInt* lhs, BigInt* rhs) {

    if (lhs->length > rhs->length) {
        rhs->values = realloc(rhs->values, sizeof(ull) * lhs->length);
        memset(rhs->values + rhs->length, 0, (lhs->length - rhs->length) * sizeof(ull));
    }
    else if (lhs->length < rhs->length) {
        lhs->values = realloc(lhs->values, sizeof(ull) * rhs->length);
        memset(lhs->values + lhs->length, 0, (rhs->length - lhs->length) * sizeof(ull));
        lhs->length = rhs->length;
    }

    ull carrySet;
    ull loops = lhs->length;
    asm ("mov $0, %%rax;"
         "mov %3, %%rcx;"
         "clc;"
         "bint_add_loop%=:"
             "mov (%2, %%rax, 8), %%rbx;"
             "adc %%rbx, (%1, %%rax, 8);"
             "inc %%rax;"
         "loop bint_add_loop%=;"
         "mov $0, %%rax;"
         "adc $0, %%rax;"
         "mov %%rax, %0;"
         : "=r" (carrySet)
         : "r"  (lhs->values), "r" (rhs->values), "r" (loops)
         : "%rax", "%rbx", "%rcx"
        );

    if (carrySet) {
        lhs->length += 1;
        lhs->values = realloc(lhs->values, sizeof(ull) * lhs->length);
        lhs->values[lhs->length - 1] = 1;
    }

    return lhs;
}

BigInt* bint_sub(BigInt* lhs, BigInt* rhs, int* neg) {
    if (lhs->length > rhs->length) {
        rhs->values = realloc(rhs->values, sizeof(ull) * lhs->length);
        memset(rhs->values + rhs->length, 0, (lhs->length - rhs->length) * sizeof(ull));
    }

    ull carrySet;
    ull loops = lhs->length;
    asm ("mov $0, %%rax;"
         "mov %3, %%rcx;"
         "clc;"
         "bint_sub_loop%=:"
             "mov (%2, %%rax, 8), %%rbx;"
             "sbb %%rbx, (%1, %%rax, 8);"
             "inc %%rax;"
         "loop bint_sub_loop%=;"
         "mov $0, %%rax;"
         "adc $0, %%rax;"
         "mov %%rax, %0;"
         : "=r" (carrySet)
         : "r"  (lhs->values), "r" (rhs->values), "r" (loops)
         : "%rax", "%rbx", "%rcx"
        );

    if (carrySet) {
        for (uint i = 0; i < lhs->length; ++i) {
            lhs->values[i] = ~(lhs->values[i]);
        }
        bint_addWord(lhs, 1, 0);
    }

    if (neg) *neg = carrySet;
    return lhs;
}

/*
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
*/

BigInt* bint_div(BigInt* _lhs, BigInt* _rhs, BigInt** rem) {
    if (_lhs->length < _rhs->length) {
        if (rem) *rem = bint_clone(_lhs);
        return bint_fromWord(0);
    }
    BigInt* lhs = bint_clone(_lhs);
    BigInt* rhs = bint_clone(_rhs);

    ull normFactor;
    normFactor = WORD_MAX / rhs->values[rhs->length - 1];
    if (rhs->values[rhs->length - 1] < HALF_SIZE) normFactor = 1;

    bint_mulWord(lhs, normFactor, 0);
    bint_mulWord(rhs, normFactor, 0);

    assert(rhs->length == _rhs->length);

    BigInt* quot = malloc(sizeof(BigInt));
    quot->length = _lhs->length - rhs->length + 1;
    quot->values = malloc(sizeof(ull) * quot->length);

    if (lhs->length == _lhs->length) {
        lhs->length += 1;
        lhs->values = realloc(lhs->values, sizeof(ull) * lhs->length);
        lhs->values[lhs->length - 1] = 0;
    }

    assert(lhs->length == _lhs->length + 1);

    for (uint i = quot->length; i != 0;) {
        --i;
        ull dividendHi = lhs->values[i + rhs->length];
        ull dividendLo = lhs->values[i + rhs->length - 1];
        ull divisor    = rhs->values[rhs->length - 1];

        ull testQuot = bigDiv(dividendHi, dividendLo, divisor, NULL, NULL);
        BigInt* testVal = bint_mulWord(bint_clone(rhs), testQuot, 0);

        BigInt pseudoShifted;
        pseudoShifted.length = rhs->length + 1;
        pseudoShifted.values = lhs->values + i;

        int neg;
        bint_sub(&pseudoShifted, testVal, &neg);

        int count = 0;
        while (neg) {
            assert(++count <= 2);
            while (pseudoShifted.length > 1 && pseudoShifted.values[pseudoShifted.length - 1] == 0) pseudoShifted.length -= 1;
            bint_sub(&pseudoShifted, rhs, &neg);
            neg = !neg;
            --testQuot;
        }

        quot->values[i] = testQuot;
        bint_destroy(testVal);
    }

    if (rem) {
        while (lhs->length > 1 && lhs->values[lhs->length - 1] == 0) lhs->length -= 1;
        *rem = bint_divWord(lhs, normFactor, NULL);
    }
    else bint_destroy(lhs);

    bint_destroy(rhs);

    while (quot->length > 1 && quot->values[quot->length - 1] == 0) quot->length -= 1;

    return quot;
}

BigInt* bint_mulClassical(BigInt* lhs, BigInt* rhs) {
    BigInt* ret = bint_fromWord(0);
    ret->values = realloc(ret->values, sizeof(ull) * (lhs->length + lhs->length));

    for (uint i = 0; i < lhs->length; ++i) {
        for (uint j = 0; j < rhs->length; ++j) {
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
    if (lhs->length < CUTOFF || rhs->length < CUTOFF) {
        return bint_mulClassical(lhs, rhs);
    }
    BigInt* ret;

    uint lhsHalfLength = lhs->length/2;
    uint rhsHalfLength = rhs->length/2;
    uint halfLength = (lhsHalfLength > rhsHalfLength) ? lhsHalfLength : rhsHalfLength;

    BigInt* lhsLo = lo(lhs, halfLength);
    BigInt* lhsHi = hi(lhs, halfLength);
    BigInt* rhsLo = lo(rhs, halfLength);
    BigInt* rhsHi = hi(rhs, halfLength);

    BigInt *p0, *p1, *p2;

    p0 = bint_mulKaratsuba(lhsLo, rhsLo);
    p2 = bint_mulKaratsuba(lhsHi, rhsHi);

    ret = bint_clone(p2);
    bint_leftWordShift(ret, 2*halfLength);
    bint_add(ret, p0);

    bint_add(lhsHi, lhsLo);
    bint_add(rhsHi, rhsLo);

    p1 = bint_mulKaratsuba(lhsHi, rhsHi);
    bint_sub(p1, p0, NULL);
    bint_sub(p1, p2, NULL);

    bint_leftWordShift(p1, halfLength);
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

BigInt* bint_leftWordShift(BigInt* lhs, uint words) {
    lhs->length += words;
    lhs->values = realloc(lhs->values, sizeof(ull) * lhs->length);

    for (uint i = lhs->length; i != words;) {
        --i;
        lhs->values[i] = lhs->values[i - words];
    }

    for (uint i = words; i != 0;) {
        --i;
        lhs->values[i] = 0;
    }

    return lhs;
}

BigInt* bint_radMulWord(BigInt* lhs, ull rhsVal, ull rhsExp, ull rad) {
}

BigInt* bint_radAddWord(BigInt* lhs, ull rhsVal, ull rhsExp, ull rad) {
}

BigInt* bint_radAdd(BigInt* lhs, BigInt* rhs, ull rad) {
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

    while (ret->length > 1 && ret->values[ret->length - 1] == 0) ret->length -= 1;
    bint_print(ret);

    bint_destroy(ret);
    free(threadPool);

    return 0;
}
