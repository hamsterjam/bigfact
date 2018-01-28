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
BigInt* bint_shrink(BigInt* num);
void bint_destroy(BigInt* num);
void bint_print(BigInt* num, uint threads);

BigInt* bint_toDecClassical(BigInt* num);
BigInt* bint_toDecDivAndConq(BigInt* num, uint threads);

bool bint_lessThan(BigInt* lhs, BigInt* rhs);

BigInt* bint_addWord(BigInt* lhs, ull rhsVal, uint rhsExp);            /* inplace */
BigInt* bint_subWord(BigInt* lhs, ull rhsVal, uint rhsExp, int* neg);  /* inplace */
BigInt* bint_mulWord(BigInt* lhs, ull rhsVal, uint rhsExp);            /* inplace */
BigInt* bint_divWord(BigInt* lhs, ull rhsVal, ull* rem);               /* inplace */
ull     bint_modWord(BigInt* lhs, ull rhsVal);

BigInt* bint_add(BigInt* lhs, BigInt* rhs);                            /* inplace */
BigInt* bint_sub(BigInt* lhs, BigInt* rhs, int* neg);                  /* inplace */
BigInt* bint_mul(BigInt* lhs, BigInt* rhs);
BigInt* bint_div(BigInt* lhs, BigInt* rhs, BigInt** rem);

BigInt* bint_mulThreaded(BigInt* lhs, BigInt* rhs, uint threads);
BigInt* bint_divThreaded(BigInt* lhs, BigInt* rhs, BigInt** rem, uint threads);

BigInt* bint_addNoLastCarry(BigInt* lhs, BigInt* rhs, bool* carry);    /* inplace */
BigInt* bint_subNoLastCarry(BigInt* lhs, BigInt* rhs, bool* carry);    /* inplace */

BigInt* bint_mulClassical(BigInt* lhs, BigInt* rhs);
BigInt* bint_mulKaratsuba(BigInt* lhs, BigInt* rhs);
BigInt* bint_divClassical(BigInt* lhs, BigInt* rhs, BigInt** rem);
BigInt* bint_divDivAndConq(BigInt* lhs, BigInt* rhs, BigInt** rem, uint threads);

BigInt* bint_leftShift( BigInt* lhs, uint bits, uint words);            /* inplace */
BigInt* bint_rightShift(BigInt* lhs, uint bits, uint words);            /* inplace */

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
    if (rhs == 0) {
        return lhsLo/rhs;
    }
    if (rhs <= lhsHi) {
        if (overflow) *overflow = true;
        if (_rem) {
            ull radix, resHi, resLo;
            bigDiv(1, 0, rhs, &radix, NULL);

            resLo = bigMul(radix, lhsHi, &resHi);
            resLo += lhsLo;
            if (resLo < lhsLo) ++resHi;

            bigDiv(resHi, resLo, rhs, _rem, NULL);
        }
        return 0;
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

ull div3by2(BigInt* u, BigInt* v, BigInt** _rem) {
    ull rem;
    bool over;
    ull q = bigDiv(u->values[2], u->values[1], v->values[1], &rem, &over);

    BigInt* R;
    if (over) {
        q = WORD_MAX;
        R = bint_leftShift(bint_fromWord(u->values[2]), 0, 1);
        R->values[0] = u->values[1];

        bint_subWord(R, v->values[1], 1, NULL);
        bint_addWord(R, v->values[1], 0);
    }
    else {
        R = bint_fromWord(rem);
    }

    BigInt* D = bint_mulWord(bint_fromWord(q), v->values[0], 0);
    bint_leftShift(R, 0, 1);
    R->values[0] = u->values[0];
    bool neg;
    bint_subNoLastCarry(R, D, &neg);

    while (neg) {
        --q;
        bint_addNoLastCarry(R, v, &neg);
        neg = !neg;
    }

    bint_destroy(D);

    if (_rem) *_rem = R;
    else bint_destroy(R);
    return q;
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

BigInt* bint_shrink(BigInt* num) {
    while (num->length > 1 && num->values[num->length - 1] == 0) num->length -= 1;
    return num;
}

void bint_destroy(BigInt* num) {
    free(num->values);
    free(num);
}

void bint_print(BigInt* num, uint threads) {
    BigInt* dec = bint_toDecDivAndConq(num, threads);

    printf("%llu", dec->values[dec->length - 1]);;
    for (uint i = dec->length - 1; i != 0;) {
        --i;
        printf("%019llu", dec->values[i]);
    }
    printf("\n");

    bint_destroy(dec);
}


BigInt* bint_toDecClassical(BigInt* num) {
    const double binToDecLengthFactor = (double) WORD_LENGTH * M_LN2 / M_LN10 / DECIMAL_LENGTH;

    BigInt* ret = malloc(sizeof(BigInt));
    ret->length = (double) num->length * binToDecLengthFactor + 1;
    ret->values = malloc(sizeof(ull) * ret->length);

    BigInt* dividend = bint_clone(num);
    for (uint i = 0; i < ret->length; ++i) {
        ull rem;
        bint_divWord(dividend, DECIMAL_SIZE, &rem);
        ret->values[i] = rem;
    }
    bint_destroy(dividend);

    return bint_shrink(ret);
}

typedef struct _info_bint_toDecDivAndConq {
    BigInt*          num;
    pthread_mutex_t* poolLock;
    pthread_t*       threadPool;
    uint*            currThreads;
    uint             maxThreads;
} _info_bint_toDecDivAndConq;

void* _thread_bint_toDecDivAndConq(void* _info) {
    _info_bint_toDecDivAndConq* info = (_info_bint_toDecDivAndConq*) _info;

    BigInt*          num         = info->num;
    pthread_mutex_t* poolLock    = info->poolLock;
    uint*            currThreads = info->currThreads;
    uint             maxThreads  = info->maxThreads;

    free(_info);

    const double binToDecLengthFactor = (double) WORD_LENGTH * M_LN2 / M_LN10 / DECIMAL_LENGTH;
    const int CUTOFF = 10;
    if (num->length < CUTOFF) {
        return (void*) bint_toDecClassical(num);
    }

    BigInt* ret = malloc(sizeof(BigInt));
    ret->length = (double) num->length * binToDecLengthFactor + 1;
    ret->values = malloc(sizeof(ull) * ret->length);
    memset(ret->values, 0, sizeof(ull) * ret->length);

    uint wordLength = ret->length / 2;
    BigInt* divisor = bint_fromWord(1e19);
    for (uint i = 1; i < wordLength; ++i) {
        bint_mulWord(divisor, 1e19, 0);
    }

    bool splitThread = false;
    pthread_mutex_lock(poolLock);
    uint threadsPerThread = maxThreads / *currThreads;
    if (*currThreads < maxThreads) {
        splitThread = true;
        *currThreads += 1;
    }
    pthread_mutex_unlock(poolLock);

    BigInt *numLo, *numHi;
    numHi = bint_divThreaded(num, divisor, &numLo, threadsPerThread);

    bint_destroy(divisor);


    BigInt *retHi, *retLo;

    _info_bint_toDecDivAndConq* infoHi = malloc(sizeof(_info_bint_toDecDivAndConq));
    infoHi->num         = numHi;
    infoHi->poolLock    = poolLock;
    infoHi->currThreads = currThreads;
    infoHi->maxThreads  = maxThreads;

    pthread_t split;
    if (splitThread) pthread_create(&split, NULL, _thread_bint_toDecDivAndConq, (void*) infoHi);
    else retHi = _thread_bint_toDecDivAndConq((void*) infoHi);

    _info_bint_toDecDivAndConq* infoLo = malloc(sizeof(_info_bint_toDecDivAndConq));
    infoLo->num         = numLo;
    infoLo->poolLock    = poolLock;
    infoLo->currThreads = currThreads;
    infoLo->maxThreads  = maxThreads;

    retLo = _thread_bint_toDecDivAndConq((void*) infoLo);

    if (splitThread) {
        pthread_join(split, (void**) &retHi);
        pthread_mutex_lock(poolLock);
        *currThreads -= 1;
        pthread_mutex_unlock(poolLock);
    }

    memcpy(ret->values,              retLo->values, sizeof(ull) * retLo->length);
    memcpy(ret->values + wordLength, retHi->values, sizeof(ull) * retHi->length);

    bint_destroy(numHi);
    bint_destroy(numLo);
    bint_destroy(retHi);
    bint_destroy(retLo);

    return (void*) bint_shrink(ret);
}

BigInt* bint_toDecDivAndConq(BigInt* num, uint threads) {
    pthread_mutex_t poolLock;
    pthread_mutex_init(&poolLock, NULL);

    uint currThreads = 1;

    _info_bint_toDecDivAndConq* info = malloc(sizeof(_info_bint_toDecDivAndConq));
    info->num         = num;
    info->poolLock    = &poolLock;
    info->currThreads = &currThreads;
    info->maxThreads  = threads;

    BigInt* ret = _thread_bint_toDecDivAndConq((void*) info);

    pthread_mutex_destroy(&poolLock);

    return ret;
}

bool bint_lessThan(BigInt* lhs, BigInt* rhs) {
    if (lhs->length != rhs->length) {
        return lhs->length < rhs->length;
    }

    for (uint i = lhs->length; i != 0;) {
        --i;
        if (lhs->values[i] < rhs->values[i]) return true;
        if (lhs->values[i] > rhs->values[i]) return false;
    }
    return false;
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
    bool carry = false;
    if (lhs->values[rhsExp] < rhsVal) {
        carry = true;
        for (uint i = rhsExp + 1; i < lhs->length; ++i) {
            lhs->values[i] += 1;
            if (lhs->values[i] != 0) {
                carry = false;
                break;
            }
        }
    }

    if (carry) {
        lhs->length += 1;
        lhs->values = realloc(lhs->values, sizeof(ull) * lhs->length);
        lhs->values[lhs->length - 1] = 1;
    }

    return bint_shrink(lhs);
}

BigInt* bint_subWord(BigInt* lhs, ull rhsVal, uint rhsExp, int* _neg) {
    uint diffLength = lhs->length;
    if (lhs->length < rhsExp + 1) {
        diffLength = rhsExp + 1;
        lhs->values = realloc(lhs->values, sizeof(ull) * diffLength);
        memset(lhs->values + lhs->length, 0, (diffLength - lhs->length) * sizeof(ull));
    }

    ull oldVal = lhs->values[rhsExp] - rhsVal;
    lhs->values[rhsExp] -= rhsVal;
    bool neg = false;
    if (lhs->values[rhsExp] > oldVal) {
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

    if (_neg) *_neg = neg;
    return bint_shrink(lhs);
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
    bint_shrink(lhs);
    bint_leftShift(lhs, 0, rhsExp);

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

    if (_rem) *_rem = rem;
    return bint_shrink(lhs);
}

ull bint_modWord(BigInt* lhs, ull rhsVal) {
    ull radMod;
    bigDiv(1, 0, rhsVal, &radMod, NULL);
    ull ret = lhs->values[lhs->length - 1] % rhsVal;
    for (uint i = lhs->length - 1; i != 0;) {
        --i;
        ull resHi, resLo;
        resLo = bigMul(ret, radMod, &resHi);

        resLo += lhs->values[i] % rhsVal;
        if (resLo < lhs->values[i] % rhsVal)
            resHi += 1;

        bigDiv(resHi, resLo, rhsVal, &ret, NULL);
    }

    return ret;
}

BigInt* bint_add(BigInt* lhs, BigInt* rhs) {
    bool carry;
    bint_addNoLastCarry(lhs, rhs, &carry);

    if (carry) {
        lhs->length += 1;
        lhs->values = realloc(lhs->values, sizeof(ull) * lhs->length);
        lhs->values[lhs->length - 1] = 1;
    }

    return lhs;
}

BigInt* bint_sub(BigInt* lhs, BigInt* rhs, bool* neg) {
    bool carry;
    bint_subNoLastCarry(lhs, rhs, &carry);

    if (carry) {
        for (uint i = 0; i < lhs->length; ++i) {
            lhs->values[i] = ~(lhs->values[i]);
        }
        bint_addWord(lhs, 1, 0);
    }

    if (neg) *neg = carry;
    return lhs;
}

BigInt* bint_mul(BigInt* lhs, BigInt* rhs) {
    return bint_mulKaratsuba(lhs, rhs);
}

BigInt* bint_div(BigInt* lhs, BigInt* rhs, BigInt** rem) {
    return bint_divThreaded(lhs, rhs, rem, 1);
}

void* _thread_bint_mulKaratsuba(void* _args) {
    BigInt** args = (BigInt**) _args;

    BigInt* lhs = args[0];
    BigInt* rhs = args[1];
    free(args);

    return (void*) bint_mulKaratsuba(lhs, rhs);
}

BigInt* bint_mulThreaded(BigInt* lhs, BigInt* rhs, uint threads) {
    if (threads == 1) return bint_mul(lhs, rhs);

    uint lengthPerThread = rhs->length / 4;

    BigInt* sliceHi = bint_clone(rhs);

    BigInt** slices = malloc(sizeof(BigInt*) * threads);
    for (int i = 0; i < threads - 1; ++i) {
        BigInt* oldSliceHi = sliceHi;
        slices[i] = lo(sliceHi, lengthPerThread);
        sliceHi   = hi(sliceHi, lengthPerThread);

        bint_destroy(oldSliceHi);
    }
    slices[threads - 1] = sliceHi;

    pthread_t* threadPool = malloc(sizeof(pthread_t) * threads);
    for (int i = 0; i < threads; ++i) {
        BigInt** args = malloc(sizeof(BigInt*) * 2);
        args[0] = lhs;
        args[1] = slices[i];
        pthread_create(&threadPool[i], NULL, _thread_bint_mulKaratsuba, (void*) args);
    }

    BigInt* ret;
    pthread_join(threadPool[0], (void**) &ret);
    bint_destroy(slices[0]);

    for (int i = 1; i < threads; ++i) {
        BigInt* part;
        pthread_join(threadPool[i], (void**) &part);
        bint_destroy(slices[i]);

        bint_leftShift(part, 0, lengthPerThread * i);
        bint_add(ret, part);

        bint_destroy(part);
    }

    free(slices);
    free(threadPool);
    return ret;
}

BigInt* bint_divThreaded(BigInt* lhs, BigInt* rhs, BigInt** _rem, uint threads) {
    BigInt *quot, *rem;
    quot = bint_divDivAndConq(lhs, rhs, &rem, threads);
    if (_rem) *_rem = rem;
    return quot;
}

BigInt* bint_addNoLastCarry(BigInt* lhs, BigInt* rhs, bool* _carry) {
    if (lhs->length > rhs->length) {
        rhs->values = realloc(rhs->values, sizeof(ull) * lhs->length);
        memset(rhs->values + rhs->length, 0, (lhs->length - rhs->length) * sizeof(ull));
    }
    else if (lhs->length < rhs->length) {
        lhs->values = realloc(lhs->values, sizeof(ull) * rhs->length);
        memset(lhs->values + lhs->length, 0, (rhs->length - lhs->length) * sizeof(ull));
        lhs->length = rhs->length;
    }

    ull carry;
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
         : "=r" (carry)
         : "r"  (lhs->values), "r" (rhs->values), "r" (loops)
         : "%rax", "%rbx", "%rcx"
        );

    if (_carry) *_carry = carry;
    return bint_shrink(lhs);
}

BigInt* bint_subNoLastCarry(BigInt* lhs, BigInt* rhs, bool* _carry) {
    if (lhs->length > rhs->length) {
        rhs->values = realloc(rhs->values, sizeof(ull) * lhs->length);
        memset(rhs->values + rhs->length, 0, (lhs->length - rhs->length) * sizeof(ull));
    }

    ull carry;
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
         : "=r" (carry)
         : "r"  (lhs->values), "r" (rhs->values), "r" (loops)
         : "%rax", "%rbx", "%rcx"
        );

    if (_carry) *_carry = carry;
    return bint_shrink(lhs);
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

    return bint_shrink(ret);
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
    bint_leftShift(ret, 0, 2*halfLength);
    bint_add(ret, p0);

    bint_add(lhsHi, lhsLo);
    bint_add(rhsHi, rhsLo);

    p1 = bint_mulKaratsuba(lhsHi, rhsHi);
    bint_sub(p1, p0, NULL);
    bint_sub(p1, p2, NULL);

    bint_leftShift(p1, 0, halfLength);
    bint_add(ret, p1);

    bint_destroy(lhsLo);
    bint_destroy(lhsHi);
    bint_destroy(rhsLo);
    bint_destroy(rhsHi);
    bint_destroy(p0);
    bint_destroy(p1);
    bint_destroy(p2);

    return bint_shrink(ret);
}

// Based on Knuth, The Art of Computer Programming, Vol II, pg 272, Algorithm D
BigInt* bint_divClassical(BigInt* lhs, BigInt* rhs, BigInt** rem) {
    if (lhs->length < rhs->length) {
        if (rem) *rem = bint_clone(lhs);
        return bint_fromWord(0);
    }

    BigInt* u = bint_clone(lhs);
    BigInt* v = bint_clone(rhs);

    uint n = rhs->length;
    uint m = lhs->length - rhs->length;

    BigInt* q = malloc(sizeof(BigInt));
    q->length = m + 1;
    q->values = malloc(sizeof(ull) * q->length);

    // Normalize (D1)
    uint d = 0;
    for (ull top = v->values[n-1]; top < (1L << 63); top <<= 1) ++d;

    bint_leftShift(u, d, 0);
    bint_leftShift(v, d, 0);

    v->values = realloc(v->values, sizeof(ull) * (n+1));
    v->values[n] = 0;

    if (u->length < n + m + 1) {
        u->length = n + m + 1;
        u->values = realloc(u->values, sizeof(ull) * (n + m + 1));
        u->values[n + m] = 0;
    }

    // Initialize j (D2) / Loop on j (D7)
    for (uint j = m + 1; j != 0;) {
        --j;
        // Calculate q_hat (D3)
        BigInt vHi;
        vHi.length = 2;
        vHi.values = v->values + n-2;

        BigInt uHi;
        uHi.length = 3;
        uHi.values = u->values + j+n-2;

        ull q_hat = div3by2(&uHi, &vHi, NULL);

        // Multiply and subtract (D4)
        BigInt uPart;
        uPart.values = u->values + j;
        uPart.length = n + 1;

        BigInt* qv = bint_clone(v);
        bint_mulWord(qv, q_hat, 0);

        bool neg;
        bint_subNoLastCarry(&uPart, qv, &neg);

        bint_destroy(qv);

        // Test remainder (D5)
        q->values[j] = q_hat;

        if (neg) {
            // Add back (D6)
            q->values[j] -= 1;
            bint_addNoLastCarry(&uPart, v, NULL);
        }
    }

    if (rem) {
        // Unnormalize (D8)
        bint_rightShift(u, d, 0);
        *rem = bint_shrink(u);
    }
    else {
        bint_destroy(u);
    }
    bint_destroy(v);

    return bint_shrink(q);
}

BigInt* bint_divDivAndConq2by1(BigInt* lhs, BigInt* rhs, BigInt** rem, uint threads);
BigInt* bint_divDivAndConq3by2(BigInt* lhs, BigInt* rhs, BigInt** rem, uint threads);

BigInt* bint_divDivAndConq2by1(BigInt* lhs, BigInt* rhs, BigInt** rem, uint threads) {
    const uint CUTOFF = 10;
    if (rhs->length < CUTOFF) {
        return bint_divClassical(lhs, rhs, rem);
    }

    if (lhs->length < rhs->length) {
        if (rem) *rem = bint_clone(lhs);
        return bint_fromWord(0);
    }

    uint wordLength = (rhs->length + 1) / 2;

    BigInt* u = bint_clone(lhs);
    BigInt* v = bint_clone(rhs);

    // Renormalize
    uint D = 0;
    if (wordLength * 2 > v->length) D = 1;

    bint_leftShift(u, 0, D);
    bint_leftShift(v, 0, D);

    BigInt* uHi3 = hi(u, wordLength);
    BigInt* uLo1 = lo(u, wordLength);

    BigInt* R1;
    BigInt* Q1 = bint_divDivAndConq3by2(uHi3, v, &R1, threads);

    bint_leftShift(R1, 0, wordLength);
    bint_add(R1, uLo1);

    BigInt* R2;
    BigInt* Q2 = bint_divDivAndConq3by2(R1, v, &R2, threads);

    bint_leftShift(Q1, 0, wordLength);
    bint_add(Q1, Q2);

    bint_destroy(uHi3);
    bint_destroy(uLo1);
    bint_destroy(R1);
    bint_destroy(Q2);

    bint_destroy(u);
    bint_destroy(v);

    if (rem) {
        bint_rightShift(R2, 0, D);
        *rem = bint_shrink(R2);
    }
    else bint_destroy(R2);
    return bint_shrink(Q1);
}

BigInt* bint_divDivAndConq3by2(BigInt* lhs, BigInt* rhs, BigInt** rem, uint threads) {
    if (lhs->length < rhs->length) {
        if (rem) *rem = bint_clone(lhs);
        return bint_fromWord(0);
    }

    uint wordLength = (rhs->length + 1) / 2;

    BigInt* lhsHi1 = hi(lhs, wordLength*2);
    BigInt* lhsHi2 = hi(lhs, wordLength);
    BigInt* rhsHi  = hi(rhs, wordLength);

    BigInt* lhsLo1 = lo(lhs, wordLength);
    BigInt* rhsLo  = lo(rhs, wordLength);

    BigInt *Q, *R;
    if (bint_lessThan(lhsHi1, rhsHi)) {
        Q = bint_divDivAndConq2by1(lhsHi2, rhsHi, &R, threads);
    }
    else {
        Q = malloc(sizeof(BigInt));
        Q->length = wordLength;
        Q->values = malloc(sizeof(ull) * wordLength);
        memset(Q->values, 0xff, sizeof(ull) * wordLength);

        R = bint_clone(lhsHi2);
        BigInt* rhsHiShifted = bint_leftShift(bint_clone(rhsHi), 0, wordLength);
        bint_sub(R, rhsHiShifted, NULL);
        bint_add(R, rhsHi);

        bint_destroy(rhsHiShifted);
    }

    BigInt* D = bint_mulThreaded(Q, rhsLo, threads);

    bint_leftShift(R, 0, wordLength);
    bint_add(R, lhsLo1);
    bool neg;
    bint_subNoLastCarry(R, D, &neg);

    while (neg) {
        bint_subWord(Q, 1, 0, NULL);
        bint_addNoLastCarry(R, rhs, &neg);
        neg = !neg;
    }

    bint_destroy(lhsHi1);
    bint_destroy(lhsHi2);
    bint_destroy(rhsHi);
    bint_destroy(lhsLo1);
    bint_destroy(rhsLo);
    bint_destroy(D);

    if (rem) *rem = bint_shrink(R);
    else bint_destroy(R);
    return bint_shrink(Q);
}

BigInt* bint_divDivAndConq(BigInt* lhs, BigInt* rhs, BigInt** rem, uint threads) {
    BigInt* u = bint_clone(lhs);
    BigInt* v = bint_clone(rhs);

    uint wordLength = (rhs->length > (lhs->length + 1) / 2) ? rhs->length : (lhs->length + 1) / 2;

    uint D = 0;
    if (wordLength > rhs->length) D = wordLength - rhs->length;

    uint d = 0;
    for (ull top = v->values[v->length - 1]; top < (1L << 63); top <<= 1) ++d;

    bint_leftShift(u, d, D);
    bint_leftShift(v, d, D);

    BigInt* R;
    BigInt* Q = bint_divDivAndConq2by1(u, v, &R, threads);

    bint_destroy(u);
    bint_destroy(v);

    if (rem) {
        bint_rightShift(R, d, D);
        *rem = R;
    }
    else bint_destroy(R);
    return Q;
}

BigInt* bint_rightShift(BigInt* lhs, uint bits, uint words) {
    if (words >= lhs->length) {
        lhs->length = 1;
        lhs->values[0] = 0;
        return lhs;
    }

    lhs->length -= words;

    uint lBits = WORD_LENGTH - bits;
    uint rBits = bits;

    ull mask = (1L << lBits) - 1;
    if (bits == 0) mask = ~mask;

    for (uint i = 0; i < lhs->length; ++i) {
        ull carry;
        carry  =  lhs->values[i + words] << lBits;
        carry &= ~mask;

        ull val;
        val  = lhs->values[i + words] >> rBits;
        val &= mask;

        lhs->values[i] = val;
        if (i != 0) lhs->values[i-1] |= carry;
    }

    return lhs;
}

BigInt* bint_leftShift(BigInt* lhs, uint bits, uint words) {
    lhs->length += words;
    lhs->length += 1;
    lhs->values = realloc(lhs->values, sizeof(ull) * lhs->length);

    uint lBits = bits;
    uint rBits = WORD_LENGTH - bits;

    ull mask = (1L << lBits) - 1;

    lhs->values[lhs->length - 1] = 0;
    for (uint i = lhs->length - 1; i != words;) {
        --i;
        ull carry;
        carry  = lhs->values[i - words] >> rBits;
        carry &= mask;

        ull val;
        val  =  lhs->values[i - words] << lBits;
        val &= ~mask;

        lhs->values[i]      = val;
        lhs->values[i + 1] |= carry;
    }

    memset(lhs->values, 0, sizeof(ull) * words);

    if (lhs->values[lhs->length - 1] == 0) lhs->length -= 1;

    return lhs;
}

BigInt* bint_radMulWord(BigInt* lhs, ull rhsVal, ull rhsExp, ull rad) {
    BigInt carry;
    carry.length = lhs->length + 1;
    carry.values = malloc(sizeof(ull) * carry.length);
    carry.values[0] = 0;

    for (uint i = 0; i < lhs->length; ++i) {
        ull prod, over;
        prod = bigMul(lhs->values[i], rhsVal, &over);
        over = bigDiv(over, prod, rad, &prod, NULL);

        lhs->values[i]      = prod;
        carry.values[i + 1] = over;
    }

    bint_radAdd(lhs, &carry, rad);
    bint_shrink(lhs);
    bint_leftShift(lhs, 0, rhsExp);

    free(carry.values);
    return lhs;
}

BigInt* bint_radAddWord(BigInt* lhs, ull rhsVal, ull rhsExp, ull rad) {
    if (lhs->length < rhsExp + 1) {
        lhs->values = realloc(lhs->values, sizeof(ull) * (rhsExp + 1));
        memset(lhs->values + lhs->length, 0 ,sizeof(ull) * (rhsExp + 1 - lhs->length));
        lhs->length = rhsExp + 1;
        lhs->values[rhsExp] = rhsVal;

        return lhs;
    }

    ull sum   = lhs->values[rhsExp] + rhsVal;
    bool over = sum < rhsVal;
    bigDiv(over, sum, rad, &sum, NULL);

    lhs->values[rhsExp] = sum;
    bool carry = false;
    if (lhs->values[rhsExp] < rhsVal) {
        carry = true;
        for (uint i = rhsExp + 1; i < lhs->length; ++i) {
            lhs->values[i] += 1;
            if (lhs->values[i] >= rad){
                lhs->values[i] = 0;
            }
            else {
                carry = false;
                break;
            }
        }
    }

    if (carry) {
        lhs->length += 1;
        lhs->values = realloc(lhs->values, sizeof(ull) * lhs->length);
        lhs->values[lhs->length - 1] = 1;
    }

    return bint_shrink(lhs);
}

BigInt* bint_radAdd(BigInt* lhs, BigInt* rhs, ull rad) {
    ull sumLength = (lhs->length > rhs->length) ? lhs->length : rhs->length;

    if (rhs->length < sumLength) {
        rhs->values = realloc(rhs->values, sizeof(ull) * sumLength);
        memset(rhs->values + rhs->length, 0, (sumLength - rhs->length) * sizeof(ull));
    }
    else if (lhs->length < sumLength) {
        lhs->values = realloc(lhs->values, sizeof(ull) * sumLength);
        memset(lhs->values + lhs->length, 0, (sumLength - lhs->length) * sizeof(ull));
        lhs->length = sumLength;
    }

    bool carry = false;
    for (uint i = 0; i < sumLength; ++i) {
        ull sum = lhs->values[i] + rhs->values[i] + carry;
        bool overflow = lhs->values[i] != 0 &&
                        rhs->values[i] != 0 &&
                        sum <= lhs->values[i] &&
                        sum <= rhs->values[i];

        carry = overflow || sum >= rad;

        bigDiv(overflow, sum, rad, &sum, NULL);
        lhs->values[i] = sum;
    }

    if (carry) {
        lhs->length = sumLength + 1;
        lhs->values = realloc(lhs->values, sizeof(ull) * lhs->length);
        lhs->values[lhs->length - 1] = 1;
    }

    return bint_shrink(lhs);
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

        BigInt* prod = bint_mulThreaded(ret, partial, threads);
        bint_destroy(ret);
        bint_destroy(partial);

        ret = prod;
    }

    bint_print(ret, threads);

    bint_destroy(ret);
    free(threadPool);

    return 0;
}
