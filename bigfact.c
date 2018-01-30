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

typedef unsigned long long bint_word_t;
typedef unsigned int bint_exp_t;

typedef unsigned int uint;

const bint_exp_t  WORD_LENGTH = 64;
const bint_word_t WORD_MAX    = ULLONG_MAX;
const bint_exp_t  HALF_LENGTH = 32;
const bint_word_t HALF_SIZE   = 1L << 32;

const bint_exp_t  DECIMAL_LENGTH = 19;
const bint_word_t DECIMAL_SIZE   = 1e19;

typedef struct bint_struct {
    bint_word_t* values;
    bint_exp_t   length;
} bint_struct;

typedef bint_struct* bint_t;

bint_t bint_fromWord(bint_word_t value);
bint_t bint_clone(bint_t num);
bint_t bint_shrink(bint_t num);
void bint_destroy(bint_t num);
void bint_print(bint_t num, uint threads);

bint_t bint_structoDecClassical(bint_t num);
bint_t bint_structoDecDivAndConq(bint_t num, uint threads);

bool bint_lessThan(bint_t lhs, bint_t rhs);

bint_t bint_addWord(bint_t lhs, bint_word_t rhsVal, bint_exp_t rhsExp);            /* inplace */
bint_t bint_subWord(bint_t lhs, bint_word_t rhsVal, bint_exp_t rhsExp, int* neg);  /* inplace */
bint_t bint_mulWord(bint_t lhs, bint_word_t rhsVal, bint_exp_t rhsExp);            /* inplace */
bint_t bint_divWord(bint_t lhs, bint_word_t rhsVal, bint_word_t* rem);             /* inplace */
bint_word_t     bint_modWord(bint_t lhs, bint_word_t rhsVal);

bint_t bint_add(bint_t lhs, bint_t rhs);           /* inplace */
bint_t bint_sub(bint_t lhs, bint_t rhs, int* neg); /* inplace */
bint_t bint_mul(bint_t lhs, bint_t rhs);
bint_t bint_div(bint_t lhs, bint_t rhs, bint_t* rem);

bint_t bint_mulThreaded(bint_t lhs, bint_t rhs, uint threads);
bint_t bint_divThreaded(bint_t lhs, bint_t rhs, bint_t* rem, uint threads);

bint_t bint_addNoLastCarry(bint_t lhs, bint_t rhs, bool* carry); /* inplace */
bint_t bint_subNoLastCarry(bint_t lhs, bint_t rhs, bool* carry); /* inplace */

bint_t bint_mulClassical(bint_t lhs, bint_t rhs);
bint_t bint_mulKaratsuba(bint_t lhs, bint_t rhs);
bint_t bint_divClassical(bint_t lhs, bint_t rhs, bint_t* rem);
bint_t bint_divDivAndConq(bint_t lhs, bint_t rhs, bint_t* rem, uint threads);

bint_t bint_leftShift( bint_t lhs, uint bits, bint_exp_t words); /* inplace */
bint_t bint_rightShift(bint_t lhs, uint bits, bint_exp_t words); /* inplace */

bint_t bint_radMulWord(bint_t lhs, bint_word_t rhsVal, bint_word_t rhsExp, bint_word_t rad); /* inplace */
bint_t bint_radAddWord(bint_t lhs, bint_word_t rhsVal, bint_word_t rhsExp, bint_word_t rad); /* inplace */

bint_t bint_radAdd(bint_t lhs, bint_t rhs, bint_word_t rad); /* inplace */

bint_t hi(bint_t num, bint_exp_t cut) {
    if (cut >= num->length) {
        return bint_fromWord(0);
    }

    bint_t ret = malloc(sizeof(bint_struct));
    bint_exp_t testSize = num->length - cut;
    ret->length = (num->length < testSize) ? num->length : testSize;
    ret->values = malloc(sizeof(bint_word_t) * ret->length);
    memcpy(ret->values, num->values + cut, sizeof(bint_word_t) * ret->length);

    return ret;
}

bint_t lo(bint_t num, bint_exp_t cut) {
    if (cut >= num->length) {
        return bint_clone(num);
    }
    bint_t ret = malloc(sizeof(bint_struct));
    ret->length = (num->length < cut) ? num->length : cut;
    ret->values = malloc(sizeof(bint_word_t) * ret->length);
    memcpy(ret->values, num->values, sizeof(bint_word_t) * ret->length);

    return ret;
}

bint_word_t bigMul(bint_word_t lhs, bint_word_t rhs, bint_word_t* _over) {
    bint_word_t ret, over;
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

bint_word_t bigDiv(bint_word_t lhsHi, bint_word_t lhsLo, bint_word_t rhs, bint_word_t* _rem, int* overflow) {
    if (rhs == 0) {
        return lhsLo/rhs;
    }
    if (rhs <= lhsHi) {
        if (overflow) *overflow = true;
        if (_rem) {
            bint_word_t radix, resHi, resLo;
            bigDiv(1, 0, rhs, &radix, NULL);

            resLo = bigMul(radix, lhsHi, &resHi);
            resLo += lhsLo;
            if (resLo < lhsLo) ++resHi;

            bigDiv(resHi, resLo, rhs, _rem, NULL);
        }
        return 0;
    }
    bint_word_t ret, rem;
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

bint_word_t div3by2(bint_t u, bint_t v, bint_t* _rem) {
    bint_word_t rem;
    bool over;
    bint_word_t q = bigDiv(u->values[2], u->values[1], v->values[1], &rem, &over);

    bint_t R;
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

    bint_t D = bint_mulWord(bint_fromWord(q), v->values[0], 0);
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

bint_t bint_fromWord(bint_word_t value) {
    bint_t ret = malloc(sizeof(bint_struct));
    ret->values = malloc(sizeof(bint_word_t));
    ret->length = 1;

    ret->values[0] = value;

    return ret;
}

bint_t bint_clone(bint_t num) {
    bint_t ret = malloc(sizeof(bint_struct));
    ret->length = num->length;
    ret->values = malloc(sizeof(bint_word_t) * ret->length);
    memcpy(ret->values, num->values, sizeof(bint_word_t) * ret->length);

    return ret;
}

bint_t bint_shrink(bint_t num) {
    while (num->length > 1 && num->values[num->length - 1] == 0) num->length -= 1;
    return num;
}

void bint_destroy(bint_t num) {
    free(num->values);
    free(num);
}

void bint_print(bint_t num, uint threads) {
    bint_t dec = bint_structoDecDivAndConq(num, threads);

    printf("%llu", dec->values[dec->length - 1]);;
    for (bint_exp_t i = dec->length - 1; i != 0;) {
        --i;
        printf("%019llu", dec->values[i]);
    }
    printf("\n");

    bint_destroy(dec);
}


bint_t bint_structoDecClassical(bint_t num) {
    const double binToDecLengthFactor = (double) WORD_LENGTH * M_LN2 / M_LN10 / DECIMAL_LENGTH;

    bint_t ret = malloc(sizeof(bint_struct));
    ret->length = (double) num->length * binToDecLengthFactor + 1;
    ret->values = malloc(sizeof(bint_word_t) * ret->length);

    bint_t dividend = bint_clone(num);
    for (bint_exp_t i = 0; i < ret->length; ++i) {
        bint_word_t rem;
        bint_divWord(dividend, DECIMAL_SIZE, &rem);
        ret->values[i] = rem;
    }
    bint_destroy(dividend);

    return bint_shrink(ret);
}

typedef struct _info_bint_structoDecDivAndConq {
    bint_t          num;
    pthread_mutex_t* poolLock;
    pthread_t*       threadPool;
    uint*            currThreads;
    uint             maxThreads;
} _info_bint_structoDecDivAndConq;

void* _thread_bint_structoDecDivAndConq(void* _info) {
    _info_bint_structoDecDivAndConq* info = (_info_bint_structoDecDivAndConq*) _info;

    bint_t          num         = info->num;
    pthread_mutex_t* poolLock    = info->poolLock;
    uint*            currThreads = info->currThreads;
    uint             maxThreads  = info->maxThreads;

    free(_info);

    const double binToDecLengthFactor = (double) WORD_LENGTH * M_LN2 / M_LN10 / DECIMAL_LENGTH;
    const int CUTOFF = 10;
    if (num->length < CUTOFF) {
        return (void*) bint_structoDecClassical(num);
    }

    bint_t ret = malloc(sizeof(bint_struct));
    ret->length = (double) num->length * binToDecLengthFactor + 1;
    ret->values = malloc(sizeof(bint_word_t) * ret->length);
    memset(ret->values, 0, sizeof(bint_word_t) * ret->length);

    uint wordLength = ret->length / 2;
    bint_t divisor = bint_fromWord(1e19);
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

    bint_t numLo, numHi;
    numHi = bint_divThreaded(num, divisor, &numLo, threadsPerThread);

    bint_destroy(divisor);

    bint_t retHi, retLo;

    _info_bint_structoDecDivAndConq* infoHi = malloc(sizeof(_info_bint_structoDecDivAndConq));
    infoHi->num         = numHi;
    infoHi->poolLock    = poolLock;
    infoHi->currThreads = currThreads;
    infoHi->maxThreads  = maxThreads;

    pthread_t split;
    if (splitThread) pthread_create(&split, NULL, _thread_bint_structoDecDivAndConq, (void*) infoHi);
    else retHi = _thread_bint_structoDecDivAndConq((void*) infoHi);

    _info_bint_structoDecDivAndConq* infoLo = malloc(sizeof(_info_bint_structoDecDivAndConq));
    infoLo->num         = numLo;
    infoLo->poolLock    = poolLock;
    infoLo->currThreads = currThreads;
    infoLo->maxThreads  = maxThreads;

    retLo = _thread_bint_structoDecDivAndConq((void*) infoLo);

    if (splitThread) {
        pthread_join(split, (void**) &retHi);
        pthread_mutex_lock(poolLock);
        *currThreads -= 1;
        pthread_mutex_unlock(poolLock);
    }

    memcpy(ret->values,              retLo->values, sizeof(bint_word_t) * retLo->length);
    memcpy(ret->values + wordLength, retHi->values, sizeof(bint_word_t) * retHi->length);

    bint_destroy(numHi);
    bint_destroy(numLo);
    bint_destroy(retHi);
    bint_destroy(retLo);

    return (void*) bint_shrink(ret);
}

bint_t bint_structoDecDivAndConq(bint_t num, uint threads) {
    pthread_mutex_t poolLock;
    pthread_mutex_init(&poolLock, NULL);

    uint currThreads = 1;

    _info_bint_structoDecDivAndConq* info = malloc(sizeof(_info_bint_structoDecDivAndConq));
    info->num         = num;
    info->poolLock    = &poolLock;
    info->currThreads = &currThreads;
    info->maxThreads  = threads;

    bint_t ret = _thread_bint_structoDecDivAndConq((void*) info);

    pthread_mutex_destroy(&poolLock);

    return ret;
}

bool bint_lessThan(bint_t lhs, bint_t rhs) {
    if (lhs->length != rhs->length) {
        return lhs->length < rhs->length;
    }

    for (bint_exp_t i = lhs->length; i != 0;) {
        --i;
        if (lhs->values[i] < rhs->values[i]) return true;
        if (lhs->values[i] > rhs->values[i]) return false;
    }
    return false;
}

bint_t bint_addWord(bint_t lhs, bint_word_t rhsVal, bint_exp_t rhsExp) {
    if (lhs->length < rhsExp + 1) {
        lhs->values = realloc(lhs->values, sizeof(bint_word_t) * (rhsExp + 1));
        memset(lhs->values + lhs->length, 0, sizeof(bint_word_t) * (rhsExp + 1 - lhs->length));
        lhs->length = rhsExp + 1;
        lhs->values[rhsExp] = rhsVal;

        return lhs;
    }

    lhs->values[rhsExp] += rhsVal;
    bool carry = false;
    if (lhs->values[rhsExp] < rhsVal) {
        carry = true;
        for (bint_exp_t i = rhsExp + 1; i < lhs->length; ++i) {
            lhs->values[i] += 1;
            if (lhs->values[i] != 0) {
                carry = false;
                break;
            }
        }
    }

    if (carry) {
        lhs->length += 1;
        lhs->values = realloc(lhs->values, sizeof(bint_word_t) * lhs->length);
        lhs->values[lhs->length - 1] = 1;
    }

    return bint_shrink(lhs);
}

bint_t bint_subWord(bint_t lhs, bint_word_t rhsVal, bint_exp_t rhsExp, int* _neg) {
    bint_exp_t diffLength = lhs->length;
    if (lhs->length < rhsExp + 1) {
        diffLength = rhsExp + 1;
        lhs->values = realloc(lhs->values, sizeof(bint_word_t) * diffLength);
        memset(lhs->values + lhs->length, 0, (diffLength - lhs->length) * sizeof(bint_word_t));
    }

    bint_word_t oldVal = lhs->values[rhsExp] - rhsVal;
    lhs->values[rhsExp] -= rhsVal;
    bool neg = false;
    if (lhs->values[rhsExp] > oldVal) {
        neg = true;
        for (bint_exp_t i = rhsExp + 1; i < diffLength; ++i) {
            lhs->values[i] -= 1;
            if (lhs->values[i] != WORD_MAX) {
                neg = false;
                break;
            }
        }
    }

    if (neg) {
        lhs->length = diffLength;
        for (bint_exp_t i = 0; i < diffLength; ++i) {
            lhs->values[i] = ~lhs->values[i];
        }
        bint_addWord(lhs, 1, 0);
    }

    if (_neg) *_neg = neg;
    return bint_shrink(lhs);
}

bint_t bint_mulWord(bint_t lhs, bint_word_t rhsVal, bint_exp_t rhsExp) {
    bint_struct carry;
    carry.length = lhs->length + 1;
    carry.values = malloc(sizeof(bint_word_t) * carry.length);
    carry.values[0] = 0;

    for (bint_exp_t i = 0; i < lhs->length; ++i) {
        bint_word_t over;
        lhs->values[i] = bigMul(lhs->values[i], rhsVal, &over);
        carry.values[i+1] = over;
    }

    bint_add(lhs, &carry);
    bint_shrink(lhs);
    bint_leftShift(lhs, 0, rhsExp);

    free(carry.values);
    return lhs;
}

bint_t bint_divWord(bint_t lhs, bint_word_t rhsVal, bint_word_t* _rem) {
    bint_word_t rem = 0;
    for (bint_exp_t i = lhs->length; i != 0;) {
        --i;

        bint_word_t divHi = rem;
        bint_word_t divLo = lhs->values[i];

        lhs->values[i] = bigDiv(divHi, divLo, rhsVal, &rem, NULL);
    }

    if (_rem) *_rem = rem;
    return bint_shrink(lhs);
}

bint_word_t bint_modWord(bint_t lhs, bint_word_t rhsVal) {
    bint_word_t radMod;
    bigDiv(1, 0, rhsVal, &radMod, NULL);
    bint_word_t ret = lhs->values[lhs->length - 1] % rhsVal;
    for (bint_exp_t i = lhs->length - 1; i != 0;) {
        --i;
        bint_word_t resHi, resLo;
        resLo = bigMul(ret, radMod, &resHi);

        resLo += lhs->values[i] % rhsVal;
        if (resLo < lhs->values[i] % rhsVal)
            resHi += 1;

        bigDiv(resHi, resLo, rhsVal, &ret, NULL);
    }

    return ret;
}

bint_t bint_add(bint_t lhs, bint_t rhs) {
    bool carry;
    bint_addNoLastCarry(lhs, rhs, &carry);

    if (carry) {
        lhs->length += 1;
        lhs->values = realloc(lhs->values, sizeof(bint_word_t) * lhs->length);
        lhs->values[lhs->length - 1] = 1;
    }

    return lhs;
}

bint_t bint_sub(bint_t lhs, bint_t rhs, bool* neg) {
    bool carry;
    bint_subNoLastCarry(lhs, rhs, &carry);

    if (carry) {
        for (bint_exp_t i = 0; i < lhs->length; ++i) {
            lhs->values[i] = ~(lhs->values[i]);
        }
        bint_addWord(lhs, 1, 0);
    }

    if (neg) *neg = carry;
    return lhs;
}

bint_t bint_mul(bint_t lhs, bint_t rhs) {
    return bint_mulKaratsuba(lhs, rhs);
}

bint_t bint_div(bint_t lhs, bint_t rhs, bint_t* rem) {
    return bint_divThreaded(lhs, rhs, rem, 1);
}

typedef struct _info_bint_mulKaratsuba {
    bint_t    lhs;
    bint_t    rhs;
    bint_exp_t offset;
    bint_exp_t stride;
} _info_bint_mulKaratsuba;

void* _thread_bint_mulKaratsuba(void* _info) {
    _info_bint_mulKaratsuba* info = (_info_bint_mulKaratsuba*) _info;

    if (info->stride + info->offset > info->rhs->length)
        info->stride = info->rhs->length - info->offset;

    bint_t lhs = info->lhs;
    bint_struct  rhs;
    rhs.values = info->rhs->values + info->offset;
    rhs.length = info->stride;

    bint_exp_t rhslen = info->rhs->length;
    bint_exp_t offset = info->offset;
    bint_exp_t stride = info->stride;

    free(info);

    return (void*) bint_mulKaratsuba(lhs, &rhs);
}

bint_t bint_mulThreaded(bint_t lhs, bint_t rhs, uint threads) {
    if (threads == 1) return bint_mul(lhs, rhs);

    uint lengthPerThread = rhs->length / 4;

    pthread_t* threadPool = malloc(sizeof(pthread_t) * threads);
    for (int i = 0; i < threads; ++i) {
        _info_bint_mulKaratsuba* info = malloc(sizeof(_info_bint_mulKaratsuba));
        info->lhs    = lhs;
        info->rhs    = rhs;
        info->offset = i * lengthPerThread;
        info->stride = lengthPerThread;

        if (i == threads - 1)
            info->stride = rhs->length - info->offset;

        pthread_create(&threadPool[i], NULL, _thread_bint_mulKaratsuba, (void*) info);
    }

    bint_t ret;
    pthread_join(threadPool[0], (void**) &ret);

    for (int i = 1; i < threads; ++i) {
        bint_t part;
        pthread_join(threadPool[i], (void**) &part);

        bint_leftShift(part, 0, lengthPerThread * i);
        bint_add(ret, part);

        bint_destroy(part);
    }

    free(threadPool);
    return ret;
}

bint_t bint_divThreaded(bint_t lhs, bint_t rhs, bint_t* _rem, uint threads) {
    bint_t quot, rem;
    quot = bint_divDivAndConq(lhs, rhs, &rem, threads);
    if (_rem) *_rem = rem;
    return quot;
}

bint_t bint_addNoLastCarry(bint_t lhs, bint_t rhs, bool* _carry) {
    if (lhs->length > rhs->length) {
        rhs->values = realloc(rhs->values, sizeof(bint_word_t) * lhs->length);
        memset(rhs->values + rhs->length, 0, (lhs->length - rhs->length) * sizeof(bint_word_t));
    }
    else if (lhs->length < rhs->length) {
        lhs->values = realloc(lhs->values, sizeof(bint_word_t) * rhs->length);
        memset(lhs->values + lhs->length, 0, (rhs->length - lhs->length) * sizeof(bint_word_t));
        lhs->length = rhs->length;
    }

    bint_word_t carry;
    bint_word_t loops = lhs->length;
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

bint_t bint_subNoLastCarry(bint_t lhs, bint_t rhs, bool* _carry) {
    if (lhs->length > rhs->length) {
        rhs->values = realloc(rhs->values, sizeof(bint_word_t) * lhs->length);
        memset(rhs->values + rhs->length, 0, (lhs->length - rhs->length) * sizeof(bint_word_t));
    }

    bint_word_t carry;
    bint_word_t loops = lhs->length;
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

bint_t bint_mulClassical(bint_t lhs, bint_t rhs) {
    bint_t ret = bint_fromWord(0);
    ret->values = realloc(ret->values, sizeof(bint_word_t) * (lhs->length + lhs->length));

    for (bint_exp_t i = 0; i < lhs->length; ++i) {
        for (bint_exp_t j = 0; j < rhs->length; ++j) {
            bint_word_t prod, over;
            prod = bigMul(lhs->values[i], rhs->values[j], &over);
            bint_addWord(ret, prod, i + j);
            if (over != 0) bint_addWord(ret, over, i + j + 1);
        }
    }

    return bint_shrink(ret);
}

bint_t bint_mulKaratsuba(bint_t lhs, bint_t rhs) {
    const bint_exp_t CUTOFF = 10;
    if (lhs->length < CUTOFF || rhs->length < CUTOFF) {
        return bint_mulClassical(lhs, rhs);
    }
    bint_t ret;

    bint_exp_t lhsHalfLength = lhs->length/2;
    bint_exp_t rhsHalfLength = rhs->length/2;
    bint_exp_t halfLength = (lhsHalfLength > rhsHalfLength) ? lhsHalfLength : rhsHalfLength;

    bint_t lhsLo = lo(lhs, halfLength);
    bint_t lhsHi = hi(lhs, halfLength);
    bint_t rhsLo = lo(rhs, halfLength);
    bint_t rhsHi = hi(rhs, halfLength);

    bint_t p0, p1, p2;

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
bint_t bint_divClassical(bint_t lhs, bint_t rhs, bint_t* rem) {
    if (lhs->length < rhs->length) {
        if (rem) *rem = bint_clone(lhs);
        return bint_fromWord(0);
    }

    bint_t u = bint_clone(lhs);
    bint_t v = bint_clone(rhs);

    bint_exp_t n = rhs->length;
    bint_exp_t m = lhs->length - rhs->length;

    bint_t q = malloc(sizeof(bint_struct));
    q->length = m + 1;
    q->values = malloc(sizeof(bint_word_t) * q->length);

    // Normalize (D1)
    uint d = 0;
    for (bint_word_t top = v->values[n-1]; top < (1L << 63); top <<= 1) ++d;

    bint_leftShift(u, d, 0);
    bint_leftShift(v, d, 0);

    v->values = realloc(v->values, sizeof(bint_word_t) * (n+1));
    v->values[n] = 0;

    if (u->length < n + m + 1) {
        u->length = n + m + 1;
        u->values = realloc(u->values, sizeof(bint_word_t) * (n + m + 1));
        u->values[n + m] = 0;
    }

    // Initialize j (D2) / Loop on j (D7)
    for (bint_exp_t j = m + 1; j != 0;) {
        --j;
        // Calculate q_hat (D3)
        bint_struct vHi;
        vHi.length = 2;
        vHi.values = v->values + n-2;

        bint_struct uHi;
        uHi.length = 3;
        uHi.values = u->values + j+n-2;

        bint_word_t q_hat = div3by2(&uHi, &vHi, NULL);

        // Multiply and subtract (D4)
        bint_struct uPart;
        uPart.values = u->values + j;
        uPart.length = n + 1;

        bint_t qv = bint_clone(v);
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

bint_t bint_divDivAndConq2by1(bint_t lhs, bint_t rhs, bint_t* rem, uint threads);
bint_t bint_divDivAndConq3by2(bint_t lhs, bint_t rhs, bint_t* rem, uint threads);

bint_t bint_divDivAndConq2by1(bint_t lhs, bint_t rhs, bint_t* rem, uint threads) {
    const bint_exp_t CUTOFF = 10;
    if (rhs->length < CUTOFF) {
        return bint_divClassical(lhs, rhs, rem);
    }

    if (lhs->length < rhs->length) {
        if (rem) *rem = bint_clone(lhs);
        return bint_fromWord(0);
    }

    bint_exp_t wordLength = (rhs->length + 1) / 2;

    bint_t u = bint_clone(lhs);
    bint_t v = bint_clone(rhs);

    // Renormalize
    bint_exp_t D = 0;
    if (wordLength * 2 > v->length) D = 1;

    bint_leftShift(u, 0, D);
    bint_leftShift(v, 0, D);

    bint_t uHi3 = hi(u, wordLength);
    bint_t uLo1 = lo(u, wordLength);

    bint_t R1;
    bint_t Q1 = bint_divDivAndConq3by2(uHi3, v, &R1, threads);

    bint_leftShift(R1, 0, wordLength);
    bint_add(R1, uLo1);

    bint_t R2;
    bint_t Q2 = bint_divDivAndConq3by2(R1, v, &R2, threads);

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

bint_t bint_divDivAndConq3by2(bint_t lhs, bint_t rhs, bint_t* rem, uint threads) {
    if (lhs->length < rhs->length) {
        if (rem) *rem = bint_clone(lhs);
        return bint_fromWord(0);
    }

    bint_exp_t wordLength = (rhs->length + 1) / 2;

    bint_t lhsHi1 = hi(lhs, wordLength*2);
    bint_t lhsHi2 = hi(lhs, wordLength);
    bint_t rhsHi  = hi(rhs, wordLength);

    bint_t lhsLo1 = lo(lhs, wordLength);
    bint_t rhsLo  = lo(rhs, wordLength);

    bint_t Q, R;
    if (bint_lessThan(lhsHi1, rhsHi)) {
        Q = bint_divDivAndConq2by1(lhsHi2, rhsHi, &R, threads);
    }
    else {
        Q = malloc(sizeof(bint_struct));
        Q->length = wordLength;
        Q->values = malloc(sizeof(bint_word_t) * wordLength);
        memset(Q->values, 0xff, sizeof(bint_word_t) * wordLength);

        R = bint_clone(lhsHi2);
        bint_t rhsHiShifted = bint_leftShift(bint_clone(rhsHi), 0, wordLength);
        bint_sub(R, rhsHiShifted, NULL);
        bint_add(R, rhsHi);

        bint_destroy(rhsHiShifted);
    }

    bint_t D = bint_mulThreaded(Q, rhsLo, threads);

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

bint_t bint_divDivAndConq(bint_t lhs, bint_t rhs, bint_t* rem, uint threads) {
    bint_t u = bint_clone(lhs);
    bint_t v = bint_clone(rhs);

    bint_exp_t wordLength = (rhs->length > (lhs->length + 1) / 2) ? rhs->length : (lhs->length + 1) / 2;

    bint_exp_t D = 0;
    if (wordLength > rhs->length) D = wordLength - rhs->length;

    uint d = 0;
    for (bint_word_t top = v->values[v->length - 1]; top < (1L << 63); top <<= 1) ++d;

    bint_leftShift(u, d, D);
    bint_leftShift(v, d, D);

    bint_t R;
    bint_t Q = bint_divDivAndConq2by1(u, v, &R, threads);

    bint_destroy(u);
    bint_destroy(v);

    if (rem) {
        bint_rightShift(R, d, D);
        *rem = R;
    }
    else bint_destroy(R);
    return Q;
}

bint_t bint_rightShift(bint_t lhs, uint bits, bint_exp_t words) {
    if (words >= lhs->length) {
        lhs->length = 1;
        lhs->values[0] = 0;
        return lhs;
    }

    lhs->length -= words;

    uint lBits = WORD_LENGTH - bits;
    uint rBits = bits;

    bint_word_t mask = (1L << lBits) - 1;
    if (bits == 0) mask = ~mask;

    for (bint_exp_t i = 0; i < lhs->length; ++i) {
        bint_word_t carry;
        carry  =  lhs->values[i + words] << lBits;
        carry &= ~mask;

        bint_word_t val;
        val  = lhs->values[i + words] >> rBits;
        val &= mask;

        lhs->values[i] = val;
        if (i != 0) lhs->values[i-1] |= carry;
    }

    return lhs;
}

bint_t bint_leftShift(bint_t lhs, uint bits, bint_exp_t words) {
    lhs->length += words;
    lhs->length += 1;
    lhs->values = realloc(lhs->values, sizeof(bint_word_t) * lhs->length);

    uint lBits = bits;
    uint rBits = WORD_LENGTH - bits;

    bint_word_t mask = (1L << lBits) - 1;

    lhs->values[lhs->length - 1] = 0;
    for (bint_exp_t i = lhs->length - 1; i != words;) {
        --i;
        bint_word_t carry;
        carry  = lhs->values[i - words] >> rBits;
        carry &= mask;

        bint_word_t val;
        val  =  lhs->values[i - words] << lBits;
        val &= ~mask;

        lhs->values[i]      = val;
        lhs->values[i + 1] |= carry;
    }

    memset(lhs->values, 0, sizeof(bint_word_t) * words);

    if (lhs->values[lhs->length - 1] == 0) lhs->length -= 1;

    return lhs;
}

bint_t bint_radMulWord(bint_t lhs, bint_word_t rhsVal, bint_word_t rhsExp, bint_word_t rad) {
    bint_struct carry;
    carry.length = lhs->length + 1;
    carry.values = malloc(sizeof(bint_word_t) * carry.length);
    carry.values[0] = 0;

    for (bint_exp_t i = 0; i < lhs->length; ++i) {
        bint_word_t prod, over;
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

bint_t bint_radAddWord(bint_t lhs, bint_word_t rhsVal, bint_word_t rhsExp, bint_word_t rad) {
    if (lhs->length < rhsExp + 1) {
        lhs->values = realloc(lhs->values, sizeof(bint_word_t) * (rhsExp + 1));
        memset(lhs->values + lhs->length, 0 ,sizeof(bint_word_t) * (rhsExp + 1 - lhs->length));
        lhs->length = rhsExp + 1;
        lhs->values[rhsExp] = rhsVal;

        return lhs;
    }

    bint_word_t sum   = lhs->values[rhsExp] + rhsVal;
    bool over = sum < rhsVal;
    bigDiv(over, sum, rad, &sum, NULL);

    lhs->values[rhsExp] = sum;
    bool carry = false;
    if (lhs->values[rhsExp] < rhsVal) {
        carry = true;
        for (bint_exp_t i = rhsExp + 1; i < lhs->length; ++i) {
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
        lhs->values = realloc(lhs->values, sizeof(bint_word_t) * lhs->length);
        lhs->values[lhs->length - 1] = 1;
    }

    return bint_shrink(lhs);
}

bint_t bint_radAdd(bint_t lhs, bint_t rhs, bint_word_t rad) {
    bint_word_t sumLength = (lhs->length > rhs->length) ? lhs->length : rhs->length;

    if (rhs->length < sumLength) {
        rhs->values = realloc(rhs->values, sizeof(bint_word_t) * sumLength);
        memset(rhs->values + rhs->length, 0, (sumLength - rhs->length) * sizeof(bint_word_t));
    }
    else if (lhs->length < sumLength) {
        lhs->values = realloc(lhs->values, sizeof(bint_word_t) * sumLength);
        memset(lhs->values + lhs->length, 0, (sumLength - lhs->length) * sizeof(bint_word_t));
        lhs->length = sumLength;
    }

    bool carry = false;
    for (bint_exp_t i = 0; i < sumLength; ++i) {
        bint_word_t sum = lhs->values[i] + rhs->values[i] + carry;
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
        lhs->values = realloc(lhs->values, sizeof(bint_word_t) * lhs->length);
        lhs->values[lhs->length - 1] = 1;
    }

    return bint_shrink(lhs);
}

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

    bint_print(ret, threads);

    bint_destroy(ret);
    free(threadPool);

    return 0;
}
