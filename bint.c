#include <stdlib.h>
#include <stdio.h>
#include <math.h>   // M_LN10, M_LN2
#include <string.h> // memcpy, memset
#include <pthread.h>
#include "bint.h"

#define bool int
#define false 0
#define true !(false)

#define LNDEC_BIN (M_LN2 * WORD_LENGTH / M_LN10 / DECIMAL_LENGTH)
#define LNBIN_DEC (1/LNDEC_BIN)

const bint_exp_t  WORD_LENGTH = 64;
const bint_word_t WORD_MAX    = ULLONG_MAX;
const bint_exp_t  HALF_LENGTH = 32;
const bint_word_t HALF_SIZE   = 1L << 32;

const bint_exp_t  DECIMAL_LENGTH = 19;
const bint_word_t DECIMAL_SIZE   = 1e19;

static bint_t hi(bint_t num, bint_exp_t cut) {
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

static bint_t lo(bint_t num, bint_exp_t cut) {
    if (cut >= num->length) {
        return bint_clone(num);
    }
    bint_t ret = malloc(sizeof(bint_struct));
    ret->length = (num->length < cut) ? num->length : cut;
    ret->values = malloc(sizeof(bint_word_t) * ret->length);
    memcpy(ret->values, num->values, sizeof(bint_word_t) * ret->length);

    return ret;
}

// The x86 MUL instruction returns a 2 word product. Only the low part is
// available (easily) in C. This function returns both words
static bint_word_t bigMul(bint_word_t lhs, bint_word_t rhs, bint_word_t* _over) {
    // [over, ret] = lhs * rhs
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

// The x86 DIV instruction can accept a 2 word dividend. Thisfeature isn't used
// by C, probably because if the quotient overflows (which is possible) it will
// trigger a #DE (Divide Error) exception so you have to be pretty careful
static bint_word_t bigDiv(bint_word_t lhsHi, bint_word_t lhsLo, bint_word_t rhs, bint_word_t* _rem, int* overflow) {
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
    // ret = [lhsHi, lhsLo] / rhs
    // rem = [lhsHi, lhsLo] % rhs
    // (Assuming the quotient doesn't overflow)
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

static bint_word_t div3by2(bint_t u, bint_t v, bint_t* _rem) {
    bint_word_t rem;
    bool over;
    bint_word_t q = bigDiv(u->values[2], u->values[1], v->values[1], &rem, &over);

    bint_t R;
    if (over) {
        q = WORD_MAX;
        // R = [u2, u1] - [v1, 0] + v1
        R = bint_leftShift(bint_fromWord(u->values[2]), 0, 1);
        R->values[0] = u->values[1];

        bint_subWord(R, v->values[1], 1, NULL);
        bint_addWord(R, v->values[1], 0);
    }
    else {
        R = bint_fromWord(rem);
    }

    // D = q * v0
    bint_t D = bint_mulWord(bint_fromWord(q), v->values[0], 0);

    // R = [R, u0] - D
    bint_leftShift(R, 0, 1);
    R->values[0] = u->values[0];
    // Intentionally don't do a borrow, It will be replaced by a carry later
    bool neg;
    bint_subNoLastCarry(R, D, &neg);

    // This loop should run at most twice (assuming normalized divisor)
    while (neg) {
        --q;
        // We didn't do a borrow before, not doing a carry here cancels that out
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

// TODO // Get an exact reference
// This is a pretty standard way to calculate powers, don't have a copy on me
// for an exact reference, but SICP discusses this. It is much faster than the
// naive method (repeated multiplication)
bint_t bint_powerOfTen(bint_exp_t exp) {
    if (exp == 0) return bint_fromWord(0);
    if (exp == 1) return bint_fromWord(1e19);

    bint_t prev = bint_powerOfTen(exp/2);
    bint_t ret  = bint_mul(prev, prev);
    if (exp % 2 == 1) bint_mulWord(ret, 1e19, 0);

    bint_destroy(prev);
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
    // Recalling that bint_t = bint_struct*
    free(num->values);
    free(num);
}

void bint_print(bint_t num) {
    bint_printThreaded(num, 1);
}

void bint_printThreaded(bint_t num, uint threads) {
    bint_t dec = bint_toDecDivAndConq(num, threads);

    printf("%llu", dec->values[dec->length - 1]);;
    for (bint_exp_t i = dec->length - 1; i != 0;) {
        --i;
        printf("%019llu", dec->values[i]);
    }
    printf("\n");

    bint_destroy(dec);
}


bint_t bint_toDecClassical(bint_t num) {
    bint_t ret = malloc(sizeof(bint_struct));
    ret->length = (double) num->length * LNDEC_BIN + 1;
    ret->values = malloc(sizeof(bint_word_t) * ret->length);

    // Repeatedly divide by 1e19 storing the remainder into the return values
    bint_t dividend = bint_clone(num);
    for (bint_exp_t i = 0; i < ret->length; ++i) {
        bint_word_t rem;
        bint_divWord(dividend, DECIMAL_SIZE, &rem);
        ret->values[i] = rem;
    }
    bint_destroy(dividend);

    return bint_shrink(ret);
}

typedef struct _info_bint_toDecDivAndConq {
    bint_t          num;
    pthread_mutex_t* poolLock;
    pthread_t*       threadPool;
    uint*            currThreads;
    uint             maxThreads;
} _info_bint_toDecDivAndConq;

static void* _thread_bint_toDecDivAndConq(void* _info) {
    _info_bint_toDecDivAndConq* info = (_info_bint_toDecDivAndConq*) _info;

    bint_t          num         = info->num;
    pthread_mutex_t* poolLock    = info->poolLock;
    uint*            currThreads = info->currThreads;
    uint             maxThreads  = info->maxThreads;

    free(_info);

    const int CUTOFF = 10;
    if (num->length < CUTOFF) {
        return (void*) bint_toDecClassical(num);
    }

    bint_t ret = malloc(sizeof(bint_struct));
    ret->length = (double) num->length * LNDEC_BIN + 1;
    ret->values = malloc(sizeof(bint_word_t) * ret->length);
    memset(ret->values, 0, sizeof(bint_word_t) * ret->length);

    // This divisor is a power of 10^19 that is roughly the square root of num
    uint wordLength = ret->length / 2;
    bint_t divisor = bint_powerOfTen(wordLength);

    // How many threads can we use to calculate the division
    pthread_mutex_lock(poolLock);
    uint threadsPerThread = maxThreads / *currThreads;
    pthread_mutex_unlock(poolLock);

    // This division is the crux, it splits num (roughly) into 2 equal pieces
    // by dividing by (1e19)^n for some large n. This gives 2 "digits" of a
    // number in base (1e19)^n. We can simply print these in order, to do so we
    // repeat the process on each of these "digits".
    bint_t numLo, numHi;
    numHi = bint_divThreaded(num, divisor, &numLo, threadsPerThread);
    bint_destroy(divisor);

    bint_t retHi, retLo;

    // If there are still more threads we can use, launch this function in a
    // new thread, otherwise, run it in this thread
    bool splitThread = false;
    pthread_mutex_lock(poolLock);
    if (*currThreads < maxThreads) {
        splitThread = true;
        *currThreads += 1;
    }
    pthread_mutex_unlock(poolLock);

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

    // If we ran the first recursive call in a new thread, we need to join the
    // thread back (obviously)
    if (splitThread) {
        pthread_join(split, (void**) &retHi);
        pthread_mutex_lock(poolLock);
        *currThreads -= 1;
        pthread_mutex_unlock(poolLock);
    }

    // TODO //
    // Copy the values into our ret here, It would be a fairly simple speed up
    // to just write them into the correct spot in the first place
    memcpy(ret->values,              retLo->values, sizeof(bint_word_t) * retLo->length);
    memcpy(ret->values + wordLength, retHi->values, sizeof(bint_word_t) * retHi->length);

    bint_destroy(numHi);
    bint_destroy(numLo);
    bint_destroy(retHi);
    bint_destroy(retLo);

    return (void*) bint_shrink(ret);
}

bint_t bint_toDecDivAndConq(bint_t num, uint threads) {
    // This just calls _thread_bint_toDecDivAndConq
    pthread_mutex_t poolLock;
    pthread_mutex_init(&poolLock, NULL);

    uint currThreads = 1;

    _info_bint_toDecDivAndConq* info = malloc(sizeof(_info_bint_toDecDivAndConq));
    info->num         = num;
    info->poolLock    = &poolLock;
    info->currThreads = &currThreads;
    info->maxThreads  = threads;

    bint_t ret = _thread_bint_toDecDivAndConq((void*) info);

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

bool bint_greaterThan(bint_t lhs, bint_t rhs) {
    return bint_lessThan(rhs, lhs);
}

bool bint_equals(bint_t lhs, bint_t rhs) {
    return !bint_lessThan(lhs, rhs) && !bint_lessThan(rhs, lhs);
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
        // If the sum is less than either summand, overflow has occurd
        carry = true;
        for (bint_exp_t i = rhsExp + 1; i < lhs->length; ++i) {
            lhs->values[i] += 1;
            if (lhs->values[i] != 0) {
                // If we add 1 and get 0, we have again overflowed
                carry = false;
                break;
            }
        }
    }

    if (carry) {
        // If carry is still true we ran out of space to flow the carry up.
        lhs->length += 1;
        lhs->values = realloc(lhs->values, sizeof(bint_word_t) * lhs->length);
        lhs->values[lhs->length - 1] = 1;
    }

    return bint_shrink(lhs);
}

// See bint_addWord above, this is much the same
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
    // It is *substantially* faster to keep track of all the high parts of the
    // multiplications inside this carry variable and add them all at once
    // rather than adding them as they come up
    bint_struct carry;
    carry.length = lhs->length + 1;
    carry.values = malloc(sizeof(bint_word_t) * carry.length);
    carry.values[0] = 0;

    // For all i, set [carry[i+1], lhs[i]] = lhs[i] * rhsVal
    for (bint_exp_t i = 0; i < lhs->length; ++i) {
        bint_word_t over;
        lhs->values[i] = bigMul(lhs->values[i], rhsVal, &over);
        carry.values[i+1] = over;
    }

    bint_add(lhs, &carry);

    // If we shrink first, leftShift has to do (slightly) less work
    bint_shrink(lhs);
    bint_leftShift(lhs, 0, rhsExp);

    free(carry.values);
    return lhs;
}

// This is classical long division. See: Elementary School
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

// This function is very, very slightly faster than divWord so I'm going to
// leave it here
//
// It relies on the fact that (a b)%m = ((a%m) (b%m))%m and similarly,
// (a + b)%m = ((a%m) + (b&m))%m. This allows you to do arithmetic within the
// bounds of a word. Here we essentially just evaluate the polynomial.
bint_word_t bint_modWord(bint_t lhs, bint_word_t rhsVal) {
    bint_word_t radMod;
    bigDiv(1, 0, rhsVal, &radMod, NULL);
    bint_word_t ret = lhs->values[lhs->length - 1] % rhsVal;
    for (bint_exp_t i = lhs->length - 1; i != 0;) {
        --i;
        bint_word_t resHi, resLo;
        // [resHi, resLo] = ret * (WORD_MAX % rhsVal) + lhs[i]
        // Note that this can't overflow as the worst case scenario doesnt:
        // (WORD_MAX - 1) * (WORD_MAX - 1) + (WORD_MAX - 1) < WORD_MAX^2
        resLo = bigMul(ret, radMod, &resHi);

        resLo += lhs->values[i] % rhsVal;
        if (resLo < lhs->values[i] % rhsVal)
            resHi += 1;

        bigDiv(resHi, resLo, rhsVal, &ret, NULL);
    }

    return ret;
}

bint_t bint_leftShift(bint_t lhs, uint bits, bint_exp_t words) {
    lhs->length += words;
    lhs->length += 1;
    lhs->values = realloc(lhs->values, sizeof(bint_word_t) * lhs->length);

    uint lBits = bits;
    uint rBits = WORD_LENGTH - bits;

    // This contains rBits 0's followed by lBits 1's
    // that might make you think those names are backwards but I'm confused
    // enough as it is OK
    bint_word_t mask = (1L << lBits) - 1;

    // You have to be *very* careful how you do this
    lhs->values[lhs->length - 1] = 0;
    for (bint_exp_t i = lhs->length - 1; i != words;) {
        --i;
        // The carry is the bits that were shifted out beyond the word, we can
        // calculate this with a *right* shift
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

    // If bits is 0 then lBits is 64 = 0 % 64, so the mask is then all zeros.
    // That is the exact oposite of what we want
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

bint_t bint_add(bint_t lhs, bint_t rhs) {
    bool carry;
    bint_addNoLastCarry(lhs, rhs, &carry);

    // If we have a carry, just allocated space for a new word and set it to 1
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

    // If we have a carry (rightly a borrow) we already have the answer is in
    // twos-complement, flip all the bits and add 1 to change it back
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

// This essentially just calls bint_mulKaratsuba in a pthread compliant manner
static void* _thread_bint_mulKaratsuba(void* _info) {
    _info_bint_mulKaratsuba* info = (_info_bint_mulKaratsuba*) _info;

    if (info->stride + info->offset > info->rhs->length)
        info->stride = info->rhs->length - info->offset;

    bint_t lhs = info->lhs;
    bint_struct  rhs;
    rhs.values = info->rhs->values + info->offset;
    rhs.length = info->stride;

    free(info);

    return (void*) bint_mulKaratsuba(lhs, &rhs);
}

bint_t bint_mulThreaded(bint_t lhs, bint_t rhs, uint threads) {
    if (threads == 1) return bint_mulKaratsuba(lhs, rhs);

    uint lengthPerThread = rhs->length / 4;

    // We split rhs into as many limbs as there are threads and multiply the rhs
    // by each limb in an individual thread.
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

    // When reading back the results, we just need to shift the naswers to
    // correspond to where the limb was taken from rhs
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

// The only multi-threading going on here is actually just calls to bint_mulThreaded
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

    // This asm snippet sets the carry flag to 0 then starting from the least
    // significant bit calculates:
    //      lhs[i] = lhs[i] + rhs[i] + carry.
    // If the calculation carries (overflows) it sets the carry flag to 1
    //
    // Upon completion it outputs the carry flag to the variable carry
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

// The name of this is a slight misnomer, it should properly be called
// bint_subNolastBorrow but I wanted it to be consistent with the add above
bint_t bint_subNoLastCarry(bint_t lhs, bint_t rhs, bool* _carry) {
    if (lhs->length > rhs->length) {
        rhs->values = realloc(rhs->values, sizeof(bint_word_t) * lhs->length);
        memset(rhs->values + rhs->length, 0, (lhs->length - rhs->length) * sizeof(bint_word_t));
    }

    // this asm snippet sets the carry flag to 0 then starting from the least
    // significant word calculates:
    //      lhs[i] = lhs[i] - (rhs[i] + carry).
    // If this calculation borrows (underflows) it sets the carry flag to 1
    //
    // Upon completion it outputs the carry flag to the variable carry
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
    // Realloc once here to make sure it doesn't get resized every time
    // bint_add is called

    for (bint_exp_t i = 0; i < lhs->length; ++i) {
        for (bint_exp_t j = 0; j < rhs->length; ++j) {
            // [ret[i+j+1], ret[i+j]] += lhs[i] * rhs[j]
            bint_word_t prod, over;
            prod = bigMul(lhs->values[i], rhs->values[j], &over);
            // We must ADD it, you can't set it. Multiple pairs map to the same
            // exponent!
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
    // p0 = lhsLo * rhsLo

    // p1 = lhsLo * rhsHi + lhsHi * rhsLo
    //    = (lhsHi + lhsLo) * (rhsHi + rhsLo) - p0 - p2

    // p2 = lhsHi * rhsHi

    // Defined as such we can write
    // lhs * rhs = p2 * RADIX^2 + p1 * RADIX + p0

    p0 = bint_mulKaratsuba(lhsLo, rhsLo);
    p2 = bint_mulKaratsuba(lhsHi, rhsHi);

    ret = bint_clone(p2);
    bint_leftShift(ret, 0, 2*halfLength);
    bint_add(ret, p0);

    // Modifying values so I don't have to allocate more space
    // These variables are really misnamed from here on. They represet the sum
    // of the hi part and the lo part
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

    // Different normalization factor that uses shifts instead of
    // multiplication for speed, the factor here is actually 2^d
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

        // A do a very different estimate for q_hat than what Knuth uses.
        // I divide the top 3 words of the dividend by the top 2 words of the
        // divisor which will produce a quotient that is at most 1 too big
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

// These functions are based on:
//      Christoph Burnikel, and Joachim Ziegler
//      Fast Recursive Division (1998)
//
// I will refer to the notation in the paper in comments

static bint_t bint_divDivAndConq2by1(bint_t lhs, bint_t rhs, bint_t* rem, uint threads);
static bint_t bint_divDivAndConq3by2(bint_t lhs, bint_t rhs, bint_t* rem, uint threads);

// Algorithm 1/ (D_{2n/1n})
static bint_t bint_divDivAndConq2by1(bint_t lhs, bint_t rhs, bint_t* rem, uint threads) {
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
    // This step isnt included in Burnikel and Ziegler, this allows me to
    // further subdivide odd word sizes by renormalizing it up one
    bint_exp_t D = 0;
    if (wordLength * 2 > v->length) D = 1;

    bint_leftShift(u, 0, D);
    bint_leftShift(v, 0, D);

    bint_t uHi3 = hi(u, wordLength); // [A1, A2, A3]
    bint_t uLo1 = lo(u, wordLength); // A4

    // Calculate [A1, A2, A3] / [B1, B2]
    bint_t R1;
    bint_t Q1 = bint_divDivAndConq3by2(uHi3, v, &R1, threads);

    // R1 = [R11, R12, A4]
    bint_leftShift(R1, 0, wordLength);
    bint_add(R1, uLo1);

    // Calculate [R11, R12, A4] / [B1, B2]
    bint_t R2;
    bint_t Q2 = bint_divDivAndConq3by2(R1, v, &R2, threads);

    // Q1 = [Q1, Q2]
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

// Algorithm 2. (D_{3n/2n})
static bint_t bint_divDivAndConq3by2(bint_t lhs, bint_t rhs, bint_t* rem, uint threads) {
    if (lhs->length < rhs->length) {
        if (rem) *rem = bint_clone(lhs);
        return bint_fromWord(0);
    }

    bint_exp_t wordLength = (rhs->length + 1) / 2;

    bint_t lhsHi1 = hi(lhs, wordLength*2); // A1
    bint_t lhsHi2 = hi(lhs, wordLength);   // [A1, A2]
    bint_t rhsHi  = hi(rhs, wordLength);   // B1

    bint_t lhsLo1 = lo(lhs, wordLength);   // A3
    bint_t rhsLo  = lo(rhs, wordLength);   // B2

    bint_t Q, R;
    if (bint_lessThan(lhsHi1, rhsHi)) {
        // Calculate [A1, A2] / B1
        Q = bint_divDivAndConq2by1(lhsHi2, rhsHi, &R, threads);
    }
    else {
        // Set Q to the \beta^n - 1
        Q = malloc(sizeof(bint_struct));
        Q->length = wordLength;
        Q->values = malloc(sizeof(bint_word_t) * wordLength);
        memset(Q->values, 0xff, sizeof(bint_word_t) * wordLength);

        // R = [A1, A2] - [B1, 0] + B1
        R = bint_clone(lhsHi2);
        bint_t rhsHiShifted = bint_leftShift(bint_clone(rhsHi), 0, wordLength);
        bint_sub(R, rhsHiShifted, NULL);
        bint_add(R, rhsHi);

        bint_destroy(rhsHiShifted);
    }

    // D = Q * B2
    bint_t D = bint_mulThreaded(Q, rhsLo, threads);

    // R = R * \beta^n + A3 - D
    bint_leftShift(R, 0, wordLength);
    bint_add(R, lhsLo1);
    bool neg;
    // Don't carry during the subtraction, leave it in twos-complement.
    // This will be rectified by ignoring the carry in the addition later.
    bint_subNoLastCarry(R, D, &neg);

    // This loop should run at most twice
    while (neg) {
        // Q = Q - 1
        // R = R + B
        bint_subWord(Q, 1, 0, NULL);
        bint_addNoLastCarry(R, rhs, &neg);
        neg = !neg;
        // As mentioned before, dont carry as it will cancel with the borrow we
        // didn't do with the subtraction.
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

    // All we are doing here is normalizing so that we have a 2 word dividend
    // and a normalized 1 word divisior (normalized as in its top bit is set)
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

// Apart from carefully enforcing the radix, this is the same as bint_wordAdd
bint_t bint_radAddWord(bint_t lhs, bint_word_t rhsVal, bint_word_t rhsExp, bint_word_t rad) {
    if (lhs->length < rhsExp + 1) {
        lhs->values = realloc(lhs->values, sizeof(bint_word_t) * (rhsExp + 1));
        memset(lhs->values + lhs->length, 0 ,sizeof(bint_word_t) * (rhsExp + 1 - lhs->length));
        lhs->length = rhsExp + 1;
        lhs->values[rhsExp] = rhsVal;

        return lhs;
    }

    // Note that doing a simple 1 word modulus is insufficent
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

// Apart from carefully enforcing the radix, this is the same as bint_mulWord
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

// Apart from carefully enforcing the radix, this is the same as bint_add
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

        // Note that doing a simple 1 word modulus is insufficent
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
