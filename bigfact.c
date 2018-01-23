#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
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
void bint_print(BigInt* num);

BigInt* bint_toDecClassical(BigInt* num);
BigInt* bint_toDecDivAndConq(BigInt* num);

bool bint_lessThan(BigInt* lhs, BigInt* rhs);

BigInt* bint_addWord(BigInt* lhs, ull rhsVal, uint rhsExp);            /* inplace */
BigInt* bint_subWord(BigInt* lhs, ull rhsVal, uint rhsExp, int* neg);  /* inplace */
BigInt* bint_mulWord(BigInt* lhs, ull rhsVal, uint rhsExp);            /* inplace */
BigInt* bint_divWord(BigInt* lhs, ull rhsVal, ull* rem);               /* inplace */
ull     bint_modWord(BigInt* lhs, ull rhsVal);

BigInt* bint_add(BigInt* lhs, BigInt* rhs);                            /* inplace */
BigInt* bint_sub(BigInt* lhs, BigInt* rhs, int* neg);                  /* inplace */
BigInt* bint_mul(BigInt* lhs, BigInt* rhs, uint threads);
BigInt* bint_div(BigInt* lhs, BigInt* rhs, BigInt** rem);

BigInt* bint_addNoLastCarry(BigInt* lhs, BigInt* rhs, bool* carry);    /* inplace */
BigInt* bint_subNoLastCarry(BigInt* lhs, BigInt* rhs, bool* carry);    /* inplace */

BigInt* bint_mulClassical(BigInt* lhs, BigInt* rhs);
BigInt* bint_mulKaratsuba(BigInt* lhs, BigInt* rhs);
BigInt* bint_divClassical(BigInt* lhs, BigInt* rhs, BigInt** rem);

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

ull div3by2(BigInt* u, BigInt* v) {
    BigInt* rn = bint_fromWord(bigDiv(1, 0, v->values[1], NULL, NULL));
    bint_leftWordShift(rn, 2);
    for (uint n = 0; n < 6; ++n) {
        BigInt* a = bint_mulClassical(rn, rn);
        BigInt* b = bint_mulClassical(a, v);
        bint_rightWordShift(b, 4);

        bint_mulWord(rn, 2, 0);
        bint_sub(rn, b, NULL);

        bint_destroy(a);
        bint_destroy(b);
    }

    BigInt* q1 = bint_mulClassical(u, rn);
    ull q2 = q1->values[4];
    bint_destroy(q1);
    bint_destroy(rn);

    return q2;
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

void bint_print(BigInt* num) {
    BigInt* dec = bint_toDecClassical(num);

    printf("%llu", dec->values[dec->length - 1]);;
    for (uint i = dec->length - 1; i != 0;) {
        --i;
        printf("%019llu", dec->values[i]);
    }
    printf("\n");

    bint_destroy(dec);
}

void bint_printDivAndConq(BigInt* num) {
    BigInt* dec = bint_toDecDivAndConq(num);

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

BigInt* bint_toDecDivAndConq(BigInt* num) {
    const double binToDecLengthFactor = (double) WORD_LENGTH * M_LN2 / M_LN10 / DECIMAL_LENGTH;
    const int CUTOFF = 10;
    if (num->length < CUTOFF) {
        return bint_toDecClassical(num);
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

    BigInt *numLo, *numHi;
    numHi = bint_div(num, divisor, &numLo);

    bint_destroy(divisor);

    BigInt* retHi = bint_toDecDivAndConq(numHi);
    BigInt* retLo = bint_toDecDivAndConq(numLo);

    memcpy(ret->values,              retLo->values, sizeof(ull) * retLo->length);
    memcpy(ret->values + wordLength, retHi->values, sizeof(ull) * retHi->length);

    bint_destroy(numHi);
    bint_destroy(numLo);
    bint_destroy(retHi);
    bint_destroy(retLo);

    return bint_shrink(ret);
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

void* _thread_bint_mulKaratsuba(void* _args) {
    BigInt** args = (BigInt**) _args;

    BigInt* lhs = args[0];
    BigInt* rhs = args[1];
    free(args);

    return (void*) bint_mulKaratsuba(lhs, rhs);
}

BigInt* bint_mul(BigInt* lhs, BigInt* rhs, uint threads) {
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

        bint_leftWordShift(part, lengthPerThread * i);
        bint_add(ret, part);

        bint_destroy(part);
    }

    free(slices);
    free(threadPool);
    return ret;
}

BigInt* bint_div(BigInt* lhs, BigInt* rhs, BigInt** _rem) {
    BigInt *quot, *rem;
    quot = bint_divClassical(lhs, rhs, &rem);
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

    return bint_shrink(ret);
}

// Based on Knuth, The Art of Computer Programming, Vol II, pg 272, Algorithm D
BigInt* bint_divClassical(BigInt* lhs, BigInt* rhs, BigInt** rem) {
    BigInt* u = bint_clone(lhs);
    BigInt* v = bint_clone(rhs);

    uint n = rhs->length;
    uint m = lhs->length - rhs->length;

    BigInt* q = malloc(sizeof(BigInt));
    q->length = m + 1;
    q->values = malloc(sizeof(ull) * q->length);

    // Normalize (D1)
    ull d;
    if (v->values[n-1] > HALF_SIZE) {
        d = 1;
    }
    else {
        d = HALF_SIZE;
    }
    bint_mulWord(u, d, 0);
    bint_mulWord(v, d, 0);

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

        ull q_hat = div3by2(&uHi, &vHi);

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
        bint_divWord(u, d, NULL);
        *rem = bint_shrink(u);
    }
    else {
        bint_destroy(u);
    }
    bint_destroy(v);

    return bint_shrink(q);
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

BigInt* bint_rightWordShift(BigInt*lhs, uint words) {
    lhs->length -= words;
    for (uint i = 0; i < lhs->length; ++i) {
        lhs->values[i] = lhs->values[i + words];
    }

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
    bint_leftWordShift(lhs, rhsExp);

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

        BigInt* prod = bint_mul(ret, partial, threads);
        bint_destroy(ret);
        bint_destroy(partial);

        ret = prod;
    }

    bint_printDivAndConq(ret);

    bint_destroy(ret);
    free(threadPool);

    return 0;
}
