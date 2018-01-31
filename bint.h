#ifndef BINT_H_DEFINED
#define BINT_H_DEFINED

#include <limits.h>

typedef unsigned long long bint_word_t;
typedef unsigned int bint_exp_t;

typedef unsigned int uint;

extern const bint_exp_t  WORD_LENGTH;
extern const bint_word_t WORD_MAX;
extern const bint_exp_t  HALF_LENGTH;
extern const bint_word_t HALF_SIZE;

extern const bint_exp_t  DECIMAL_LENGTH;
extern const bint_word_t DECIMAL_SIZE;

typedef struct bint_struct {
    bint_word_t* values;
    bint_exp_t   length;
} bint_struct;

typedef bint_struct* bint_t;

// Any function labeled as inplace will store its result in lhs without
// allocating any new memory.

/****************
 * Miscelaneous *
 ****************/

// Returns a bint_t representing value.
bint_t bint_fromWord(bint_word_t value);

// Returns (1e19)^exp
bint_t bint_powerOfTen(bint_exp_t exp);

// Returns a new bint_t with value equal to num.
bint_t bint_clone(bint_t num);

// Removes any leading zeros present in num.
bint_t bint_shrink(bint_t num);

// Frees all memory used by num.
void bint_destroy(bint_t num);

// Prints the value of num in base 10 to stdout.
void bint_print(bint_t num);

// Prints the value of num in base 10 to stdout using threads threads.
void bint_printThreaded(bint_t num, uint threads);

// Calculates an equivalent representation of num with radix 1e19.
bint_t bint_toDecClassical(bint_t num);

// Calculates an equivalent representation of num with radix 1e19.
bint_t bint_toDecDivAndConq(bint_t num, uint threads);

/***************
 * Comparisons *
 ***************/

// Returns (lhs < rhs)
int bint_lessThan(bint_t lhs, bint_t rhs);

// Returns (lhs > rhs)
int bint_greaterThan(bint_t lhs, bint_t rhs);

// Returns (lhs == rhs)
int bint_equals(bint_t lhs, bint_t rhs);

/********************************************
 * Arithmetic on a bint_t and a bint_word_t *
 ********************************************/

// Returns lhs + rhsVal * 2^(WORD_LENGTH * rhsExp).
bint_t bint_addWord(bint_t lhs, bint_word_t rhsVal, bint_exp_t rhsExp); /* inplace */

// Returns  lhs - rhsVal * 2^(WORD_LENGTH * rhsExp).
// If neg is not NULL, the value it points to will be set to 0 if the result is
// positive and set to 1 if the result is negative.
bint_t bint_subWord(bint_t lhs, bint_word_t rhsVal, bint_exp_t rhsExp, int* neg); /* inplace */

// Returns lhs * rhsVal * 2^(WORD_LENGTH * rhsExp).
bint_t bint_mulWord(bint_t lhs, bint_word_t rhsVal, bint_exp_t rhsExp); /* inplace */

// Returns floor(lhs / rhsVal).
// If rem is not NULL, the value it points to will be set to the remainder
bint_t bint_divWord(bint_t lhs, bint_word_t rhsVal, bint_word_t* rem); /* inplace */

// Returns the remainder of lhs / rhsVal.
bint_word_t bint_modWord(bint_t lhs, bint_word_t rhsVal);

// Returns lhs shifted left by (bits + WORD_LENGTH * words) bits.
bint_t bint_leftShift(bint_t lhs, uint bits, bint_exp_t words); /* inplace */

// Returns rhs shifted right by (bits + WORD_LENGTH * words) bits.
bint_t bint_rightShift(bint_t lhs, uint bits, bint_exp_t words); /* inplace */

/*************************************
 * Arithmetic on two bint_t operands *
 *************************************/

// Returns lhs + rhs.
bint_t bint_add(bint_t lhs, bint_t rhs); /* inplace */

// Returns lhs - rhs.
// If neg is not NULL, the value it points to will be set to 0 if the result is
// positive and set to 1 if the result is negative.
bint_t bint_sub(bint_t lhs, bint_t rhs, int* neg); /* inplace */

// Returns lhs * rhs.
bint_t bint_mul(bint_t lhs, bint_t rhs);

// Returns floor(lhs / rhs).
// If rem is not NULL, the value it points to will be set to the remainder.
bint_t bint_div(bint_t lhs, bint_t rhs, bint_t* rem);

/*****************************
 * Multi-Threaded Arithmetic *
 *****************************/

// Returns lhs * rhs using threads threads.
bint_t bint_mulThreaded(bint_t lhs, bint_t rhs, uint threads);

// Returns floor(lhs / rhs) using threads threads.
// If rem is not NULL, the value it points to will be set to the remainder
bint_t bint_divThreaded(bint_t lhs, bint_t rhs, bint_t* rem, uint threads);

/***********************
 * Specific Algorithms *
 ***********************/

// Returns the value of lhs + rhs without performing the final carry.
// If carry is not NULL, the value it points to will be set to 1 if a carry is
// required (but wasn't performed) and 0 otherwise.
bint_t bint_addNoLastCarry(bint_t lhs, bint_t rhs, int* carry); /* inplace */

// Returns the value of lhs - rhs without performing the final borrow.
// If carry is not NULL, the value it points to will be set to 1 if a borrow is
// required (but wasn't performed) and 0 otherwise.
bint_t bint_subNoLastCarry(bint_t lhs, bint_t rhs, int* carry); /* inplace */

// Returns lhs * rhs using a the classic long multiplication algorithm.
// Technically this is multiplication by parts, but it is equivalent.
//
// Runs in O(n^2) time.
// (with n = max(lhs->length, rhs->length))
bint_t bint_mulClassical(bint_t lhs, bint_t rhs);

// Returns lhs * rhs using the Karatsuba algorithm
//
// Runs in O(n^(log_2(3))) time, that's roughly O(n^1.585).
// (with n = max(lhs->length, rhs->length))
bint_t bint_mulKaratsuba(bint_t lhs, bint_t rhs);

// Returns floor(lhs / rhs) using classic long division.
// If rem is not NULL, the value it points to will be set to the remainder.
// This implementation is based on:
//      Donald Knuth
//      The Art of Computer Programming, Vol II,
//      pg 272, Algorithm D
//
// Runs in O((m-n)*n) time.
// (with m = lhs->length, n = rhs->length)
bint_t bint_divClassical(bint_t lhs, bint_t rhs, bint_t* rem);

// Returns floor(lhs / rhs) using recursive division.
// If rem is not NULL, the value it points to will be set to the remainder.
// This implementation is based on:
//      Christoph Burnikel, and Joachim Ziegler
//      Fast Recursive Division (1998)
//
// Runs in O(n log(n)) time.
// (with n = rhs->length);
bint_t bint_divDivAndConq(bint_t lhs, bint_t rhs, bint_t* rem, uint threads);

/******************************
 * Abritrary Radix Arithmetic *
 ******************************/

// Returns lhs + rhsVal * rad^rhsExp.
// This function enforces (and expects) a radix of rad.
bint_t bint_radAddWord(bint_t lhs, bint_word_t rhsVal, bint_word_t rhsExp, bint_word_t rad); /* inplace */

// Returns lhs * rhsVal * rad^rhsExp.
// This function enforces (and expects) a radix of rad.
bint_t bint_radMulWord(bint_t lhs, bint_word_t rhsVal, bint_word_t rhsExp, bint_word_t rad); /* inplace */

// Returns lhs + rhs.
// This function enforces (and expects) a radix of rad.
bint_t bint_radAdd(bint_t lhs, bint_t rhs, bint_word_t rad); /* inplace */

#endif
