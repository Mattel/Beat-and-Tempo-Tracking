#ifndef BTT_FIXED_MATH_H
#define BTT_FIXED_MATH_H

#include <limits.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef int32_t q31_t;

static inline q31_t q31_from_float(float v)
{
    if (v >= 0.9999999995343387f) {
        return INT32_MAX;
    }
    if (v <= -1.0f) {
        return INT32_MIN;
    }
    return (q31_t)(v * 2147483648.0f);
}

static inline float q31_to_float(q31_t v)
{
    return (float)v / 2147483648.0f;
}

static inline q31_t q31_saturating_add(q31_t a, q31_t b)
{
    int64_t res = (int64_t)a + (int64_t)b;
    if (res > INT32_MAX) return INT32_MAX;
    if (res < INT32_MIN) return INT32_MIN;
    return (q31_t)res;
}

static inline q31_t q31_saturating_sub(q31_t a, q31_t b)
{
    int64_t res = (int64_t)a - (int64_t)b;
    if (res > INT32_MAX) return INT32_MAX;
    if (res < INT32_MIN) return INT32_MIN;
    return (q31_t)res;
}

static inline q31_t q31_mul(q31_t a, q31_t b)
{
    int64_t prod = (int64_t)a * (int64_t)b;
    prod >>= 31;
    if (prod > INT32_MAX) return INT32_MAX;
    if (prod < INT32_MIN) return INT32_MIN;
    return (q31_t)prod;
}

static inline q31_t q31_mul_shr(q31_t a, q31_t b, unsigned shift)
{
    int64_t prod = (int64_t)a * (int64_t)b;
    if (shift < 31) {
        prod >>= shift;
    } else {
        prod >>= 31;
    }
    if (prod > INT32_MAX) return INT32_MAX;
    if (prod < INT32_MIN) return INT32_MIN;
    return (q31_t)prod;
}

static inline q31_t q31_shift(q31_t v, int shift)
{
    if (shift > 0) {
        int64_t res = ((int64_t)v) << shift;
        if (res > INT32_MAX) return INT32_MAX;
        if (res < INT32_MIN) return INT32_MIN;
        return (q31_t)res;
    }
    return (q31_t)(v >> (-shift));
}

#ifdef __cplusplus
}
#endif

#endif // BTT_FIXED_MATH_H
