#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

uint64_t asuint(double t) {
  union {
    double f;
    uint64_t k;
  } u = {t};
  return u.k;
}

uint32_t asuint_32(double f) {
  union {
    float f;
    uint32_t i;
  } u = {static_cast<float>(f)};
  return u.i;
}

double asdouble(uint64_t t) {
  union {
    uint64_t k;
    double f;
  } u = {t};
  return u.f;
}

typedef union
{
  double value;
  struct
  {
    __uint32_t lsw;
    __uint32_t msw;
  } parts;
} ieee_double_shape_type;

#define EXTRACT_WORDS(ix0,ix1,d)				\
do {								\
  ieee_double_shape_type ew_u;					\
  ew_u.value = (d);						\
  (ix0) = ew_u.parts.msw;					\
  (ix1) = ew_u.parts.lsw;					\
} while (0)

int main(int argc, char **argv) {
  char* pEnd;
  if (argc == 2)
  {
    double value = strtod(argv[1], &pEnd);
    uint64_t long_word = *reinterpret_cast<uint64_t*>(&value);
    int32_t long_word_32 = *reinterpret_cast<int32_t*>(&value);
    printf("size of uint64 = %d, size of uint32 = %d, size of double = %d, size of float = %d\n", sizeof(uint64_t), sizeof(uint32_t), sizeof(double), sizeof(float));
    uint32_t msw = static_cast<uint32_t>(long_word >> (32 * 1));
    uint32_t lsw = static_cast<uint32_t>(long_word >> (32 * 0));
    printf("msw = %d, lsw = %d\n", msw, lsw);
    printf("result of asuint = %lu\t%lx\t%lu\n", asuint(value), asuint(value), asuint_32(value));
    printf("result of asdouble = %f\n", asdouble(value));
    uint32_t hx, lx;
    EXTRACT_WORDS(hx, lx, value);
    printf("result of GET_HIGH_WORD = %d\n", hx);
    printf("result of GET_LOW_WORD = %d\n", lx);

    printf("value: %f\n", value);
    printf("static cast: %lu\n", static_cast<uint64_t>(value));
    printf("long_word_32 static cast: %lu\n", static_cast<double>(long_word_32));
    printf("32bits static cast: %d\n", static_cast<uint32_t>(long_word>>32));
    uint64_t u64Num;
    memcpy(&u64Num, &value, sizeof(double));
    printf("memcpy: %lu\t%lx\t%lf\n", u64Num, u64Num, u64Num);
    double dbNum;
    memcpy(&dbNum, &u64Num, sizeof(uint64_t));
    printf("memcpy back: %f\n", dbNum);
    float float_value = (float)value;
    uint64_t float_u64Num = 0;
    memcpy(&float_u64Num, &float_value, sizeof(float));
    printf("32bits memcpy: %lu\n", float_u64Num);
    float fpNum;
    memcpy(&fpNum, &float_u64Num, sizeof(uint64_t));
    printf("32bits memcpy back: %f\n", fpNum);
    printf("%f\n", static_cast<double>(float_u64Num));
    printf("%lx\n", *(unsigned long*)(&value));
    printf("%lx\n", &value);
    printf("%lu\n", (unsigned long)(value));
  }
  return 0;
}

