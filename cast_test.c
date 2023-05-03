/*
 * Test with Cast Instruction to see if `__builtin_assume` support it.
 * */

#include <stdio.h>

#define HI(x) *(1+(int*)&x)

int main() {
	double x;
	scanf("%lf", &x);
	__builtin_assume(x > -16 && x < 16);
//	if (x <= -16 || x >= 16) {
//		__builtin_unreachable();
//	}
	int int_x = (int)x;
	printf("int_x = %d\n", int_x);
	int a = int_x + 3;
	printf("before branch, a = %d\n", a);
	if (a > 20) {
		a = a + 100;
	} else {
		a = a * 100;
	}
	printf("a = %d\n", a);
	return 0;
}
