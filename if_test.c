/*
 * Run with command: clang if_test.c -S -emit-llvm -O3 -o if_test.ll
 * */

#include <stdio.h>

int main() {
	int x;
	scanf("%d", &x);
	__builtin_assume(x > -16 && x < 16);
//	if (x <= -16 || x >= 16) {
//		__builtin_unreachable();
//	}
	int a = x + 3;
	if (a > 20) {
		a = a + 100;
	} else {
		a = a * 100;
	}
	printf("a = %d\n", a);
	return 0;
}
