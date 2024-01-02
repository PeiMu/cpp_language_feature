/*
 * Run with command: clang function_test.c -S -emit-llvm -O3 -o function_test.ll
 * */

#include <stdio.h>

__attribute__((noinline))
int testFunc(int x) {
	int a = x + 3;
	if (a > 20) {
		a = a + 100;
	} else {
		a = a * 100;
	}
	return a;
}

int main() {
	int x;
	scanf("%d", &x);
	__builtin_assume(x > -16 && x < 16);
//	if (x <= -16 || x >= 16) {
//		__builtin_unreachable();
//	}
	int a = testFunc(x);
	printf("a = %d\n", a);
	return 0;
}
