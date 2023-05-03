/*
 * Run with command: clang main.c -S -emit-llvm -O3 -o main.ll
 * */

#include <stdio.h>

int main() {
	int x;
	scanf("%d", &x);
	__builtin_assume(x > -16 && x < 16);
	int a = x + 3;
	if (a > 20) {
		a = a + 100;
	} else {
		a = a * 100;
	}
	printf("a = %d\n", a);
	return 0;
}
