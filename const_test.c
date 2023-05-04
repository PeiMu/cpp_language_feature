/*
 * Check if the `__builtin_assume` can broadcast constant value
 * */

#include <stdio.h>

int main() {
	int x;
	scanf("%d", &x);
	__builtin_assume(x == 10);
//	if (x != 10) {
//		__builtin_unreachable();
//	}
	int a = x + 3;
	printf("a = %d\n", a);
	return 0;
}
