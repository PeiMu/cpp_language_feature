#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

/*
 * arg0: lhs of shift
 * arg1: rhs of shift
 * arg2: shift operation name
 * */
int main(int argc, char** argv)
{
    int parameters[2];
    char* pEnd;
    if (argc == 4) {
        for (size_t idx = 0; idx < argc-2; idx++) {
            parameters[idx] = atoi(argv[idx+1]);
        }
    }
    else {
	printf("error input! arg0: lhs of shift, arg1: rhs of shift, arg2: shift operation name, e.g. << or >>");
    }


    if (argv[3] == "<<") {
	printf("%d\n", parameters[0] << parameters[1]);
    } else {
	printf("%d\n", parameters[0] >> parameters[1]);
    }

    return 0;
}
