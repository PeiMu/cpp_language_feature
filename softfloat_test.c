#include <stdio.h>

unsigned long long convertDouble(double my_double)
{
    unsigned long long my_ulonglong;

    my_ulonglong = *(unsigned long long*)(&my_double);

    return my_ulonglong;
}

double convertBack(unsigned long long my_ulonglong)
{
    double my_double;

    my_double = *(double*)(&my_ulonglong);

    return my_double;
}


int main()
{
    double my_double = 1.5;
    unsigned long long my_ulonglong = convertDouble(my_double);

    printf("my_double: %lf\n", convertBack(my_ulonglong));
    printf("my_ulonglong: 0x%llx\n", my_ulonglong);

    return 0;
}