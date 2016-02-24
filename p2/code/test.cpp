#include<iostream>
#include<fstream>


double* f(void )
{
    double* x = new double[1];
    x[0] = 10213123;
    return x;
}

void g(double** x)
{
    (*x) = new double[1];
    (*x)[0] = 123123213;
    return ;
}
int main()
{
    FILE* origin;
    origin = fopen("./at.txt", "w");
    fprintf(origin, "ASDFASDFDSA\n");
    fclose(origin);
    return 0;
}
