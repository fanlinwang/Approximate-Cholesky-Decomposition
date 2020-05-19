#include <iostream>
int main ()
{
    double a[10000];
    for (int i = 0; i < 10000; i++)
    {
        a[i] = double(i);
        a[i] *= a[i];
        std::cout << a[i] << std::endl;
    }
    return 0;
}