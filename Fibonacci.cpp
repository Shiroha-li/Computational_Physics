#include <iostream>
using namespace std;
int main()
{
    int index [3] ={7,15,31};
    int Fibonacci [31] = {0};
    Fibonacci[0] = 1;
    Fibonacci[1] = 1;
    for (int i = 2;i < 31;i++){
        Fibonacci[i] = Fibonacci[i-1] + Fibonacci[i-2];
    }
    for (int i = 0;i < 3;i++){
        cout << "F_" << index[i] << " = " << Fibonacci[index[i]-1] << "\n";
    }
    return 0;
}