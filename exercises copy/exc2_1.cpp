using namespace std;
#include <iostream>
#include <math.h>


int main()
  {
    int i=0, ;
    double number, fract, integer;
    cout << "Please enter a floating point number:";
    cin >> number;
    fract = modf(number, &integer);

    for(i ; i<32 ; i++){}
    return 0;
  }
  
