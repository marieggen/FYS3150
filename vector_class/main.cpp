#include <iostream>
#include "vec.h"
#include <vector>

using namespace std;

int main()
{
//    vec myVector(3);
//    myVector[0] = 3;
//    cout << "Dimension:" << myVector.dimension << endl;
//    cout << "Hello World!" << endl;

//    double* myVector = new double[10]; example vetors
//    vector<double> myAwesomeVector;
//    myAwesomeVector.resize(20);
    vec myVector(3);
    cout << "Dimension:" << myVector.dimension() << endl;
    myVector[0] = 18;
    myVector[1] = 2;
    myVector[2] = 5;
    cout << myVector[0] << endl;

    vec myVector1(3);
    cout << "Dimension:" << myVector.dimension() << endl;
    myVector1[0] = 1;
    myVector1[1] = 2;
    myVector1[2] = 5;
    cout << myVector1[0] << endl;

    vec myVector2 = myVector+myVector1;
    cout << myVector2[0] << "," << myVector2[1] << "," << myVector2[2] << endl;

    cout << "LengthSquare: " << myVector2.lengthSquared() << endl;


    return 0;
}

