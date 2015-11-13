#include "vec.h"
#include <iostream>
#include <cmath>
using std::cout; using std::endl;

vec::vec(int dim)
{
    components.resize(dim, 0);

}

vec vec::operator+(vec rhs)
{
    if(this->dimension() != rhs.dimension()){
        cout << "Tried to add two vectors of different dimensions" << endl;
        //when the object is an pointer it is rached to with ->
        exit(1);
        //now we know that the diensions are the same(if message not appears)
    }

    vec newVector(this->dimension());
    for(int i = 0; i<this->dimension(); i++){
        newVector[i] = components[i] + rhs[i];
    }
    return newVector;
}


double vec::length(){
    return sqrt(lengthSquared());
}

double vec::lengthSquared()
{
  double sum = 0;
//    for(int i=0; i<dimension(); i++){
//        sum+= components[i]*components[i];
//    }


    //same as above:

    for(double value : components){
        sum += value*value; //this for-loop is preferred.
        //loop through all elements in components, and call them value
    }

    return sum;
}












