#ifndef VEC_H//if not defined
#define VEC_H//then define it (nothing defined twice or more)
#include <vector>
//using namespace std;
using std::vector;

class vec
{
public:
    vec(int dim);

    //int dimension = 0; //always set a value
    int dimension() {return components.size();}

    vector<double> components;

    double &operator[](int index){return components[index]; }//gives copies without &
    //If you want to send in the actual value, use &operator
    double &operator()(int index){return components.at(index); }
    //gives error if there is segmentation fault (vector not that big)
    vec operator+(vec rhs); //rhs=right hand side. plus operator(see vec.cpp for spez.)
    //vec operator+(double rhs);//add a number to every element in vec

    double lengthSquared();
    double length();
};

#endif // VEC_H









