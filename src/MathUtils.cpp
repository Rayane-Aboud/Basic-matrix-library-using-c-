#include "../headers/MathUtils.hpp"

#include<math.h>

int MathUtils::sign(double d)
{
    if (d<0)
        return -1;
    else 
        return 1;
}

double MathUtils::abs(double d){
    if(d<0) return -1*d;
    else return d;
}
