#include <iostream>
#include <vector>
#include <math.h>
#include "../headers/Vector.hpp"
#include "../headers/Matrix.hpp"
Vector::Vector(int length){
    this->length=length;
    for(int i=0;i<this->length;i++){
        this->vector.push_back(0);
    }
}

Vector::Vector(std::vector<double>& vector){
    this->length=vector.size();
    for(int i=0;i<this->length;i++){
        this->vector.push_back(vector[i]);
    }
}

Vector::Vector(Matrix& Matrix_,int row){
    if(row>Matrix_.getNumRows()){
        std::cout<<"the row does not exist";
        exit(0);
    }
    this->length=Matrix_.getNumCols();
    for(int j=0;j<this->length;j++){
        this->vector.push_back(Matrix_.getElement(row,j));
    }
}

Vector::~Vector(){}

int Vector::getLength(){
    return this->length;
}

double Vector::getElement(int index){
    if (index<0||index>=this->length)
    {
        std::cout<<"couldn't get because index";
        exit(0);
    }
    return this->vector[index];

}

bool Vector::setElement(int index, double value){
    if(index<0 || index>=this->length)
        return false;
    this->vector[index]=value;
    return true;
}

//show

void Vector::printVector(){
    std::cout<<"[ ";
    for(int i=0;i<this->length;i++){
        std::cout<<this->vector[i]<<" ";
    }
    std::cout<<" ]";
}

//fill
void Vector::fill(double value){
    for (int i=0;i<this->length;i++){
        this->vector[i]=value;
    }
}
//boolean operation
bool Vector::isNullVector(){
    for(int i=0;i<this->length;i++){
        if(this->vector[i]!=0)
            return false;
    }
    return true;
}


bool Vector::swapElements(int index1,int index2){
    if(index1<0 ||index2<0)
        return false;
    double tmp;
    tmp=this->vector[index1];
    this->vector[index1]=this->vector[index2];
    this->vector[index2]=tmp;
    return true;
}
//operation
bool Vector::operator ==(const Vector& vector){
    for(int i=0;i<this->length;i++){
        if(this->vector[i]!=vector.vector[i])
            return false;
    }
    return true;
}

double Vector::scalarProduct(const Vector& vector){
    if(this->length!=vector.length)
    {
        std::cout<<"cannot do scalar product not same";
        exit(0);
    }
    double scalar_product=0;
    for(int i=0;i<this->length;i++){
        scalar_product+=this->vector[i]*vector.vector[i];
    }
    return scalar_product;
}

Vector Vector::operator*(double scalar){
    Vector retVect(*this);
    for(int i=0;i<this->length;i++)
        retVect.vector[i]=this->vector[i]*scalar;
    return retVect;
}



void Vector::thisScalarMultiplication(const double scalar){
    for(int i=0;i<this->length;i++){
        this->vector[i]*=scalar;
    }
}

Matrix Vector::convertToColMatrix(){
    Matrix retMatrix(this->length,1);
    for(int i=0;i<this->length;i++){
        retMatrix.matrix[i][1]=this->vector[i];
    }
    return retMatrix;
}

Matrix Vector::convertToRowMatrix(){
    Matrix retMatrix(1,this->length);
    for(int j=0;j<this->length;j++){
        retMatrix.matrix[1][j]=this->vector[j];
    }
    return retMatrix;
}