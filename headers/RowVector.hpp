#ifndef ROWVECTOR_HPP
#define ROWVECTOR_HPP

#include <vector>
#include "Matrix.hpp"

class RowVector
{
private:
    int length;
    std::vector<double>Vector;
public:
    RowVector(int length);
    RowVector(int length,std::vector<double>& vector);
    RowVector(const RowVector& rowVector);
    RowVector(const Matrix& Matrix_,int row);

    ~RowVector();

    //getters and setters
    int getLength();
    double getElement(int row,int column);
    bool setElement(int row,int column,double value);
    //show
    void printVector();
    //basic operation
    void fill(double value);
    void randomize();
    bool isVectorNull();
    //swapper
    bool swapElements(int index1,int index2);
    //operation
    bool equals(RowVector rowVector);
    RowVector add(const RowVector&);
    /*multiplication with a col vector*/
    double multiply(double scalar);
    RowVector multiply(Matrix& matrix);
    RowVector operator*(double scalar);
    RowVector operator*(Matrix& matrix);
    bool operator==(const RowVector rowVector);
    void thisScalarMultiplication();
    //convertion
    Matrix convertToMatrix();

};



#endif