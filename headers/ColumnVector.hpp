#ifndef COLUMNVECTOR_HPP
#define COLUMNVECTOR_HPP

#include <vector>
#include "Matrix.hpp"

class ColumnVector
{
private:
    int length;
    std::vector<double>Vector;
public:
    ColumnVector(int length);
    ColumnVector(int length,std::vector<double>& vector);
    ColumnVector(const RowVector& rowVector);
    ColumnVector(const Matrix& Matrix_,int row);

    ~ColumnVector();

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
    bool equals(ColumnVector columnVector);
    ColumnVector add(const ColumnVector&);
    /*multiplication with a col vector*/
    double multiply(double scalar);
    ColumnVector multiply(Matrix& matrix);
    ColumnVector operator*(double scalar);
    ColumnVector operator*(Matrix& matrix);
    bool operator==(const ColumnVector columnVector);
    void thisScalarMultiplication();
    //convertion
    Matrix convertToMatrix();

};



#endif