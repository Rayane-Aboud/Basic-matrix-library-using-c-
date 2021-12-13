#ifndef VECTOR_HPP
#define VECTOR_HPP
#include <vector>
class Matrix;
class Vector
{
private:
    int length;
    std::vector<double>vector;
    friend class Matrix;
public:
    Vector(int length);
    Vector(std::vector<double>& vector);
    //ColumnVector(const RowVector& rowVector);
    Vector(Matrix& Matrix_,int row);

    ~Vector();

    //getters and setters
    int getLength();
    double getElement(int index);
    bool setElement(int index,double value);
    //show
    void printVector();
    //basic operation
    void fill(double value);
    void randomize();
    bool isNullVector();
    //swapper
    bool swapElements(int index1,int index2);
    
    /*multiplication with a col vector*/
    double scalarProduct(const Vector&);
    Vector operator*(double scalar);
    Vector operator+(const Vector&);
    Vector operator-(const Vector&);
    Vector pseudoMultiplication(const Vector&);
    bool operator==(const Vector& vector);
    void thisScalarMultiplication(const double scalar);
    double euclidianNorm();
    //convertion
    Matrix convertToColMatrix();
    Matrix convertToRowMatrix();
    

};



#endif