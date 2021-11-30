#ifndef MATRIX_HPP
#define MATRIX_HPP
#include <vector>
#include "Vector.hpp"
class Matrix
{
private:
    int rows,columns;
    std::vector<std::vector<double>> matrix;
    friend class Vector;
public:
    
    /*Constructors*/

    /*1*/
    Matrix(int rows,int columns);//reserve causes memory leak
    /*2*/
    Matrix(int rows,int columns,std::vector<std::vector<double>>& matrix);
    /*3*/
    Matrix(const Matrix& Matrix_);

    /*Destructor*/
    ~Matrix();

    /*Getters and setters*/
    /*1*/
    int getNumRows();
    /*2*/
    int getNumCols();
    /*3*/
    double getElement(int row, int column);
    /*4*/
    bool setElement(int row,int column,double valueElement );

    /*Display */
    void printMatrix();

    /*Data filling operation*/
    /*1*/
    void id();
    /*2*/
    void fill(double value);
    /*3*/
    void randomize();
    
    /*Boolean Methodes*/
    /*1*/
    bool isRowNull(int row);
    /*2*/
    bool isColumnNull(int col);
    /*3*/
    bool isNullMatrix();
    /*4*/
    bool isIdMatrix();
    /*5*/
    bool isSymetricMatrix();
    /*6*/
    bool isRowVector();
    /*7*/
    bool isColumnVector();
    /*8*/
    bool isVector();
    /*9*/
    bool isDouble();
    
    /*Swappers*/
    /*1*/
    void swapRows(int row1,int row2);
    /*2*/
    void swapColumns(int column1,int column2);
    
    /*Operations*/
    /*1*/
    Matrix operator*(const Matrix& Matrix_);
    /*2*/
    Matrix operator*(const double scalar);
    /*3*/
    Matrix operator+(const Matrix& Matrix_);
    /*4*/
    bool operator==(const Matrix& Matrix_);

    /*internal methods*/
    /*1*/
    void thisScalarMultiplication(const double scalar);
    /*2*/
    bool multiplyRowScalar(int row,double scalar);
    /*3*/
    bool multiplyColScalar(int col,double scalar);
    /*4*/
    double getMinMatrix();
    /*5*/
    double getMaxMatrix();

    /*Numerical methods*/
    /*1*/
    Matrix getTranspose();
    /*2*/
    std::vector<Matrix> PLU();
    /*3*/
    double determinant();
    /*4*/
    Matrix Choleskey();
    /*5*/
    Matrix getInverse();
    /*Vector transformation*/
    Vector vectorTransform();
    /*Static methods*/
    static Matrix HilbertMatrix(int n);
    /*Solving Equations System*/
};



#endif


