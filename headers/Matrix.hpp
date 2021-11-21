#ifndef MATRIX_HPP
#define MATRIX_HPP
#include <vector>
class Matrix
{
private:
    //add a isMatrix boolean
    int rows,columns;
    std::vector<std::vector<double>> matrix;
public:
    //Matrix();
    Matrix(int rows,int columns);//reserve causes memory leak
    Matrix(int rows,int columns,std::vector<std::vector<double>>& matrix);
    Matrix(const Matrix& Matrix_);
    
    ~Matrix();
    //getters and setters
    int getNumRows();
    int getNumCols();
    double getElement(int row, int column);
    bool setElement(int row,int column,double valueElement );
    //show
    void printMatrix();
    //basic operation
    void id();
    void fill(double value);
    void randomize();
    bool isRowNull(int row);
    bool isColumnNull(int col);
    bool isNullMatrix();
    //
    //Swappers
    void swapRows(int row1,int row2);
    void swapColumns(int column1,int column2);
    
    //operations
    bool equals(Matrix );
    Matrix add(Matrix matrix);
    Matrix multiply(const Matrix& Matrix_);
    Matrix multiply(double scalar);
    Matrix operator*(const Matrix& Matrix_);
    Matrix operator*(const double scalar);
    Matrix operator+(const Matrix& Matrix_);
    bool operator==(const Matrix& Matrix_);
    //special
    void thisScalarMultiplication(const double scalar);
    bool multiplyRowScalar(int row,double scalar);
    bool multiplyColScalar(int col,double scalar);
    //numerical methods
    Matrix getTranspose();
    Matrix getInverse();
    std::vector<Matrix> PLU();
    std::vector<Matrix> Choleskey();
    double determinant();
    
    
};



#endif


