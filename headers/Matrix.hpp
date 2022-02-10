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
    Matrix(std::vector<std::vector<double>>& matrix);
    /*3*/
    Matrix(const Matrix& Matrix_);
    /*4*/
    Matrix(int row,int col,int size,const Matrix& Matrix_);
    /*5*/
    Matrix(int beginRow,int beginCol,int endRow,int endCol,int size,const Matrix& Matrix_);
    /*6*/
    Matrix(int size);
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
    bool isNullRow(int row);
    /*2*/
    bool isNullColumn(int col);
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
    /*10*/
    bool isTriangularSup();
    /*11*/
    bool isTriangularInf(); 

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
    /*6*/
    void thisMatrixMultiplicationLeft(const Matrix& Matrix);
    /*7*/
    void thisMatrixMultiplicationRight(const Matrix& Matrix);
    /*Numerical methods*/
    /*1*/
    Matrix getTranspose();
    /*2*/
    std::vector<Matrix> plu();
    /*3*/
    double determinant();
    /*4*/
    Matrix Choleskey();
    /*5*/
    Matrix getInverse();//ordre
    /*6*/
    std::vector<Matrix>qr();
    /*Extracting elements*/
    /*1*/
    Vector convertVector();
    /*2*/
    Vector extractRowVector(int row);
    /*3*/
    Matrix extractRowMatrix(int row);
    /*4*/
    Vector extractColumnVector(int column);
    /*5*/
    Matrix extractColumnMatrix(int column);
    /*Static methods*/
    static Matrix HilbertMatrix(int n);
    static Matrix cannonicRowMat(int range,int dimension);
    static Matrix cannonicColMat(int range,int dimension);
    static Matrix getId(int dimension);
    void placeMatrix(int i,int j,const Matrix& matrix);
    /*Solving Equations System*/
    bool isSquare();//added lately
    bool isDiagonallyDominant();//added lately
    bool isColinearWith(const Matrix& matrix);
    double getEuclidianNorm();
    Matrix extractMatrix(int rowBegin,int colBegin,int rowEnd,int colEnd);
    Vector solve(Vector y);
    Vector solve_GS(Vector y,Vector x0,int it_max,float error);
    double power_iteration(Vector x0,Vector* eigenVector,int it_max,float error);
    std::vector< Matrix> getDiagonalEquivalent_Gauss();
    Matrix getDiagonalEquivalent_LR(int);
    Matrix getDiagonalEquivalent_QR(int);
};



#endif


