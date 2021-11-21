#include <iostream>
#include <vector>
#include "../headers/Includes.hpp"
/*Matrix::Matrix()
{   
    this->isMatrix=false;
    this->rows=1;
    this->columns=1;
    std::vector<double> vector;
    vector.push_back(0);
    this->matrix.push_back(vector);
    this->matrix[0][0]=0;
}
*/
Matrix::Matrix(int rows,int columns)
{
    this->rows=rows;
    this->columns=columns;
    if(this->rows<1 || this->columns<1)
    {
        std::cout<<"an error has occured! the matrix input is not valid";
        exit(0);

    }
    /*this->matrix.reserve(this->rows);
    for(int i=0;i<this->rows;i++){
        this->matrix[i].reserve(this->columns);
    }*/
    //push_back vs emplace_back
    for(int i=0;i<this->rows;i++){
        std::vector<double> vector;
        for(int j=0;j<this->columns;j++){
            vector.push_back(0);
        }
        this->matrix.push_back(vector);
    }

}
Matrix::Matrix(int rows,int columns,std::vector<std::vector<double>>& matrix){
    this->rows=rows;
    this->columns=columns;
    for(int i=0;i<this->rows;i++){
        std::vector<double> vector;
        for(int j=0;j<this->columns;j++){
            vector.push_back(matrix[i][j]);
        }
        this->matrix.push_back(vector);
    }
}
Matrix::Matrix(const Matrix& Matrix_){
    this->rows=Matrix_.rows;
    this->columns=Matrix_.columns;
    for(int i=0;i<this->rows;i++){
        std::vector<double> vector;
        for(int j=0;j<this->columns;j++){
            vector.push_back(Matrix_.matrix[i][j]);
        }
        this->matrix.push_back(vector);
    }
}
//getters and setters-----------------------
int Matrix::getNumRows(){
    return this->rows;
}

int Matrix::getNumCols(){
    return this->columns;
}

double Matrix::getElement(int row, int column){
    return this->matrix[row][column];
}

bool Matrix::setElement(int row, int column,double valueElement){
    if(!(row<this->rows && column<this->columns))
        return false;
    this->matrix[row][column]=valueElement;
    return true;
}
//show---------------------------------------------
void Matrix::printMatrix(){
        
    for(int i=0;i<this->rows;i++){
        for (int j=0;j<this->columns;j++){
            std::cout<<this->matrix[i][j]<<"  ";
        }
        std::cout<<std::endl;
    }
}

bool Matrix::isRowNull(int row){
    if (this->rows<row)
        return false;
    for(int j=0;j<this->columns;j++){
        if (this->matrix[row][j]!=0) return false;
    }
    return true;
}

bool Matrix::isColumnNull(int col){
    if (this->columns<col)
        return false;
    for(int i=0;i<this->rows;i++){
        if(this->matrix[i][col]!=0) return false;
    }
    return true;
}

bool Matrix::isNullMatrix(){
    for(int i=0;i<this->rows;i++){
        for (int j=0;j<this->columns;j++){
            if(this->matrix[i][j]!=0)   return false;
        }
    }
    return true;
}
//Swapper----------------------------------------------
void Matrix::swapRows(int row1,int row2){
    int tmp;
    for(int j=0;j<this->columns;j++){
        tmp=this->matrix[row1][j];
        this->matrix[row1][j]=this->matrix[row2][j];
        this->matrix[row2][j]=tmp;
    }
}

void Matrix::swapColumns(int col1,int col2){
    int tmp;
    for(int i=0;i<this->rows;i++){
        tmp=this->matrix[i][col1];
        this->matrix[i][col1]=this->matrix[i][col2];
        this->matrix[i][col2]=tmp;

    }
}



//basic operation----------------------------
void Matrix::id(){
    for(int i=0;i<this->rows;i++){
        for(int j=0;j<this->columns;j++){
            if(i==j)
                this->matrix[i][j]=1;
            else
                this->matrix[i][j]=0;
        }
    }
}
void Matrix::fill(double value){
    for(int i=0;i<this->rows;i++)
        for(int j=0;j<this->columns;j++)
            this->matrix[i][j]=value;
}

//void Matrix::randomize(){}
/*Matrix Matrix::multiply(const Matrix& Matrix_){
    Matrix returnMatrix(this->rows,Matrix_.columns);
    if (this->columns!=Matrix_.rows)
    {
        std::cout<<"cannot multiply the two matrix:Math Error";
        exit(0);
    }
    
    double retElement;
    for(int i=0;i<this->rows;i++){
        for(int j=0;j<Matrix_.columns;j++){
            retElement=0;
            for(int k=0;k<this->columns;k++){
                retElement+=this->matrix[i][k]*Matrix_.matrix[k][j];
            }
            returnMatrix.setElement(i,j,retElement);
        }
    }
    return returnMatrix;
}
*/
Matrix Matrix::operator*(const Matrix& Matrix_){
    if(this->columns!=Matrix_.rows)
    {
            std::cout<<"cannot multiply the two matrix:Math Error";
            exit(0);
    }
    Matrix retMatrix(this->rows,Matrix_.columns);
    for(int i=0;i<retMatrix.rows;i++){
        for(int j=0;j<retMatrix.columns;j++){
            for(int k=0;k<this->columns;k++){
                retMatrix.matrix[i][j]+=this->matrix[i][k]*Matrix_.matrix[k][j];
            }
        }
    }
    return retMatrix;
}

Matrix Matrix::operator*(const double scalar){
    Matrix retMatrix(*this);
    for(int i=0;i<this->rows;i++){
        for(int j=0;j<this->columns;j++){
            retMatrix.matrix[i][j]*=scalar;
        }
    }
    return retMatrix;
}

void Matrix::thisScalarMultiplication(const double scalar){
    for(int i=0;i<this->rows;i++){
        for(int j=0;j<this->columns;j++)
            this->matrix[i][j]*=scalar;
    }
}





Matrix Matrix::operator+(const Matrix& Matrix_){
    if(this->rows!=Matrix_.rows||this->columns!=Matrix_.columns )
    {
            std::cout<<"cannot add two matrixes of different size:Math Error";
            exit(0);
    }
    Matrix sumMatrix(this->rows,this->columns);
    for(int i=0;i<this->rows;i++){
        for(int j=0;j<this->columns;j++){
            sumMatrix.matrix[i][j]=this->matrix[i][j]+Matrix_.matrix[i][j];
        }
    }
    return sumMatrix;
}


Matrix Matrix::getTranspose(){
    Matrix retMatrix(this->columns,this->rows);
    for(int i=0;i<this->rows;i++){
        for(int j=0;j<this->columns;j++){
            retMatrix.matrix[j][i]=this->matrix[i][j];
        }
    }
    return retMatrix;
}
/*
Matrix Matrix::getInverse(){
    if (this->rows!=this->columns)
    {
        std::cout<<"cannot get the inverse! it is not a square matrix";
        exit(0);
    }
    Matrix thisMatrix(*this);
    Matrix inverseMatrix(this->rows,this->columns);
}
*/

std::vector<Matrix> Matrix::PLU(){
    std::vector<Matrix> result;
    if(this->rows!=this->columns){
        std::cout<<"not yet buddy";
        exit(0);
    }
    int pivot;
    Matrix U(*this);
    Matrix L(this->rows,this->columns);L.id();
    Matrix P(this->rows,this->columns);P.id();
    for(int i=0;i<this->rows;i++){
        
        //regler le probleme du pivot
        pivot=i;
        
        if(this->matrix[pivot][pivot]==0){
            int pivotrow=pivot+1;
            for(pivotrow;pivotrow<this->rows;pivotrow++){
                if (this->matrix[pivotrow][pivot]!=0) break;
            }
            if (pivotrow==this->rows)continue;
            P.swapRows(pivot,pivotrow);
            P.thisScalarMultiplication(-1);
            L.swapRows(pivot,pivotrow);
            L.swapColumns(pivot,pivotrow);
            U.swapRows(pivot,pivotrow);
        }
        //effectuer les operations
        for(int j=pivot+1;j<this->rows;j++){
            L.matrix[j][pivot]=U.matrix[j][pivot]/U.matrix[pivot][pivot];
            //printf("%f \n",L.matrix[j][pivot]);
            for(int k=0;k<this->columns;k++){
                U.matrix[j][k]+=-1*L.matrix[j][pivot]*U.matrix[pivot][k];
                
            }
            if(U.isRowNull(j)){
                std::cout<<"matrice non diagonalisable ! ligne:"<<j<<" est null";
                //return results that makes the determinant equals to 0
                P.fill(0);
                result.push_back(P);
                result.push_back(L);
                result.push_back(U);
                return result;
            }
        }
    }
    result.push_back(P);
    result.push_back(L);
    result.push_back(U);
    return result;
}


double Matrix::determinant(){
    if(this->rows!=this->columns){
        return 0;
    }
    std::vector<Matrix>pluList=this->PLU();
    if(pluList[0].isNullMatrix())
        return 0;
    int sign;
    for(int i=0;i<this->rows;i++){
        if(pluList[0].matrix[i][0]!=0){
            sign=pluList[0].matrix[i][0];
            break;
        }
    }
    double determinant=1;
    for(int i=0;i<pluList[0].rows;i++){
        determinant=determinant*pluList[1].matrix[i][i];
    }
    for(int i=0;i<pluList[0].columns;i++){
        determinant=determinant*pluList[2].matrix[i][i];
    }
    return determinant*sign;
}




Matrix::~Matrix()
{
}
