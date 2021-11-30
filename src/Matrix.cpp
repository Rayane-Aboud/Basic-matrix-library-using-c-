#include <iostream>
#include <vector>
#include <math.h>
#include "../headers/Matrix.hpp"

/*Constructors*/

/*1*/
Matrix::Matrix(int rows,int columns)
{
    this->rows=rows;
    this->columns=columns;
    if(this->rows<1 || this->columns<1)
    {
        std::cout<<"an error has occured! the matrix input is not valid";
        exit(0);

    }
    
    for(int i=0;i<this->rows;i++){
        std::vector<double> vector;
        for(int j=0;j<this->columns;j++){
            vector.push_back(0);
        }
        this->matrix.push_back(vector);
    }

}

/*2*/
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

/*3*/
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

/*Destructor*/
Matrix::~Matrix()
{
}


/*Getters and setters*/
/*1*/
int Matrix::getNumRows(){
    return this->rows;
}
/*2*/
int Matrix::getNumCols(){
    return this->columns;
}
/*3*/
double Matrix::getElement(int row, int column){
    return this->matrix[row][column];
}
/*4*/
bool Matrix::setElement(int row, int column,double valueElement){
    if(!(row<this->rows && column<this->columns))
        return false;
    this->matrix[row][column]=valueElement;
    return true;
}
/*Display */
void Matrix::printMatrix(){
        
    for(int i=0;i<this->rows;i++){
        for (int j=0;j<this->columns;j++){
            std::cout<<this->matrix[i][j]<<std::endl;
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}



/*Data filling operations*/
/*1*/
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
/*2*/
void Matrix::fill(double value){
    for(int i=0;i<this->rows;i++)
        for(int j=0;j<this->columns;j++)
            this->matrix[i][j]=value;
}
/*3*/
//void Matrix::randomize(){}


/*Boolean Methods*/
/*1*/
bool Matrix::isRowNull(int row){
    if (this->rows<row)
        return false;
    for(int j=0;j<this->columns;j++){
        if (this->matrix[row][j]!=0) return false;
    }
    return true;
}
/*2*/
bool Matrix::isColumnNull(int col){
    if (this->columns<col)
        return false;
    for(int i=0;i<this->rows;i++){
        if(this->matrix[i][col]!=0) return false;
    }
    return true;
}
/*3*/
bool Matrix::isNullMatrix(){
    for(int i=0;i<this->rows;i++){
        for (int j=0;j<this->columns;j++){
            if(this->matrix[i][j]!=0)   return false;
        }
    }
    return true;
}
/*4*/
bool Matrix::isIdMatrix(){
    if(this->rows!=this->columns)
        return false;
    for(int i=0;i<this->rows;i++){
        for(int j=0;j<this->columns;j++){
            if((i==j && this->matrix[i][j]!=1)||(i!=j && this->matrix[i][j]!=0))
                return false;
        }
    }
    return true;
}
/*5*/
bool Matrix::isSymetricMatrix(){
    if(this->rows!=this->columns){
        return false;
    }
    for(int i=0;i<this->rows;i++){
        for(int j=0;j<this->columns;j++){
            if(this->matrix[i][j]!=this->matrix[j][i])
                return false;
        }
    }
    return true;
}

/*6*/
bool Matrix::isRowVector(){
    return(this->rows==1);
}
/*7*/
bool Matrix::isColumnVector(){
    return(this->columns==1);
}
/*8*/
bool Matrix::isVector(){
    return(this->rows==1||this->columns==1);
}
/*9*/
bool Matrix::isDouble(){
    return(this->rows==1 && this->columns==1);
}
/*Swappers*/
/*1*/
void Matrix::swapRows(int row1,int row2){
    int tmp;
    for(int j=0;j<this->columns;j++){
        tmp=this->matrix[row1][j];
        this->matrix[row1][j]=this->matrix[row2][j];
        this->matrix[row2][j]=tmp;
    }
}
/*2*/
void Matrix::swapColumns(int col1,int col2){
    int tmp;
    for(int i=0;i<this->rows;i++){
        tmp=this->matrix[i][col1];
        this->matrix[i][col1]=this->matrix[i][col2];
        this->matrix[i][col2]=tmp;

    }
}



/*Operations*/
/*1*/
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
/*2*/
Matrix Matrix::operator*(const double scalar){
    Matrix retMatrix(*this);
    for(int i=0;i<this->rows;i++){
        for(int j=0;j<this->columns;j++){
            retMatrix.matrix[i][j]*=scalar;
        }
    }
    return retMatrix;
}
/*3*/
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
/*4*/
bool Matrix::operator==(const Matrix& Matrix_){
    if (this->rows!=Matrix_.rows || this->columns!=Matrix_.columns)
        return false;
    for(int i=0;i<this->rows;i++){
        for(int j=0;j<this->columns;j++){
            if (this->matrix[i][j]!=Matrix_.matrix[i][j])
                return false;
        }
    }
    return true;
}
/*internal methods*/
/*1*/
void Matrix::thisScalarMultiplication(const double scalar){
    
    for(int i=0;i<this->rows;i++){
        for(int j=0;j<this->columns;j++)
            this->matrix[i][j]=this->matrix[i][j]*scalar;
    }
    
}
/*2*/
bool Matrix::multiplyRowScalar(int row,double scalar){
    if (row>=this->rows)
        return false;
    for(int j=0;j<this->columns;j++){
        this->matrix[row][j]=this->matrix[row][j]*scalar;
    }
    return true;
}
/*3*/
bool Matrix::multiplyColScalar(int col,double scalar){
    if (col>=this->columns)
        return false;
    for(int i=0;i<this->rows;i++){
        this->matrix[i][col]=this->matrix[i][col]*scalar;
    }
    return true;
}
/*4*/
double Matrix::getMinMatrix(){
    double min=this->matrix[0][0];
    for(int i=0;i<this->rows;i++){
        for(int j=0;j<this->columns;j++){
            if(this->matrix[i][j]<min)
                min=this->matrix[i][j];
        }
    }
    return min;
}
/*5*/
double Matrix::getMaxMatrix(){
    double max=this->matrix[0][0];
    for(int i=0;i<this->rows;i++){
        for(int j=0;j<this->columns;j++){
            if(this->matrix[i][j]>max)
                max=this->matrix[i][j];
        }
    }
    return max;
}


/*Numerical methods*/
/*1*/
Matrix Matrix::getTranspose(){
    Matrix retMatrix(this->columns,this->rows);
    for(int i=0;i<this->rows;i++){
        for(int j=0;j<this->columns;j++){
            retMatrix.matrix[j][i]=this->matrix[i][j];
        }
    }
    return retMatrix;
}

/*2*/
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
            P.swapColumns(pivot,pivotrow);
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

/*3*/
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
/*4*/
Matrix Matrix::Choleskey(){
    std::vector<Matrix>matrixList;
    if(this->isSymetricMatrix()==false){
        std::cout<<"error couldn't extract choleskey:non symetrical matrix";
        exit(0);
    }
    matrixList=this->PLU();
    if(matrixList[0].isIdMatrix()==false){
        std::cout<<"error couldn't extract choleskey:P is not equal to ID";
        exit(0);
    }
    
    for(int i=0;i<this->rows;i++){
        if(matrixList[2].matrix[i][i]<0){
            std::cout<<"no garantee the matrix is positive";
            exit(0);
        }
        for(int j=0;j<this->columns;j++){
            matrixList[1].matrix[i][j]*=sqrt(matrixList[2].matrix[j][j]);
        }
    }
    return matrixList[1];
}
/*5*/
Matrix Matrix::getInverse(){
    if (this->rows!=this->columns)
    {
        std::cout<<"cannot get the inverse! it is not a square matrix";
        exit(0);
    }
    Matrix thisMatrix(*this);
    Matrix inverseMatrix(this->rows,this->columns);inverseMatrix.id();
    //double regulator=this->getMinMatrix();
    //thisMatrix.thisScalarMultiplication(1/regulator);
    int pivot;
    for(int i=0;i<this->rows;i++){
        pivot=i;
        //probleme de la diagonale nulle
        if(this->matrix[pivot][pivot]==0){
            int pivotrow=pivot+1;
            for(pivotrow;pivotrow<this->rows;pivotrow++){
                if (this->matrix[pivotrow][pivot]!=0) break;
            }
            if (pivotrow==this->rows)continue;
            thisMatrix.swapRows(pivot,pivotrow);
            inverseMatrix.swapRows(pivot,pivotrow);
        }
        inverseMatrix.multiplyRowScalar(pivot,1/thisMatrix.matrix[pivot][pivot]);
        thisMatrix.multiplyRowScalar(pivot,1/thisMatrix.matrix[pivot][pivot]);
        
        for(int j=pivot+1;j<this->rows;j++){
            
            double coeff=thisMatrix.matrix[j][pivot];
            for(int k=0;k<this->columns;k++){
                thisMatrix.matrix[j][k]+=-1*coeff*thisMatrix.matrix[pivot][k];
                inverseMatrix.matrix[j][k]+=-1*coeff*inverseMatrix.matrix[pivot][k];
            }
            if(thisMatrix.isRowNull(j)){
                std::cout<<"matrice non diagonalisable ! ligne:"<<j<<" est null";
                //return results that makes the determinant equals to 0
                return thisMatrix;
            }
        }
    }

    for(int i=this->rows-1;i>-1;i--){
        pivot=i;
        inverseMatrix.multiplyRowScalar(pivot,1/thisMatrix.matrix[pivot][pivot]);
        thisMatrix.multiplyRowScalar(pivot,1/thisMatrix.matrix[pivot][pivot]);
        for(int j=pivot-1;j>-1;j--){
            double coeff=thisMatrix.matrix[j][pivot];
            for(int k=this->columns-1;k>-1;k--){
                thisMatrix.matrix[j][k]+=-1*coeff*thisMatrix.matrix[pivot][k];
                inverseMatrix.matrix[j][k]+=-1*coeff*inverseMatrix.matrix[pivot][k];
            }
        }
    }
   /* if(thisMatrix.isIdMatrix()==false){
        std::cout<<"something wrong I can feel it";
        exit(0);
    }*/
    //inverseMatrix.thisScalarMultiplication(regulator);
    return inverseMatrix;
}
/*Vector transformation*/
Vector Matrix::vectorTransform(){
    if(!this->isVector())
    {
        std::cout<<"this is not a vector";
        exit(0);
    }
    if(this->isRowVector()){
        Vector v(this->columns);
        for(int i=0;i<this->columns;i++){
            v.vector[i]=this->matrix[0][i];
        }
    }
    if(this->isColumnVector()){
        Vector v(this->rows);
        for(int i=0;i<this->rows;i++){
            v.vector[i]=this->matrix[i][0];
        }
    }
    

}
/*Static function*/
Matrix Matrix::HilbertMatrix(int n){
    Matrix H(n,n);
    for(double i=0;i<H.rows;i++){
        for(double j=0;j<H.columns;j++){
            H.matrix[i][j]=1/(i+j+1);
        }
    }
    return H;
}




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
