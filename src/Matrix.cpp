#include <iostream>
#include <vector>
#include "../headers/MathUtils.hpp"
#include <math.h>
#include "../headers/MathUtils.hpp"
#include "../headers/Matrix.hpp"
/*Constructors*/

/*1*/
//create null matrix
Matrix::Matrix(int rows,int columns)
{
    try{  
        this->rows=rows;
        this->columns=columns;
        if(this->rows<1 || this->columns<1) throw MATRIX_SIZE_ERROR;
        
        for(int i=0;i<this->rows;i++){
            std::vector<double> vector;
            for(int j=0;j<this->columns;j++){
                vector.push_back(0);
            }
            this->matrix.push_back(vector);
        }
    }  
    catch(...){
        std::cout<<"size of the matrix < 0";
        exit(MATRIX_SIZE_ERROR);
    }

}

/*2*/
//create Matrix object and copy matrix to it
Matrix::Matrix(std::vector<std::vector<double>>& matrix){
    this->rows=matrix.size();
    this->columns=matrix.size();
    for(int i=0;i<this->rows;i++){
        std::vector<double> vector;
        for(int j=0;j<this->columns;j++){
            vector.push_back(matrix[i][j]);
        }
        this->matrix.push_back(vector);
    }
}

/*3*/
//create Matrix object and copy Matrix_ to it
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
/*4*/
//create square matrix with size "size" from a specified "row" and "col"
Matrix::Matrix(int row,int col,int size,const Matrix& Matrix_){
    if(row<0||row>=Matrix_.rows||col<0||col>Matrix_.columns)
        exit(0);
    if(size<0||row+size>Matrix_.rows||col+size>Matrix_.columns)
        exit(0);
    this->rows=size;
    this->columns=size;
    std::vector<double>vector;
    for(int i=0;i<size;i++){
        for(int j=0;j<size;j++){
            vector.push_back(Matrix_.matrix[row+i][col+j]);
        }
        this->matrix.push_back(vector);
    }
}

/*5*/
//create a matrix by copying a part of Matrix_
Matrix::Matrix(int beginRow,int beginCol,int endRow,int endCol,int size,const Matrix& Matrix_){

}


/*6*/
//square matrix
Matrix::Matrix(int size)
{
    try{  
        this->rows=size;
        this->columns=size;
        if(this->rows<1 || this->columns<1) throw MATRIX_SIZE_ERROR;
        
        for(int i=0;i<this->rows;i++){
            std::vector<double> vector;
            for(int j=0;j<this->columns;j++){
                vector.push_back(0);
            }
            this->matrix.push_back(vector);
        }
    }  
    catch(...){
        std::cout<<"size of the matrix < 0";
        exit(MATRIX_SIZE_ERROR);
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
            printf("%.3lf  ",this->matrix[i][j]);
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


/*Boolean Methodes*/
/*1*/
bool Matrix::isNullRow(int row){
    if (this->rows<row || row<0)
        return false;
    for(int j=0;j<this->columns;j++){
        if (this->matrix[row][j]!=0) return false;
    }
    return true;
}
/*2*/
bool Matrix::isNullColumn(int col){
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
/*10*/
bool Matrix::isTriangularSup(){
    for(int i=0;i<this->rows;i++){
        for(int j=0;j<this->columns;j++)
        {
            if(i>j && this->matrix[i][j]!=0)
                return false;
        }
        
    }
    return true;
}
/*11*/
bool Matrix::isTriangularInf(){
    for(int i=0;i<this->rows;i++){
        for(int j=0;j<this->columns;j++)
        {
            if(j>i && this->matrix[i][j]!=0)
                return false;
        }
    }
    return true;
}
/*Swappers*/
/*1*/
void Matrix::swapRows(int row1,int row2){
    int tmp;
    try{
        for(int j=0;j<this->columns;j++){
            tmp=this->matrix[row1][j];
            this->matrix[row1][j]=this->matrix[row2][j];
            this->matrix[row2][j]=tmp;
        }
    }
    catch(...){
        std::cout<<"one or both of the rows is out of boundary";
    }
}
/*2*/
void Matrix::swapColumns(int col1,int col2){
    int tmp;
    try{
        for(int i=0;i<this->rows;i++){
            tmp=this->matrix[i][col1];
            this->matrix[i][col1]=this->matrix[i][col2];
            this->matrix[i][col2]=tmp;

        }
    }
    catch(...){
        std::cout<<"one or both of the columns is out of boundary";
    }
}



/*Operations*/
/*1*/
Matrix Matrix::operator*(const Matrix& Matrix_){
    try{
        if(this->columns!=Matrix_.rows)throw MATRIX_MULTIPLICATION_ERROR;
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
    catch(...){
        std::cout<<"multiplication not possible !";
        exit(MATRIX_MULTIPLICATION_ERROR);
    }
}
/*2*/
Matrix Matrix::operator*(const double scalar){
    try{
        Matrix retMatrix(*this);
        for(int i=0;i<this->rows;i++){
            for(int j=0;j<this->columns;j++){
                retMatrix.matrix[i][j]*=scalar;
            }
        }
        return retMatrix;
    }
    catch(...){
        std::cout<<"error at 'operator*(const double scalar)'";
        exit(0);
    }
}
/*3*/
Matrix Matrix::operator+(const Matrix& Matrix_){
    try{
        if(this->rows!=Matrix_.rows||this->columns!=Matrix_.columns )throw MATRIX_ADDITION_ERROR;
        Matrix sumMatrix(this->rows,this->columns);
        for(int i=0;i<this->rows;i++){
            for(int j=0;j<this->columns;j++){
                sumMatrix.matrix[i][j]=this->matrix[i][j]+Matrix_.matrix[i][j];
            }
        }
    return sumMatrix;
    }
    catch(...)
    {
            std::cout<<"cannot add two matrixes of different size:Math Error";
            exit(MATRIX_ADDITION_ERROR);
    }
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
    try{
        for(int i=0;i<this->rows;i++){
            for(int j=0;j<this->columns;j++)
                this->matrix[i][j]=this->matrix[i][j]*scalar;
        }
    }
    catch(...){
        std::cout<<"ERROR at 'thisScalarMultiplication(const double scalar)'";
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

/*6*/
void Matrix::thisMatrixMultiplicationLeft(const Matrix& Matrix_){
    try{
        if(this->columns!=Matrix_.rows)throw MATRIX_MULTIPLICATION_ERROR;
        //must find a way to resize the matrix
        if(this->columns!=Matrix_.columns ||this->rows!=Matrix_.rows) throw MATRIX_MULTIPLICATION_ERROR;
        Matrix tmpMatrix(this->rows,Matrix_.columns);
        for(int i=0;i<tmpMatrix.rows;i++){
            for(int j=0;j<tmpMatrix.columns;j++){
                for(int k=0;k<this->columns;k++){
                    tmpMatrix.matrix[i][j]+=this->matrix[i][k]*Matrix_.matrix[k][j];
                }
            }
        }
        for(int i=0;i<this->rows;i++){
            for(int j=0;j<this->columns;j++){
                this->matrix[i][j]=tmpMatrix.matrix[i][j];
            }
        }
       
    }
    catch(...)
        {
            std::cout<<"Cannot multiply Matrix";
            exit(0);
    }
    
}

void Matrix::thisMatrixMultiplicationRight(const Matrix& Matrix_){
try{
        if(this->columns!=Matrix_.rows)throw MATRIX_MULTIPLICATION_ERROR;
        //must find a way to resize the matrix
        if(this->columns!=Matrix_.columns ||this->rows!=Matrix_.rows) throw MATRIX_MULTIPLICATION_ERROR;
        Matrix tmpMatrix(this->rows,Matrix_.columns);
        for(int i=0;i<tmpMatrix.rows;i++){
            for(int j=0;j<tmpMatrix.columns;j++){
                for(int k=0;k<this->columns;k++){
                    tmpMatrix.matrix[i][j]+=Matrix_.matrix[i][k]*this->matrix[k][j];
                }
            }
        }
        for(int i=0;i<this->rows;i++){
            for(int j=0;j<this->columns;j++){
                this->matrix[i][j]=tmpMatrix.matrix[i][j];
            }
        }
       
    }
    catch(...)
        {
            std::cout<<"Cannot multiply Matrix";
            exit(0);
    }
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
std::vector<Matrix> Matrix::plu(){
    try{
        if(this->rows!=this->columns)throw PLU_EXCEPTION;
        std::vector<Matrix> result;
        int pivot;
        Matrix U(*this);
        Matrix L(this->rows,this->columns);L.id();
        Matrix P(this->rows,this->columns);P.id();
        for(int i=0;i<this->rows;i++){
            
            //changing the pivot if it is null
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
            
            for(int j=pivot+1;j<this->rows;j++){
                L.matrix[j][pivot]=U.matrix[j][pivot]/U.matrix[pivot][pivot];
                for(int k=0;k<this->columns;k++){
                    U.matrix[j][k]+=-1*L.matrix[j][pivot]*U.matrix[pivot][k];
                    
                }
                if(U.isNullRow(j)){
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
    catch(...){
            std::cout<<"Non square matrix cannot apply plu correctly";
            exit(PLU_EXCEPTION);
        }
}

/*3*/
double Matrix::determinant(){
    if(this->rows!=this->columns){
        return 0;
    }
    std::vector<Matrix>pluList=this->plu();
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
    try{
        if(this->isSymetricMatrix()==false)throw CHOL_ERR1;
        std::vector<Matrix>matrixList;
        matrixList=this->plu();
        if(matrixList[0].isIdMatrix()==false)throw CHOL_ERR2;
        
        for(int i=0;i<this->rows;i++){
            if(matrixList[2].matrix[i][i]<0)throw CHOL_ERR3;
            for(int j=0;j<this->columns;j++){
                matrixList[1].matrix[i][j]*=sqrt(matrixList[2].matrix[j][j]);
            }
        }
        return matrixList[1];
    }

    catch(int e){
        if(e==CHOL_ERR1)
        {
            std::cout<<"error couldn't extract choleskey:non symetrical matrix";
            exit(CHOL_ERR1);
        }
        else if(e==CHOL_ERR2){
            std::cout<<"error couldn't extract choleskey:P is not equal to ID";
            exit(CHOL_ERR2);
        }
        else 
        {
            std::cout<<"no garantee the matrix is positive";
            exit(CHOL_ERR3);
        }
    }
}
/*5*/
Matrix Matrix::getInverse(){
    try{
        if (this->rows!=this->columns) throw MATRIX_INVERSE_EXCEPTION;
        Matrix thisMatrix(*this);
        Matrix inverseMatrix(this->rows,this->columns);inverseMatrix.id();
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
                if(thisMatrix.isNullRow(j)){
                    std::cout<<"matrice non diagonalisable ! ligne:"<<j<<" est null";
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
        return inverseMatrix;
    }
    catch(...){
        std::cout<<"cannot get the inverse! it is not a square matrix";
        exit(0);
    }
}

/*Extracting elements*/
/*1*/
Vector Matrix::convertVector(){
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
        return v;
    }
    else{//it is column vector
        Vector v(this->rows);
        for(int i=0;i<this->rows;i++){
            v.vector[i]=this->matrix[i][0];
        }
        return v;
    }
    
}
/*2*/
Vector Matrix::extractRowVector(int row){
    try{
        if(row>=this->rows)throw EXTR_ROW_VEC_EXCP;
        Vector v(this->columns);
        for(int j=0;j<v.length;j++)
            v.vector[j]=this->matrix[row][j];
        return v;
    }
    catch(...){
    
        std::cout<<"couldn't extract the row"<<std::endl;
        exit(EXTR_ROW_VEC_EXCP);
    }
}
/*3*/
Matrix Matrix::extractRowMatrix(int row){
    try{
        if(row>=this->rows)throw EXTR_ROW_MAT_EXCP; 
        Matrix retMat(1,this->columns);
        for(int j=0;j<retMat.columns;j++)
            retMat.matrix[0][j]=this->matrix[row][j];
        return retMat;
    }
    catch(...)
    {
        std::cout<<"couldn't extract the row"<<std::endl;
        exit(0);
    }
}
/*4*/
Vector Matrix::extractColumnVector(int column){
    try{
        if(column>=this->columns)throw EXTR_COL_VEC_EXCP;
        Vector v(this->rows);
        for(int i=0;i<v.length;i++)
            v.vector[i]=this->matrix[i][column];
        return v;
    }
    catch(...){
        std::cout<<"couldn't extract the column"<<std::endl;
        exit(EXTR_COL_VEC_EXCP);
    }
}
/*5*/
Matrix Matrix::extractColumnMatrix(int column){
    try{
        if(column>=this->columns)throw EXTR_COL_MAT_EXCP;
        Matrix retMat(this->rows,1);
        for(int i=0;i<retMat.rows;i++)
            retMat.matrix[i][0]=this->matrix[i][column];
        return retMat;
    }
    catch(...){
        std::cout<<"couldn't extract the column"<<std::endl;
        exit(0);
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


//---------------the QR project-----------------
double Matrix::getEuclidianNorm(){
    try{
        double euclNorm=0;
        if(this->isColumnVector()){
            for(int i=0;i<this->rows;i++){
                euclNorm+=this->matrix[i][0]*this->matrix[i][0];
            }
            
        }
        else if(this->isRowVector()){
            for(int j=0;j<this->columns;j++){
                euclNorm+=this->matrix[0][j];
            }
        }
        else throw EUCL_NORM_EXCEP;
        euclNorm=sqrt(euclNorm);
        return euclNorm;
    }
    catch(...){
        std::cout<<"cannot get euclidian norm for now";
        exit(EUCL_NORM_EXCEP);
    }
}
Matrix Matrix::extractMatrix(int rowBegin,int colBegin,int rowEnd,int colEnd){
    try{
        if(rowBegin>rowEnd||colBegin>colEnd) throw EXTR_MAT_EXCP;
            
        if(rowBegin<0||rowBegin>this->rows||rowEnd<0||rowEnd>this->rows)throw EXTR_MAT_EXCP;
    
        Matrix extractedMatrix(rowEnd-rowBegin+1,colEnd-colBegin+1);
        for(int i=rowBegin;i<=rowEnd;i++){
            for(int j=colBegin;j<=colEnd;j++){
                extractedMatrix.matrix[i-rowBegin][j-colBegin]=this->matrix[i][j];
            }
        }
        return extractedMatrix;
    }
    catch(...){
        exit (EXTR_MAT_EXCP);
    }
}
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------


Matrix Matrix::cannonicRowMat(int index,int dimension){
    try{
        if(index>=dimension||index<0)
            throw CAN_ROW_MAT_EXCP;
    Matrix v =Matrix(1,dimension);
    v.matrix[0][index]=1;
    return v;
    }
    catch(...){
        printf("creation of row");
        exit(CAN_ROW_MAT_EXCP);
    }
    
}

Matrix Matrix::cannonicColMat(int index,int dimension){
    try{
        if(index>=dimension||index<0)
        exit(0);
        Matrix v =Matrix(dimension,1);
        v.matrix[index][0]=1;
        return v;
    }
    catch(...){
        printf("error in cannonic column Mat");
        exit(CAN_ROW_MAT_EXCP);
    }
    
}
Matrix Matrix::getId(int dimension){
    Matrix retMatrix(dimension,dimension);
    for(int i=0;i<dimension;i++){
        retMatrix.matrix[i][i]=1;
    }
    return retMatrix;
}

void Matrix::placeMatrix(int rowBegin,int colBegin,const Matrix& matrix){
    for(int i=0;i<matrix.rows;i++){
        for(int j=0;j<matrix.columns;j++){
            this->matrix[i+rowBegin][j+colBegin]=matrix.matrix[i][j];
        }
    }
}
bool Matrix::isColinearWith(const Matrix& matrix){
    double scalar;
    bool scalarFound=false;
    if(this->rows!=matrix.rows || this->columns!=matrix.columns) exit(0);
    for(int i=0;i<matrix.rows ;i++){
        for(int j=0;j<matrix.columns ;j++){
            if (scalarFound==true){
                if(matrix.matrix[i][j]*scalar!=this->matrix[i][j])
                    return false;
            }
            else{
                if((matrix.matrix[i][j]!=0 && this->matrix[i][j]==0) ||(matrix.matrix[i][j]!=0 && this->matrix[i][j]==0)){
                    return false;
                }
                else if(matrix.matrix[i][j]!=0 && this->matrix[i][j]!=0){
                    scalarFound=true;
                    scalar=this->matrix[i][j]/matrix.matrix[i][j];
                }
            }
        }
    }
    return true;
}
std::vector<Matrix> Matrix::qr(){
    if(this->rows!=this->columns)
        exit(0);
    
    int size=this->rows;
    
    std::vector<Matrix>Hvector;
    Matrix Q=getId(size);
    Matrix A(*this);
    for(int i=0;i<size;i++){
        Matrix Hsize=getId(size);
        Matrix a=A.extractMatrix(i,i,size-1,i);
        Matrix e=cannonicColMat(0,size-i);
        e.thisScalarMultiplication(MathUtils::sign(a.matrix[0][0])*a.getEuclidianNorm());
        Matrix u=(a+e);
        u.thisScalarMultiplication(1/u.getEuclidianNorm());
        Matrix tmpMatrix=u*u.getTranspose();
        tmpMatrix.thisScalarMultiplication(-2);
        Matrix H=Matrix::getId(size-i)+tmpMatrix;
        Hsize.placeMatrix(i,i,H);
        A.thisMatrixMultiplicationRight(Hsize);
        Hvector.push_back(Hsize);
    }
    
   for(int i=0;i<size;i++){
        Q=Q*Hvector[i];
    }
    std::vector<Matrix> retVector;
    retVector.push_back(Q);
    retVector.push_back(A);
    return retVector;
}

Vector Matrix::solve(Vector y){
    try{
        
        if(this->columns!=y.length) throw SOLV_EQ_EXCP;
        Matrix yMat=y.convertToColMatrix();
        std::vector<Matrix>thisQR=this->qr();
        Matrix tQ=thisQR[0].getTranspose();
        Matrix y2=tQ*yMat;
        Vector solution(y.length);
        for(int i=solution.length-1;i>=0;i--){
            solution.vector[i]=y2.matrix[i][0];
            for(int j=thisQR[1].columns-1;j>i;j--){
                solution.vector[i]-=thisQR[1].matrix[i][j]*solution.vector[j];
            }
            solution.vector[i]*=1.0/thisQR[1].matrix[i][i];
        }
        return solution;
    }
    catch(...){
        std::cout<<"couldn't solve it";
        exit(SOLV_EQ_EXCP);
    }
    
}
bool Matrix::isSquare(){
    return this->rows==this->columns;
}

bool Matrix::isDiagonallyDominant(){
    if(!this->isSquare())
        return false;
    double diag_elt;
    double other_elt=0;
    for(int i=0;i<this->rows;i++){
        diag_elt=this->matrix[i][i];
        for(int j=0;j<this->rows;j++){
            if(i!=j){
                other_elt+=MathUtils::abs(this->matrix[i][j]);
            }
            if(diag_elt<other_elt)
                return false;
        }
    }
    return true;
}



Vector Matrix::solve_GS(Vector y,Vector x0,int it_max,float error){
    try{
        if(this->columns!=y.getLength())throw GS_EXCP;
        Vector xRet(y);
        
        double sum_temp;
        int nb_it=0;
        if(this->isDiagonallyDominant()){
            std::cout<<"iterates until convergence"<<std::endl;
            do
            {
                if(nb_it!=0)
                    x0=xRet;
                
                for(int i=0;i<x0.length;i++){
                    sum_temp=y.vector[i];
                    for(int j=0;j<i;j++){
                        sum_temp-=this->matrix[i][j]*xRet.vector[j];
                    }
                    for(int j=i+1;j<this->columns;j++){
                        sum_temp-=this->matrix[i][j]*x0.vector[j];
                        
                    }
                    xRet.vector[i]=(1.0/this->matrix[i][i])*sum_temp;
                    
                }   
                
                nb_it++;
            }while (MathUtils::abs((xRet-x0).getEuclidianNorm())>error);
            
            printf("iteration number is:%d\n",nb_it);
        }
        else{
            printf("iterates until it_max:\n");
            for(int it=0;it<it_max;it++){
                
                if(nb_it!=0)
                    x0=xRet;
                
                for(int i=0;i<x0.length;i++){
                    sum_temp=y.vector[i];
                    for(int j=0;j<i;j++){
                        sum_temp-=this->matrix[i][j]*xRet.vector[j];
                    }
                    for(int j=i+1;j<this->columns;j++){
                        sum_temp-=this->matrix[i][j]*x0.vector[j];
                        
                    }
                    xRet.vector[i]=(1.0/this->matrix[i][i])*sum_temp;
                    
                }   
                
                nb_it++; 
            }
        }
        return xRet;
    }
    catch(...){
        exit(GS_EXCP);
    }
    
}

double Matrix::power_iteration(Vector x0,Vector* eigenVector,int it_max,float error){
    double err=1;

    Matrix xMat=x0.convertToColMatrix();
    Matrix qMat(x0.getLength(),1);
    
    for(int i=0;i<it_max && err>error;i++){
        qMat=xMat*(1.0/xMat.getEuclidianNorm());
        xMat=(*this)*qMat;
        err=MathUtils::abs((xMat.convertVector()-qMat.convertVector()).getEuclidianNorm());
    }
    
    (*eigenVector)=qMat.convertVector();
    double eigenValue=(xMat.getTranspose()*qMat).getElement(0,0);
    return eigenValue;
}

std::vector< Matrix> Matrix::getDiagonalEquivalent_Gauss(){
    std::vector<Matrix>list;
    if(this->rows!=this->columns)exit(0);
    if(!this->isSymetricMatrix())exit(0);
    Matrix retMat(*this);
    Matrix qRotation(this->rows,this->columns);
    Matrix eigenVectors(this->rows,this->columns);
    int nb_iteration=0;
    qRotation.id();
    eigenVectors.id();
    int index_i,index_j;
    double max;
    double epsilon=0.0001;
    do{
        index_i=0;
        index_j=0;
        max=0;
        for(int i=0;i<this->rows;i++){
            for(int j=i+1;j<this->columns;j++){
                if(MathUtils::abs(this->matrix[i][j])>max && MathUtils::abs(MathUtils::abs(retMat.matrix[i][j])-max)>epsilon){
                    max=MathUtils::abs(this->matrix[i][j]);
                    
                    index_i=i;
                    index_j=j;
                }
            }
        }
        if(index_i!=0 || index_j!=0){
            double beta=(retMat.matrix[index_i][index_i]-retMat.matrix[index_j][index_j])/(2*retMat.matrix[index_i][index_j]);
            double t,cos_teta,sin_teta;
            if(beta<0){
                t=-1*beta-sqrt(beta*beta+1);
            }
            else{
                t=-1*beta+sqrt(beta*beta+1);
            }
            cos_teta=1.0/sqrt(1+t*t);
            sin_teta=t/sqrt(1+t*t);
            qRotation.setElement(index_i,index_i,cos_teta);
            qRotation.setElement(index_i,index_j,sin_teta);
            qRotation.setElement(index_j,index_i,-1*sin_teta);
            qRotation.setElement(index_j,index_j,cos_teta);
            retMat=qRotation.getTranspose()*retMat*qRotation;
            eigenVectors=eigenVectors*qRotation;
            qRotation.id();
        }
        nb_iteration++;
    }while(index_i!=0 || index_j!=0);
    list.push_back(retMat);
    list.push_back(eigenVectors);
    return list;
}