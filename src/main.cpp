#include <iostream>
#include "../headers/Includes.hpp"

using namespace std;
  
  
int main()
{   
    std::vector<Matrix>L;
    std::vector<std::vector<double>>matrix;
    for(int i=0;i<3;i++){
        std::vector<double> vector;
        for(int j=0;j<3;j++){
            vector.push_back(0);
        }
        matrix.push_back(vector);
    }
    matrix[0][0]=0;
    matrix[0][1]=1;
    matrix[0][2]=1;
    matrix[1][0]=1;
    matrix[1][1]=5;
    matrix[1][2]=5;
    matrix[2][0]=1;
    matrix[2][1]=5;
    matrix[2][2]=14;

    Matrix A(3,3,matrix);
    //A.printMatrix();
    /*L=A.PLU();
    L[0].printMatrix();
    L[1].printMatrix();
    L[2].printMatrix();
    
    Matrix C(3,3);
    C=A.Choleskey();
    C.printMatrix();
    (C*C.getTranspose()).printMatrix();*/
    //Matrix B=Matrix::HilbertMatrix(100);
    //B.printMatrix();
    (Matrix::HilbertMatrix(3)*Matrix::HilbertMatrix(3).getInverse()).printMatrix();/*non stable à partir */
    //std::cout<<A.determinant();
    //std::cout<<A.determinant();
}
/*// Initializing the vector of vectors
    vector<vector<int> > vec;
  
    // Elements to insert in column
    int num = 10;
  
    // Inserting elements into vector
    for (int i = 0; i < ROW; i++) {
        // Vector to store column elements
        vector<int> v1;
  
        for (int j = 0; j < COL; j++) {
            v1.push_back(num);
            num += 5;
        }
  
        // Pushing back above 1D vector
        // to create the 2D vector
        vec.push_back(v1);
    }
  
    // Displaying the 2D vector
    for (int i = 0; i < vec.size(); i++) {
        for (int j = 0; j < vec[i].size(); j++)
            cout << vec[i][j] << " ";
        cout << endl;
    }
    return 0;*/
    
  /*  A.getNumRows();
    A.getNumCols();
    A.getElement(1,1);
    A.setElement(1,1,1);
    A.fill(5);
    Matrix B(2,3);
    B.fill(7);
    Matrix C=A*B;
    C.printMatrix();*/