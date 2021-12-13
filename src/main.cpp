#include <iostream>
#include "../headers/Includes.hpp"

using namespace std;
  
  
int main()
{   
   std::vector<std::vector<double>>matrix;
   std::vector<double> v;
   v.push_back(2);
   v.push_back(-2);
   v.push_back(18);
   matrix.push_back(v);
   v[0]=2;
   v[1]=1;
   v[2]=0;
   matrix.push_back(v);
   v[0]=1;
   v[1]=2;
   v[2]=0;
   matrix.push_back(v);
   Matrix A(3,3,matrix);
   std::vector<Matrix> listMat;
   listMat=A.QR();
   listMat[0].printMatrix();
   listMat[1].printMatrix();
}
/*
Matrix A(3,3);
    Matrix B(2,3);
    A=B;
    A.printMatrix();
    Matrix C=A*B;
    C=B*A;
    */
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