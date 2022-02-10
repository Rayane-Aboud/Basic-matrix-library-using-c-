#ifndef MATHUTILS_HPP
#define MATHUTILS_HPP

#define MATRIX_SIZE_ERROR 1
#define MATRIX_MULTIPLICATION_ERROR 2
#define MATRIX_ADDITION_ERROR 3
#define PLU_EXCEPTION 4
#define CHOL_ERR1 5
#define CHOL_ERR2 6
#define CHOL_ERR3 7
#define MATRIX_INVERSE_EXCEPTION 8
#define EXTR_ROW_VEC_EXCP 9
#define EXTR_ROW_MAT_EXCP 10
#define EXTR_COL_VEC_EXCP 11
#define EXTR_COL_MAT_EXCP 12
#define EUCL_NORM_EXCEP 13
#define EXTR_MAT_EXCP 14
#define CAN_ROW_MAT_EXCP 15
#define CAN_COL_MAT_EXCP 16
#define SOLV_EQ_EXCP 17
#define GS_EXCP 18

namespace MathUtils{
int sign(double d);
double abs(double d);
}


#endif