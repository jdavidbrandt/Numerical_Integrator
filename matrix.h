#ifndef MATRIX_H
#define MATRIX_H
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <iostream>

//signum function (A.K.A. Sign function)  Returns 1, -1, or 0
template <typename T>
int sgn(T val){
    return (T(0) < val) - (val < T(0));
}

bool bIsEven(int n){
    if(n % 2 == 0){
        return true;
    }
    else{
        return false;
    }
}


bool bEpsilonIsEvenPermutation(int N, std::vector<int> eps_indices){
    //create a temporary vector, and sort it to look for duplicates
    std::vector<int> temp_eps_indices = eps_indices;
    std::sort(temp_eps_indices.begin(), temp_eps_indices.end());
    bool bHasDuplicates = std::adjacent_find(temp_eps_indices.begin(), temp_eps_indices.end()) != temp_eps_indices.end();
    //if it has duplicate indices, return false
    if(bHasDuplicates)
    {
        return false;
    }
    //initialize final sign for product equation, then calculate PI(0<= i < j <= N) * sgn(a_j - a_i), where a_j and a_i are the COLUMN (j) and ROW (i) of epsilon indices that refer to levi-civita
    int final_sign = 1;
        for(int j = 1; j <= N-1; ++j){
            for(int i = 0; i < N; ++i){
                if(i >= j){
                    continue;
                }
                else{
                    //perform calc only if i < j
                    final_sign*=sgn(eps_indices[j] - eps_indices[i]);
                } 
            }
        }
        if(final_sign == 1){   
            //the order of indices is either cyclic or even
            return true;
        }
        else{
            //the order of indices is either anticyclic or odd
            return false;
        }
}

bool bEpsilonIsOddPermutation(int N, std::vector<int> eps_indices){
    //create a temporary vector, and sort it to search for duplicates
    std::vector<int> temp_eps_indices = eps_indices;
    std::sort(temp_eps_indices.begin(), temp_eps_indices.end());
    bool bHasDuplicateIndices = std::adjacent_find(temp_eps_indices.begin(), temp_eps_indices.end()) != temp_eps_indices.end();
    //if it has duplicates, return false
    if(bHasDuplicateIndices){
        return false;
    }
    //initialize final sign for product equation, then calculate PI(0<= i < j <= N) * sgn(a_j - a_i), where a_j and a_i are the COLUMN (j) and ROW (i) of epsilon indices that refer to levi-civita
    int final_sign = 1;
    for(int j = 1; j <= N-1; ++j){
        for(int i = 0; i < N; ++i){
            if(i >= j){
                continue;
            }
            else{
                //perform calc only if i < j
                final_sign*=sgn(eps_indices[j] - eps_indices[i]);
            }
        }
    }
    if(final_sign == -1){
        //the order of indices is either anticyclic or odd
        return true;
    }
    else{
        //the order of indices is either cyclic or even
        return false;
    }
}

int LeviCivita(int n, std::vector<int> epsilon_indices){
    double epsilon;
    if(bEpsilonIsEvenPermutation(n, epsilon_indices)){
        epsilon = 1;
    }
    else if(bEpsilonIsOddPermutation(n, epsilon_indices)){
        epsilon = -1;
    }
    else{
        epsilon = 0;
    }
    return epsilon;
}

class Matrix
{
public:
    //attributes
    std::vector<std::vector<double>> elements;
    int rows;
    int cols;
    double determinant;
    double trace = 0;
    bool bIsHermitian;
    bool bIsSymmetric;
    bool bIsSquare;
    
    
    
    //constructor that requires elements as parameters
    Matrix(std::vector<std::vector<double>> passed_elements)
    {
        //there is only one initial property about a matrix, and it's that it's a row and col of numbers
        elements = passed_elements;
        rows = passed_elements.size();
        cols = passed_elements[0].size();
        AssignIsSquare();
    }

    //creates null matrix if bIsIdentity = false,
    //creates identity matrix if bIsIdentity = true;
    Matrix(int passed_rows, int passed_cols, bool bIsIdentity){
        elements.resize(passed_rows);
        for(int j = 0; j < passed_rows; ++j){
            elements[j].resize(passed_cols, 0);
        }
        rows = passed_rows;
        cols = passed_cols;
        if(bIsIdentity){
            //reassign diagonal values to 1 if bIsIdentity is true
            for(int it = 0; it < rows; ++it){
                elements[it][it] = 1;
            }
        }
        AssignIsSquare();
    }
    

    double GetDeterminant()
    {
        //define indices, i being rows, j being columns, and dimensions
        determinant = 0;
        int n = elements.size();
        int initial_depth = 0;
        std::vector<int> subindices(n);
        DescendDimensions(initial_depth, n, subindices);
        double result = determinant;
        return result;
    }

    double GetTrace(){
        if(bIsSquare){
        int j = 0;
        //sum diagonal element one-dimensionally
            for(int i = 0; i < rows; ++i){
                trace += elements[i][j];
                ++j;
            }
        return trace;
        }
        else{
            //if not a square matrix, assign trace as -1 for now
            std::cout<<"Not a square matrix"<<std::endl;
            trace = -1;
            return trace;
        }
        
    }


    void AssignIsSquare(){
        //will assign bIsSquare attribute AND return true or false
        bool result;
        if(rows == cols){
            bIsSquare = true;
        }
        else{
            bIsSquare = false;
        }    
    }

    void process_calc(const std::vector<int>& indices){
        
        //create temp vector to sort, so I can throw out duplicates, saving computation time
        std::vector<int> temp_indices;
        std::vector<int> epsilon_indices = indices;
        int n = elements[0].size();

        //initiate and keep track of multiplier, and current sum for each loop
        double multiplier = 1; 
        for (int index : indices){
            temp_indices.push_back(index);
        }
        int letter = 0;
        // iterate through constants 0, 1, 2, ... n for einstein notation equation
        for(int einstein_eqn_const=0; einstein_eqn_const < n; ++einstein_eqn_const){
            multiplier *= elements[einstein_eqn_const][epsilon_indices[letter]];
            letter+=1;
        }
        //get Levi_civita for each epsilon indices permutation
        int epsilon = LeviCivita(n, epsilon_indices);

        //add epsilon*multiplier to determinant (many times because of recursion)
        determinant += epsilon*multiplier;
        return;
}

//Use recursion to create nested loops
void DescendDimensions(int current_depth, int max_depth, std::vector<int>& current_indices){
    //if current_depth == max _depth, begin calculations
    if(current_depth == max_depth){
        process_calc(current_indices);
        return;
    }
    //if current_depth != max_depth, descend further into loops
    for (int i = 0; i<max_depth; ++i){
        current_indices[current_depth] = i;
        DescendDimensions(current_depth +1, max_depth,  current_indices);
    }
    
}

};

void PrintMatrix(Matrix matrix){
    std::cout<<std::endl;
    for(int i = 0; i < matrix.rows; ++i){
        for(int j = 0; j < matrix.cols; ++j){
            if(j == matrix.cols - 1){
                std::cout<<matrix.elements[i][j];
            }
            else{
                std::cout<<matrix.elements[i][j]<<", ";
            }
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

Matrix AddMatrices(Matrix matrix2, Matrix matrix1){
    //create temp null matrix
    Matrix temp_matrix = Matrix(matrix1.rows, matrix2.cols, false);
    if(matrix2.rows==matrix1.rows && matrix2.cols==matrix1.cols){
        for(int i = 0; i < matrix1.rows; ++i){
            for(int j = 0; j < matrix1.cols; ++j){
                //add elements of the same index
                temp_matrix.elements[i][j] += matrix2.elements[i][j] + matrix1.elements[i][j];
            }
        }
    Matrix result = temp_matrix;
    return result;
    }
    else{
        throw "Matrices not the same dimensions!";
    }
}

Matrix ScalarMultiply(double scalar, Matrix matrix){
    //iterate through each element, and multiply by scalar

    //create and return a copy
    Matrix result = matrix;
    for(int i = 0; i < matrix.rows; ++i){
        for(int j = 0; j < matrix.cols; ++j){
            result.elements[i][j] *= scalar;
        }
    }
    return result;
}


//The subtraction order is matrix2 - matrix1
Matrix SubtractMatrices(Matrix matrix2, Matrix matrix1){

    Matrix temp_matrix = Matrix(matrix1.rows, matrix2.cols, false);
    if(matrix2.rows==matrix1.rows && matrix2.cols==matrix1.cols){
        for(int i = 0; i < matrix1.rows; ++i){
            for(int j = 0; j < matrix1.cols; ++j){
                //subtract elements of the same index
                temp_matrix.elements[i][j] += matrix2.elements[i][j] - matrix1.elements[i][j];
            }
        }
    Matrix result = temp_matrix;
    return result;
    }
    else{
        throw "Matrices not the same dimensions";
    }
    
}

//vector2 dot vector1 but dot product is commutative so it doesn't matter
double DotProduct(std::vector<double> vector1, std::vector<double> vector2){
    double dot_product = 0;
    //row dot column
    if(vector1.size() == vector2.size()) {
        //define when its row * column vector or column * row vector
        for(int i = 0; i < vector1.size(); ++i){  
            dot_product += vector1[i]*vector2[i]; 
            }
        }

    else{
        throw "Incorrect dimensions";
    }
    return dot_product;
}

    //matrix1 cross matrix 2 in 3 dimensions
std::vector<double> CrossProduct(std::vector<double> vector1, std::vector<double> vector2){
    std::vector<double> result = {0,0,0};
    int dim = 3;

   //from einstein notation equation:   u cross v =  eps_{ijk}*u_{j}*v_{k}e_{i}
    for(int i = 0; i < vector1.size(); ++i){
        for(int j = 0; j < vector1.size(); ++j){
            for(int k = 0; k < vector1.size(); ++k){
                std::vector<int> epsilon_indices = {i,j,k};
                int epsilon = LeviCivita(dim, epsilon_indices);
                result[i] += epsilon*vector1[j]*vector2[k];
            }
        }
    }
    return result;

}

//Matrix * Vector
std::vector<double> MatrixVectorProduct(Matrix matrix, std::vector<double> vector){
    if(matrix.cols == vector.size()){

        //creates null vector to store results
        std::vector<double> result(matrix.cols, 0);

        //Using equation: u_{i} = A_{ij}*v_{j}
        for(int i = 0; i < matrix.rows; ++i){
            for(int j = 0; j < matrix.cols; ++j){
                result[i] += matrix.elements[i][j]*vector[j];
            }
        }
        return result;
    }
    else{
        throw "Incorrect dimensions!";
    }
}

//matrix1 times matrix2
Matrix MatrixProduct(Matrix matrix1, Matrix matrix2){
    if(matrix1.cols == matrix2.rows){
        //creates null matrix to store results
        Matrix result = Matrix(matrix1.rows, matrix2.cols, false);
        //Using equation: C{ik} = A_{ij}*B_{ik}
        for(int i = 0; i < matrix1.cols; ++i){
            for(int j = 0; j < matrix2.cols; ++j){
                for(int k = 0; k < matrix1.cols; ++k){
                    result.elements[i][k] += matrix1.elements[i][j]*matrix2.elements[j][k];
                }

            }
        }
        return result;
    }
    else{
        throw "Incorrect dimensions!";
    }
}

//matrix to the nth power
Matrix RaiseMatrixToPower(Matrix matrix, int n){
    Matrix result = matrix;
    int temp_n;
    if(bIsEven(n)){
        temp_n = n/2;
    }
    else{
        temp_n = (n+1)/2;
    }
    while(temp_n > 0){
        if((temp_n == 1) && (n % 2 == 1)){
            //if odd, take the result of the matrix to an even power, and multiply it by matrix 1 more time and return
            result = MatrixProduct(result, matrix);
            --temp_n;
            return result;
        }
        //group by powers of 2 to help computation time
        result = MatrixProduct(result, result);
        --temp_n;
    }
    if(n == 1){
        return matrix;
    }
    else{
        return result;
    }
}

Matrix Transpose(Matrix matrix){
    //create null matrix, store values, return
    Matrix transpose = Matrix(matrix.cols, matrix.rows, false);
    for(int i = 0; i < transpose.rows; ++i){
        for(int j = 0; j < transpose.cols; ++j){
            transpose.elements[i][j] = matrix.elements[j][i];
        }
    }
    return transpose;
}


#endif // MATRIX_H