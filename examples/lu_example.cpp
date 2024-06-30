#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include "math_utilities.h"
#include "lu_decomposition.h"

#ifndef TEST_DIR
#error "TEST_DIR is not defined"
#endif

// Function to read a matrix from a file
DNR::MatrixDouble read_matrix(const std::string& filename, int size) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    DNR::MatrixDouble matrix(size, size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (!(file >> matrix[i][j])) {
                throw std::runtime_error("Error reading matrix element");
            }
        }
    }
    return matrix;
}

// Function to read a vector from a file
DNR::VectorDouble read_vector(const std::string& filename, int size) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    DNR::VectorDouble vector(size);
    for (int i = 0; i < size; ++i) {
        if (!(file >> vector[i])) {
            throw std::runtime_error("Error reading vector element");
        }
    }
    return vector;
}

int main() {
    using namespace DNR;

    try {
        // Define the size of the matrix
        int matrix_size = 1000;  // Make sure this matches the size used in Python

        // Read the matrix and vector from files
        std::string matrix_file = std::string(TEST_DIR) + "/matrix.txt";
        std::string vector_file = std::string(TEST_DIR) + "/vector.txt";
        MatrixDouble A = read_matrix(matrix_file, matrix_size);
        VectorDouble b = read_vector(vector_file, matrix_size);
        VectorDouble x(matrix_size);

        // Debug: Print first few elements to verify correct loading
        std::cout << "Matrix A (first 3x3):" << std::endl;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                std::cout << A[i][j] << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "Vector b (first 3 elements):" << std::endl;
        for (int i = 0; i < 3; ++i) {
            std::cout << b[i] << " ";
        }
        std::cout << std::endl;

        // Measure the time taken to perform LU decomposition and solve the system
        auto start = std::chrono::high_resolution_clock::now();

        LU alu(A);
        alu.solve(b, x);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Time taken to solve a " << matrix_size << "x" << matrix_size << " system: " << elapsed.count() << " seconds" << std::endl;

    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}

