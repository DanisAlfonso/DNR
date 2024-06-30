import numpy as np
from scipy.linalg import lu, lu_solve, lu_factor
import time

# Function to generate a random matrix and vector
def generate_random_system(size):
    A = np.random.rand(size, size)
    b = np.random.rand(size)
    return A, b

# Save the matrix and vector to files
def save_to_files(A, b, matrix_file, vector_file):
    np.savetxt(matrix_file, A)
    np.savetxt(vector_file, b)

# Main function to generate and save the system
def main():
    # Define the size of the matrix
    matrix_size = 1000  # You can change this to test larger sizes
    
    # Generate a random system
    A, b = generate_random_system(matrix_size)
    
    # Save the matrix and vector to files
    save_to_files(A, b, 'matrix.txt', 'vector.txt')
    
    print(f"Generated {matrix_size}x{matrix_size} matrix and vector")

if __name__ == "__main__":
    main()

