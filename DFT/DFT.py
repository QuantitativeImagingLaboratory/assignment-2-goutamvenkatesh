# For this part of the assignment, please implement your own code for all computations,
# Do not use inbuilt functions like fft from either numpy, opencv or other libraries
import numpy as np
import math

class DFT:

    def forward_transform(self, matrix):
        """Computes the forward Fourier transform of the input matrix
        takes as input:
        matrix: a 2d matrix
        returns a complex matrix representing fourier transform"""
        
        
        output=np.zeros((15,15), dtype=complex)
        
        
        

        for u in range(0,15):
            for v in range(0,15):
                total_sum = 0
                for i in range(0,15):
                    for j in range(0,15):
                        t=u*i
                        k=v*j
                        y=t + k
                        x=2.0 * math.pi * y
                        angle = x/ 15
                        expcal=math.cos(angle) - 1j *math.sin(angle)
                        total_sum += matrix[i,j] * expcal
                    

            
                output[u,v]=total_sum

            

        return output




        

    def inverse_transform(self, matrix):
        """Computes the inverse Fourier transform of the input matrix
        matrix: a 2d matrix (DFT) usually complex
        takes as input:
        returns a complex matrix representing the inverse fourier transform"""



        return matrix


    def discrete_cosine_tranform(self, matrix):
        """Computes the discrete cosine transform of the input matrix
        takes as input:
        matrix: a 2d matrix
        returns a matrix representing discrete cosine transform"""



        return matrix


    def magnitude(self, matrix):
        """Computes the magnitude of the DFT
        takes as input:
        matrix: a 2d matrix
        returns a matrix representing magnitude of the dft"""

        return matrix
