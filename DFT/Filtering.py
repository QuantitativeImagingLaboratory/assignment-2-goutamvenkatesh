# For this part of the assignment, You can use inbuilt functions to compute the fourier transform
# You are welcome to use fft that are available in numpy and opencv

import numpy as np
import math


class Filtering:
    image = None
    filter = None
    cutoff = None
    order = None
    
    

    def __init__(self, image, filter_name, cutoff, order = 0):
        """initializes the variables frequency filtering on an input image
        takes as input:
        image: the input image
        filter_name: the name of the mask to use
        cutoff: the cutoff frequency of the filter
        order: the order of the filter (only for butterworth
        returns"""
        self.image = image
        if filter_name == 'ideal_l':
            self.filter = self.get_ideal_low_pass_filter
        elif filter_name == 'ideal_h':
            self.filter = self.get_ideal_high_pass_filter
        elif filter_name == 'butterworth_l':
            self.filter = self.get_butterworth_low_pass_filter
        elif filter_name == 'butterworth_h':
            self.filter = self.get_butterworth_high_pass_filter
        elif filter_name == 'gaussian_l':
            self.filter = self.get_gaussian_low_pass_filter
        elif filter_name == 'gaussian_h':
            self.filter = self.get_gaussian_high_pass_filter

        self.cutoff = cutoff
        self.order = order


    def get_ideal_low_pass_filter(self, shape, cutoff):
        """Computes a Ideal low pass mask
        takes as input:
        shape: the shape of the mask to be generated
        cutoff: the cutoff frequency of the ideal filter
        returns a ideal low pass mask"""
        
        rows, cols = shape
        mask = np.zeros((rows, cols),np.uint8)
        for i in range(rows):
            for j in range(cols):
                value = math.sqrt((i - (rows / 2)) ** 2 + (j - (cols / 2)) ** 2)
                if (value <= cutoff):
                    mask[i, j] = 1

                else:
                    mask[i, j] = 0
                    
        return mask


    def get_ideal_high_pass_filter(self, shape, cutoff):
        """Computes a Ideal high pass mask
        takes as input:
        shape: the shape of the mask to be generated
        cutoff: the cutoff frequency of the ideal filter
        returns a ideal high pass mask"""

        #Hint: May be one can use the low pass filter function to get a high pass mask
        
        rows, cols = shape
        mask = np.zeros((rows, cols), np.uint8)
        for i in range(rows):
            for j in range(cols):
                value = math.sqrt((i - (rows / 2)) ** 2 + (j - (cols / 2)) ** 2)
                if (value <= cutoff):
                    mask[i, j] = 0

                else:
                    mask[i, j] = 1

        
        return mask



    def get_butterworth_low_pass_filter(self, shape, cutoff, order):
        """Computes a butterworth low pass mask
        takes as input:
        shape: the shape of the mask to be generated
        cutoff: the cutoff frequency of the butterworth filter
        order: the order of the butterworth filter
        returns a butterworth low pass mask"""
        
        def dist( x, y):
            return np.sqrt(np.sum(((x[0] - y[0]) ** 2) + ((x[1] - y[1]) ** 2)))

        mask=np.zeros(shape)
        center=[shape[0]/2,shape[1]/2]
        for i in range(shape[0]):
            for j in range(shape[1]):

                mask[i,j]=1/(1+np.power((dist([i,j],center)/cutoff),2*order))


        return mask

    def get_butterworth_high_pass_filter(self, shape, cutoff, order):
        """Computes a butterworth high pass mask
        takes as input:
        shape: the shape of the mask to be generated
        cutoff: the cutoff frequency of the butterworth filter
        order: the order of the butterworth filter
        returns a butterworth high pass mask"""

        #Hint: May be one can use the low pass filter function to get a high pass mask
        
        def dist( x, y):
            return np.sqrt(np.sum(((x[0] - y[0]) ** 2) + ((x[1] - y[1]) ** 2)))

        mask=np.zeros(shape)
        center=[shape[0]/2,shape[1]/2]
        for i in range(shape[0]):
            for j in range(shape[1]):

                mask[i,j]=1/(1+(np.power(cutoff/(dist([i,j],center)),2*order)))


        return mask


    def get_gaussian_low_pass_filter(self, shape, cutoff):
        """Computes a gaussian low pass mask
        takes as input:
        shape: the shape of the mask to be generated
        cutoff: the cutoff frequency of the gaussian filter (sigma)
        returns a gaussian low pass mask"""
        
        def dist( x, y):
            return np.sqrt(np.sum(((x[0] - y[0]) ** 2) + ((x[1] - y[1]) ** 2)))

        
        mask = np.zeros(shape)
        center = [shape[0] / 2, shape[1] / 2]
        for i in range(shape[0]):
            for j in range(shape[1]):
                mask[i, j] = np.exp(-(dist([i, j], center) ** 2) / (2 * (cutoff ** 2)))

        return mask
        

    def get_gaussian_high_pass_filter(self, shape, cutoff):
        """Computes a gaussian high pass mask
        takes as input:
        shape: the shape of the mask to be generated
        cutoff: the cutoff frequency of the gaussian filter (sigma)
        returns a gaussian high pass mask"""

        #Hint: May be one can use the low pass filter function to get a high pass mask

        def dist( x, y):
            return np.sqrt(np.sum(((x[0] - y[0]) ** 2) + ((x[1] - y[1]) ** 2)))
        print("this is checckpoint 1 ")
        mask = np.zeros(shape)
        mask1=self.get_gaussian_low_pass_filter(shape,cutoff)
        center = [shape[0] / 2, shape[1] / 2]
        for i in range(shape[0]):
            for j in range(shape[1]):

                mask[i, j] = 1-np.exp((-(dist([i, j], center)) ** 2) / (2 * (cutoff ** 2)))


        return mask

    def post_process_image(self, image):
        """Post process the image to create a full contrast stretch of the image
        takes as input:
        image: the image obtained from the inverse fourier transform
        return an image with full contrast stretch
        -----------------------------------------------------
        1. Full contrast stretch (fsimage)
        2. take negative (255 - fsimage)
        """
        
        #Written this as part of the filtering() method

        return image


    def filtering(self):
        """Performs frequency filtering on an input image
        returns a filtered image, magnitude of DFT, magnitude of filtered DFT        
        ----------------------------------------------------------
        You are allowed to used inbuilt functions to compute fft
        There are packages available in numpy as well as in opencv
        Steps:
        1. Compute the fft of the image
        2. shift the fft to center the low frequencies
        3. get the mask (write your code in functions provided above) the functions can be called by self.filter(shape, cutoff, order)
        4. filter the image frequency based on the mask (Convolution theorem)
        5. compute the inverse shift
        6. compute the inverse fourier transform
        7. compute the magnitude
        8. You will need to do a full contrast stretch on the magnitude and depending on the algorithm you may also need to
        take negative of the image to be able to view it (use post_process_image to write this code)
        Note: You do not have to do zero padding as discussed in class, the inbuilt functions takes care of that
        filtered image, magnitude of DFT, magnitude of filtered DFT: Make sure all images being returned have grey scale full contrast stretch and dtype=uint8 
        """
        
        dft_image = np.fft.fft2(self.image)
        image_shift = np.fft.fftshift(dft_image)
        mask = self.filter(np.shape(image_shift), self.cutoff) # there is a problem here. check once
        filtered_image = mask * image_shift
        inverse_shift = np.fft.ifftshift(filtered_image)
        filtered_image2 = np.fft.ifft2(inverse_shift)
        ImgRows, ImgCols = filtered_image2.shape

        outMat = np.zeros((ImgRows, ImgCols), np.uint8)
        finalImg = np.zeros((ImgRows, ImgCols), np.uint8)

        for u in range(0, ImgRows):
            for v in range(0, ImgCols):
                outMat[u, v] = np.absolute(filtered_image2[u, v])

        B = np.amax(outMat)
        A = np.amin(outMat)

        difInt = B - A
        for u in range(0, ImgRows):
            for v in range(0, ImgCols):
                firstMul = 255 / difInt
                finalVal = firstMul * (outMat[u, v] - A)
                finalImg[u, v] = int(np.round(finalVal))

        for i in range(ImgRows):
            for j in range(ImgCols):
                if(mask[i,j]!=0):
                    mask[i,j]=255





        return [mask, finalImg,finalImg]
    
        
     
