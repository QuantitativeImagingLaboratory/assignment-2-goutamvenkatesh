# Report

1. DFT-

  1. Forward Transform -
      
      This was a straight-forward application of the formula explained in the class. 
      First, we define a matrix full of zeros, in a grid of 15 rows and 15 columns. Each of these cells are filled with random 
      numbers as defined in the main function. Each of this value is taken and multiplied with the negative polar form of the 
      complex number that is represented by the angle that is calculated in the method.
      
      This is summed over first under two inner loops and is summed again, that is represented by the outer loops, that sums it over 
      every cell twice. This represents the double summation.
      
 2. Inverse transform - 
    
    The same method as forward transform is followed, except that instead of multiplying the matrix's randomly generated number,
    we use the DFT result, which is the grid of complex numbers and then multiply that with the polar form of the complex number
    that is represented by the angle that is calculated in the method.
    Again, this is summed over all the values twice and the final result is returned.
    
3. Cosine transform -
   
   Similar to the forward transform, we perform the forward transform, but only with the real part of the polar form, which is 
   represented by the cosine part (the real part) of the whole imaginary representation.
   
   
4. Magnitude - 
   
   here, we have looped over every cell in the input matrix and taken the square values of the cosine and sine part, and added it first.
   Then, we compute the square root pf this value, which gives us the magnitude of the polar representation. This is usually used to 
   represent the image as a whole.
   
  
  
  
2. Filtering-

 a. Ideal Low pass filter: The idea here is to set a cutoff frequency such that if the frequency of the image is lower than the cutoff
  frequency, then the output is 1 and 0 if the frequency is greater than the cutoff frequency. So first, convert the image in frequency
  domain and then compare it with the cutoff frequency.
  A lot of fringes are present in the image since we immediately cut off all frequencies that go beyond this radius of the mask.
 

 b. Ideal High pass filter: Similarly like the ideal low pass filter, the process is same just the difference comes in the final
  condition. In high pass filter, the output is 0 if the image frequency is less than the cutoff frequency and 1 if it is higher.

 c. Butterworth low pass filter: To obtain the results of a filter, we either convolute the image and the filter mask in the spatial
  domain or multiply the image and mask in the frequency domain. The mask for butterworth low pass filter is given by 1/1+[D(u,v)/D0]^2N
  where D0 is a positive constant and D(u,v) is the distance between a point (u,v) in the frequency domain and the center of the
  frequency rectangle and N is the order of the filter.
  Here, we observe that the fringes that existed in the filtered image when we applied an ideal low pass filter are all gone, since we
  don't abruptly cut off all frequencies beyonf the radius, instead, we build the filter in such a way that it allows the image's 
  frequencies beyond the radius of the filter.

 d. Butterworth High pass filter: In this filter, the steps are same like butterworth low pass. The only difference is in the mask. The
  mask is given by 1/1+[D0/D(u,v)]^2N. This mask allows only high frequencies to pass and blocks the low frequency.

 e. Gaussian low pass filter: For this filter, the mask is given as e^-[D(u,v)^2]/2D0^2. This mask is multiplied by the dft image
  obtained and the output is the filtered image.
  A similar effect is seen as in the butterworth low pass filter where it eliminates fringes in the image, which occurs as a result
  of applying a filter.

 f. Gaussian high pass filter: For this filter, the mask is given as 1 - e^-[D(u,v)^2]/2D0^2. This is simple 1 - the mask for gaussian
  low pass filter as low pass + high pass = 1. So by subtracting the low frequencies from 1 gives output as high frequencies.


3. Post Processing-

  I have performed this function as part of the filtering() method that follows this, as I found it easier that way.
  
  
4. Filtering-

  All the steps provided in the method was followed, which included -
  
  -Applying a fourier transform on the image which results in the image being transformed to the frequency domain.
  -Shifting the lower frequencies to the top corners by using fftshift() method.
  -Performing convolution by using a desired mask.
  -Post processing is done logarithmically to display the DFT of the filtered image.
  -The lower frequencies are shifted back to the center.
  -Inverse transform is performed on this to get the DFT of the filtered image.
  -This is then brought back into the spatial domain and displayed in the output folder.
  
