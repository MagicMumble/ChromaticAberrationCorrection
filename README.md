# ChromaticAberrationCorrection
trying to implement an algorithm which is dedicated to eliminate chromatic aberrations from image


# About this project

The program takes one testing image in PGM format with large amount of black circles on the white (or light-colored) background with the same radius size. Implemented 
algorithm makes measurements according to the testing image and is able to work with any other image exposed with the same camera settings. 

We suppose that camera wasn't set perpendicular to the image plane so the circles could become ellipses. We need to build affine transformation 
to match green and red (similarly green and blue) channels. Green channel was chosen to be the main channel because it prevails in the Bayer filter.

Program defines the centers and radii of circles by checking the luminance of each pixel (check the structure of pgm files). Supposingly luminance value is 
decreasing by moving away from the circle center to its boundary.

The main idea of the algorithm is to build polinomials from two arguments (coords x and y - circle center) which move red and blue channels (of the testing image)
to the green one. By firstly matching centers of all three channels and then every pixel of the red and blue channel images to the green one we'll get an image free 
of chromatic aberration and be able to use these polynomials for any image with the same camera settings.




