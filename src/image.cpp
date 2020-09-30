#include "image.h"

void horizontalProjectionHistogram(png::image<png::gray_pixel> const &img, double* hist, bool normalize){
    if(hist){

        // Img info
        unsigned height = img.get_height();
        unsigned width = img.get_width();

        unsigned histLeng = *(&hist + 1) - hist;
        if(height <= histLeng){

            // For normilization
            unsigned totalSum = 0;

            // For each row count pixels that are below the threshold
            unsigned thresh = 64;
            for(unsigned r = 0; r < height; r++){

                // Read row once
                const std::vector<png::byte, std::allocator<png::byte>> row = img[r];

                unsigned sum = 0;
                for(unsigned c = 0; c < width; c++){
                    if(row[c] < thresh){
                        sum++;
                        totalSum++;
                    }
                }

                hist[r] = sum;
            }

            // Normalize the histogram if requested
            if(normalize){
                double normScale = 1.0 / (double)totalSum;
                for(unsigned r = 0; r < height; r++)
                    hist[r] *= normScale;
            }
        }
        else{
            std::cout << "\tERROR! horizontalProjectionHistogram() recieved a histogram smaller["<<histLeng<<"] than image height["<<height<<"]. Failed to make a histogram." << std::endl;
            std::cout << sizeof(hist) << std::endl;
            std::cout << sizeof(hist[0]) << std::endl;
        }
    }
    else{
        std::cout << "\tERROR! horizontalProjectionHistogram() recieved a null histogram. Failed to make a histogram." << std::endl;
    }
}

void verticalProjectionHistogram(png::image<png::gray_pixel> const &img, double* hist, bool normalize){
    if(hist){

        // Img info
        unsigned height = img.get_height();
        unsigned width = img.get_width();

        unsigned histLeng = sizeof(hist)/sizeof(hist[0]);
        if(width <= histLeng){


            // For normilization
            unsigned totalSum = 0;

            // For each col count pixels that are below the threshold
            unsigned thresh = 127;
            for(unsigned c = 0; c < width; c++){

                unsigned sum = 0;
                for(unsigned r = 0; r < height; r++){
                    if(img[r][c] < thresh){
                        sum++;
                        totalSum++;
                    }
                }

                hist[c] = sum;
            }

            // Normalize the histogram if requested
            if(normalize){
                double normScale = 1.0 / (double)totalSum;
                for(unsigned r = 0; r < height; r++)
                    hist[r] *= normScale;
            }
        }
        else{
            std::cout << "\tERROR! verticalProjectionHistogram() recieved a histogram smaller["<<histLeng<<"] than image width["<<width<<"]. Failed to make a histogram." << std::endl;
        }
    }
    else{
        std::cout << "\tERROR! verticalProjectionHistogram() recieved a null histogram. Failed to make a histogram." << std::endl;
    }
}

// Function to create Gaussian kernel on length n
template<class T>
T** gaussiankernel(unsigned const masklength, T const sigma){

    T** kernel;

    // Mask length must be odd
    if(masklength & 1 == 1){

        const int radius = (masklength - 1) / 2;

        // intialising standard deviation
        T var = 2.0 * sigma * sigma;
        T denom = M_PI * var;
    
        // Allocating kernel 
        kernel = (T**)malloc(masklength * sizeof(T*));
        for (unsigned i = 0; i < masklength; i++){
            kernel[i] = (T*)malloc(masklength * sizeof(T));
        }

        T sum = 0.0; // sum is for normalization 
        for (int x = -radius; x <= radius; x++) {
            T rx = x * x;
            for (int y = -radius; y <= radius; y++) {
                T Gxy = exp(-(rx + y * y) / var) / denom;
                kernel[x + radius][y + radius] = Gxy;
                sum += Gxy; 
            } 
        }

        // normalising the Kernel 
        for (int i = 0; i < masklength; ++i){
            for (int j = 0; j < masklength; ++j){
                kernel[i][j] /= sum;
            }
        }
    }
    
    return kernel;
} 

// Perform a 2d convolution on an image object
template<typename T>
png::image<png::rgb_pixel> convd2D(png::image<png::rgb_pixel> &image, unsigned const maskHeight, unsigned const maskWidth, T** kernel){

    // Image dimm
    unsigned height = image.get_height();
    unsigned width = image.get_width();

    // Make a new image of same dimmension
    png::image<png::rgb_pixel> imageTgt(height, width);

    // Verify odd mask dimmensions
    if(maskHeight & 1 == 1 && maskWidth & 1 == 1){

        // Get kernel radii
        unsigned maskHeightRadius = (maskHeight-1) / 2;
        unsigned maskWidthRadius = (maskWidth-1) / 2;

        // Loop through th image pixels
        for (unsigned i = 0; i < height; i++){
            for (unsigned j = 0; j < width; j++){
                
                //convolve with kernel
                T sR = 0.0f;
                T sG = 0.0f;
                T sB = 0.0f;

                for(unsigned m = 0; m < maskHeight; m++){
                    int ii = i - maskHeightRadius + m;
                    
                    // Mirror off edges
                    if(ii < 0){
                        ii = -ii;
                    }
                    if(height <= ii){
                        ii = height - 2 - (ii - height);
                    }

                    for(unsigned n = 0; n < maskWidth; n++){
                        int jj = j - maskWidthRadius + n;
                        // Mirror off edges
                        if(jj < 0){
                            jj = -jj;
                        }
                        if(width <= ii){
                            jj = width - 2 - (jj - width);
                        }

                        png::rgb_pixel pixel = image[ii][jj];

                        T r = (T)pixel.red;
                        T g = (T)pixel.green;
                        T b = (T)pixel.blue;

                        T a = kernel[m][n];

                        sR += r * a;
                        sG += g * a;
                        sB += b * a;
                    }
                }

                imageTgt[i][j] = png::rgb_pixel(sR, sG, sB);
            }
        }
    }
    else{
        std::cout << "\tERROR! 2Dconv cannot handle even dim kernel.";
    }

    return imageTgt;
}

void drawASmile(char const* addrOut){
    
    unsigned const srcRes = 8;
    unsigned const tgtRes = 1024;
    unsigned const scale = tgtRes/srcRes;

    // Lower resolution image template
    std::vector<std::vector<unsigned>> imgTemplate {
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 1, 0, 0},
        {0, 0, 1, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 1, 0},
        {0, 1, 0, 0, 0, 0, 1, 0},
        {0, 0, 1, 1, 1, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
    };
    
    // Create the image
    png::image< png::rgb_pixel > image(tgtRes, tgtRes);
    for (unsigned y = 0; y < tgtRes; y++){
        unsigned ySrc = y/scale;
        for (unsigned x = 0; x < tgtRes; x++){
            unsigned xSrc = x/scale;
            
            if(imgTemplate[ySrc][xSrc]) image[y][x] = png::rgb_pixel(0, 0, 0);
            else image[y][x] = png::rgb_pixel(127, 200, 0);
        }
    }

    // Gaussian blur
    unsigned const masklength = 25u;
    float const sigma = 54.0f;
    
    // Gaussian kernel
    float** gKernel = gaussiankernel<float>(masklength, sigma);

    image = convd2D<float>(image, masklength, masklength, gKernel);

    // Write image out
    image.write(addrOut);
}
