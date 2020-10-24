#include "imagesegmentation.h"
#include "imageprocessing.h"

/*
 * Segmentation
 */

png::image<png::gray_pixel> preProcessDocumentImage(png::image<png::gray_pixel> const& imgDoc){
    
    png::image<png::gray_pixel> imgEdges = edgeMapImg(imgDoc);
    png::image<png::gray_pixel> imgEroded = erodeImg(imgEdges, 3);
    imgEroded.write("../eroded.png");
    
    return imgEroded;
}

std::vector<std::vector<png::image<png::gray_pixel>>> wordSegmentation(png::image<png::gray_pixel> const &imgDoc){

    unsigned const histThreshold = 127u;                    // Minimum intensity for the histogram to count a pixel
    unsigned const minBinCountY = 5u;                        // minBinCount before the histogram supresses that bin
    unsigned const minBinCountX = 1u;                        // minBinCount before the histogram supresses that bin
    unsigned const edgeMapThreshold = 60u;                  // Threshold for the edge map to use
    unsigned const yPadding = 25u;                          // Number of pixels to pad the output with in the y direction
    unsigned const xPadding = 25u;                          // Number of pixels to pad the output with in the x direction
    unsigned const numColumnsForHist = 750u;                // Number of lefthand columns to use for the histogram projection
    // unsigned const histBinWidth = 25u;                      // Helps with forgiving smaller gaps in a word
    // unsigned const histBinWidthHalf = histBinWidth / 2;     // for centering the midpoint after hist bin change
    
    // Image dimensions
    unsigned const height = imgDoc.get_height();
    unsigned const width = imgDoc.get_width();

    // ---------------
    // Pre-Processing
    // ---------------

    png::image<png::gray_pixel> imgDocPreProc = preProcessDocumentImage(imgDoc);

    // -------------------------
    // Determine Line Midpoints
    // -------------------------

    // Generate a horizontal histogram of the image to find the segments
    // NOTE: We supress bins with too few counts to clean the results a little
    // NOTE: We only use the first X00 columns of the image (helps isolate peaks)
    std::vector<unsigned> histHorizontal = horizontalProjectionHistogram(
        crop(imgDocPreProc, 0, numColumnsForHist, 0, height-1), 
        histThreshold
    );
        
    // Determine the first and last non-zero bins
    // NOTE: the padding is applied here
    int yFirstBin = 0;
    int yLastBin  = height-1;
    for(unsigned i = 0; i < height-1; i++){
        if(0 < histHorizontal[i]){
            yFirstBin = i - yPadding;
            if(yFirstBin < 0){
                yFirstBin = 0;
            }
            break;
        }
    }
    for(unsigned i = height-1; 0 < i; i--){
        if(0 < histHorizontal[i]){
            yLastBin = i + yPadding;
            if(height <= yLastBin){
                yLastBin = height-1;
            }
            break;
        }
    }

    // Determine the midpoints between each line
    std::vector<unsigned> midPointsY = getMidPoints(histHorizontal, minBinCountY);
    unsigned numLines = midPointsY.size() + 1;

    // -----------------------------
    // Trace through Line Midpoints
    // -----------------------------

    // For storing staggered Y values that segment each line
    std::vector<std::vector<unsigned>> yArrays(numLines+1);

    // Initialize the list of lines
    // One will retain preprocessing effects, the other: the original image
    std::vector<png::image<png::gray_pixel>> imgLinesPreProc(numLines);
    std::vector<png::image<png::gray_pixel>> imgLines(numLines);

    // Determine the y values to cut along, cut, and save where we cut
    // Initialize the current yarray to be a straight line at the first non-zero bin
    unsigned iMPY = 1;
    yArrays[0] = std::vector<unsigned>(width, yFirstBin);
    for(unsigned const& midPoint : midPointsY){

        // Determine a line through each midpoint that seperates the lines
        // Used that line to crop out the next line
        // NOTE: Opeation is repeated on both the preprocessed doc and the og
        yArrays[iMPY] = spaceTraceHorizontal(imgDocPreProc, midPoint);
        imgLinesPreProc[iMPY-1] = staggeredYCrop(imgDocPreProc, yArrays[iMPY-1], yArrays[iMPY], yPadding);
        imgLines[iMPY-1] = staggeredYCrop(imgDoc, yArrays[iMPY-1], yArrays[iMPY], yPadding);
        
        iMPY++;
    }

    // Set the last yarray to be a straight line at the last non-zero bin
    // NOTE: iMPY would have incremented to the last index
    yArrays[iMPY] = std::vector<unsigned>(width, yLastBin);
    imgLinesPreProc[iMPY-1] = staggeredYCrop(imgDocPreProc, yArrays[iMPY-1], yArrays[iMPY], yPadding);
    imgLines[iMPY-1] = staggeredYCrop(imgDoc, yArrays[iMPY-1], yArrays[iMPY], yPadding);

    // ------------------
    // Word Segmentation
    // ------------------
    
    // Loop through each line image to segment the words
    // Store the words in this list
    std::vector<std::vector<png::image<png::gray_pixel>>> linesImgs(numLines);
    unsigned iLine = 0;
    for(png::image<png::gray_pixel> const& imgLinePreProc: imgLinesPreProc){
        png::image<png::gray_pixel> imgLine = imgLines[iLine];

        std::vector<unsigned> histVertical = verticalProjectionHistogram(
            imgLinePreProc,
            histThreshold
        );
        unsigned histVerticalLength = histVertical.size();

        // for(unsigned i = 0; i < histVertical.size(); i++)
        //     std::cout << i << "," << histVertical[i] << std::endl;
        // std::cout << std::endl << std::endl;

        // Determine the first and last non-zero bins
        // NOTE: the padding is applied here
        int xFirstBin = 0;
        int xLastBin = histVerticalLength-1;
        for(unsigned i = 0; i < histVerticalLength-1; i++){
            if(0 < histVertical[i]){
                xFirstBin = i - xPadding;
                if(xFirstBin < 0){
                    xFirstBin = 0;
                }
                break;
            }
        }
        for(unsigned i = histVerticalLength-1; 0 < i; i--){
            if(0 < histVertical[i]){
                xLastBin = i + xPadding;
                if(width <= xLastBin){
                    xLastBin = width-1;
                }
                break;
            }
        }

        // Determine the midpoints between each word
        std::vector<unsigned> midPointsX = getMidPoints(histVertical, minBinCountX, 25);
        unsigned numWords = midPointsX.size() + 1;

        // Adjust midpoints based on histogram bin width
        for(unsigned i = 0; i < midPointsX.size(); i++){
            unsigned midPoint = midPointsX[i];
            midPointsX[i] = (width <= midPoint)? width-1 : midPoint;
        }
        
        // For storing staggered X values that segment each line
        std::vector<std::vector<unsigned>> xArrays(numWords+1);

        // Initialize the list of lines
        std::vector<png::image<png::gray_pixel>> imgWords(numWords);

        // Determine the y values to cut along, cut, and save where we cut
        // Initialize the current yarray to be a straight line at the first non-zero bin
        unsigned iMPX = 1;
        xArrays[0] = std::vector<unsigned>(height, xFirstBin);
        for(unsigned const& midPoint : midPointsX){

            // Determine a line through each midpoint that seperates the lines
            // Used that line to crop out the next line
            xArrays[iMPX] = spaceTraceVertical(imgLinePreProc, midPoint);
            imgWords[iMPX-1] = staggeredXCrop(imgLine, xArrays[iMPX-1], xArrays[iMPX], yPadding);
            iMPX++;
        }

        // Set the last x array to be a straight line at the last non-zero bin
        // NOTE: iMPX would have incremented to the last index
        xArrays[iMPX] = std::vector<unsigned>(width, xLastBin);
        imgWords[iMPX-1] = staggeredXCrop(imgLine, xArrays[iMPX-1], xArrays[iMPX], yPadding);

        linesImgs[iLine] = imgWords;

        iLine++;
    }

    return linesImgs;
}

/*
 * Tracing
 */

std::vector<unsigned> spaceTraceHorizontal(png::image<png::gray_pixel> const& img, unsigned const midPoint){

    // Img info
    unsigned height = img.get_height();
    unsigned width = img.get_width();

    std::vector<unsigned> yArray(width);

    // TODO: tracing
    // DEBUG
    for(unsigned i = 0; i < width; i++){
        yArray[i] = midPoint;
    }

    return yArray;
}

std::vector<unsigned> spaceTraceVertical(png::image<png::gray_pixel> const& img, unsigned const midPoint){

    // Img info
    unsigned height = img.get_height();
    unsigned width = img.get_width();

    std::vector<unsigned> xArray(height);

    // TODO: tracing
    // DEBUG
    for(unsigned i = 0; i < height; i++){
        xArray[i] = midPoint;
    }

    return xArray;
}
