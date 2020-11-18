#include "imagesegmentation.h"
#include "imageprocessing.h"

/*
 * IAM Utilities
 */

std::map<std::string, std::array<unsigned, 7>> loadIAMDocumentsMetadata(std::string metaDataPath){
    
    // Open the metadata file
    std::ifstream fp(metaDataPath);
    if(fp){

        std::map<std::string, std::array<unsigned, 7>> metaData;

        // Parse each line, skip comments #---, read the fields into the array
        std::string line;
        while (std::getline(fp, line)) {
            if(line[0] == '#') continue;
            else{
                unsigned const i0 = line.find(' ');
                unsigned const i1 = line.find(' ', i0+1);
                unsigned const i2 = line.find(' ', i1+1);
                unsigned const i3 = line.find(' ', i2+1);
                unsigned const i4 = line.find(' ', i3+1);
                unsigned const i5 = line.find(' ', i4+1);
                unsigned const i6 = line.find(' ', i5+1);
                unsigned const i7 = line.size();
                
                std::string const docID = line.substr(0, i0);
                unsigned const writerID = std::stoi(line.substr(i0+1, i1 - (i0+1)));
                unsigned const numSent = std::stoi(line.substr(i1+1, i2 - (i1+1)));
                std::string const segQualStr = line.substr(i2+1, i3 - (i2+1));
                unsigned const segQual = (segQualStr == "all")? ALL : PRT;
                unsigned const numLines = std::stoi(line.substr(i3+1, i4 - (i3+1)));
                unsigned const numLinesGood = std::stoi(line.substr(i4+1, i5 - (i4+1)));
                unsigned const numWords = std::stoi(line.substr(i5+1, i6 - (i5+1)));
                unsigned const numWordsGood = std::stoi(line.substr(i6+1, i7 + (i6+1)));

                metaData[docID] = {writerID, numSent, segQual, numLines, numLinesGood, numWords, numWordsGood};
            }
        }

        return metaData;
    }
    else{
        std::cerr << "\tWARNING! Cannot access \""<<metaDataPath<<"\". Could not load IAM document metadata." << std::endl; 
        throw std::exception();
    }
}

/*
 * Segmentation
 */

png::image<png::gray_pixel> preProcessDocumentImage(png::image<png::gray_pixel> const& imgDoc){
    
    png::image<png::gray_pixel> imgEdges = edgeMapImg(imgDoc, SOBEL, 30);
    // imgEdges.write("../edges.png");
    png::image<png::gray_pixel> imgEroded = erodeImg(imgEdges, 3);
    // imgEroded.write("../eroded.png");
    
    return imgEroded;
}

png::image<png::gray_pixel> preProcessWordImage(png::image<png::gray_pixel> const& img){
    
    unsigned t = determineKittlerThreshold(img);
    png::image<png::gray_pixel> imgThresh = segmentImageThreshold(img, t);
    imgThresh.write("../thresh.png");
    png::image<png::gray_pixel> imgSkeleton = skeletonizeImg(imgThresh);
    imgSkeleton.write("../skeleton.png");
    // png::image<png::gray_pixel> imgClean = noiseReduxSPImg(imgSkeleton);
    // imgClean.write("../clean.png");
    
    return imgSkeleton;
}

std::vector<std::vector<png::image<png::gray_pixel>>> wordSegmentation(png::image<png::gray_pixel> const &imgDoc){

    unsigned const histThreshold = 127u;                    // Minimum intensity for the histogram to count a pixel
    unsigned const peakThresholdY = 55u;                      // Minumum threshold for a hist peak to be considered a line
    unsigned const peakThresholdX = 1u;                       // Minumum threshold for a hist peak to be considered a word
    unsigned const minNoiseFloorY = 1u;                     // Floor that is considered the bounds of a line
    unsigned const minNoiseFloorX = 0u;                     // Floor that is considered the bounds of a word
    unsigned const edgeMapThreshold = 60u;                  // Threshold for the edge map to use
    unsigned const yPadding = 5u;                          // Number of pixels to pad the output with in the y direction
    unsigned const xPadding = 5u;                          // Number of pixels to pad the output with in the x direction
    unsigned const numColumnsForHist = 450u;                // Number of lefthand columns to use for the histogram projection
    unsigned const minWordSpacing = 18u;                    // Min spacing that will be called a seperate word
    
    // Image dimensions
    unsigned const height = imgDoc.get_height();
    unsigned const width = imgDoc.get_width();

    // ---------------
    // Pre-Processing
    // ---------------

    png::image<png::gray_pixel> imgDocPreProc = preProcessDocumentImage(imgDoc);
    // crop(imgDocPreProc, 0, numColumnsForHist, 0, height-1).write("/home/seabass/extra/pprocimg.png");

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
        bool isEdgeInRow = false;
        std::vector<png::gray_pixel> row = imgDocPreProc[i];
        for(unsigned j = 0; j < width; j++){
            if(row[j] < histThreshold){
                isEdgeInRow = true;
                break;
            }
        }
        if(isEdgeInRow){
            yFirstBin = i - yPadding;
            if(yFirstBin < 0){
                yFirstBin = 0;
            }
            break;
        }
    }
    for(unsigned i = height-1; 0 < i; i--){
        bool isEdgeInRow = false;
        std::vector<png::gray_pixel> row = imgDocPreProc[i];
        for(unsigned j = 0; j < width; j++){
            if(row[j] < histThreshold){
                isEdgeInRow = true;
                break;
            }
        }
        if(isEdgeInRow){
            yLastBin = i + yPadding;
            if(height <= yLastBin){
                yLastBin = height-1;
            }
            break;
        }
    }

    // Determine the midpoints between each line
    std::vector<tuple2<int>> midPointBoundsY = getMidPointBounds(histHorizontal, peakThresholdY, minNoiseFloorY);
    unsigned numLines = midPointBoundsY.size() + 1;

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
    for(tuple2<int> const& midPointBounds : midPointBoundsY){

        // Determine a line through each midpoint that seperates the lines
        // Used that line to crop out the next line
        // NOTE: Opeation is repeated on both the preprocessed doc and the og
        yArrays[iMPY] = spaceTraceHorizontal(imgDocPreProc, midPointBounds._0, midPointBounds._1);
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
        std::vector<tuple2<int>> midPointBoundsX = getMidPointBounds(histVertical, peakThresholdX, minNoiseFloorX, minWordSpacing);
        unsigned numWords = midPointBoundsX.size() + 1;
        
        // For storing staggered X values that segment each line
        std::vector<std::vector<unsigned>> xArrays(numWords+1);

        // Initialize the list of lines
        std::vector<png::image<png::gray_pixel>> imgWords(numWords);

        // Determine the y values to cut along, cut, and save where we cut
        // Initialize the current yarray to be a straight line at the first non-zero bin
        unsigned iMPX = 1;
        xArrays[0] = std::vector<unsigned>(height, xFirstBin);
        for(tuple2<int> const& midPointBounds : midPointBoundsX){

            // Determine a line through each midpoint that seperates the lines
            // Used that line to crop out the next line
            xArrays[iMPX] = spaceTraceVertical(imgLinePreProc, midPointBounds._0, midPointBounds._1);
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

std::vector<png::image<png::gray_pixel>> charSegmentation(png::image<png::gray_pixel> const &imgword){
    // std::cout << "TODO: Check to see if you can measure and adjust skew in a word." << std::endl;
    
    // Image dimensions
    unsigned const height = imgword.get_height();
    unsigned const width = imgword.get_width();

    // ---------------
    // Pre-Processing
    // ---------------

    png::image<png::gray_pixel> imgDocPreProc = preProcessWordImage(imgword);
    // crop(imgDocPreProc, 0, numColumnsForHist, 0, height-1).write("/home/seabass/extra/pprocimg.png");

    return std::vector<png::image<png::gray_pixel>>();
}

/*
 * Tracing
 */

std::vector<unsigned> spaceTraceHorizontal(png::image<png::gray_pixel> const& img, unsigned const yMin, unsigned const yMax){

    unsigned const minHistThreshold = 3;    // Minimum count for a bin to be considered populated

    // Img info
    unsigned const height = img.get_height();
    unsigned const width = img.get_width();

    // Histogram will have 3 bins in between the yMin and yMax
    // (yMin)|       |(y0)|       |(y1)|       |(yMax)
    float const dy = (yMax - yMin) / 3.0f;
    unsigned const y0 = (unsigned)(yMin + dy);
    unsigned const y1 = (unsigned)(y0 + dy);
    unsigned hist[3] = {0};

    // Precompute the midpoints inside the bins
    unsigned const midPointHigh = (yMin + y0) / 2;
    unsigned const midPointMid = (y0 + y1) / 2;
    unsigned const midPointLow = (y1 + yMax) / 2;

    // Saftey code
    if(yMin < yMax){
        if(yMax < height){

            // Sweep from left to right
            //   on finding an edge
            //   begin a histogram until edges stop
            //   draw a mipoint based on the histogram
            //   
            std::vector<unsigned> yArray(width);
            unsigned x = 0;
            while(x < width){

                // yArray[x] = midPointMid;
                // x++;
                // continue;

                // Look for an edge (pixel == 0)
                bool edgeFound = false;
                for(unsigned y = yMin; y <= yMax; y++){
                    if(img[y][x] < 127){
                        edgeFound = true;
                        break;
                    }
                }

                // Move forward building a 3bin histogram along the way
                // build until no edges found
                // Bisect the lines based on histogram cases
                if(edgeFound){
                    hist[0] = 0;
                    hist[1] = 0;
                    hist[2] = 0;

                    // Break when we hit the end of image of stop seeing edges
                    edgeFound = true;
                    unsigned xx = x; // NOTE we need a seperate xx var so we can update the yArray after from x -> xx
                    while(xx < width && edgeFound){

                        // Fill hist
                        edgeFound = false;
                        for(unsigned y = yMin; y < y0; y++){
                            if(img[y][xx] < 127){
                                hist[0]++;
                                edgeFound = true;
                            }
                        }
                        for(unsigned y = y0; y < y1; y++){
                            if(img[y][xx] < 127){
                                hist[1]++;
                                edgeFound = true;
                            }
                        }
                        for(unsigned y = y1; y <= yMax; y++){
                            if(img[y][xx] < 127){
                                hist[2]++;
                                edgeFound = true;
                            }
                        }

                        xx++;
                    }

                    // Determine the midpoint for each case in the histogram
                    unsigned midPoint = midPointMid;
                    if(minHistThreshold < hist[0]){
                        if(minHistThreshold < hist[1]){
                            if(minHistThreshold < hist[2]){
                                midPoint = yMax;
                            }
                            else{
                                midPoint = midPointLow;
                            }
                        }
                        else{
                            midPoint = midPointMid;
                        }
                    }
                    else{
                        midPoint = midPointHigh;
                    }

                    // Update yArray
                    // NOTE: this also updates the main x variable 
                    // NOTE: we dont write to xx, as it is past the edges
                    for(x; x < xx; x++){
                        yArray[x] = midPoint;
                    }
                }
                else{
                    yArray[x] = midPointMid;
                    x++;
                }
            }

            return yArray;
        }
        else{
            std::cerr << "\tERROR! In spaceTraceHorizontal() yMax: "<<yMax<<" goes outside height bounds of image with height: "<<height<<"." << std::endl;
            throw std::exception();
        }
    }
    else{
        std::cerr << "\tERROR! In spaceTraceHorizontal() yMin: "<<yMin<<" needs to be less than than yMax: "<<yMax<<"." << std::endl;
        throw std::exception();
    }
}

std::vector<unsigned> spaceTraceVertical(png::image<png::gray_pixel> const& img, unsigned const xMin, unsigned const xMax){

    // Img info
    unsigned const height = img.get_height();
    unsigned const width = img.get_width();
    unsigned const midPoint = (xMin + xMax) / 2;

    // Saftey code
    if(xMin < xMax){
        if(xMax < width){

            std::vector<unsigned> xArray(height);

            // TODO: tracing
            // DEBUG
            for(unsigned i = 0; i < height; i++){
                xArray[i] = midPoint;
            }

            return xArray;

        }
        else{
            std::cerr << "\tERROR! In spaceTraceVertical() xMax: "<<xMax<<" goes outside width bounds of image with width: "<<width<<"." << std::endl;
            throw std::exception();
        }
    }
    else{
        std::cerr << "\tERROR! In spaceTraceVertical() xMin: "<<xMin<<" needs to be less than than xMax: "<<xMax<<"." << std::endl;
        throw std::exception();
    }
}
