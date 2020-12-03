#include "../src/imageprocessing.h"
#include "../src/utils.h"

#include <iostream>
#include <string>
#include <chrono> 

// TODO: Automate test a little bit
void edgeDetection_test(std::string testDir, std::string testImg){

    // Address for test image edgemaps
    std::string testPathNoExt = testDir + getPathNoExtension(testImg);
    std::string pathEdgeSobel = testPathNoExt+"_edge-sobel.png";
    std::string pathEdgeMedian = testPathNoExt+"_edge-median.png";

    // Timers
    std::chrono::_V2::system_clock::time_point start, stop;
    std::chrono::microseconds duration;

    // Inform user
    std::cout << "Using image \""<<testImg<<"\" to test edge detectors." << std::endl;

    png::image<png::gray_pixel> img(testImg);

    // Sobel Edge Map
    std::cout << "Running Sobel edge detection." << std::endl;
    start = std::chrono::high_resolution_clock::now();
    png::image<png::gray_pixel> edgeMapImage = edgeMapImg(img, SOBEL, 30);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "\tTime Elapsed: " << duration.count() / 1000.0f << "ms." << std::endl;
    edgeMapImage.write(pathEdgeSobel);

    // Median Edge Map
    std::cout << "Running JR-Median edge detection." << std::endl;
    start = std::chrono::high_resolution_clock::now();
    edgeMapImage = edgeMapImg(img, MEDIAN, 40);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "\tTime Elapsed: " << duration.count() / 1000.0f << "ms." << std::endl;
    edgeMapImage.write(pathEdgeMedian);
}

void cleanIAMDocument_test(std::string testDir, std::string testImg){
    if(fs::exists(testImg)){
        std::cout << "Cleaning test document \""<<testImg<<"\"." << std::endl;

        // Timers
        std::chrono::_V2::system_clock::time_point start, stop;
        std::chrono::microseconds duration;

        png::image<png::gray_pixel> imgDoc(testImg);

        start = std::chrono::high_resolution_clock::now();
        png::image<png::gray_pixel> imgDocClean = cleanIAMDocument(imgDoc);
        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::cout << "\tTime Elapsed: " << duration.count() / 1000.0f << "ms." << std::endl;

        std::string fileOut = testDir+"document_clean.png";
        imgDocClean.write(fileOut);
    }
    else{
        std::cerr << "\tERROR! preProcessing_test() cannot find \""<<testImg<<"\" needed for test." << std::endl;
    }
}

void erosion_test(std::string testDir, std::string testDocumentPath){
    if(fs::exists(testDocumentPath)){
        std::cout << "Edge-mapping and Eroding image \""<<testDocumentPath<<"\"." << std::endl;

        // Timers
        std::chrono::_V2::system_clock::time_point start, stop;
        std::chrono::microseconds duration;

        std::string testDocumentPathNoExt = testDir + getPathNoExtension(testDocumentPath);

        png::image<png::gray_pixel> imgDoc(testDocumentPath);
        
        // Edge detection step
        std::cout << "Edge-Mapping Image..." << std::endl;
        start = std::chrono::high_resolution_clock::now();
        png::image<png::gray_pixel> imgEdge = edgeMapImg(imgDoc, SOBEL, 60);
        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::cout << "\tTime Elapsed: " << duration.count() / 1000.0f << "ms." << std::endl;

        // Image erosion step
        std::cout << "Eroding Image..." << std::endl;
        start = std::chrono::high_resolution_clock::now();
        png::image<png::gray_pixel> imgEroded = erodeImg(imgEdge, 3);
        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::cout << "\tTime Elapsed: " << duration.count() / 1000.0f << "ms." << std::endl;

        imgEdge.write(testDocumentPathNoExt+"_edge.png");
        imgEroded.write(testDocumentPathNoExt+"_eroded.png");
    }
    else{
        std::cerr << "\tERROR! preProcessing_test() cannot find \""<<testDocumentPath<<"\" needed for test." << std::endl;
    }
}

void kittler_test(std::string testDir, std::string testImgPath){
    if(fs::exists(testImgPath)){
        std::cout << "Thresholding image \""<<testImgPath<<"\" with kittler's method." << std::endl;

        // Timers
        std::chrono::_V2::system_clock::time_point start, stop;
        std::chrono::microseconds duration;

        png::image<png::gray_pixel> img(testImgPath);
        
        // Determine threshold
        std::cout << "Finding Kittler threshold..." << std::endl;
        start = std::chrono::high_resolution_clock::now();
        unsigned t = determineKittlerThreshold(img);
        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::cout << "\tTime Elapsed: " << duration.count() / 1000.0f << "ms." << std::endl;

        // Threshold image
        std::cout << "Thresholding Image..." << std::endl;
        start = std::chrono::high_resolution_clock::now();
        png::image<png::gray_pixel> imgKittler = segmentImageThreshold(img, t);
        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::cout << "\tTime Elapsed: " << duration.count() / 1000.0f << "ms." << std::endl;

        std::string outputImgPath = testDir + getPathNoExtension(testImgPath)+"_kittlerThresh.png";
        imgKittler.write(outputImgPath);
    }
    else{
        std::cerr << "\tERROR! kittler_test() cannot find \""<<testImgPath<<"\" needed for test." << std::endl;
    }
}

int main(int argc, char const *argv[]){

    // Create a directory for doing tests
    std::string testDir = "imageprocessing_test/";
    if(!fs::exists(testDir)) fs::create_directory(testDir);

    // Draw a smile
    std::string pathSmile = testDir+"smile.png";
    std::cout << "Drawing a smile at \""<<pathSmile<<"\"." << std::endl;
    drawASmile(pathSmile);

    std::string testDir0 = testDir+"edgeDetection_test/";
    if(!fs::exists(testDir0)) fs::create_directory(testDir0);
    edgeDetection_test(testDir0, "document.png");

    std::string testDir1 = testDir+"cleanIAMDocument_test/";
    std::string testImg1 = "../data/IAM/documents/a02-000.png";
    if(!fs::exists(testDir1)) fs::create_directory(testDir1);
    cleanIAMDocument_test(testDir1, testImg1);

    std::string testDir2 = testDir+"erosion_test/";
    if(!fs::exists(testDir2)) fs::create_directory(testDir2);
    erosion_test(testDir2, "document_cleaned.png");

    std::string testDir3 = testDir+"kittler_test/";
    if(!fs::exists(testDir3)) fs::create_directory(testDir3);
    kittler_test(testDir3, "document.png");

    return 0;
}