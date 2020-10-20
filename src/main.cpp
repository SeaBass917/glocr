#include "imageprocessing.h"
#include "utils.h"

#include <iostream>

void cleanDocuments(char const* documentsDirectory, char const* documentsCleanedDirectory){

    // Check that the data documents are on file, and
    // make an ouput directory if we need it
    if(fs::exists(documentsDirectory)){
        if(!fs::exists(documentsCleanedDirectory)) fs::create_directories(documentsCleanedDirectory);

        // Trackers for communicating with user
        unsigned numFiles = 0;
        for(const auto& documentFile : fs::directory_iterator(documentsDirectory))
            numFiles++;
        
        // Loop through each document in the directory, 
        // clean the document, 
        // and store it in the clean directory
        unsigned iFile = 1;
        for(const auto& documentFile : fs::directory_iterator(documentsDirectory)){
            std::string documentPath = documentFile.path();
            std::string ext = documentFile.path().extension();
            std::string filename = documentFile.path().filename();
            std::string pathOut = documentsCleanedDirectory+filename;

            // Skip files alread processed and non-png files
            if(!fs::exists(pathOut) && ext == ".png"){
                std::cout << "Cleaning file ("<<iFile<<"/"<<numFiles<<")\""<<filename<<"\"." << std::endl;
                png::image<png::gray_pixel> imgDoc(documentPath);

                png::image<png::gray_pixel> imgDocClean = preProcessDocument(imgDoc);

                // Use dimension ratio as a sanity check, if its too wide something went wrong
                float const ratio = (float)imgDocClean.get_width() / (float)imgDocClean.get_height();
                if (ratio < 3.3f){
                    imgDocClean.write(pathOut);
                }
                else{
                    std::cout << "\tWarning. Error detected in file cleaning. Won't save result." << std::endl;
                }
            }
            else{
                std::cout << "Skipping file ("<<iFile<<"/"<<numFiles<<")\""<<filename<<"\"." << std::endl;
            }
            iFile++;
        }
    }
    else{
        std::cerr << "\tERROR! cleanDocuments() recieved bad address to documents directory. Cannot reach \"" << documentsDirectory << "\"." << std::endl;
        throw std::exception();
    }
}

int main(int argc, char const *argv[]){

    char const* documentsDirectory = "data/IAM/documents/";
    char const* documentsCleanedDirectory = "data/IAM/documentsCleaned/";

    cleanDocuments(documentsDirectory, documentsCleanedDirectory);

    return 0;
}
