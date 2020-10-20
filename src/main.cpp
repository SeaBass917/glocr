#include "imageprocessing.h"
#include "utils.h"

#include <iostream>

void cleanDocuments(char const* documentsDirectory, char const* documentsCleanedDirectory){

    // Check that the data documents are on file, and
    // make an ouput directory if we need it
    if(fs::exists(documentsDirectory)){
        if(!fs::exists(documentsCleanedDirectory)) fs::create_directories(documentsCleanedDirectory);

        // Loop through each document in the directory, 
        // clean the document, 
        // and store it in the clean directory
        for(const auto& documentFile : fs::directory_iterator(documentsDirectory)){
            std::string documentPath = documentFile.path();
            std::string ext = documentFile.path().extension();
            std::string filename = documentFile.path().filename();

            if(ext == ".png"){
                png::image<png::gray_pixel> imgDoc(documentPath);

                png::image<png::gray_pixel> imgDocClean = preProcessDocument(imgDoc);

                imgDocClean.write(documentsCleanedDirectory+filename);
            }
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
