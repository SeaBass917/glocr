#include "imageprocessing.h"
#include "utils.h"

#include <iostream>

void cleanDocuments(char const* documentsDirectory, char const* documentsCleanedDirectory){

    // Check that the data documents are on file, and
    // make an ouput directory if we need it
    if(fs::exists(documentsDirectory)){
        if(!fs::exists(documentsDirectory)) fs::create_directories(documentsCleanedDirectory);

        fs::directory_iterator(documentsDirectory)
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
