#include "utils.h"

#include <string>
#include <iostream>

/*
 * Filesystem Helpers
 */

char const* getPathNoExtension(char const* path){
    std::string sPath = std::string(path);
    return getPathNoExtension(sPath).c_str();
}

std::string getPathNoExtension(std::string sPath){
    if(!fs::is_directory(sPath)){

        // Trim off the leading dirs in the path to get just the filename
        int i = sPath.find_last_of('/');
        bool iValid = 0 <= i && i+1 < sPath.length();
        std::string sParentPath = (iValid)? sPath.substr(0,i+1) : "";
        std::string sFilename = (iValid)? sPath.substr(i+1) : sPath;

        i = sFilename.find_last_of('.');
        if(0 <= i){
            return sParentPath + sFilename.substr(0, i);
        }
    }
    return sPath;
}

char const* getFileExtension(char const* path){
    std::string sPath = std::string(path);
    return getFileExtension(sPath).c_str();
}

std::string getFileExtension(std::string sPath){
    if(!fs::is_directory(sPath)){

        // Trim off the leading dirs in the path
        int i = sPath.find_last_of('/');
        std::string sFilename = (0 <= i && i+1 < sPath.length())? sPath.substr(i+1) : sPath;

        // Find the extension after the last '.'
        i = sFilename.find_last_of('.');
        if(0 <= i && i+1 < sFilename.length()) return sFilename.substr(i+1).c_str();
    }
    return "";
}

/*
 * Strings
 */

std::string toLower(std::string s){
    unsigned len = s.length();
    char* cs = (char*)malloc(sizeof(char)*len);
    if(cs){
        for(unsigned i = 0; i < len; i++){
            cs[i] = std::tolower(s[i]);
        }

        std::string sLower(cs);
        free(cs);
        return sLower;
    }
    else{
        std::cerr << "ERROR! Failed to allocate memory for toLower()." << std::endl;
        throw std::exception();
    }
}

/*
 * Maths
 */

unsigned median(unsigned* a, unsigned len){
    std::sort(a, a+len);
    unsigned iHalf = len / 2;
    return a[iHalf];
}
int median(int* a, unsigned len){
    std::sort(a, a+len);
    unsigned iHalf = len / 2;
    return a[iHalf];
}
float median(float* a, unsigned len){
    std::sort(a, a+len);
    unsigned iHalf = len / 2;
    return a[iHalf];
}
double median(double* a, unsigned len){
    std::sort(a, a+len);
    unsigned iHalf = len / 2;
    return a[iHalf];
}
char median(char* a, unsigned len){
    std::sort(a, a+len);
    unsigned iHalf = len / 2;
    return a[iHalf];
}