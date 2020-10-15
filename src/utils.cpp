#include "utils.h"

#include <string>

char const* getPathNoExtension(char const* path){
    if(!fs::is_directory(path)){
        std::string sPath = std::string(path);

        // Trim off the leading dirs in the path to get just the filename
        int i = sPath.find_last_of('/');
        bool iValid = 0 <= i && i+1 < sPath.length();
        std::string sParentPath = (iValid)? sPath.substr(0,i+1) : "";
        std::string sFilename = (iValid)? sPath.substr(i+1) : sPath;

        i = sFilename.find_last_of('.');
        if(0 <= i){
            return (sParentPath + sFilename.substr(0, i)).c_str();
        }
    }
    return path;
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
    if(!fs::is_directory(path)){
        std::string sPath = std::string(path);

        // Trim off the leading dirs in the path
        int i = sPath.find_last_of('/');
        std::string sFilename = (0 <= i && i+1 < sPath.length())? sPath.substr(i+1) : sPath;

        // Find the extension after the last '.'
        i = sFilename.find_last_of('.');
        if(0 <= i && i+1 < sFilename.length()) return sFilename.substr(i+1).c_str();
    }
    return "";
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
