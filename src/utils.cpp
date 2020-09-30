#include "utils.h"

bool pathExists(char const* path){
    struct stat buffer;   
    return (stat (path, &buffer) == 0); 
}