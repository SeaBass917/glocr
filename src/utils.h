#ifndef UTILS_H_DEFINED
#define UTILS_H_DEFINED
/**************************************
 * Generic functions/classes/structs 
 * for QoL improvments.
 * 
 * - Tuples
 * - OS checks
 * - filesystem namespace/library in fs
 */

#include <filesystem>
namespace fs = std::filesystem;

template<typename T>
struct tuple2{
    T _0;
    T _1;
};

template<typename T>
struct tuple3{
    T _0;
    T _1;
    T _2;
};

template<typename T>
struct tuple4{
    T _0;
    T _1;
    T _2;
    T _3;
};

// Returns the path given with the extension removed
//  NOTE: Returns "" If path is a directory or has no extension
char const* getPathNoExtension(char const* path);
std::string getPathNoExtension(std::string path);

// Returns the extension of the given path
char const* getFileExtension(char const* path);
std::string getFileExtension(std::string path);

#endif