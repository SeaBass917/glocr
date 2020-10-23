#ifndef UTILS_H_DEFINED
#define UTILS_H_DEFINED
/**************************************
 * Generic functions/classes/structs 
 * for QoL improvments.
 * 
 * - Tuples
 * - OS checks
 * - filesystem namespace/library in fs
 * 
 */

#include <filesystem>
namespace fs = std::filesystem;

/*
 * Tuples
 */

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

/*
 * Timer code
 */

// #include <chrono> 
// auto start = std::chrono::high_resolution_clock::now();
// auto stop = std::chrono::high_resolution_clock::now();
// auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
// std::cout << duration.count() << std::endl;

/*
 * Filesystem Helpers
 */

// Returns the path given with the extension removed
//  NOTE: Returns "" If path is a directory or has no extension
char const* getPathNoExtension(char const* path);
std::string getPathNoExtension(std::string path);

// Returns the extension of the given path
char const* getFileExtension(char const* path);
std::string getFileExtension(std::string path);

/*
 * Strings
 */

// Converts the string to lowercase
std::string toLower(std::string s);

/*
 * Maths
 */

// Returns the median value of the given list
unsigned median(unsigned* a, unsigned len);
int median(int* a, unsigned len);
float median(float* a, unsigned len);
double median(double* a, unsigned len);
char median(char* a, unsigned len);

#endif