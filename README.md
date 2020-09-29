# glocr

## Requirements
 - [zlib](https://zlib.net/)
 - [libpng](http://www.libpng.org/pub/png/libpng.html)
 - [png++](https://www.nongnu.org/pngpp/)

## Simple Example Setup
 - In linux env
 - install packages in your env
    * zlib (sudo apt-get install zlib1g-dev)
    * libpng (sudo apt-get install libpng-dev)
 - Grab requirements (zlib, libpng, png++)
 - build them (NOTE: png++ is header only library, and libpng needs zlib)
 - headers in the INCLUDE_DIR, .a, .so in the LIB_DIR
 - g++ -g src/main.cpp src/image.cpp -I ${INCLUDE_DIR} -L ${LIB_DIR} -lz -lpng -o glocr
 - ./glocr