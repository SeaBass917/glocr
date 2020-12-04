add_library(utils utils.cpp utils.h)
add_library(imageprocessing imageprocessing.cpp imageprocessing.h utils.cpp utils.h)
add_library(imagesegmentation imagesegmentation.cpp imagesegmentation.h imageprocessing.cpp imageprocessing.h utils.cpp utils.h)