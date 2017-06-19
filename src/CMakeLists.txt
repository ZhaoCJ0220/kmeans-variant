cmake_minimum_required(VERSION 3.5)

find_package(OpenCV REQUIRED)
find_package(Boost REQUIRED COMPONENTS regex filesystem)
find_package(Qt5 REQUIRED COMPONENTS Widgets Core Concurrent Test Gui)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(preprocess preprocess.cc utils.cc)
target_link_libraries(preprocess ${OpenCV_LIBS} ${Qt5_LIBS} ${Boost_LIBRARIES})

add_executable(drawPalette drawPalette.cc utils.cc)
target_link_libraries(drawPalette ${OpenCV_LIBS} ${Qt5_LIBS} ${Boost_LIBRARIES})

add_executable(extractPalette extractPalette.cc utils.cc)
target_link_libraries(extractPalette ${OpenCV_LIBS} ${Qt5_LIBS} ${Boost_LIBRARIES})

add_executable(kmeans grid_kmeans.cc utils.cc)
target_link_libraries(kmeans ${OpenCV_LIBS} ${Qt5_LIBS} ${Boost_LIBRARIES})

add_executable(recolor recolor.cc palette_based_recolor.cc utils.cc)
target_link_libraries(recolor ${OpenCV_LIBS} ${Qt5_LIBS} ${Boost_LIBRARIES})

add_executable(recolor_one4one recolor_one4one.cc palette_based_recolor.cc utils.cc)
target_link_libraries(recolor_one4one ${OpenCV_LIBS} ${Qt5_LIBS} ${Boost_LIBRARIES})

add_executable(test_opencv test_opencv.cc)
target_link_libraries(test_opencv ${OpenCV_LIBS} ${Qt5_LIBS} ${Boost_LIBRARIES})