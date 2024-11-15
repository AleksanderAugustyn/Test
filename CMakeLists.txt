cmake_minimum_required(VERSION 3.30)
project(Test)

set(CMAKE_CXX_STANDARD 17)

# Add ROOT
list(APPEND CMAKE_PREFIX_PATH /usr/share/root/cmake)
find_package(ROOT REQUIRED)
# Add ROOT include directories
include(${ROOT_USE_FILE})
# Add ROOT libraries
link_libraries(${ROOT_LIBRARIES})
# Add ROOT definitions
add_definitions(${ROOT_DEFINITIONS})
# Add ROOT flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
# Print ROOT version
message(STATUS "ROOT version: ${ROOT_VERSION}")

add_executable(Test main.cpp
        TraceGraphs.cpp
        EventSelection.cpp
        main.h
        RootInput.cpp
        FitAnalysis.cpp
        AnalyseTraces.cpp)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -O2")

#[[if(MSVC)
    target_compile_options(ROOT_Test PRIVATE /W4 /WX /O2)
else()
    target_compile_options(ROOT_Test PRIVATE -Wall -Wextra -Wpedantic -Werror -O2)
endif()]]