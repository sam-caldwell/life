cmake_minimum_required(VERSION 3.16)
project(life VERSION 0.2.0 LANGUAGES CXX)

# Use a widely available standard for portability
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(Curses REQUIRED)
find_package(Threads REQUIRED)

add_executable(life
    cmd/life/main.cpp
    src/Board.cpp
    src/Automaton.cpp
    src/Logger.cpp
)

target_include_directories(life PRIVATE include)

target_compile_options(life PRIVATE
    -Wall -Wextra -Wpedantic
)

target_link_libraries(life PRIVATE
    ${CURSES_LIBRARIES}
    Threads::Threads
)

target_compile_definitions(life PRIVATE
    NCURSES_WIDECHAR=1
)

if (CURSES_INCLUDE_DIR)
    target_include_directories(life PRIVATE ${CURSES_INCLUDE_DIR})
endif()

# Physics variant (Newtonian particles)
add_executable(physics
    cmd/physics/main.cpp
    src/PhysicsWorld.cpp
    src/Logger.cpp
)

target_include_directories(physics PRIVATE include)

target_compile_options(physics PRIVATE
    -Wall -Wextra -Wpedantic
)

target_link_libraries(physics PRIVATE
    ${CURSES_LIBRARIES}
)

target_compile_definitions(physics PRIVATE
    NCURSES_WIDECHAR=1
)

if (CURSES_INCLUDE_DIR)
    target_include_directories(physics PRIVATE ${CURSES_INCLUDE_DIR})
endif()
