add_executable(main main.c)

target_compile_options(main PRIVATE -Wall -O3 -march=native -ffast-math)

add_custom_target(run
    COMMAND main
    DEPENDS main
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)