add_executable(main main.c)

target_compile_options(main PRIVATE -Wall -O3 -march=native -ffast-math)

target_link_libraries(main compphy-lib blas)

add_custom_target(run
    COMMAND main
    DEPENDS main
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)