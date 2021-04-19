add_executable(main main.c functions.c)

target_compile_options(main PRIVATE -Wall -O3 -march=native -ffast-math)

target_link_libraries(main compphy-lib blas gsl)

add_custom_target(run
    COMMAND main
    DEPENDS main
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/frequency.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(plot
    COMMAND python frequency.py
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/frequency.csv
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)