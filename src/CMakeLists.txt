add_executable(ehrenfest ehrenfest.c functions.c)
target_compile_options(ehrenfest PRIVATE -Wall -O3 -march=native -ffast-math)
target_link_libraries(ehrenfest compphy-lib blas gsl)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/ehrenfest.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
) 

add_custom_target(run-ehrenfest
    COMMAND ehrenfest && python ehrenfest.py 2>&1 | grep -v 'libGL'
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)
add_custom_target(plot-ehrenfest
    COMMAND python ehrenfest.py 2>&1 | grep -v 'libGL'
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)
