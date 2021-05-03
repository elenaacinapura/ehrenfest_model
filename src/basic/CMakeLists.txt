#############################################################
# BASIC EHRENFEST MODEL
#############################################################
add_executable(basic.out ehrenfest_basic.c ../functions.c)
target_compile_options(basic.out PRIVATE -Wall -O3 -march=native -ffast-math)
target_link_libraries(basic.out compphy-lib blas gsl)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plots_basic.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
) 

add_custom_target(basic
    COMMAND basic.out && python plots_basic.py 2>&1 | grep -v 'libGL'
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)
add_custom_target(run-basic
    COMMAND basic.out
    DEPENDS basic.out
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)
add_custom_target(plot-basic
    COMMAND python plots_basic.py 2>&1 | grep -v 'libGL'
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)