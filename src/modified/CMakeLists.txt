#############################################################
# MODIFIED EHRENFEST MODEL
#############################################################
add_executable(modified.out ehrenfest_modified.c ../functions.c)
target_compile_options(modified.out PRIVATE -Wall -O3 -march=native -ffast-math)
target_link_libraries(modified.out compphy-lib blas gsl)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plots_modified.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)

add_custom_target(modified
    COMMAND modified.out && python plots_modified.py 2>&1 | grep -v 'libGL'
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)
add_custom_target(run-modified
    COMMAND modified.out
    DEPENDS modified.out
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)
add_custom_target(plot-modified
    COMMAND python plots_modified.py 2>&1 | grep -v 'libGL'
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

