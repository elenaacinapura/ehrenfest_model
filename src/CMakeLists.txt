add_executable(basic ehrenfest_basic.c functions.c)
add_executable(modified ehrenfest_modified.c functions.c)

target_compile_options(basic PRIVATE -Wall -O3 -march=native -ffast-math)
target_compile_options(modified PRIVATE -Wall -O3 -march=native -ffast-math)

target_link_libraries(basic compphy-lib blas gsl)
target_link_libraries(modified compphy-lib blas gsl)

add_custom_target(run-basic
    COMMAND basic
    DEPENDS basic
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)
add_custom_target(run-modified
    COMMAND modified
    DEPENDS modified
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plots_basic.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plots_modified.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(plot-basic
    COMMAND python plots_basic.py
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)
add_custom_target(plot-modified
    COMMAND python plots_modified.py
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

add_custom_target(run-plot-modified
    COMMAND modified && python plots_modified.py
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

add_custom_target(run-plot-basic
    COMMAND basic && python plots_basic.py
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)
