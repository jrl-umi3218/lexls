macro(add_lexls_test NAME)
  add_executable(${NAME} ${NAME}.cpp)
  target_link_libraries(${NAME} PUBLIC lexls)
  add_test(NAME ${NAME} COMMAND ${NAME} "${CMAKE_CURRENT_SOURCE_DIR}/${NAME}.dat")
endmacro()

add_lexls_test(test_01)
add_lexls_test(test_numerical_error)
