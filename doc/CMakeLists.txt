add_custom_target(doc ALL)

if(INSTALL_HTML_DOCUMENTATION)

  find_program(DOXYGEN doxygen)

  if(DOXYGEN)
    configure_file(Doxyfile.in "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile")
    add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/doxygen-html/index.html"
      COMMAND ${DOXYGEN} Doxyfile
      DEPENDS ${LEXLS_HEADERS} "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile" "${CMAKE_CURRENT_SOURCE_DIR}/mainpage.dox" "${CMAKE_CURRENT_SOURCE_DIR}/newcommand.js"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )
    add_custom_target(doc-html DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/doxygen-html/index.html")
    add_dependencies(doc doc-html)
    install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/doxygen-html" DESTINATION share/doc/lexls)
  endif()

endif()

if(INSTALL_PDF_DOCUMENTATION)

  find_program(PDFLATEX pdflatex)

  if(PDFLATEX)
    set(PDF_SRC "${CMAKE_CURRENT_SOURCE_DIR}/hot_start.tex")
    set(PDF_OUT "${CMAKE_CURRENT_BINARY_DIR}/hot_start.pdf")
    add_custom_command(OUTPUT ${PDF_OUT}
      COMMAND ${PDFLATEX} ${PDF_SRC}
      COMMAND ${PDFLATEX} ${PDF_SRC}
      DEPENDS ${PDF_SRC}
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )
    add_custom_target(doc-pdf DEPENDS ${PDF_OUT})
    add_dependencies(doc doc-pdf)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hot_start.pdf" DESTINATION share/doc/lexls)
  endif()

endif()
