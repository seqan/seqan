find_package (PythonInterp)

set(SEQAN_DOCS_PATH "${CMAKE_SOURCE_DIR}/docs")

add_custom_target(docs
    WORKING_DIRECTORY ${SEQAN_DOCS_PATH}
    COMMAND ${PYTHON_EXECUTABLE} ${SEQAN_DOCS_PATH}/main.py ${CMAKE_SOURCE_DIR}/core/include ${CMAKE_SOURCE_DIR}/extras/include -d ${SEQAN_DOCS_PATH}/concepts -d ${SEQAN_DOCS_PATH}/pages -o ${CMAKE_BINARY_DIR}/docs/html)

# You have to run "make docs" before calling make install.  We cannot
# add dependencies to the install target at the moment.  Also see:
# http://public.kitware.com/Bug/view.php?id=8438

install(DIRECTORY ${CMAKE_BINARY_DIR}/docs/html
        DESTINATION share/seqan/docs)
