cmake_minimum_required(VERSION 3.10)

project(KarameloDoc)

set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

project(KarameloDoc)

if(CMAKE_COMPILER_IS_GNUCXX)
  if(CMAKE_SYSTEM_NAME MATCHES Linux)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1 -g1")
  endif()
endif()

#option(KARAMELO_INTERNAL_DOCUMENTATION "Build internal documentation" OFF)


find_package(Git)
if(Git_FOUND)
  message("Git found: ${GIT_EXECUTABLE}")
endif()
# the commit's SHA1, and whether the building workspace was dirty or not
execute_process(COMMAND
  "${GIT_EXECUTABLE}" describe --match=NeVeRmAtCh --always --abbrev=8
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_SHA1
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)


# Set some Doxygen flags
set(KARAMELO_VERSION                       "${GIT_SHA1}")
message(STATUS "Version: ${KARAMELO_VERSION}")
set(Karamelo_SOURCE_DIR ".")
set(Karamelo_BINARY_DIR ".")
set(KARAMELO_DOXY_PROJECT_NAME             "Karamelo")
set(KARAMELO_DOXY_OUTPUT_DIRECTORY_SUFFIX  "")
set(KARAMELO_DOXY_INPUT                    "\"${Karamelo_SOURCE_DIR}/src\" \"${Karamelo_SOURCE_DIR}/doc\"")
set(KARAMELO_DOXY_HTML_COLORSTYLE_HUE      "220")
set(KARAMELO_DOXY_TAGFILES                 "")
if(KARAMELO_INTERNAL_DOCUMENTATION)
  set(KARAMELO_DOXY_INTERNAL                 "YES")
else()
  set(KARAMELO_DOXY_INTERNAL                 "NO")
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
  ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
)

#set(KARAMELO_DOXY_PROJECT_NAME             "Karamelo-unsupported")
#set(KARAMELO_DOXY_OUTPUT_DIRECTORY_SUFFIX  "/unsupported")
#set(KARAMELO_DOXY_INPUT                    "\"${Karamelo_SOURCE_DIR}/unsupported/Karamelo\" \"${Karamelo_SOURCE_DIR}/unsupported/doc\"")
#set(KARAMELO_DOXY_HTML_COLORSTYLE_HUE      "0")
#set(KARAMELO_DOXY_TAGFILES                 "\"${Karamelo_BINARY_DIR}/doc/Karamelo.doxytags=..\"")
#set(KARAMELO_DOXY_TAGFILES                 "")

#configure_file(
#  ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
#  ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile-unsupported
#)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/karamelodoxy_header.html.in
  ${CMAKE_CURRENT_BINARY_DIR}/karamelodoxy_header.html
)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/karamelodoxy_footer.html.in
  ${CMAKE_CURRENT_BINARY_DIR}/karamelodoxy_footer.html
)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/karamelodoxy_layout.xml.in
  ${CMAKE_CURRENT_BINARY_DIR}/karamelodoxy_layout.xml
)


set(examples_targets "")
set(snippets_targets "")

add_definitions("-DKARAMELO_MAKING_DOCS")
add_custom_target(all_examples)

#add_subdirectory(examples)
#add_subdirectory(special_examples)
#add_subdirectory(snippets)

add_custom_target(
  doc-karamelo-prerequisites
  ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/html/
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/karamelo_navtree_hacks.js           ${CMAKE_CURRENT_BINARY_DIR}/html/
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/Karamelo_Silly_Professor_64x64.png  ${CMAKE_CURRENT_BINARY_DIR}/html/
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/ftv2pnode.png                    ${CMAKE_CURRENT_BINARY_DIR}/html/
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/ftv2node.png                     ${CMAKE_CURRENT_BINARY_DIR}/html/
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/AsciiQuickReference.txt          ${CMAKE_CURRENT_BINARY_DIR}/html/
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_custom_target(
  doc-unsupported-prerequisites
  ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory ${Karamelo_BINARY_DIR}/doc/html/unsupported
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/karamelo_navtree_hacks.js           ${CMAKE_CURRENT_BINARY_DIR}/html/unsupported/
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/Karamelo_Silly_Professor_64x64.png  ${CMAKE_CURRENT_BINARY_DIR}/html/unsupported/
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/ftv2pnode.png                    ${CMAKE_CURRENT_BINARY_DIR}/html/unsupported/
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/ftv2node.png                     ${CMAKE_CURRENT_BINARY_DIR}/html/unsupported/
  WORKING_DIRECTORY ${Karamelo_BINARY_DIR}/doc
)

#add_dependencies(doc-karamelo-prerequisites all_snippets all_examples)
#add_dependencies(doc-unsupported-prerequisites unsupported_snippets unsupported_examples)

add_custom_target(doc ALL
  COMMAND doxygen
  COMMAND doxygen Doxyfile-unsupported
  COMMAND ${CMAKE_COMMAND} -E copy ${Karamelo_BINARY_DIR}/doc/html/group__TopicUnalignedArrayAssert.html ${Karamelo_BINARY_DIR}/doc/html/TopicUnalignedArrayAssert.html
  COMMAND ${CMAKE_COMMAND} -E rename html karamelo-doc
  COMMAND ${CMAKE_COMMAND} -E remove karamelo-doc/karamelo-doc.tgz karamelo-doc/unsupported/_formulas.log karamelo-doc/_formulas.log
  COMMAND ${CMAKE_COMMAND} -E tar cfz karamelo-doc.tgz karamelo-doc
  COMMAND ${CMAKE_COMMAND} -E rename karamelo-doc.tgz karamelo-doc/karamelo-doc.tgz
  COMMAND ${CMAKE_COMMAND} -E rename karamelo-doc html
  WORKING_DIRECTORY ${Karamelo_BINARY_DIR}/doc)

add_dependencies(doc doc-karamelo-prerequisites doc-unsupported-prerequisites)
