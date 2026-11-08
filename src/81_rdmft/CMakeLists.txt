#
# 81_rdmft
#
add_library(81_rdmft STATIC
  m_rdmft.F90
  )

set_target_properties(81_rdmft
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

target_include_directories(81_rdmft
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_BINARY_DIR}/modules
  )

target_compile_definitions(81_rdmft
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(81_rdmft
  PUBLIC
  abinit::10_defs
  abinit::14_hidewrite
  )

add_library(abinit::81_rdmft ALIAS 81_rdmft)
