PROJECT(skel2graph)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


SET(Required_Packages
  VTK
  ITK
)

# FOREACH(Package ${Required_Packages})
#   LOADPACKAGE(${Package})
# ENDFOREACH(Package)

FIND_PACKAGE ( ITK )
IF ( ITK_FOUND )
      INCLUDE( ${USE_ITK_FILE} )
ENDIF( ITK_FOUND )
FIND_PACKAGE ( VTK )
IF ( VTK_FOUND )
      INCLUDE( ${USE_VTK_FILE} )
ENDIF( VTK_FOUND )

SET(VTKGraphConversionsDir  "VTK-GraphConversions/")

INCLUDE_DIRECTORIES(
  "."
  ${VTKGraphConversionsDir}
  )

SET(Libraries
  # vtkFiltering
  # vtkGraphics
  # vtkImaging
  # vtkIO
  ${VTK_LIBRARIES} #this replaces naming of ITK-libs in itk >=4.0.0
  ## using ITK-libs only where necessary
  #${ITK_LIBRARIES} #this replaces naming of ITK-libs in itk >=4.0.0
)

SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

SET(CurrentExe "poly2graph")
ADD_EXECUTABLE(${CurrentExe} ${CurrentExe}.cxx 
  ${VTKGraphConversionsDir}/vtkPolyDataToGraph.cxx)
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "vo-skel2poly")
ADD_EXECUTABLE(${CurrentExe} ${CurrentExe}.cxx)
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries} ${ITK_LIBRARIES})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "shortest-loop-path")
ADD_EXECUTABLE(${CurrentExe} ${CurrentExe}.cxx)
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries} ${ITK_LIBRARIES})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)

SET(CurrentExe "fill_loops")
ADD_EXECUTABLE(${CurrentExe} ${CurrentExe}.cxx  ${VTKGraphConversionsDir}/vtkPolyDataToGraph.cxx)
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries} ${ITK_LIBRARIES})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)



# ##ToDo
# -add tests for each file in e.g. tests/vo-skel2poly/auto/
# -- testing if vo-skel2poly results in the same EPC as tests/scripts/euler3D.m
# --- add test for correct execution of tests/scripts/euler3D.m at the beginning

