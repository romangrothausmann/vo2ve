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

SET(CurrentExe "skel2graph")
ADD_EXECUTABLE(${CurrentExe} ${CurrentExe}.cxx 
  ${VTKGraphConversionsDir}/vtkPolyDataToGraph.cxx)
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
install (TARGETS ${CurrentExe} RUNTIME DESTINATION bin)
