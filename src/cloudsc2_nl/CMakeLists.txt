# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC2_NL DEFAULT ON
    DESCRIPTION "Build the Fortran CLOUDSC2 non-linear dwarf"
)

if( HAVE_CLOUDSC2_NL )
   ecbuild_enable_fortran(MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
function(PRINT_VAR VARNAME)
  message(STATUS "${VARNAME}: ${${VARNAME}}")
endfunction()
PRINT_VAR("CMAKE_CURRENT_SOURCE_DIR")
PRINT_VAR("CMAKE_CURRENT_BINARY_DIR")

    # Define the binary build target for this variant
   #ecbuild_add_executable( TARGET dwarf-cloudsc2-nl
   #    SOURCES
   #        dwarf_cloudsc.F90
   #        cloudsc_driver_mod.F90
   #        satur.F90
   #        cuadjtqs.F90
   #        cloudsc2.F90
   #        yomphyder.F90
   #    DEFINITIONS ${CLOUDSC_DEFINITIONS}
   #)

    ecbuild_add_library( TARGET dwarf-cloudsc2-nl-lib
 #TYPE STATIC
        SOURCES
            cloudsc_driver_mod.F90
            satur.F90
            cuadjtqs.F90
            cloudsc2.F90
            yomphyder.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )

    target_link_libraries( dwarf-cloudsc2-nl-lib PRIVATE cloudsc-common-lib HDF5::HDF5 )
    ecbuild_add_executable( TARGET dwarf-cloudsc2-nl
        SOURCES
            dwarf_cloudsc.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    target_link_libraries( dwarf-cloudsc2-nl PRIVATE cloudsc-common-lib dwarf-cloudsc2-nl-lib )

    if( HAVE_OMP AND TARGET OpenMP::OpenMP_Fortran )
        target_link_libraries( dwarf-cloudsc2-nl PRIVATE OpenMP::OpenMP_Fortran )
    endif()

    if( HAVE_MPI AND TARGET MPI::MPI_Fortran )
        target_link_libraries( dwarf-cloudsc2-nl PRIVATE MPI::MPI_Fortran )
    endif()

    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/input.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../input.h5 )
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/reference.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../reference.h5 )

    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/input.h5 ${CMAKE_CURRENT_SOURCE_DIR}/input.h5 )
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/reference.h5 ${CMAKE_CURRENT_SOURCE_DIR}/reference.h5 )
     set(F90WRAP_SOURCES
             yomphyder.F90
             parkind1.F90
             cloudsc_driver_mod.F90
        )
     set(PYTHON_MODN Exampledwarf)
   add_custom_command(TARGET dwarf-cloudsc2-nl POST_BUILD
  # add_custom_command(OUTPUT ${F90WRAP_SOURCES} TARGET dwarf-cloudsc2-nl POST_BUILD
         COMMAND pwd
         COMMAND ls 
         COMMAND rm -f ${CMAKE_CURRENT_SOURCE_DIR}/libcloudsc-common-lib.dylib 
         COMMAND rm -f ${CMAKE_CURRENT_SOURCE_DIR}/libdwarf-cloudsc2-nl-lib.dylib 
         COMMAND ln -s 
${CMAKE_CURRENT_SOURCE_DIR}/../../build/lib/libcloudsc-common-lib.dylib 
${CMAKE_CURRENT_SOURCE_DIR}/libcloudsc-common-lib.dylib 
         COMMAND ln -s
${CMAKE_CURRENT_SOURCE_DIR}/../../build/lib/libdwarf-cloudsc2-nl-lib.dylib 
${CMAKE_CURRENT_SOURCE_DIR}/libdwarf-cloudsc2-nl-lib.dylib 
         COMMAND  f90wrap -m${PYTHON_MODN} ${F90WRAP_SOURCES} -k kind_map -P
         COMMAND  f2py-f90wrap   -c -m _${PYTHON_MODN}  -L. -lcloudsc-common-lib  -ldwarf-cloudsc2-nl-lib f90wrap_cloudsc_driver_mod.f90 f90wrap_yomphyder.f90  
 #f90wrap_parkind1.f90
         WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
         COMMENT "Running F90wrap"
         DEPENDS dwarf-cloudsc2-nl-lib 
         VERBATIM
   )    
#
endif()
