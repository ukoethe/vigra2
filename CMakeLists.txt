cmake_minimum_required(VERSION 2.6)

##################################################
#
#     basic initializations
#
##################################################

PROJECT(muli)

# read the current version from config_version.hxx
file(READ ${muli_SOURCE_DIR}/include/muli/config_version.hxx MULI_VERSION_FILE)
string(REGEX MATCH "MULI_VERSION_MAJOR[ \t\n]+[^ \t\n]+" MULI_VERSION_MAJOR ${MULI_VERSION_FILE})
string(REGEX REPLACE "MULI_VERSION_MAJOR[ \t\n]" "" MULI_VERSION_MAJOR ${MULI_VERSION_MAJOR})
string(REGEX MATCH "MULI_VERSION_MINOR[ \t\n]+[^ \t\n]+" MULI_VERSION_MINOR ${MULI_VERSION_FILE})
string(REGEX REPLACE "MULI_VERSION_MINOR[ \t\n]" "" MULI_VERSION_MINOR ${MULI_VERSION_MINOR})
string(REGEX MATCH "MULI_VERSION_PATCH[ \t\n]+[^ \t\n]+" MULI_VERSION_PATCH ${MULI_VERSION_FILE})
string(REGEX REPLACE "MULI_VERSION_PATCH[ \t\n]" "" MULI_VERSION_PATCH ${MULI_VERSION_PATCH})

SET(muli_version_short ${MULI_VERSION_MAJOR}.${MULI_VERSION_MINOR})
SET(muli_version ${MULI_VERSION_MAJOR}.${MULI_VERSION_MINOR}.${MULI_VERSION_PATCH})
message(STATUS "Configuring MULI version ${muli_version}")

# SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${muli_SOURCE_DIR}/config)

# include(VigraSetDefaults)
# include(VigraCMakeUtils)

INCLUDE_DIRECTORIES(${muli_SOURCE_DIR}/include)

IF(MULI_STATIC_LIB)
    SET(LIBTYPE STATIC)
    ADD_DEFINITIONS(-DMULI_STATIC_LIB)
ELSE(MULI_STATIC_LIB)
    SET(LIBTYPE SHARED)
ENDIF(MULI_STATIC_LIB)

IF (MSVC)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE -D_SCL_SECURE_NO_DEPRECATE)
ENDIF ()

# IF("${CMAKE_CXX_COMPILER_ID}" MATCHES "(GNU|Clang)")
    # set(SUFFICIENT_TEMPLATE_DEPTH FALSE)
    # TRY_COMPILE(SUFFICIENT_TEMPLATE_DEPTH
                # ${CMAKE_BINARY_DIR} ${PROJECT_SOURCE_DIR}/config/checkTemplateDepth.cxx
                # COMPILE_DEFINITIONS "-DDEPTH=900") 
            
    # IF(SUFFICIENT_TEMPLATE_DEPTH)
        # MESSAGE(STATUS "Checking template recursion depth: ok")
    # ELSE()
        # if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
            # # As of cmake 2.8.10, there is a variable CMAKE_CXX_COMPILER_VERSION.
            # # Emulate it by a compiler call if not present.
            # IF(NOT CMAKE_CXX_COMPILER_VERSION)
                # exec_program(${CMAKE_CXX_COMPILER}
                    # ARGS ${CMAKE_CXX_COMPILER_ARG1} -dumpversion
                    # OUTPUT_VARIABLE CMAKE_CXX_COMPILER_VERSION)
                # string(REGEX REPLACE "([0-9])\\.([0-9])(\\.[0-9])?" "\\1.\\2"
                      # CMAKE_CXX_COMPILER_VERSION ${CMAKE_CXX_COMPILER_VERSION})
            # endif()

            # message ("Detected GCC version ${CMAKE_CXX_COMPILER_VERSION}")
            # if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.4.0")
                # MESSAGE(WARNING "GNU c++ < 4.4 cannot build MULINumPy; disabling (found GCC ${CMAKE_CXX_COMPILER_VERSION})")
                # SET(WITH_MULINUMPY 0)
            # elseif(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.5.0")
                # # gcc 4.4.x syntax is -ftemplate-depth-900 (not =900)
                # SET(CMAKE_CXX_FLAGS "-ftemplate-depth-900 ${CMAKE_CXX_FLAGS}" CACHE STRING  "Flags used by the compiler during all build types" FORCE)        
                # MESSAGE(STATUS "Checking template recursion depth: using -ftemplate-depth-900")
            # else()
                # # new syntax is =900, not -900
                # SET(CMAKE_CXX_FLAGS "-ftemplate-depth=900 ${CMAKE_CXX_FLAGS}" CACHE STRING  "Flags used by the compiler during all build types" FORCE)        
                # MESSAGE(STATUS "Checking template recursion depth: using -ftemplate-depth=900")
            # endif()
        # else()
            # # clang uses newer gcc syntax
            # SET(CMAKE_CXX_FLAGS "-ftemplate-depth=900 ${CMAKE_CXX_FLAGS}" CACHE STRING  "Flags used by the compiler during all build types" FORCE)        
            # MESSAGE(STATUS "Checking template recursion depth: using -ftemplate-depth=900")
        # endif()
    # ENDIF()
# ENDIF()

# ##################################################
# #
# #     search for dependencies
# #
# ##################################################

# INCLUDE(VigraFindPackage)
# MULI_FIND_PACKAGE(ZLIB)
# MULI_FIND_PACKAGE(TIFF NAMES libtiff_i libtiff) # prefer DLL on Windows
# MULI_FIND_PACKAGE(JPEG NAMES libjpeg)
# MULI_FIND_PACKAGE(PNG)
# MULI_FIND_PACKAGE(FFTW3 NAMES libfftw3-3 libfftw-3.3)
# MULI_FIND_PACKAGE(FFTW3F NAMES libfftw3f-3 libfftwf-3.3)

# IF(WITH_OPENEXR)
    # MULI_FIND_PACKAGE(OpenEXR)
# ENDIF()

# IF(WITH_HDF5)
    # MULI_FIND_PACKAGE(HDF5)
# ENDIF()

# set(WITH_BOOST "OFF")
# set(WITH_BOOST_COMPONENTS "")

# IF(WITH_BOOST_GRAPH)
    # set(WITH_BOOST "ON")
# ENDIF()

# IF(WITH_MULINUMPY)
    # set(WITH_BOOST "ON")
    # set(WITH_BOOST_COMPONENTS ${WITH_BOOST_COMPONENTS} python)
# ENDIF()

# IF(WITH_BOOST_THREAD)
    # set(WITH_BOOST "ON")
    # set(WITH_BOOST_COMPONENTS ${WITH_BOOST_COMPONENTS} thread system date_time chrono)
    # ADD_DEFINITIONS(-DUSE_BOOST_THREAD)
# ENDIF()

# IF(WITH_BOOST)
    # MULI_FIND_PACKAGE( Boost 1.40.0 COMPONENTS ${WITH_BOOST_COMPONENTS})
    
    # if(Boost_FOUND AND WITH_BOOST_COMPONENTS)
        # # configure boost's autolink magic to use the right library name
        # # (default on Windows is a mangled name like 'boost_python-vc110-mt-1_51.lib')
        # if(("${Boost_PYTHON_LIBRARY}" MATCHES "boost_python\\.lib") OR
           # ("${Boost_SYSTEM_LIBRARY}" MATCHES "boost_system\\.lib"))
            # ADD_DEFINITIONS(-DBOOST_AUTO_LINK_NOMANGLE)
        # endif()
    # endif()
# ENDIF()

# IF(WITH_LEMON)
    # MULI_FIND_PACKAGE(LEMON)
# ENDIF()

# SET(DOXYGEN_SKIP_DOT TRUE)
# FIND_PACKAGE(Doxygen)
# FIND_PACKAGE(PythonInterp 2)

# IF(WITH_MULINUMPY)
    # FIND_PACKAGE( MULINUMPY_DEPENDENCIES )
# ENDIF()

# IF(WITH_VALGRIND)
    # FIND_PROGRAM(VALGRIND_EXECUTABLE valgrind)
# ENDIF()

# # Must be included AFTER boost is found.
# include(VigraDetectThreading)
# include(VigraConfigureThreading)

# include(VigraDetectCppVersion)
# MULI_DETECT_CPP_VERSION()

# ##################################################
# #
# #     setup testing environment
# #
# ##################################################

# ADD_CUSTOM_TARGET(check)
# ADD_CUSTOM_TARGET(ctest COMMAND ${CMAKE_CTEST_COMMAND})

# ADD_CUSTOM_TARGET(check_cpp)
# ADD_DEPENDENCIES(check check_cpp)


# OPTION(CREATE_CTEST_TARGETS "Create special targets for CTest" OFF)
# IF(CREATE_CTEST_TARGETS)
  # INCLUDE(CTest)
# ENDIF()

# IF(WITH_MULINUMPY)
    # ADD_CUSTOM_TARGET(check_python)
    # ADD_DEPENDENCIES(check check_python)
# ENDIF()

# ENABLE_TESTING()

# ##################################################
# #
# #     setup documentation generator
# #
# ##################################################

# ADD_CUSTOM_TARGET(doc)

# # Automatically push the latest documentation to github's gh-pages.
# # This only works when ${DOCDIR} refers to a MULI repository that is 
# # checked out in branch 'gh-pages', and the present repository is 
# # in branch 'master' (these conditions are checked in the script).
# if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

    # configure_file(${PROJECT_SOURCE_DIR}/config/update-gh-pages.sh.in 
                   # ${PROJECT_BINARY_DIR}/update-gh-pages.sh @ONLY)
    # EXECUTE_PROCESS(COMMAND chmod u+x ${PROJECT_BINARY_DIR}/update-gh-pages.sh OUTPUT_QUIET ERROR_QUIET)
       
    # ADD_CUSTOM_TARGET(gh-pages 
          # ${PROJECT_BINARY_DIR}/update-gh-pages.sh
          # DEPENDS doc
          # WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
          # COMMENT "Pushing documentation to github pages.")
# endif()

# ##################################################
# #
# #     setup private experimentation target
# #
# ##################################################

# ADD_CUSTOM_TARGET(experiments)

##################################################
#
#     subdirectories
#
##################################################

# ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(test)
# ADD_SUBDIRECTORY(docsrc)

# IF(WITH_MULINUMPY)
    # ADD_SUBDIRECTORY(mulinumpy)
# ENDIF()

# IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/private/CMakeLists.txt)
    # ADD_SUBDIRECTORY(private EXCLUDE_FROM_ALL)
# ENDIF()


# set(LIB_SUFFIX "" CACHE STRING "Define suffix of lib directory name (32/64)" )

# ##################################################
# #
# #     global installation commands
# #
# ##################################################

# INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/include/muli
        # DESTINATION include)

# ##################################################
# #
# #     enable CPack package generation
# #
# ##################################################
# INCLUDE(CPackConfig)

# ##################################################
# #
# #     add target muli-$version-src.tar.gz
# #
# ##################################################

# find_program(GIT_EXECUTABLE 
              # NAMES git git.exe git.cmd
              # HINTS $ENV{ProgramFiles}/Git/bin
              # DOC "git command line client")

# CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/config/package-src.cmake.in
                # ${PROJECT_BINARY_DIR}/package-src.cmake
                # @ONLY IMMEDIATE)

# add_custom_target(PACKAGE_SRC_TAR
                   # COMMAND ${CMAKE_COMMAND} -P package-src.cmake
                   # WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
                   # COMMENT "Creating ${PROJECT_BINARY_DIR}/muli-${muli_version}-src.tar.gz")        

# ADD_DEPENDENCIES(PACKAGE_SRC_TAR check)
# ADD_DEPENDENCIES(PACKAGE_SRC_TAR doc_cpp)
# IF(WITH_MULINUMPY AND PYTHON_SPHINX)
    # ADD_DEPENDENCIES(PACKAGE_SRC_TAR doc_python)
# ENDIF()

# ##################################################
# #
# #  enable other non-CMake projects to find Vigra
# #
# ##################################################

# IF(NOT WIN32)
# # muli-config script for (non-CMake-based) third party packages
# CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/config/muli-config.in
        # ${PROJECT_BINARY_DIR}/bin/muli-config
        # @ONLY IMMEDIATE)
# INSTALL(FILES ${PROJECT_BINARY_DIR}/bin/muli-config
        # DESTINATION bin
        # PERMISSIONS
            # OWNER_READ OWNER_WRITE OWNER_EXECUTE
            # GROUP_READ GROUP_EXECUTE
            # WORLD_READ WORLD_EXECUTE)
# ENDIF(NOT WIN32)

# ##################################################
# #
# #     enable other CMake projects to find Vigra
# #
# ##################################################

# # config file for CMake FIND_PACKAGE command (for using the build tree)
# CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/config/VigraConfigBuild.cmake.in
        # ${PROJECT_BINARY_DIR}/VigraConfig.cmake
        # @ONLY IMMEDIATE)
# CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/config/VigraConfigVersion.cmake.in
        # ${PROJECT_BINARY_DIR}/VigraConfigVersion.cmake
        # @ONLY IMMEDIATE)

# # config file for CMake FIND_PACKAGE command
# CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/config/VigraConfig.cmake.in
        # ${PROJECT_BINARY_DIR}/lib/muli/CMake/VigraConfig.cmake
        # @ONLY IMMEDIATE)
# CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/config/VigraConfigVersion.cmake.in
        # ${PROJECT_BINARY_DIR}/lib/muli/CMake/VigraConfigVersion.cmake
        # @ONLY IMMEDIATE)

# # export targets:
# INSTALL(EXPORT muli-targets
        # DESTINATION lib${LIB_SUFFIX}/muli)
# INSTALL(FILES ${PROJECT_BINARY_DIR}/lib/muli/CMake/VigraConfig.cmake
            # ${PROJECT_BINARY_DIR}/lib/muli/CMake/VigraConfigVersion.cmake
        # DESTINATION lib${LIB_SUFFIX}/muli)
# EXPORT(TARGETS muliimpex FILE muli-targets.cmake)

# ##################################################
# #
# #     status output
# #
# ##################################################

# MESSAGE( STATUS "---------------------------------------------------------" )
# MESSAGE( STATUS "MULI configuration information:" )
# MESSAGE( STATUS "---------------------------------------------------------" )

# IF(ZLIB_FOUND)
    # MESSAGE( STATUS "  Using ZLIB  libraries: ${ZLIB_LIBRARIES}" )
# ELSE()
    # MESSAGE( STATUS "  ZLIB libraries not found (ZLIB support disabled)" )
# ENDIF()

# IF(PNG_FOUND)
    # MESSAGE( STATUS "  Using PNG  libraries: ${PNG_LIBRARIES}" )
# ELSE()
    # MESSAGE( STATUS "  PNG libraries not found (PNG support disabled)" )
# ENDIF()

# IF(TIFF_FOUND)
    # MESSAGE( STATUS "  Using TIFF libraries: ${TIFF_LIBRARIES}" )
# ELSE()
    # MESSAGE( STATUS "  TIFF libraries not found (TIFF support disabled)" )
# ENDIF()

# IF(JPEG_FOUND)
    # MESSAGE( STATUS "  Using JPEG libraries: ${JPEG_LIBRARIES}" )
# ELSE()
    # MESSAGE( STATUS "  JPEG libraries not found (JPEG support disabled)" )
# ENDIF()

# IF(OPENEXR_FOUND)
    # MESSAGE( STATUS "  Using OpenEXR  libraries: ${OPENEXR_LIBRARIES}" )
# ELSEIF(NOT WITH_OPENEXR)
    # MESSAGE( STATUS "  OpenEXR disabled by user (WITH_OPENEXR=0)" )
# ELSE()
    # MESSAGE( STATUS "  OpenEXR libraries not found (OpenEXR support disabled)" )
# ENDIF()

# IF(FFTW3_FOUND)
    # MESSAGE( STATUS "  Using FFTW libraries: ${FFTW3_LIBRARIES}" )
# ELSE()
    # MESSAGE( STATUS "  FFTW libraries not found (FFTW support disabled)" )
# ENDIF()

# IF(HDF5_FOUND)
    # MESSAGE( STATUS "  Using HDF5 libraries: ${HDF5_LIBRARIES}" )
# ELSEIF(NOT WITH_HDF5)
    # MESSAGE( STATUS "  HDF5 disabled by user (WITH_HDF5=0)" )
# ELSE()
    # MESSAGE( STATUS "  HDF5 libraries not found (HDF5 support disabled)" )
# ENDIF()

# IF(Boost_FOUND)
    # MESSAGE( STATUS "  Using Boost Graph Library: ${Boost_INCLUDE_DIR}/boost/graph" )
# ELSEIF(NOT WITH_BOOST_GRAPH)
    # MESSAGE( STATUS "  Boost Graph Library disabled by user (WITH_BOOST_GRAPH=0)" )
# ELSE()
    # MESSAGE( STATUS "  Boost Graph Library not found (support disabled)" )
# ENDIF()

# IF(LEMON_FOUND)
    # MESSAGE( STATUS "  Using LEMON graph library: ${LEMON_LIBRARY}" )
# ELSEIF(NOT WITH_LEMON)
    # MESSAGE( STATUS "  LEMON graph library disabled by user (WITH_LEMON=0)" )
# ELSE()
    # MESSAGE( STATUS "  LEMON graph library not found (support disabled)" )
# ENDIF()

# IF(WITH_MULINUMPY)
    # IF(MULINUMPY_DEPENDENCIES_FOUND)
        # MESSAGE( STATUS "  Using Python libraries: ${MULINUMPY_LIBRARIES}" )
        # MESSAGE( STATUS "  Using Numpy includes:   ${PYTHON_NUMPY_INCLUDE_DIR}" )
    # ELSE()
        # MESSAGE( STATUS "  Vigranumpy dependencies not found (mulinumpy disabled)" )
    # ENDIF()
# ELSE()
  # MESSAGE( STATUS "  Vigranumpy disabled by user (WITH_MULINUMPY=0)" )
# ENDIF()

# MESSAGE( STATUS "---------------------------------------------------------" )

# IF(MULI_STATIC_LIB)
    # MESSAGE( STATUS "  building static lib" )
# ELSE()
    # MESSAGE( STATUS "  building shared lib" )
# ENDIF()

# MESSAGE( STATUS "  binaries will be generated in: ${CMAKE_CURRENT_BINARY_DIR}")
# IF(DOXYGEN_FOUND AND PYTHONINTERP_FOUND)
    # MESSAGE( STATUS "  manuals will be generated in: ${DOCDIR}")
# ENDIF()

# MESSAGE( STATUS "---------------------------------------------------------" )

# MESSAGE( STATUS "  includes will be installed at: ${CMAKE_INSTALL_PREFIX}/include")
# MESSAGE( STATUS "  libraries will be installed at: ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}")

# IF (MSVC)
  # MESSAGE( STATUS "  binaries and DLLs will be installed at: ${CMAKE_INSTALL_PREFIX}/bin")
# else()
  # MESSAGE( STATUS "  binaries will be installed at: ${CMAKE_INSTALL_PREFIX}/bin")
# endif()

# IF(DOXYGEN_FOUND AND PYTHONINTERP_FOUND)
  # MESSAGE( STATUS "  muli manuals will be installed at: "
    # "${CMAKE_INSTALL_PREFIX}/${DOCINSTALL}/muli/index.html")
# ELSE()
  # MESSAGE( STATUS "  muli manuals cannot be generated (Doxygen or Python not found)")
# ENDIF()

# IF(WITH_MULINUMPY AND MULINUMPY_DEPENDENCIES_FOUND)
  # MESSAGE( STATUS "  mulinumpy will be installed at ${CMAKE_INSTALL_PREFIX}/${MULINUMPY_INSTALL_DIR}" )
  # IF(PYTHON_SPHINX)
    # MESSAGE( STATUS "  mulinumpy manuals will be installed at: "
      # "${CMAKE_INSTALL_PREFIX}/${DOCINSTALL}/mulinumpy/html/index.html")
  # ELSE(PYTHON_SPHINX)
    # MESSAGE( STATUS "  mulinumpy manuals cannot be generated "
      # "(Sphinx missing)")
  # ENDIF(PYTHON_SPHINX)
# ELSE()
  # MESSAGE( STATUS "  mulinumpy will NOT be installed" )
# ENDIF()

# IF(NOT TIFF_FOUND OR NOT JPEG_FOUND OR NOT PNG_FOUND OR NOT FFTW3_FOUND OR
   # NOT HDF5_FOUND OR NOT MULINUMPY_DEPENDENCIES_FOUND)
  # MESSAGE( STATUS "" )
  # MESSAGE( STATUS "Consider setting DEPENDENCY_SEARCH_PREFIX to find missing libraries" )
# ENDIF()

# MESSAGE( STATUS "---------------------------------------------------------" )
