# Install script for directory: /home/louisapawusch/NumSim/Git/NumSim/test/louisa/SettingsParser/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/louisapawusch/NumSim/Git/NumSim/test/louisa/SettingsParser/src/../build/numsim" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/louisapawusch/NumSim/Git/NumSim/test/louisa/SettingsParser/src/../build/numsim")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/louisapawusch/NumSim/Git/NumSim/test/louisa/SettingsParser/src/../build/numsim"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/louisapawusch/NumSim/Git/NumSim/test/louisa/SettingsParser/src/../build/numsim")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/louisapawusch/NumSim/Git/NumSim/test/louisa/SettingsParser/src/../build" TYPE EXECUTABLE FILES "/home/louisapawusch/NumSim/Git/NumSim/test/louisa/SettingsParser/build/src/numsim")
  if(EXISTS "$ENV{DESTDIR}/home/louisapawusch/NumSim/Git/NumSim/test/louisa/SettingsParser/src/../build/numsim" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/louisapawusch/NumSim/Git/NumSim/test/louisa/SettingsParser/src/../build/numsim")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/louisapawusch/NumSim/Git/NumSim/test/louisa/SettingsParser/src/../build/numsim"
         OLD_RPATH "/usr/lib/x86_64-linux-gnu/hdf5/openmpi:/usr/lib/x86_64-linux-gnu/openmpi/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/louisapawusch/NumSim/Git/NumSim/test/louisa/SettingsParser/src/../build/numsim")
    endif()
  endif()
endif()

