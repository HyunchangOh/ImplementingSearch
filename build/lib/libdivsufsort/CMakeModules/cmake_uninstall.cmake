IF(NOT EXISTS "/mnt/c/Users/hyunc/simonhelpme/build/lib/libdivsufsort/install_manifest.txt")
  MESSAGE(FATAL_ERROR "Cannot find install manifest: \"/mnt/c/Users/hyunc/simonhelpme/build/lib/libdivsufsort/install_manifest.txt\"")
ENDIF(NOT EXISTS "/mnt/c/Users/hyunc/simonhelpme/build/lib/libdivsufsort/install_manifest.txt")

FILE(READ "/mnt/c/Users/hyunc/simonhelpme/build/lib/libdivsufsort/install_manifest.txt" files)
STRING(REGEX REPLACE "\n" ";" files "${files}")

SET(NUM 0)
FOREACH(file ${files})
  IF(EXISTS "$ENV{DESTDIR}${file}")
    MESSAGE(STATUS "Looking for \"$ENV{DESTDIR}${file}\" - found")
    SET(UNINSTALL_CHECK_${NUM} 1)
  ELSE(EXISTS "$ENV{DESTDIR}${file}")
    MESSAGE(STATUS "Looking for \"$ENV{DESTDIR}${file}\" - not found")
    SET(UNINSTALL_CHECK_${NUM} 0)
  ENDIF(EXISTS "$ENV{DESTDIR}${file}")
  MATH(EXPR NUM "1 + ${NUM}")
ENDFOREACH(file)

SET(NUM 0)
FOREACH(file ${files})
  IF(${UNINSTALL_CHECK_${NUM}})
    MESSAGE(STATUS "Uninstalling \"$ENV{DESTDIR}${file}\"")
    EXEC_PROGRAM(
      "/usr/bin/cmake" ARGS "-E remove \"$ENV{DESTDIR}${file}\""
      OUTPUT_VARIABLE rm_out
      RETURN_VALUE rm_retval
      )
    IF(NOT "${rm_retval}" STREQUAL 0)
      MESSAGE(FATAL_ERROR "Problem when removing \"$ENV{DESTDIR}${file}\"")
    ENDIF(NOT "${rm_retval}" STREQUAL 0)
  ENDIF(${UNINSTALL_CHECK_${NUM}})
  MATH(EXPR NUM "1 + ${NUM}")
ENDFOREACH(file)

FILE(REMOVE "/mnt/c/Users/hyunc/simonhelpme/build/lib/libdivsufsort/install_manifest.txt")
