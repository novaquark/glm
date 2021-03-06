# ---------------------------------------- #
#                                          #
#                   glm                    #
#                                          #
# ---------------------------------------- #
cmake_minimum_required(VERSION 2.8)
project(glm)

# Note:
#  This cmake file sets a specific GLM configuration
# -------------------------------------------------------------------------------------------------

# Source files helpers
# -------------------------------------------------------------------------------------------------
function(globrecurse_c_cpp_sources result folder)
	file(GLOB_RECURSE res
		"${folder}/*.c"
		"${folder}/*.cpp"
		"${folder}/*.cxx"
		"${folder}/*.h"
		"${folder}/*.hpp"
		"${folder}/*.hxx"
		"${folder}/*.inl")
	set(${result} ${res} PARENT_SCOPE)
endfunction()

function(group_sources RelDir)
	get_filename_component(RelDir ${RelDir}/void REALPATH)
	get_filename_component(RelDir ${RelDir} PATH)
	foreach(src ${ARGN})
		set(grname ${src})
		get_filename_component(grname ${grname} REALPATH)
		get_filename_component(grname ${grname} PATH)
		string(REPLACE "${RelDir}/" "" grname "${grname}")
		string(REPLACE "${RelDir}" "" grname "${grname}")
		string(REPLACE "/" "\\\\" grname "${grname}")
		source_group("${grname}" FILES ${src})
	endforeach()
endfunction()

# Folders
# -------------------------------------------------------------------------------------------------
set(GLM_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}")

# Sources
# -------------------------------------------------------------------------------------------------
globrecurse_c_cpp_sources(GLM_SRC_FILES ${GLM_ROOT_DIR}/glm)
group_sources(${GLM_ROOT_DIR}/glm ${GLM_SRC_FILES})

# Target to build
# -------------------------------------------------------------------------------------------------
add_library(glm INTERFACE)
target_include_directories(glm INTERFACE ${GLM_ROOT_DIR})
target_compile_definitions(glm
	INTERFACE GLM_FORCE_RADIANS
	INTERFACE GLM_FORCE_FLOAT_DETERMINISM)

# Per compiler config
# -------------------------------------------------------------------------------------------------
if(WIN32)
	# Enable SSE2
	target_compile_options(glm INTERFACE "/arch:SSE2")
endif()

if(MINGW OR COMPILER_IS_GNUCXX OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
	# Enable C++11 if gcc or clang
	target_compile_options(glm INTERFACE "-std=c++11")

	# Enable SSE2
	target_compile_options(glm INTERFACE "-msse2")
endif()
