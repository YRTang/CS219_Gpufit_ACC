
function( add_example module name )
	# include_directories(${CMAKE_CURRENT_SOURCE_DIR})
	# file(GLOB SOURCES "*.cpp" "*.h") #added
	# add_executable( ${name} ${SOURCES}) # add_executable( ${name} ${name}.cpp) #original
	add_executable( ${name} ${name}.cpp )
	target_link_libraries( ${name} ${module} )
	
	set_property( TARGET ${name}
		PROPERTY RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}" )
	set_property( TARGET ${name} PROPERTY FOLDER GpufitExamples )
endfunction()

function( add_cuda_example module name )
	cuda_add_executable( ${name} ${name}.cu )
	target_link_libraries( ${name} ${module} )
	set_property( TARGET ${name}
		PROPERTY RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}" )
	set_property( TARGET ${name} PROPERTY FOLDER GpufitExamples )
endfunction()

# Examples
# add_example( Gpufit mp_solver_example)
add_example( Gpufit mp_solver_example_timed)
