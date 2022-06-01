add_library(CalDict SHARED)
target_include_directories(CalDict
	PUBLIC ${CMAKE_CURRENT_SRC_DIR}
	SYSTEM PUBLIC ${ROOT_INCLUDE_DIRS}
	)
ROOT_GENERATE_DICTIONARY(cal_dict DataStructs.h
	LINKDEF LinkDef_CalDict.h
	MODULE CalDict
	)
target_link_libraries(CalDict
	${ROOT_LIBRARIES}
	)
target_sources(CalDict PRIVATE DataStructs.cpp DataStructs.h)
set_target_properties(CalDict PROPERTIES
	LIBRARY_OUTPUT_DIRECTORY ${SABRERECON_LIBRARY_DIR}
	)
add_custom_command(TARGET CalDict POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy
	${CMAKE_CURRENT_BINARY_DIR}/libCalDict_rdict.pcm
	${SABRERECON_LIBRARY_DIR}/libCalDict_rdict.pcm
	)
