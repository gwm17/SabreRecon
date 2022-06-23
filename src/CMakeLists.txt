add_subdirectory(CalDict)

add_executable(SabreRecon)
target_include_directories(SabreRecon
	PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
	PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/CalDict
	SYSTEM PUBLIC ${ROOT_INCLUDE_DIRS}
	)
target_sources(SabreRecon PUBLIC
	CutHandler.h
	CutHandler.cpp
	Histogrammer.h
	Histogrammer.cpp
	Reconstructor.h
	Reconstructor.cpp
	MassLookup.h
	MassLookup.cpp
	RandomGenerator.h
	RandomGenerator.cpp
	Detectors/SabreDetector.h
	Detectors/SabreDetector.cpp
	Detectors/FocalPlaneDetector.h
	Detectors/FocalPlaneDetector.cpp
	EnergyLoss/EnergyLossConstants.h
	EnergyLoss/EnergyLoss.h
	EnergyLoss/EnergyLoss.cpp
	EnergyLoss/Target.h
	EnergyLoss/Target.cpp
	EnergyLoss/CubicSpline.h
	EnergyLoss/CubicSpline.cpp
	EnergyLoss/PunchTable.h
	EnergyLoss/PunchTable.cpp
	EnergyLoss/ElossTable.h
	EnergyLoss/ElossTable.cpp
	main.cpp
	)
target_link_libraries(SabreRecon
	CalDict
	catima
	${ROOT_LIBRARIES}
	)
set_target_properties(SabreRecon PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${SABRERECON_BINARY_DIR}
	)
	
