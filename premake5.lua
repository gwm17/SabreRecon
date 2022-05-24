workspace "SabreRecon"
	architecture "x64"
	configurations {
		"Release",
		"Debug"
	}

ROOTIncludeDir = "/usr/include/root/"
ROOTLibDir = "/usr/lib64/root/"

project "CalDict"
	kind "SharedLib"
	language "C++"
	cppdialect "c++11"
	targetdir "./lib/"
	objdir "./objs/"

	prebuildcommands {
		"rootcint -f src/CalDict/cal_dict.cxx src/CalDict/DataStructs.h src/CalDict/LinkDef_CalDict.h",
		"{COPY} src/CalDict/*.pcm ./lib/"
	}

	postbuildcommands {
		"{COPY} src/CalDict/DataStructs.h ./include/"
	}

	files {
		"src/CalDict/DataStructs.h",
		"src/CalDict/*.cpp",
		"src/CalDict/*.cxx"
	}

	includedirs {
		"./",
		"src/CalDict",
	}

	sysincludedirs {
		ROOTIncludeDir
	}

	libdirs {
		ROOTLibDir
	}

	links {
		"Gui", "Core", "Imt", "RIO", "Net", "Hist", 
		"Graf", "Graf3d", "Gpad", "ROOTDataFrame", "ROOTVecOps",
		"Tree", "TreePlayer", "Rint", "Postscript", "Matrix",
		"Physics", "MathCore", "Thread", "MultiProc", "m", "dl"
	}

	filter "system:macosx or linux"
		linkoptions {
			"-pthread",
			"-rdynamic"
		}

	filter "configurations:Debug"
		symbols "On"

	filter "configurations:Release"
		optimize "On"

project "SabreRecon"
	kind "ConsoleApp"
	language "C++"
	cppdialect "c++11"
	targetdir "./bin/"
	objdir "./objs/"

	files {
		"src/*.cpp",
		"src/*.h",
		"src/Detectors/*.cpp",
		"src/Detectors/*.h",
		"src/EnergyLoss/*.cpp",
		"src/EnergyLoss/*.h",
		"src/CalDict/*.h"
	}

	includedirs {
		"src/",
		"src/CalDict",
		"src/Detectors",
		"src/EnergyLoss"
	}

	sysincludedirs {
		ROOTIncludeDir
	}

	libdirs {
		ROOTLibDir,
	}

	links {
		"CalDict", "Gui", "Core", "Imt", "RIO", "Net", "Hist", 
		"Graf", "Graf3d", "Gpad", "ROOTDataFrame", "ROOTVecOps",
		"Tree", "TreePlayer", "Rint", "Postscript", "Matrix",
		"Physics", "MathCore", "Thread", "MultiProc", "m", "dl"
	}

	filter "system:macosx or linux"
		linkoptions {
			"-pthread",
			"-rdynamic"
		}

	filter "configurations:Debug"
		symbols "On"

	filter "configurations:Release"
		optimize "On"