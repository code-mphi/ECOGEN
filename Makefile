#desactivation of the implicite rules
.SUFFIXES:

#Definitions
EXECUTABLE = ECOGEN
CXX = mpicxx

# -Wundef: Warn about undefined macros in #if directives
# -Wshadow: Warn when a variable shadows an outer scope variable
# -Wuninitialized: Warn about uninitialized variables
# -Wnull-dereference
CXXFLAGS += -fPIC -std=c++11 -DOMPI_SKIP_MPICXX -Wall -Wextra -Wpedantic -Wundef -Wshadow -Wuninitialized

release: CXXFLAGS += -O3 -DNDEBUG

relwithassert: CXXFLAGS += -O3

debug: CXXFLAGS += -O0 -g

coverage: CXXFLAGS += -O3 --coverage
coverage: LDFLAGS += --coverage

profile: CXXFLAGS += -O3 -pg

LIBSOURCES = $(shell find ./src -type f -name "*.cpp" -not -name 'main.cpp' -not -name '*OutputVTKHDF.cpp')
LIBOBJETS = $(LIBSOURCES:.cpp=.o)

GCOV_OBJ = $(LIBSOURCES:.cpp=.gcno) $(LIBSOURCES:.cpp=.gcda) main.gcno main.gcda

## Manage dependencies to header files. See: https://makefiletutorial.com/#makefile-cookbook
DEPS := $(LIBOBJETS:.o=.d) main.d
# Every folder in ./src will need to be passed to GCC so that it can find header files
INC_DIRS := $(shell find ./src -type d)
# Add a prefix to INC_DIRS. So moduleA would become -ImoduleA. GCC understands this -I flag
INC_FLAGS := $(addprefix -I,$(INC_DIRS))
# The -MMD and -MP flags: MMD creates dep files containing user headers, MP creates a phony target for each file
# These files will have .d instead of .o as the output.
CXXFLAGS += $(INC_FLAGS) -MMD -MP

all: release # Default is release build

# Use build from sources (because simpliest) by default
all debug release relwithassert coverage profile: exe_from_sources

# Command line to build executable from sources
exe_from_sources: src/main.o $(LIBOBJETS)
		$(CXX) $(LDFLAGS) $^ -o $(EXECUTABLE)

# Command line to build executable using shared lib (not used if not called explicitely: 'make exe_using_lib'
exe_using_lib: src/main.o lib$(EXECUTABLE).so
		$(CXX) $(LDFLAGS) $^ -o $(EXECUTABLE) -Wl,-rpath,'$$ORIGIN'

# Shared library compilation (used only if 'exe_using_lib' is called or if called explicitely: 'make libECOGEN.so')
lib$(EXECUTABLE).so: $(LIBOBJETS)
		$(CXX) $(LDFLAGS) -shared  $^ -o $@

%.o: %.cpp
		$(CXX) -c $< -o $@ $(CXXFLAGS) $(COVFLAGS)

###
# Without .PHONY, 'clean', 'cleanres'... targets will not run if there is a file whose name is the target name in the directory
.PHONY: clean distclean cleanres nonreg cleannonreg

clean:
		find ./src -type f -name "*.o" -delete -o -name "*.gcda" -delete -o -name "*.gcno" -delete -o -name "*.d" -delete

distclean: clean
		rm -f ECOGEN libECOGEN.so

cleanres:
		rm -rf ./results/*

#Use command: time make nonreg
nonreg:
	@./nonreg/nonreg.sh || true

cleannonreg:
		rm -rf ./nonreg/reports/*

# Include the .d makefiles. The - at the front suppresses the errors of missing
# Makefiles. Initially, all the .d files will be missing, and we don't want those
# errors to show up.
-include $(DEPS)
