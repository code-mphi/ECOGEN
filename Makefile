#desactivation of the implicite rules
.SUFFIXES:
.PHONY: nonreg

#Definitions
EXECUTABLE = ECOGEN
CXX = mpicxx
CXXFLAGS = -O3 -std=c++11 -Wall -Wextra -Wpedantic #release
#CXXFLAGS = -g -std=c++11 -Wall -Wextra -Wpedantic #debug
#CXXFLAGS = -O3 -std=c++11 -Wall -Wextra -Wpedantic -fprofile-arcs -ftest-coverage #code coverage

dirs = $(shell find . -type d)
SOURCES = $(foreach dir,$(dirs),$(wildcard $(dir)/*.cpp))
OBJETS = $(SOURCES:.cpp=.o)

all: $(OBJETS)
		$(CXX) $^ -o $(EXECUTABLE) $(CXXFLAGS)

%o: %cpp
		$(CXX) -c $< -o $@ $(CXXFLAGS)


###

depend:
		makedepend $(SOURCES)

clean:
		rm -rf $(OBJETS)

cleanres:
		rm -rf ./results/*

#Use command: time make nonreg
nonreg:
	@./nonreg/nonreg.sh || true

cleannonreg:
		rm -rf ./nonreg/reports/*