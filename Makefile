#desactivation of the implicite rules
.SUFFIXES:
.PHONY: nonreg

#Definitions
EXECUTABLE = ECOGEN
CXX = mpicxx


CXXFLAGS = -std=c++11 -Wall -Wextra -Wpedantic -Wno-nonnull $(MYCXXFLAGS)

release: CXXFLAGS += -O3

debug: CXXFLAGS += -O0 -g

coverage: CXXFLAGS += -O3 -fprofile-arcs -ftest-coverage

profile: CXXFLAGS += -O3 -pg

SOURCES = $(shell find ./src -type f -name "*.cpp")
OBJETS = $(SOURCES:.cpp=.o)
GCOV_OBJ = $(SOURCES:.cpp=.gcno) $(SOURCES:.cpp=.gcda)

all: release # Default is release build

all debug release coverage profile: exec

exec: $(OBJETS)
		$(CXX) $^ -o $(EXECUTABLE) $(CXXFLAGS)

%o: %cpp
		$(CXX) -c $< -o $@ $(CXXFLAGS)


###

depend:
		makedepend $(SOURCES)

clean:
		rm -rf $(OBJETS) $(GCOV_OBJ)

cleanres:
		rm -rf ./results/*

#Use command: time make nonreg
nonreg:
	@./nonreg/nonreg.sh || true

cleannonreg:
		rm -rf ./nonreg/reports/*
