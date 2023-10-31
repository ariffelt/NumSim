CXX=g++         # C++ compiler
CXXFLAGS=-Wall  # compiler flags

# source files
SRC=LouisaSettings.cpp LouisaMain.cpp


all: debug release

debug:
  $(CXX) $(SRC) $(CXXFLAGS) -g -o louisa5

release:
  $(CXX) $(SRC) $(CXXFLAGS) -Ofast -DNDEBUG -o louisa5_release

clean:
  rm louisa5 louisa5_release