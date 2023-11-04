# all: debug release

# debug:
# 	g++ MaureenSettings.cpp MaureenMain.cpp -g -o test2

# release:
# 	g++ MaureenSettings.cpp MaureenMain.cpp -Ofast -DNDEBUG -o test2_release

# clean:
# 	rm test2 test2_release

CXX=g++			# C++ compiler
CXXFLAGS=-Wall	# compiler flags

# source files
SRC=MaureenSettings.cpp MaureenMain.cpp


all: debug release

debug:
	$(CXX) $(SRC) $(CXXFLAGS) -g -o test2

release:
	$(CXX) $(SRC) $(CXXFLAGS) -Ofast -DNDEBUG -o test2_release

clean:
	rm test2 test2_release