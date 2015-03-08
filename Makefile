
CXX=clang++

all:
	$(CXX) -std=c++11 -Wall -O3 -o Lylios -DGLM_SWIZZLE `pkg-config sdl --cflags --libs` main.cpp
