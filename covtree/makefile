SRC_d = ./src
SRC_c = $(SRC_d)/covtree.cpp
SRC_h =
SRC = $(SRC_h) $(SRC_c)
EXE_d = ./bin
EXE = $(EXE_d)/covtree

CC=g++

all: $(EXE)

$(EXE): $(SRC)
	$(CC)  -pthread -std=c++14 -Wall -pedantic -Ofast -march=native -finline-functions -o $(EXE) $(SRC_c) -I$(SRC_d) -I$(SRC_d)/butils  -static-libgcc   -static-libstdc++ -Wl,--as-needed

clean:
	rm -f $(EXE)

