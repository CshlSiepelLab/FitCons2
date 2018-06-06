SRC_d = ./src
SRC_c = $(SRC_d)/Insight2.cpp
SRC_h =
SRC = $(SRC_h) $(SRC_c)
EXE_d = ./bin
EXE = $(EXE_d)/Insight2

CC=g++

all: $(EXE)

$(EXE): $(SRC)
	$(CC)  -pthread -std=c++11 -Wall -pedantic -Ofast -march=native -finline-functions -o $(EXE) $(SRC_c) -I$(SRC_d) -static-libgcc   -static-libstdc++ -Wl,--as-needed


clean:
	rm -f $(EXE_d)/Insight2

