CC=g++
DIR_INC=-I ../program/src
CFLAGS=-Wall -g -O2 $(DIR_INC) 
INC=-I ../program/eigen3
DIR_SRC= ../program/src
#ARCH = -arch x86_64
CPP_FILES=$(wildcard ${DIR_SRC}/*.cpp)
SRC=$(CPP_FILES)
OBJ=$(SRC:.cpp=.o)
#EXECUTABLE=exe

TARGET=PolymerDFT

defaut: $(TARGET)
	-rm $(OBJ)

$(TARGET): $(OBJ)
	$(CC) -o $@ $(OBJ)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

.PHONY : clean
clean: 
	rm $(TARGET) $(OBJ)

