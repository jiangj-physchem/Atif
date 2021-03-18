CC=g++
DIR_INC=-I ../source
CFLAGS=-Wall -g -O2 $(DIR_INC) 
DIR_SRC= ../source
CPP_FILES=$(wildcard ${DIR_SRC}/*.cpp)
SRC=$(CPP_FILES)
OBJ=$(SRC:.cpp=.o)

TARGET=ATIFexe

defaut: $(TARGET)
	-rm $(OBJ)

$(TARGET): $(OBJ)
	$(CC) -o $@ $(OBJ)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY : clean
clean: 
	rm $(TARGET) $(OBJ)

