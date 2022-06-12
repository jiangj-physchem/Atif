CC=cl.exe
DIR_INC=-I ../source
CFLAGS=-Wall -g -O2 $(DIR_INC) 
DIR_SRC= ../source
CPP_FILES=$(wildcard ${DIR_SRC}/*.cpp)
SRC=$(CPP_FILES)
OBJ=$(SRC:.cpp=.o)

TARGET=AtifExe

defaut: $(TARGET)
	-rm $(OBJ)

$(TARGET): $(OBJ)
	$(CC) /d $@ $(OBJ)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY : clean
clean: 
	rm $(TARGET) $(OBJ)

