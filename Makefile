# Compiler
CC = g++

# Compiler flags
CFLAGS = -lgsl -lgslcblas -std=c++17 -lboost_program_options

# Target executable
TARGET = g-plot

# Object files
OBJS = g-plot.o integral.o flag_parser.o error_handler.o

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS)

g-plot.o: g-plot.cpp ./src/integral.h
	$(CC) $(CFLAGS) -c g-plot.cpp

integral.o: ./src/integral.cpp ./src/integral.h
	$(CC) $(CFLAGS) -c ./src/integral.cpp

flag_parser.o: ./src/flag_parser.cpp ./src/flag_parser.h
	$(CC) $(CFLAGS) -c ./src/flag_parser.cpp

error_handler.o: ./src/error_handler.cpp ./src/error_handler.h
	$(CC) $(CFLAGS) -c ./src/error_handler.cpp

clean:
	rm -f *.o $(TARGET)