# https://habr.com/ru/post/155201/

CC=g++
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=main.cpp Approx.cpp Engine.cpp Equation.cpp Interval.cpp Lattice.cpp Random.cpp TabFunc.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=app

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf app *.o
