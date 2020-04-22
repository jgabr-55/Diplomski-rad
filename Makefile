CC=g++
LDFLAGS = $(shell root-config --libs)
CFLAGS=-I. $(shell root-config --cflags)
DEPS = Analyzer.h

%.o: %.cpp $(DEPS)
	$(CC) $(LDFLAGS) -c -o $@ $< $(CFLAGS)

analyse: analyse.o Analyzer.o
	$(CC) $(LDFLAGS) -o analyse analyse.o Analyzer.o


clean:
	rm -rf *.o analyse



