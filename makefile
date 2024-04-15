CC = g++
CPPFLAGS = -I. -I/usr/local/include  
LDFLAGS = -g -fopenmp -L/usr/local/lib  
LDLIBS = -lgsl -lgslcblas -lm -lgomp  # Link against GSL libraries and libgomp

DEPS = bigraph.h utility.h ldp-nb.h
OBJ = bigraph.o main.o utility.o ldp-nb.o

%.o: %.cpp $(DEPS)
	$(CC) -std=c++1y -c -O3 -o $@ $< $(CPPFLAGS) $(LDFLAGS)

ldp-nb: $(OBJ)
	$(CC) -std=c++1y -O3 -pthread -o $@ $^ $(CPPFLAGS) $(LDFLAGS) $(LDLIBS)

clean:
	-rm -f ldp-nb *.o

