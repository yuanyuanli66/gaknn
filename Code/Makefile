CC = gcc
CFLAGS = -O3 -Wall
CFLAGS += -mtune=native
#CFLAGS += -g -pg -fno-inline-functions -fprofile-arcs
#CFLAGS += -ftree-vectorizer-verbose=1

CFLAGS += -fopenmp -funroll-loops
EXE = gaknn
SRCS = $(wildcard *.c)
OBJS = $(SRCS:%.c=%.o)
LIBS = -lm

# librt is for accessing the real-time clock using clock_gettime()
LIBS += -lrt

#-------------------------

all: $(EXE)
	@echo  "all done"

$(EXE):	$(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

$(OBJS): %.o: %.c ga_knn.h
	$(CC) $(CFLAGS) -c $< 

clean:
	$(RM) $(OBJS) $(EXE) *.txt
	$(RM) *.gcda *.gcov

