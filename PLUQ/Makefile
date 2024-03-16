CC=gcc
CFLAGS=-Iinclude -O3 -mfma -mavx2
FLINT_LIB_PATH=/usr/local/lib
FLINT_FLAGS=-I$(FLINT_LIB_PATH)/include
FLINT_LIBS=-L$(FLINT_LIB_PATH) -lflint -lgmp

SRCDIR=src
INCDIR=include
BINDIR=bin

SRCS=$(wildcard $(SRCDIR)/*.c)
OBJS=$(SRCS:$(SRCDIR)/%.c=$(BINDIR)/%.o)

EXECUTABLE_TEST=$(BINDIR)/test
EXECUTABLE_BENCH=$(BINDIR)/bench

all: $(EXECUTABLE_TEST) $(EXECUTABLE_BENCH)

$(EXECUTABLE_TEST): $(OBJS)
	$(CC) $(CFLAGS) $(FLINT_FLAGS)  $(filter-out bin/main_bench.o, $^) -o $@ $(FLINT_LIBS)

$(EXECUTABLE_BENCH): $(OBJS)
	$(CC) $(CFLAGS) $(FLINT_FLAGS) $(filter-out bin/main.o, $^) -o $@ $(FLINT_LIBS)

$(BINDIR)/%.o: $(SRCDIR)/%.c | $(BINDIR)
	$(CC) $(CFLAGS) $(FLINT_FLAGS) -c $< -o $@ $(FLINT_LIBS)

$(BINDIR):
	mkdir -p $(BINDIR)

clean:
	rm -rf $(BINDIR)

.PHONY: all clean
