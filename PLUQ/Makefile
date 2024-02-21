CC=gcc
CFLAGS=-Iinclude -O3 -mfma -mavx2

SRCDIR=src
INCDIR=include
BINDIR=bin

SRCS=$(wildcard $(SRCDIR)/*.c)
OBJS=$(SRCS:$(SRCDIR)/%.c=$(BINDIR)/%.o)

EXECUTABLE_TEST=$(BINDIR)/test
EXECUTABLE_BENCH=$(BINDIR)/bench

all: $(EXECUTABLE_TEST) $(EXECUTABLE_BENCH)

$(EXECUTABLE_TEST): $(OBJS)
	$(CC) $(CFLAGS) $(filter-out bin/main_bench.o, $^) -o $@

$(EXECUTABLE_BENCH): $(OBJS)
	$(CC) $(CFLAGS) $(filter-out bin/main.o, $^) -o $@

$(BINDIR)/%.o: $(SRCDIR)/%.c | $(BINDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BINDIR):
	mkdir -p $(BINDIR)

clean:
	rm -rf $(BINDIR)

.PHONY: all clean
