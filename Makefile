# A simple Makefile for compiling small SDL projects

# set the compiler
CC := gcc

# set the compiler flags
CFLAGS := `sdl2-config --libs --cflags` -ggdb3 -O3 --std=c11 -Wall -lSDL2_image -lm 
# add header files here
HDRS :=

# add source files here
SRCS := backup2.c

# generate names of object files
OBJS := $(SRCS:.c=.o)

# name of executable
EXEC := outp-t

# default recipe
all: $(EXEC)
 
showfont: showfont.c Makefile
	$(CC) -o $@ $@.c $(CFLAGS) $(LIBS)

glfont: glfont.c Makefile
	$(CC) -o $@ $@.c $(CFLAGS) $(LIBS)

# recipe for building the final executable
$(EXEC): $(OBJS) $(HDRS) Makefile
	$(CC) -o $@ $(OBJS) $(CFLAGS)

# recipe for building object files
#$(OBJS): $(@:.o=.c) $(HDRS) Makefile
#	$(CC) -o $@ $(@:.o=.c) -c $(CFLAGS)

# recipe to clean the workspace
clean:
	rm -f $(EXEC) $(OBJS)

.PHONY: all clean
