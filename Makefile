# Makefile for liboptim - Optimization Library

# Compiler and flags
CC = gcc
CFLAGS = -Wall -Wextra -O2 -fPIC
AR = ar
ARFLAGS = rcs

# Library name
LIBNAME = optim

# Source files
SRCS = brent_min.c brent_root.c nelder_mead.c

# Object files
OBJS = $(SRCS:.c=.o)

# Static library
STATIC_LIB = lib$(LIBNAME).a

# Shared library
SHARED_LIB = lib$(LIBNAME).so

# Example executable
EXAMPLE = example

# Default target
all: $(STATIC_LIB) $(SHARED_LIB)

# Rule to create object files
%.o: %.c opt.h
	$(CC) $(CFLAGS) -c $< -o $@

# Rule to create static library
$(STATIC_LIB): $(OBJS)
	$(AR) $(ARFLAGS) $@ $^

# Rule to create shared library
$(SHARED_LIB): $(OBJS)
	$(CC) -shared -o $@ $^

# Rule to build the example
$(EXAMPLE): $(EXAMPLE).c $(STATIC_LIB)
	$(CC) $(CFLAGS) -o $@ $< -L. -l$(LIBNAME) -lm

# Clean target
clean:
	rm -f $(OBJS) $(STATIC_LIB) $(SHARED_LIB) $(EXAMPLE)

# Phony targets
.PHONY: all clean example
