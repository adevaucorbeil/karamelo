#
# 'make depend' uses makedepend to automatically generate dependencies 
#               (dependencies are added to end of Makefile)
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

SHELL = /bin/bash
# define the C compiler to use
CC = g++

# define any compile-time flags
CFLAGS = -g -std=c++11
#-O3

# define any directories containing header files other than /usr/include
#
INCLUDES = -I/usr/include/eigen3/
#-I/home/newhall/include  -I../include

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LFLAGS =
#-L/home/newhall/lib  -L../lib

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
LIBS =
#-lmylib -lm

# define the C source files
SRC :=	$(wildcard *.cpp)
INC :=	$(wildcard *.h)
OBJ := 	$(SRC:.cpp=.o)


# define the executable file 
MAIN = mpm

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean

all:    $(MAIN)
	@echo  Simple compiler named mpm has been compiled

$(MAIN): $(OBJ)
	$(CC) $(CFLAGS) $(INCLUDES) $(OBJ) $(LFLAGS) $(LIBS) -o $(MAIN) 

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
%.o:%.cpp %.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRC)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
