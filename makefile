CC = g++
CPPFLAGS = -O3 -Wall -Wshadow -ansi -pedantic -DDEBUG $(INCLUDES)
LIBS =
INCLUDES = -I./Colsamm
TARGET = mdsim

SRC = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp, %.o, $(SRC))

$(TARGET): $(OBJS)
	$(CC) $(LFLAGS) $(INCLUDES) -o $(TARGET) $(OBJS) $(LIBS)

.PHONY : clean depend

clean: 
	@/bin/rm -f $(OBJS) 
	@/bin/rm -f $(TARGET)

depend: 
	@makedepend -- $(CPPFLAGS) -- $(SRC)