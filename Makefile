CC = clang++
CFLAGS = -std=c++17 -Wall -Wextra -msse
# DEBUG mode
CFLAGS += -g
# RELEASE mode
#CFLAGS += -O3 -DNDEBUG
LDFLAGS =

APP_NAME = grasshopper
APP_SOURCES = $(wildcard *.cpp)
APP_OBJECTS = $(APP_SOURCES:%.cpp=%.o)

all: $(APP_OBJECTS)
	$(CC) $(LDFLAGS) $(APP_OBJECTS) -o $(APP_NAME)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

profile: $(APP_NAME)
	rm -f callgrind.out
	valgrind --tool=callgrind --callgrind-out-file=callgrind.out ./$(APP_NAME)
	kcachegrind callgrind.out

clean:
	rm -f *.o $(APP_NAME) callgrind.out
