TARGET = gridop1

CC = g++
#CC = icpc
RM = rm -rf 

OPTS = -std=c++14 -I ../../include -I ${BOOST_ROOT} -Wno-deprecated-declarations

SRC = gridop1.cpp

OBJS = $(SRC:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) -s $(OBJS)
# $(BOOST_LIBS)

%.o: %.cpp
	$(CC) $(OPTS) -c $< -o $@

clean:
	-$(RM) $(TARGET) $(OBJS)
