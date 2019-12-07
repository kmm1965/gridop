TARGET = gridop2

#CC = g++ -openmp
CC = icpc -qopenmp
RM = rm -rf 

BOOST_LIBS = -L ${BOOST_LIBRARY_PATH} \
	-lboost_timer -lboost_chrono
 
 OPTS = -std=c++11 -I ~/.local/include -I ${BOOST_INCLUDE_PATH}

SRC = gridop2.cpp

OBJS = $(SRC:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) -s $(OBJS) $(BOOST_LIBS)

%.o: %.cpp
	$(CC) $(OPTS) -c $< -o $@

clean:
	-$(RM) $(TARGET) $(OBJS)
