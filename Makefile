#   make        -- compiles your project into program.exe
#   make clean  -- removes compiled item
#   make handin -- creates a project Zip file for hand in
#
# All .cpp flles are included.

CC = mpicc
SRCS = $(wildcard *.cpp)
HDRS = $(wildcard *.h)
OBJS = $(SRCS:.cpp=.o)
DIRS = $(subst /, ,$(CURDIR))
PROJ = mpi_dbscan

APP = $(PROJ)
CFLAGS= -c -O3 -I/share/apps/pnetcdf-1.12.0/include
LDFLAGS=
LIBS= -lstdc++ -L/share/apps/pnetcdf-1.12.0/lib -lpnetcdf -lm -O3  

all: $(APP)

$(APP): $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o $(APP) $(LIBS)

%.o: %.cpp $(HDRS) $(MF)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o $(APP)

