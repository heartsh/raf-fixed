######################################################################
# Makefile
######################################################################

######################################################################
# Compilation flags
######################################################################

CXX = g++

CXX_FLAGS = -Wall -Wundef -O3 -DNDEBUG -fomit-frame-pointer \
	-ffast-math -funroll-all-loops -funsafe-math-optimizations \
	-fpeel-loops -Winline --param large-function-growth=100000 \
	--param max-inline-insns-single=100000 \
	--param inline-unit-growth=100000 -fpermissive

OTHER_FLAGS = 

LINK_FLAGS = -lm

######################################################################
# Compilation rules
######################################################################

RAF_SRCS = \
	AlignAndFold.cpp \
	AlignmentShell.cpp \
	IO.cpp \
	Progressive.cpp \
	RAF.cpp \
	ScoringScheme.cpp \
	Utilities.cpp \
	Tree.cpp 

RAF_OBJS = $(RAF_SRCS:%.cpp=%.o)

.PHONY: all clean

all: raf

%.o: %.cpp *.hpp *.ipp Makefile
	$(CXX) $(CXX_FLAGS) $(OTHER_FLAGS) -c $<

raf: $(RAF_OBJS)
	$(CXX) $(CXX_FLAGS) $(OTHER_FLAGS) $(RAF_OBJS) $(LINK_FLAGS) -o raf

clean:
	rm -f raf *.o

######################################################################
# Machine-specific rules
######################################################################

native:
	make all OTHER_FLAGS="-mtune=native"

# default

profile:
	make all OTHER_FLAGS="-pg -g"

multi:
	make all CXX="mpiCC" OTHER_FLAGS="-DMULTI"

# debugging

debug:
	make all CXX_FLAGS="-g -fno-inline -W -Wall"

debugmulti:
	make all CXX="mpiCC" OTHER_FLAGS="-DMULTI" CXX_FLAGS="-g -fno-inline -W -Wall"

assembly:
	make all OTHER_FLAGS="-Wa,-a,-ad"

# pentium 4

gccp4:
	make all OTHER_FLAGS="-march=pentium4 -mtune=pentium4"

gccp4profile:
	make all OTHER_FLAGS="-march=pentium4 -mtune=pentium4 -pg -g"

gccp4multi:
	make all CXX="mpiCC" OTHER_FLAGS="-DMULTI -march=pentium4 -mtune=pentium4"

# athlon64

defaultgccathlon64:
	make all OTHER_FLAGS="-march=athlon64 -mtune=athlon64"

gccathlon64profile:
	make all OTHER_FLAGS="-march=athlon64 -mtune=athlon64 -pg -g"

gccathlon64multi:
	make all CXX="mpiCC" OTHER_FLAGS="-DMULTI -march=athlon64 -mtune=athlon64"

# intel

intel:
	make all CXX="icpc" OTHER_FLAGS="-xN -no-ipo -static"

intelmulti:
	make all LAMHCP="icpc" CXX="mpiCC" OTHER_FLAGS="-DMULTI -xN -no-ipo"
