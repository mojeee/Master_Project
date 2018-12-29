# The compiler
FC = ifort
NAME = go
# flags for debugging or for maximum performance, comment as necessary
#FCFLAGS = -g -fbounds-check
FCFLAGS = -openmp
FCFLAGS += -O2
# flags forall (e.g. look for system .mod files, required in gfortran)
#FCFLAGS += -I/usr/include

# libraries needed for linking, unused in the examples
#LDFLAGS = -li_need_this_lib

# List of executables to be built within the package
all: bmp11.f
	${FC} ${FCFLAGS} bmp11.f -o ${NAME}

clean:
	rm -f *~
	rm -f case*
	rm -f PDF*
	rm -f fort*
	rm -f test.csv
	rm -f sum.csv
	rm -rf lempdf*
	rm -f log