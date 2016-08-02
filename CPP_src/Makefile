#***************************************************************************/
#**    Author: Joey Azofeifa                                               */
#**    Date:   1/20/2016	                                           */
#**    Description: Makefile for SE Calculator                             */
#***************************************************************************/

CXX             = g++
CXXFLAGS        = -Wno-unused-variable -Wno-non-virtual-dtor -std=c++11 -fopenmp -Wno-write-strings -Wno-literal-suffix 
EXEC            = ${PWD}/CORR
ARCH  		= getconf LONG_BIT

CPP_FLAGS_32  	= -D32_BIT
CPP_FLAGS_64  	= -D64_BIT

GCCVERSION 		= $(shell ${CXX} -dumpversion)
COO: main.o    load.o split.o read_in_parameters.o 
	@printf "linking           : "
	@${CXX} ${CXXFLAGS}    ${PWD}/main.o ${PWD}/load.o   ${PWD}/split.o ${PWD}/read_in_parameters.o   -o ${EXEC}
	@printf "done\n"
	@echo "========================================="
	@printf "finished :)\n"




main.o: 
	@echo "========================================="
	@echo "GCC version: "${GCCVERSION}
	@printf "main              : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/main.cpp 
	@printf "done\n"


load.o: 
	@printf "load              : "	
	@${CXX} -c ${CXXFLAGS} ${PWD}/load.cpp 
	@printf "done\n"
split.o:
	@printf "split             : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/split.cpp 
	@printf "done\n"


read_in_parameters.o:
	@printf "read_in_parameters: "
	@${CXX} -c ${CXXFLAGS} ${PWD}/read_in_parameters.cpp 
	@printf "done\n"





clean:
	@rm -f *.o