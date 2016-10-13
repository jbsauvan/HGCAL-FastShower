#makefile 


CC   =   g++

UCFLAGS = -O0 -g3 -Wall -gstabs+  
#UCFLAGS = -O3 -Wall -gstabs+

RUCFLAGS := $(shell root-config --cflags) -I./include/ -I./external/jsoncpp/
LIBS :=  $(shell root-config --libs)  

vpath %.cpp ./external/jsoncpp
vpath %.cpp ./src

SRCPP = main.cpp\
	Cell.cpp\
	Geometry.cpp\
	Generator.cpp\
	Constants.cpp\
	jsoncpp.cpp

         
#OBJCPP = $(SRCPP:.cpp=.o)
OBJCPP = $(patsubst %.cpp,lib/%.o,$(SRCPP))


all : shower_simulation.exe 
	#obj/libDictionary_C.so

lib/%.o : %.cpp
	@echo "> compiling $*"
	@mkdir -p lib/
	@$(CC) -c $< $(UCFLAGS) $(RUCFLAGS) -o $@

shower_simulation.exe : $(OBJCPP)
	@echo "> linking"
	$(CC) $^ $(LIBS) -o $@

clean:
	@echo "> Cleaning object files"
	@rm  -f lib/*.o
        
cleanall: clean
	@echo "> Cleaning executable"
	@rm -f shower_simulation.exe
