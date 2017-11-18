# Makefile for the 2D-Pred-Prey Simulation

SRC_DIR=src
HEADERS=$(wildcard $(SRC_DIR)/*.hpp)
SRC_FILES=$(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES=$(patsubst $(SRC_DIR)/%.cpp, %.o, $(SRC_FILES))


TEST_DIR=test
TEST_HEADERS=$(wildcard $(TEST_DIR)/*.hpp)
TEST_SRC_FILES=$(wildcard $(TEST_DIR)/*.cpp)
TEST_OBJ_FILES=$(patsubst $(TEST_DIR)/%.cpp, %.o, $(TEST_SRC_FILES))

CXX=c++
CPPSTD=-std=c++11 
DEBUG=-g
CPPUNITLDFLAGS=-lcppunit -lboost_program_options
INC=-I$(SRC_DIR) -I$(TEST_DIR) -I$(HOME)/include 

EXE_FILE=hmc
EXE_TEST=hmc-cppunittests


$(EXE_FILE): $(OBJ_FILES) 
	$(CXX) $(CPPSTD) $(DEBUG) -o $@  $^ $(CPPUNITLDFLAGS) 

$(EXE_TEST): $(TEST_OBJ_FILES) 
	$(CXX) $(CPPSTD) $(DEBUG) -o $@  $^ $(INC) $(CPPUNITLDFLAGS) 

## objs      : create object files
.PHONY : objs
objs : $(OBJ_FILES) $(TEST_OBJ_FILES)

%.o : $(SRC_DIR)/%.cpp $(HEADERS)
	$(CXX) $(CPPSTD) $(DEBUG) -c $< -o $@ $(INC)

%.o : $(TEST_DIR)/%.cpp $(TEST_HEADERS)
	$(CXX) $(CPPSTD) $(DEBUG) -c $< -o $@ $(INC) 


## test    : create and run unit tests
.PHONY : test
test : $(EXE_TEST)
	./$<


## clean     : remove auto generated files
.PHONY : clean
clean :
	rm -f $(OBJ_FILES)
	rm -f $(TEST_OBJ_FILES)

	rm -f $(EXE_FILE)
	rm -f $(EXE_TEST)

	
	rm -f TestResults.xml

## variables : Print variables
.PHONY :variables
variables:
	@echo SRC_DIR:        $(SRC_DIR)
	@echo SRC_FILES:      $(SRC_FILES)
	@echo OBJ_FILES:      $(OBJ_FILES)

	@echo TEST_DIR:       $(TEST_DIR)
	@echo TEST_SRC_FILES: $(TEST_SRC_FILES)
	@echo TEST_OBJ_FILES: $(TEST_OBJ_FILES)



## help      : Print help
.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<
