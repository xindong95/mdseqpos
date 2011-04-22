SWIG := swig
SWIGINCLUDE := -Igsl -Iswig
SWIGFLAGS := -c++ -python $(SWIGINCLUDE)
SWIGSRC := $(filter-out $(SWIGIGN), $(wildcard swig/*.i)) 
SWIGWRAP := $(patsubst swig%, extension%, $(patsubst %.i, %_wrap.cpp, $(SWIGSRC)))

all: $(SWIGWRAP)
	cp python/* extension

$(SWIGWRAP): extension/%_wrap.cpp: swig/%.i
	$(SWIG) $(SWIGFLAGS) -o $@ $<

variables:
	@echo $(SWIGSRC)
	@echo $(SWIGWRAP)

clean: 
	rm -f extension/*
