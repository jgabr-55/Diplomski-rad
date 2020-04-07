SHELL=/usr/bin/env bash
-include Makefile.inc

CXX_COMMON:=-I$(PREFIX_INCLUDE) $(CXX_COMMON)
CXX_COMMON+= -L$(PREFIX_LIB) -Wl,-rpath,$(PREFIX_LIB) -lpythia8 -ldl 

###########################

.SECONDEXPANSION:

Analyzer: $$@.cc $(PREFIX_LIB)/libpythia8.a
ifeq ($(FASTJET3_USE),true)
	$(CXX) $< -o $@ -I$(FASTJET3_INCLUDE) $(CXX_COMMON)\
	 -L$(FASTJET3_LIB) -Wl,-rpath,$(FASTJET3_LIB) -lfastjet\
	  $(GZIP_INC) $(GZIP_FLAGS)
else
	@echo "Error: $@ requires FASTJET3"
 endif



# Clean.
clean:
	@rm -f Analyzer; rm -f out[0-9][0-9];\
