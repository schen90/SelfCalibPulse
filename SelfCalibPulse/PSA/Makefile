TARGET = AnaData_LibTrap

ROOTCFLAGS  = $(shell root-config --cflags)
ROOTLIBS    = $(shell root-config --libs) -lMinuit
CFLAGS = -I./ 

GXX = g++

all: $(TARGET)

$(TARGET) : $(TARGET).cc
	$(GXX) -o $@ $(TARGET).cc $(CFLAGS) $(ROOTCFLAGS) $(ROOTLIBS)

clean:
	@rm -f $(TARGET) *.o
