TARGET = SelfCalib
Objects = TreeReaderPulse.o

ROOTCFLAGS  = $(shell root-config --cflags)
ROOTLIBS    = $(shell root-config --libs) -lMinuit
CFLAGS = -I./ 

GXX = g++

all: $(TARGET)

$(TARGET) : $(TARGET).cc $(Objects)
	$(GXX) -o $@ $(TARGET).cc $(Objects) $(CFLAGS) $(ROOTCFLAGS) $(ROOTLIBS)

TreeReaderPulse.o : TreeReaderPulse.cc TreeReaderPulse.hh
	$(GXX) -o $@ TreeReaderPulse.cc -c $(ROOTCFLAGS)

clean:
	@rm -f $(TARGET) *.o
