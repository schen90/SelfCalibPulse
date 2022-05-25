TARGET = MakePSBase
#TARGET = MakePSBaseChi2
#TARGET = MakePSBaseRealPS
#TARGET = MakePSBasePSC

ROOTCFLAGS  = $(shell root-config --cflags)
ROOTLIBS    = $(shell root-config --libs) -lMinuit
CFLAGS = -I./ 

GXX = g++

all: $(TARGET)

$(TARGET) : $(TARGET).C
	$(GXX) -o $@ $(TARGET).C $(CFLAGS) $(ROOTCFLAGS) $(ROOTLIBS)

clean:
	@rm -f $(TARGET) *.o
