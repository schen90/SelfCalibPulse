TARGET = MakeDataG4
#TARGET = AnaDataG4

ROOTCFLAGS  = $(shell root-config --cflags)
ROOTLIBS    = $(shell root-config --libs) -lMinuit
CFLAGS = -I./ 

GXX = g++

all: $(TARGET)

$(TARGET) : $(TARGET).C
	$(GXX) -o $@ $(TARGET).C $(CFLAGS) $(ROOTCFLAGS) $(ROOTLIBS)

clean:
	@rm -f $(TARGET) *.o
