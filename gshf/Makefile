PROGNAME = gshf
SOURCEFILES = gshf.cc
OBJS    = $(patsubst %.cc, %.o, $(SOURCEFILES))

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

LDFLAGS       = -O
LIBS         += $(ROOTLIBS)
CFLAGS       += $(ROOTCFLAGS)

%.o: %.cc
	g++ ${CFLAGS} -c  -g -o $@ $<

$(PROGNAME):    $(OBJS)
	g++ -o $@ $(OBJS) $(LDFLAGS) $(LIBS)

test:
	@echo $(ROOTCFLAGS)

clean:	
	-rm -f ${PROGNAME} ${OBJS}

