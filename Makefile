CC=gcc
CXX=g++
RM=rm -rfv
CPPFLAGS=-g $(shell root-config --cflags)
LDFLAGS=-g $(shell root-config --ldflags)
LDLIBS=$(shell root-config --libs) -lRooFitCore
EXECUTABLE=clusterReconstr
DOXY=doxygen
DOCDIR=doc


SRCS=main.cpp CREntry.cpp CREvent.cpp langaufun.cpp langaufit.cpp langaupro.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(CXX) $(LDFLAGS) -o $(EXECUTABLE) $(OBJS) $(LDLIBS) 

depend: .depend

.depend: $(SRCS)
	rm -f ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(OBJS) $(DOCDIR) *~ .depend $(EXECUTABLE)

install:
	alias $(EXECUTABLE)=.$(shell pwd)/$(EXECUTABLE)

doc:
	$(DOXY) ./Doxyfile;

include .depend
