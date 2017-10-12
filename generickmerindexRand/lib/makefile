 
CC = g++                        #compiler
DEBUG = -g -std=c++0x 		# -g    adds debugging information to the executable file, 
				# -std=c++0x debugginh flags enable unseroted set/map
CFLAGS = -Wall -c $(DEBUG)      # -c create object for every cpp file & prevernt linking 
LFLAGS = -Wall $(DEBUG)         # -Wall turns on most, but not all, compiler warnings

OBJDIR  = objdir
BINDIR  =../bin
OBJS    = parsers.o basic.o db.o kmerselection.o index.o query.o auxiliary.o
EXE     = makedb makeindex runquery runauxiliary

BASIC_HEADERS=parsers.h basic.h
IDX_HEADERS=kmerselection.h index.h

all: $(OBJDIR) $(BINDIR) $(OBJS) $(EXE)

$(OBJDIR):
	mkdir -p $(OBJDIR)
$(BINDIR):
	mkdir -p $(BINDIR)

makedb: $(OBJDIR)/parsers.o $(OBJDIR)/basic.o $(OBJDIR)/db.o UI_db.cpp
	$(CC) $(LFLAGS) $^ -o $(BINDIR)/makedb

makeindex:$(OBJDIR)/parsers.o $(OBJDIR)/basic.o $(OBJDIR)/db.o \
	$(OBJDIR)/kmerselection.o $(OBJDIR)/index.o UI_index.cpp
	$(CC) $(LFLAGS) $^ -o $(BINDIR)/makeindex

runquery:$(OBJDIR)/parsers.o $(OBJDIR)/basic.o $(OBJDIR)/db.o \
	$(OBJDIR)/kmerselection.o $(OBJDIR)/index.o $(OBJDIR)/query.o UI_query.cpp
	$(CC) $(LFLAGS) $^ -o $(BINDIR)/runquery

runauxiliary:$(OBJDIR)/parsers.o $(OBJDIR)/basic.o $(OBJDIR)/auxiliary.o \
	UI_auxiliary.cpp
	$(CC) $(LFLAGS) $^ -o $(BINDIR)/runauxiliary


parsers.o: parsers.cpp parsers.h
	$(CC) $(CFLAGS) parsers.cpp -o $(OBJDIR)/parsers.o

basic.o: basic.cpp basic.h
	$(CC) $(CFLAGS) basic.cpp  -o $(OBJDIR)/basic.o

db.o: db.cpp db.h $(BASIC_HEADERS)
	$(CC) $(CFLAGS) db.cpp -o $(OBJDIR)/db.o

kmerselection.o: kmerselection.cpp kmerselection.h $(BASIC_HEADERS)
	$(CC) $(CFLAGS) kmerselection.cpp -o $(OBJDIR)/kmerselection.o

index.o: index.cpp index.h kmerselection.h $(BASIC_HEADERS)
	$(CC) $(CFLAGS) index.cpp -o $(OBJDIR)/index.o

query.o: query.cpp query.h $(BASIC_HEADERS) $(IDX_HEADERS)
	$(CC) $(CFLAGS) query.cpp -o $(OBJDIR)/query.o

auxiliary.o: auxiliary.cpp auxiliary.h basic.h parsers.h
	$(CC) $(CFLAGS) auxiliary.cpp  -o $(OBJDIR)/auxiliary.o

clean:
	-rm -r $(OBJDIR)
	-rm -r $(BINDIR) 
