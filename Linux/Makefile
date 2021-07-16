CC = g++
CFLAGS = -g

all: Program objdir

objdir:
	mkdir -p obj

OBJFILES = obj/Application.o obj/Element.o obj/Elements.o obj/FiniteElements.o obj/Node.o obj/Nodes.o


Program: $(OBJFILES)
	$(CC) $(CFLAGS) $(OBJFILES) -o Program

obj/Application.o: Application.cpp include/FiniteElements.h objdir
	$(CC) $(CFLAGS) -c Application.cpp -o obj/Application.o

obj/Element.o: src/Element.cpp include/Element.h include/Node.h objdir
	$(CC) $(CFLAGS) -c src/Element.cpp -o obj/Element.o

obj/Elements.o: src/Elements.cpp include/Elements.h include/Nodes.h include/Element.h objdir
	$(CC) $(CFLAGS) -c src/Elements.cpp -o obj/Elements.o

obj/FiniteElements.o: src/FiniteElements.cpp include/FiniteElements.h include/Elements.h objdir
	$(CC) $(CFLAGS) -c src/FiniteElements.cpp -o obj/FiniteElements.o

obj/Node.o: src/Node.cpp include/Node.h objdir
	$(CC) $(CFLAGS) -c src/Node.cpp -o obj/Node.o

obj/Nodes.o: src/Nodes.cpp include/Nodes.h include/Node.h objdir
	$(CC) $(CFLAGS) -c src/Nodes.cpp -o obj/Nodes.o

clean:
	rm -rf ./obj/* Program