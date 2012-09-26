CC = g++
CFLAGS = -O0 -g3 -W -Wall -c -fmessage-length=0
LIB_SEQAN = -IF:\Development\workspace_cpp_test\seqan\projects\library
DIR_BIN = binaries/
HEADERS = JournalIterator.h JournalString.h Node.h STNode.h Tree.h TreeIterator.h Utility.h


####################################################
#Target CPPKurs
journal: main.o
	${CC} main.o -o ${DIR_BIN}journal.exe 
	@echo 'ready'
	
main.o: main.cpp ${HEADERS}
	${CC} ${CFLAGS} $(LIB_SEQAN) main.cpp
	
####################################################
#Target clean

clean:
	rm -rf *o ${DIR_BIN}*.exe *.d