#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <string>
#include <iostream>

using namespace std;

typedef unsigned long ulong;

int main(int argc, char **argv) {
	 
	if( argc != 3 )
    {
      std::cerr << "Error usage: " << argv[0] << " fastaFile numSeq" << std::endl;
      exit(1);
    }
	
	string fastaFile = argv[1];
    ulong numSeq;
    sscanf(argv[2], "%lu", &numSeq);
    
    //DOCS
    string fndocs;
    fndocs = fastaFile + ".docs\0";
    
    FILE *InFileD = fopen(fndocs.c_str(), "wb");
    if (InFileD==NULL) {
        std::cerr << "Error opening "  << fndocs <<  "." << std::endl;
        exit (EXIT_FAILURE);
    }
    fwrite(&numSeq, sizeof(ulong), 1, InFileD);
    std::cout << "file docs: " << fndocs <<  "." << std::endl;
    fclose(InFileD);
    
	return 0;
}

