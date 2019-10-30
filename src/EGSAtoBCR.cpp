#include "Tools.h"

typedef unsigned int int_text;
typedef unsigned int int_suff; 	//-2^31 to 2^31
typedef unsigned int int_lcp;
typedef unsigned char int8; //0 to 2^8

	
typedef struct{
	dataTypeNChar	text;
	dataTypelenSeq	suff;
	dataTypelenSeq 	lcp;
		
	uchar		bwt;	
} t_GSA;
	
int main(int argc, char **argv) {
	 
	if( argc != 3 )
    {
      std::cerr << "Error usage: " << argv[0] << " fastaFile numSeq" << std::endl;
      exit(1);
    }
	
	string fastaFile = argv[1];
	string numSeq = argv[2];
	

	//EGSA
	string fnEGSA;
	fnEGSA = fastaFile + "." + numSeq + ".gesa\0";
	
	FILE *f_ESA = fopen(fnEGSA.c_str(), "rb");
	if (f_ESA == NULL) {
		std::cerr << "Error opening " << fnEGSA << std::endl;
		printf("fopen failed, errno = %d\n", errno);
        exit (EXIT_FAILURE);
	}	
	fseek(f_ESA , 0, SEEK_SET);	

	std::cerr << "file EGSA: "  << fnEGSA <<  "." << std::endl;	
	
	t_GSA GSA;
	
	//BCR 
	string fnBWT, fnDA, fnLCP;
  	fnBWT = fastaFile + ".ebwt\0";
	fnLCP = fastaFile + ".lcp\0";
	fnDA = fastaFile + ".da\0";
	
	FILE *InFileBWT = fopen(fnBWT.c_str(), "wb");
	if (InFileBWT==NULL) {
		std::cerr << "Error opening "  << fnBWT <<  "." << std::endl;
		exit (EXIT_FAILURE);
	}
	FILE *InFileLCP = fopen(fnLCP.c_str(), "wb");
	if (InFileLCP==NULL) {
		std::cerr << "Error opening "  << fnLCP << "." <<std::endl;
		exit (EXIT_FAILURE);
	}			
	FILE *InFileDA = fopen(fnDA.c_str(), "wb");
	if (InFileDA==NULL) {
		std::cerr << "Error opening " << fnDA <<  "." << std::endl;
		exit (EXIT_FAILURE);
	}
	
	std::cerr << "file ebwt : "  << fnBWT <<  "." << std::endl;
	std::cerr << "file lcp: "  << fnLCP <<  "." << std::endl;
	std::cerr << "file da: "  << fnDA <<  "." << std::endl;
	
	//read egsa
	fread(&GSA.text, sizeof(int_text), 1, f_ESA);	
	fread(&GSA.suff, sizeof(int_suff), 1, f_ESA);	
	fread(&GSA.lcp, sizeof(int_lcp), 1, f_ESA);	
	fread(&GSA.bwt, sizeof(int8), 1, f_ESA);
		
	dataTypeNChar  numEle=0;

	while ((!feof(f_ESA)) )  
	{
		fwrite(&GSA.bwt, sizeof(uchar), 1, InFileBWT);	
		fwrite(&GSA.lcp, sizeof(dataTypelenSeq), 1, InFileLCP);	
		fwrite(&GSA.text, sizeof(dataTypeNSeq), 1, InFileDA);	
		
		fread(&GSA.text, sizeof(int_text), 1, f_ESA);	
		fread(&GSA.suff, sizeof(int_suff), 1, f_ESA);	
		fread(&GSA.lcp, sizeof(int_lcp), 1, f_ESA);	
		fread(&GSA.bwt, sizeof(int8), 1, f_ESA);
		
		numEle++;
	}
	
	fclose(f_ESA);
	fclose(InFileBWT);
	fclose(InFileLCP);
	fclose(InFileDA);
	

	std::cerr <<  "The total number of elements is " << numEle << "\n";
	
	return 0;
}

