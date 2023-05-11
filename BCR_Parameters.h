/**
 ** BCR is part of:
 ** BEETL: Burrows-Wheeler Extended Tool Library
 ** Documentation in: doc/BEETL.md
 **
 ** Copyright (c) 2011-2014 Illumina, Inc. **
 ** BEETL software package is
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** Citations: 
 ** 
 ** Markus J. Bauer, Anthony J. Cox and Giovanna Rosone
 ** Lightweight BWT Construction for Very Large String Collections.
 ** Proceedings of CPM 2011, pp.219-231
 ** 
 ** Markus J. Bauer, Anthony J. Cox, Giovanna Rosone and Marinella Sciortino
 ** Lightweight LCP Construction for Next-Generation Sequencing Datasets. 
 ** Proceedings of WABI 2012, pp 326-337, 2012
 
 ** Markus J. Bauer, Anthony J. Cox, Giovanna Rosone 
 ** Lightweight algorithms for constructing and inverting the BWT of string collections. 
 ** Theoretical Computer Science 483: 134-148 (2013)
 **  
 ** Anthony J. Cox, Fabio Garofalo, Giovanna Rosone, Marinella Sciortino
 ** Lightweight LCP construction for very large collections of strings. 
 ** Journal of Discrete Algorithms (2016) 37: 17-33
 **
 ** By Giovanna Rosone
 **
 **/
 
 /*
 * Setting
 */

 
#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <iostream>
#include <fstream>

#define SIZEBUFFER 1024     //Size of the buffer for I/O partial ebwt/LCP/DA/SA

#define ext  ".aux"

#define SIZE_ALPHA 256  


typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;

#define dataTypedimAlpha uchar  //size of the alphabet (in biologic case 6 ($,A,C,G,N,T))
//note that you have to reserve two characters:
#define TERMINATE_CHAR '$'     //it is the symbol used as "end of strings", it must be lexicographically smaller than all the letters of the alphabet
#define TERMINATE_CHAR_LEN uchar(255)      //it is stored in cyc files, it is ignored by the algorithm, so it must not belong to the alphabet

/* For dataTypeLengthSequences USE: 
 *      0 (unsigned char)  - for read length <= (255-1) 
 *	1 (unsigned short) - for read length between 255 and (65.536-1)
 *	2 (unsigned int)   - for read length <= between 65.535 and (4.294.967.296-1) 
 *	3 (unsigned long)  - for read length <= otherwise 
 */

// Type size for Sequences Length (in biologic case 100)
#define dataTypeLengthSequences 0		

// Type size for Number of Sequences in the input file
#define dataTypeNumSeq 2		

// Type size for Number of Character in the input file (length of the BWT)
#define dataTypeNumChar 3		


//Set the types
#if dataTypeLengthSequences == 0
	#define dataTypelenSeq uchar
#elif dataTypeLengthSequences == 1
	#define dataTypelenSeq ushort
#elif dataTypeLengthSequences == 2
	#define dataTypelenSeq uint
#elif dataTypeLengthSequences == 3
	#define dataTypelenSeq ulong	
#endif

#if dataTypeNumSeq == 0
	#define dataTypeNSeq uchar
#elif dataTypeNumSeq == 1
	#define dataTypeNSeq ushort
#elif dataTypeNumSeq == 2
	#define dataTypeNSeq uint
#elif dataTypeNumSeq == 3
	#define dataTypeNSeq ulong	
#endif


#if dataTypeNumChar == 0
	#define dataTypeNChar uchar
#elif dataTypeNumChar == 1
	#define dataTypeNChar ushort
#elif dataTypeNumChar == 2
	#define dataTypeNChar uint
#elif dataTypeNumChar == 3
	#define dataTypeNChar ulong	
#endif

////////////

//Print of the output (BWT/DA/SA/LCP)
//Store a txt file containing (BWT/DA/SA/LCP)
#define printFinalOutput 0

//Verbose
#define verboseEncode 0
#define verboseDecode 0

//if you want to delete the partial files, please set it to 1
#define deletePartialBWT 1 
#define deletePartialLCP 1 
#define deletePartialGSA 1 
#define deletePartialQS 1
#define deletePartialSAP 1 
#define deleteCycFiles 1

#if FASTQ==1
    //Compute the quality score permutation associated to BWT   permutation
	#define USE_QS 1

    //if you want to store the titles of each read in fastQ into a file with extension .title, please set it to 1
    #define STORE_TITLE_FASTQ 1
#else
    #define USE_QS 0             //if you want to build QS permutation from fasta and qs file, please set USE_QS to 1
#endif

#if SAP==1
    //Compute the SAP-array associated with BWT permutation
	#define BUILD_SAP 1

    #define BCR_SET_ALN_RH 1
#else
    #define BUILD_SAP 0  //if you want to build the SAP array, please set BUILD_SAP to 1

	//if BCR_SET_ALN_RH=0 then BCR computes the EBWT (the input is a set) aligning strings left 
	//if BCR_SET_ALN_RH=1 then BCR computes the EBWT (the input is a set) aligning strings right 
	#define BCR_SET_ALN_RH 0
#endif

//Use kseq to read sequences
#define KSEQ_PARSER 1

//if you want to compute the LCP array, please set it to 1
#define BUILD_LCP 0

//The pair (da[i], sa[i]) is the gsa[i]
//if you want to compute the SA array, please set it to 1
#define BUILD_SA 0
//if you want to compute the DA array \in [0..nSeq], please set it to 1
#define BUILD_DA 1
//if you want to compute the bit DA array, please set it to 1
//it works in internal memory
#define BUILD_DA_bit 0



//if you want to Store the 'end positions' of the end-markers (one for each sequence), please set it to 1
#define STORE_ENDMARKER_POS 0

//if BUILD_BCR_ALTERNATE=0 then BCR computes the eBWT/SA/DA/LCP in straightforward order of the sequences (lexicographical order)
//if BUILD_BCR_ALTERNATE=1 then BCR computes the eBWT/SA/DA/LCP in alternating order of the sequences (alternating lexicographical order) See paper...
#define BUILD_BCR_ALTERNATE 0 //0 --> else we compute alternate order the BWT of the sequences

//if BCR_SET=1 then BCR computes the EBWT (the input is a set) (one can have strings of different length, so BCR uses the symbol TERMINATE_CHAR_LEN) 
//if BCR_SET=0 then BCR computes the BWT (the input is a single sequence)
#define BCR_SET 1				

//if BCR_INPUT_IN_MEMORY==1, BCR loads the input file in a string and compute the BWT of the string (it computes the BWT of the reverse string),  
//if BCR_INPUT_IN_MEMORY==0, BCR reads from files (cyc files)
#define BCR_INPUT_IN_MEMORY 0	

//if KEEP_eBWT_IN_EXT_MEMORY==1, BCR uses files for partials ebwts
//if KEEP_eBWT_IN_EXT_MEMORY==0, BCR uses strings for partials ebwts
//In both cases, SA, DA, LCP are stored in files.
#define KEEP_eBWT_IN_EXT_MEMORY  1

//Save lengths of each sequence in a file .len
#define STORE_LENGTH_IN_FILE  0

//if OUTPUT_FORMAT == 0, the output format of BCR is at most 5 files - built one after the other
//if OUTPUT_FORMAT == 1, the output format of BCR is as the output of EGSA (.gesa file). BUILD_LCP, BUILD_DA and BUILD_SA must be set to 1. Please, set the types as in eGSA
//if OUTPUT_FORMAT == 2, the output format of BCR is a unique file .egsa. BUILD_LCP must be set to 1 (we do not use a struct), BUILD_DA and BUILD_SA could be set to a either 0 or 1.  Order: ebwt, lcp, da, sa
//if OUTPUT_FORMAT == 3, the output format of BCR is at most 5 files at the same time
//if OUTPUT_FORMAT == 4, the output format of BCR is at most 3 files (ebwt, da), lcp, sa
//if OUTPUT_FORMAT == 5, the output format of BCR is at most 3 files ebwt, (lcp, da), sa
//if OUTPUT_FORMAT == 6, the output format of BCR is at most 3 files ebwt, lcp, (sa, da)
#define OUTPUT_FORMAT 0

//if OUTPUT_linear_SuffixArray == 1, BCR also computes the SA of the concatenated strings 
//if OUTPUT_linear_SuffixArray == 0, BCR does not compute the SA of the concatenated strings 
#define OUTPUT_linear_SuffixArray 0

//Computes the EGSA by starting to another EGSA data structure
//See output EGSA of Felipe Louze
//if BUILD_BCR_FROM_BCRpartials == 1, BCR takes in input the BCR partial files and adds the symbols of new sequences.
//if BUILD_BCR_FROM_BCRpartials == 0, BCR starts from scratch
#define BUILD_BCR_FROM_BCRpartials 0


//BCR reads the cyc files already computed by transpose.cpp in a previous execution.
//if BCR_FROMCYC=0 then BCR builds the cyc files before, otherwise BCR does not build the cyc files.
#define BCR_FROMCYC	0			




#endif
