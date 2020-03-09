#include "Tools.h"

/*
Step 2. analyzes alpha-clusters in order to evaluate a degree of similarity between any read and any genome in the collection.
If EBWT=1, the similarity scores are computed by using both the document array and the ebwt. 
If EBWT=0, the similarity scores are computed by using only the document array.

As preprocessing, it need to have the two data structures fileFasta.da and fileFasta.ebwt computed.

Input: fileFasta, read length, minimum similarity value beta, threads
 
Output: A .txt/.bin file containing the similarity scores of those reads whose maximum score is greater than beta.
 
*/

#ifndef EBWT
	#define EBWT 1
#endif

#ifndef SMALL
	#define SMALL 0
#endif

typedef std::pair<vector<dataTypeNSeq>, vector<dataTypeNSeq>> PairNReadRef;
typedef vector<dataTypeNSeq> vect_id;

unordered_map<uchar, uchar> umapIUPAC;
std::vector<vector<bool>> umapIUPACcorr;

struct entryDAeBWT{
    
    dataTypeNSeq da;
    dataTypedimAlpha bwt;
    
    // partial order comparison
    bool operator< (const entryDAeBWT& a) const
    {
        if (da == a.da)
            return (umapIUPAC[bwt] < umapIUPAC[a.bwt]);
        
        return  da < a.da;
    }
    
    static bool comp_bwt (const entryDAeBWT& a, const char& s)
    {
        return ( umapIUPAC[a.bwt] < umapIUPAC[s]);
    }
    
    static bool comp_da (const entryDAeBWT& a, const dataTypeNSeq& nr)
    {
        return a.da < nr;
    }
    
    operator dataTypeNSeq() const    {
        return da;
    }
    
    operator uchar() const    {
        return bwt;
    }
    
    friend ostream& operator<<(ostream& os, const entryDAeBWT& op);
    
};

ostream& operator<<(ostream& os, const entryDAeBWT& op)
{
    os << "[" << (int)op.bwt << ", " << (int)op.da << "]";
    return os;
}

bool compareByBWT(entryDAeBWT i, entryDAeBWT j)
{
    if (i.da != j.da)
        return (i.da < j.da);
    else
        return (umapIUPAC[i.bwt] < umapIUPAC[j.bwt]);
}

#if EBWT
vect_id Update_ref_symb(std::vector<entryDAeBWT>::iterator initial_pos, std::vector<entryDAeBWT>::iterator end_pos, dataTypeSim **ref_symb, dataTypeNSeq *nRead){
    
    vect_id list;
    dataTypeNSeq da, da_next;
    dataTypedimAlpha s;
    //Update ref_symb
    while( initial_pos < end_pos ){
        da = initial_pos->da;
        s =  initial_pos->bwt;
        da_next = initial_pos->da;
        for (int i=0; i<ALF; i++)
            ref_symb[da - *nRead][i]=0;
        list.push_back(da - *nRead);
        while ( ( da == da_next ) && ( initial_pos < end_pos ) )
        {
            if(ref_symb[da - *nRead][umapIUPAC[s]]<USim_MAX)
                ref_symb[da - *nRead][umapIUPAC[s]]++;
            initial_pos++;
            da = da_next;
            da_next = initial_pos->da;
            s=initial_pos->bwt;
        }
    }
    return list; //list stores ref IDs in the range [0,nRef)
}

void Analysis_and_updating(vect_id list, std::vector<entryDAeBWT>::iterator initial_pos, std::vector<entryDAeBWT>::iterator end_pos, dataTypeSim **ref_symb, dataTypeSim **SimArray_) {
    vector<dataTypeSim> symb_read(ALF,0);
    vector<dataTypeSim> symb_tmp(ALF,0);
	vector<dataTypeSim> symb_tmp_read(ALF,0);
    
    dataTypeNSeq da, da_next;
    dataTypedimAlpha s;
    //only reads in cluster
    while ( initial_pos < end_pos ){
        da = initial_pos->da;
        da_next = initial_pos->da;
        for (int i=0; i<ALF; i++)
            symb_read[i]=0;
        //Update symb_read
        while ( ( da == da_next ) && ( initial_pos < end_pos ) ) {
            s=initial_pos->bwt;
            symb_read[umapIUPAC[s]]++;
            initial_pos++;
            da = da_next;
            da_next = initial_pos->da;
        }
        //either da != da_next or initial_pos == end_pos
        for (dataTypeNSeq j=0; j<list.size();j++) //range over refs
        {
            dataTypeSim t=0;
            //Count matches between equal symbols
            for (int i=0; i<ALF; i++){
				symb_tmp_read[i]=0;
				symb_tmp[i]=0;
                if(symb_read[i]>ref_symb[list[j]][i]){
                    t+=ref_symb[list[j]][i];
                    symb_tmp_read[i]=symb_read[i]-ref_symb[list[j]][i];
                }
                else{ //symb_read[i]<=ref_symb[j][i]
                    t+=symb_read[i];
                    symb_tmp[i]=ref_symb[list[j]][i]-symb_read[i];
                }
            }
            //symb_tmp_read stores the remaining read symbols (not matched), symb_tmp stores the remaining ref symbols (not matched)
            for(int i=0; i<4;i++){
                for (int a=4; a<ALF-1; a++){
                    if(umapIUPACcorr[a][i]){
                        //Combine reads' A,C,G,Ts with refs' IUPAC symbols
                        if(symb_tmp[a]>0){
                            if(symb_tmp[a]>symb_tmp_read[i]){
                                t+=symb_tmp_read[i];
                                symb_tmp_read[i]=0;
                                symb_tmp[a]-=symb_tmp_read[i];
                            }
                            else{ //symb_tmp[a]<=symb_tmp_read[i]
                                t+=symb_tmp[a];
                                symb_tmp[a]=0;
                                symb_tmp_read[i]-=symb_tmp[a];
                            }
                        }
                        //Combine reads' IUPAC symbols with refs' A,C,G,Ts
                        if(symb_tmp_read[a]>0){
                            if(symb_tmp_read[a]>symb_tmp[i]){
                                t+=symb_tmp[i];
                                symb_tmp_read[a]-=symb_tmp[i];
                                symb_tmp[i]=0;
                            }
                            else{ //symb_tmp_read[a]<=symb_tmp[i]
                                t+=symb_tmp_read[a];
                                symb_tmp[i]-=symb_tmp_read[a];
                                symb_tmp_read[a]=0;
                            }
                        }//end-if
                    }//end-if
                }//end-for
            }//end-for
            if(t>0){
                #if OMP
                #pragma omp atomic
                #endif
                //Updating
                SimArray_[da][list[j]]+=t;
            }
        }//end-for
    }//end-while
	symb_read.shrink_to_fit();
	symb_tmp.shrink_to_fit();
	symb_tmp_read.shrink_to_fit();
}
#else
void Analysis_and_updating(std::vector<dataTypeNSeq>::iterator &lowBounds, std::vector<dataTypeNSeq>::iterator &upBounds, std::vector<dataTypeNSeq>::iterator &pos_nReadInCluster, vector<dataTypeSim> &read_symb, vector<dataTypeSim> &ref_symb, dataTypeNSeq nRead, dataTypeSim **SimArray_){
{
    PairNReadRef Pairs;
    
    dataTypeNSeq da, da_next;
	//Reads
    while ( lowBounds < pos_nReadInCluster ){
		
		da = *lowBounds;
        da_next = *lowBounds;
        read_symb[ da ]=0;
        (Pairs.first).push_back(da);
        while ( ( da == da_next ) && ( lowBounds < pos_nReadInCluster ) )
        {
            read_symb[da] ++;
            lowBounds++;
            da = da_next;
			da_next = *lowBounds;
        }
    }
    //Refs      
    while ( lowBounds < upBounds ){
		
		da = *lowBounds;
        da_next = *lowBounds;

        ref_symb[ da -nRead ]=0;
        (Pairs.second).push_back(da-nRead);
        while ( ( da == da_next ) && ( lowBounds < upBounds ))
        {
            if (ref_symb[ da-nRead ] <USim_MAX)
                ref_symb[ da-nRead ]++;
            lowBounds++;
            da = da_next;
			da_next = *lowBounds;
        }
    }
    
    dataTypeSim t=0;
	
	for (dataTypeNSeq i=0; i<Pairs.first.size();i++) //range over reads
    {
        for (dataTypeNSeq j=0; j<Pairs.second.size();j++)//range over references
        {
            //Take the minimum
            if (read_symb[Pairs.first[i]] > ref_symb[Pairs.second[j]] )
                t= ref_symb[Pairs.second[j]];
            else
                t= read_symb[Pairs.first[i]];

            //Update SimArray_
            if(t>0){
				#if OMP
					#pragma omp atomic
				#endif
                    SimArray_[Pairs.first[i]][Pairs.second[j]]+=t;
            }
        }//end-for
    }//end-for
	return;
}
#endif

#if EBWT 
dataTypeNChar clusterAnalyze(dataTypeNChar chk, FILE *InFileCluster, FILE *InDA, FILE *InBWT, dataTypeNSeq nRead, dataTypeNSeq nRef, dataTypeSim **SimArray_)
#else
dataTypeNChar clusterAnalyze(dataTypeNChar chk, FILE *InFileCluster, FILE *InDA, dataTypeNSeq nRead, dataTypeNSeq nRef, dataTypeSim **SimArray_)
#endif
{
	dataTypeNChar counter=0;
    dataTypeNChar numcharCluster, numcharDA;
	ElementCluster* clusterbuffer= new ElementCluster[BUFFERCLUSTER];
    dataTypeNSeq* DAbuffer= new dataTypeNSeq[sizeMaxBuf];
	
	#if EBWT
		dataTypeNChar numcharBWT;
		dataTypedimAlpha* BWTbuffer = new dataTypedimAlpha[sizeMaxBuf];
		vector<entryDAeBWT> elebuffer;
		std::vector<entryDAeBWT>::iterator lowBounds, upBounds;
	#else
		vector<dataTypeNSeq> elebuffer;
		std::vector<dataTypeNSeq>::iterator lowBounds, upBounds, pos_nReadInCluster;
	#endif
	
	elebuffer.resize(sizeMaxBuf);
    
	#if EBWT
		dataTypeSim **ref_symb = (dataTypeSim**) malloc(nRef * sizeof(dataTypeSim*)); 
		for(dataTypeNSeq i=0; i<nRef; ++i){
		  ref_symb[i]  = (dataTypeSim*) malloc(ALF * sizeof(dataTypeSim));
		  for(dataTypeNSeq j=0; j<ALF; ++j)
			ref_symb[i][j] = 0;
		}
	#else
		vector<dataTypeSim> read_symb;
		vector<dataTypeSim> ref_symb;
		read_symb.resize(nRead);
		ref_symb.resize(nRef);
	#endif
    
    //Read InFileCluster
    if(BUFFERCLUSTER<chk)
		numcharCluster=fread(clusterbuffer,sizeof(ElementCluster),BUFFERCLUSTER,InFileCluster);
    else
        numcharCluster=fread(clusterbuffer,sizeof(ElementCluster),chk,InFileCluster);

        
    while(numcharCluster>0)
    {
        for (dataTypeNChar c=0;c<numcharCluster; c++)
		{
            //Read cluster DA and EBWT
			fseek(InDA, (clusterbuffer[c].pStart)*sizeof(dataTypeNSeq), SEEK_SET);
            #if EBWT
				numcharDA=fread(DAbuffer,sizeof(dataTypeNSeq),clusterbuffer[c].len,InDA);
				fseek(InBWT, (clusterbuffer[c].pStart)*sizeof(dataTypedimAlpha), SEEK_SET);
				numcharBWT=fread(BWTbuffer,sizeof(dataTypedimAlpha),clusterbuffer[c].len,InBWT);
			#else
				numcharDA=fread(&elebuffer[0],sizeof(dataTypeNSeq),clusterbuffer[c].len, InDA);
			#endif

			//Sort
			#if EBWT
                //Merge DA and EBWT
				for(dataTypeNChar e=0;e<numcharBWT;e++) {
					elebuffer[e].bwt=BWTbuffer[e];
					elebuffer[e].da=DAbuffer[e];
				}
                std::sort (elebuffer.begin(), elebuffer.begin()+ numcharBWT, compareByBWT);//w.r.t. da, and then ebwt
			#else
				std::sort (elebuffer.begin(), elebuffer.begin()+ numcharDA);//w.r.t. da
			#endif
                
			//Set lowBounds and upBounds
            #if EBWT
                lowBounds=elebuffer.begin();
                //upBounds points the first genome occurrence in the cluster
                upBounds= std::lower_bound (lowBounds, elebuffer.begin()+ numcharBWT, dataTypeNSeq(entryDAeBWT{nRead}), *entryDAeBWT::comp_da);
                vect_id id_Refs=Update_ref_symb(upBounds, elebuffer.begin()+numcharBWT, ref_symb, &nRead);
                Analysis_and_updating(id_Refs, lowBounds, upBounds, ref_symb,SimArray_);
			#else
				lowBounds=elebuffer.begin();
				upBounds=elebuffer.begin()+ numcharDA;
				//pos_nReadInCluster points the first genome occurrence in the cluster
				pos_nReadInCluster=std::lower_bound (lowBounds, upBounds, nRead);
                Analysis_and_updating(lowBounds, upBounds, pos_nReadInCluster, read_symb, ref_symb, nRead, SimArray_);
            #endif
            
        }//end-for clusters
        counter+=numcharCluster;
        chk-=numcharCluster;
        if(chk>0 && BUFFERCLUSTER<chk)
            numcharCluster=fread(clusterbuffer,sizeof(ElementCluster),BUFFERCLUSTER,InFileCluster);
        else
            numcharCluster=fread(clusterbuffer,sizeof(ElementCluster),chk,InFileCluster);
    }//end-while
    elebuffer.shrink_to_fit();
	#if EBWT
	for(dataTypeNSeq i=0; i<nRef; ++i)
        free(ref_symb[i]);
        free(ref_symb);
	#else
		ref_symb.shrink_to_fit();
		read_symb.shrink_to_fit();
	#endif
	return counter;
}

#if BIN
void clusterChoose(FILE* OutFileClustersBin, FILE* OutFileClustersPos, dataTypelenSeq &norm, float &beta, dataTypeNSeq nRead, dataTypeNSeq nRef, dataTypeSim **SimArray_)
#else
void clusterChoose(FILE* OutFileClusters, dataTypelenSeq &norm, float &beta, dataTypeNSeq nRead, dataTypeNSeq nRef, dataTypeSim **SimArray_)
#endif
{
	dataTypeSim maxSim=0; //maximum similarity value

	vector<type_cluster> outputPairs;
	type_cluster pair;
	dataTypeNSeq index=0;
	float normSim;
	#if BIN
      size_t total=0, curr=0;
      size_t null=0;
    
      //Sentinel position
      pair_sim outPair;
      outPair.sim = 0.0;
      outPair.idRef = 0;//size of the list
      fwrite(&outPair, sizeof(pair_sim), 1, OutFileClustersBin);
      total++;
      curr=total;
    #endif
  
	while (index < nRead) //Stop reading the matrix SimArray_ when all reads are processed
    {
		maxSim=0;
		normSim=0;
		outputPairs.clear();
		for(dataTypeNSeq j=0;j<nRef;j++)
        {
			pair.sim= SimArray_[index][j];
			if (pair.sim>maxSim) //Update the maximum similarity
				maxSim=pair.sim;
      
			if (pair.sim>0) //Keep the pair (idRef,sim) only if the similarity is not 0
			{
				pair.idRef=j;
				outputPairs.push_back(pair);
			}
          
		}
		//write file only if the maximum similarity value is greater than beta
		normSim=static_cast<float> (maxSim)/norm; //normalized maximum similarity
		if (normSim>beta)
		{
			#if BIN
                pair_sim outPair;
                outPair.sim = normSim;
                outPair.idRef = outputPairs.size();//size of the list
                fwrite(&outPair, sizeof(pair_sim), 1, OutFileClustersBin);
                total++;
			#else 
                fprintf(OutFileClusters, "%.5f", normSim);
			#endif
			for (std::vector<type_cluster>::iterator it = outputPairs.begin() ; it != outputPairs.end(); ++it){
				normSim=static_cast<float> ((*it).sim)/norm;
				#if BIN
				outPair.sim = normSim;
				outPair.idRef = (*it).idRef;   
                fwrite(&outPair, sizeof(pair_sim), 1, OutFileClustersBin);
                total++;
				#else
				fprintf(OutFileClusters, "\t%u\t%.5f",(*it).idRef, normSim);
				#endif
			}
			#if BIN
            fwrite(&curr, sizeof(size_t), 1, OutFileClustersPos);
            curr=total;
			#endif
		}
		else{
			#if BIN
			fwrite(&null, sizeof(size_t), 1, OutFileClustersPos);
			#endif 
		}
    
		//Read processed ---> new line
		#if BIN==0
            fprintf(OutFileClusters,"\n");
		#endif
            index++;
	}//end-while
	
	#if BIN
        cout<<"Output (in bytes): "<<total<<"\t"<<sizeof(pair_sim)<<endl;
	#endif
	
	return;
}


int main(int argc, char **argv) {
    
    umapIUPAC['A'] = 0;
    umapIUPAC['C'] = 1;
    umapIUPAC['G'] = 2;
    umapIUPAC['T'] = 3;
    umapIUPAC['R'] = 4;
    umapIUPAC['Y'] = 5;
    umapIUPAC['S'] = 6;
    umapIUPAC['W'] = 7;
    umapIUPAC['K'] = 8;
    umapIUPAC['M'] = 9;
    umapIUPAC['B'] = 10;
    umapIUPAC['D'] = 11;
    umapIUPAC['H'] = 12;
    umapIUPAC['V'] = 13;
    umapIUPAC['N'] = 14;
    umapIUPAC['\0'] = 15;
    
    umapIUPACcorr.resize(ALF);
    umapIUPACcorr[0] = {1,0,0,0};
    umapIUPACcorr[1] = {0,1,0,0};
    umapIUPACcorr[2] = {0,0,1,0};
    umapIUPACcorr[3] = {0,0,0,1};
    umapIUPACcorr[4] = {1,0,1,0};  //A, G
    umapIUPACcorr[5] = {0,1,0,1};  //C, T
    umapIUPACcorr[6] = {0,1,1,0};  //G, C
    umapIUPACcorr[7] = {1,0,0,1};  //A, T
    umapIUPACcorr[8] = {0,0,1,1};  //G, T
    umapIUPACcorr[9] = {1,1,0,0};  //A, C
    umapIUPACcorr[10] = {0,1,1,1};  //C, G, T
    umapIUPACcorr[11] = {1,0,1,1};  //A, G, T
    umapIUPACcorr[12] = {1,1,0,1};  //A, C, T
    umapIUPACcorr[13] = {1,1,1,0};  //A, C, G
    umapIUPACcorr[14] = {1,1,1,1};  //A, C, G, T
    
	#if OMP
    double d_total, d_refine;
	#else
    time_t t_refine=0, t_total=0;
    clock_t c_refine=0, c_total=0;
	#endif
	 
	if( argc != 5)
	{
		std::cerr << "Error usage " << argv[0] << " fileFasta readLen beta threads"  << std::endl;
		exit(1);
	}

	int num_threads=1;
	sscanf(argv[4], "%d", &num_threads);

	#if OMP
		omp_set_num_threads(num_threads);

    int threads=0;
		#pragma omp parallel
		{
			threads = omp_get_num_threads();
			if(omp_get_thread_num()==0)
				printf("Number of threads: %d\n", threads);
		}
		printf("Number of processors: %d\n", omp_get_num_procs());
	#endif

	string fileFasta=argv[1];
	dataTypeSim readLen;
	float beta;
	sscanf(argv[2], "%hhu", &readLen);
	sscanf(argv[3], "%f", &beta);
	
	//Check read length
	if ( (readLen>USim_MAX) && (dataTypeNumSim==0))
	{
		std::cerr << "Error Usage: readLen <= USim_MAX, please change settings of dataTypeNumSim in Tools.h" << std::endl;
		exit(1);
	}
	
	//Read auxiliary file
	dataTypeNSeq numRead;
	dataTypeNSeq numRef;
	dataTypelenSeq minLCP;
	dataTypeNChar maxLen, nClusters;
	string fileaux=fileFasta.substr(0,fileFasta.find(".fasta"))+".out";
    
	FILE * outAux = fopen(fileaux.c_str(), "rb");
	if (outAux==NULL) {
		std::cerr << "Error opening " << fileaux << "." << std::endl;
		printf("fopen failed, errno = %d\n", errno);
		exit (EXIT_FAILURE);
	}
	
	fread(&numRead,sizeof(dataTypeNSeq),1,outAux);
	fread(&numRef,sizeof(dataTypeNSeq),1,outAux);
	fread(&minLCP,sizeof(dataTypelenSeq),1,outAux);
	fread(&maxLen,sizeof(dataTypeNChar),1,outAux);
  	fread(&nClusters,sizeof(dataTypeNChar),1,outAux);
   
	fclose(outAux);
    
	cout << "numRead: " << numRead << ", numRef: " << numRef << ", minLCP: " << minLCP << ", nClusters: " << nClusters << endl;

	dataTypelenSeq norm=readLen+1-minLCP; //to normalize similarity values
	
	//Check maximum clester length
	if (maxLen>sizeMaxBuf)
	{
    	cerr << "Error Usage: maximum cluster size is " << maxLen << " greater than sizeMaxBuf, please increase sizeMaxBuf in Tools.h" << endl;
		exit(1);
	}
	
	
    //Files .clrs and .da (and .ebwt #if EBWT)
	std::string fnCluster, fnDA;
    fnDA=fileFasta+".da";
	#if EBWT
		std::string fnBWT=fileFasta+".ebwt";
		FILE *InEBWT[num_threads];
	#endif
	std::stringstream ssin;
	ssin << fileFasta << "." << minLCP << ".clrs\0";
	fnCluster=ssin.str();
	
	FILE *InFileCluster[num_threads];
	FILE *InDA[num_threads];
	int t_id=0;
	#if OMP
    for(;t_id<num_threads; t_id++)
    #endif
    {
        InFileCluster[t_id] = fopen(fnCluster.c_str(), "rb");
        if ((InFileCluster[t_id]==NULL)){
			std::cerr << "Error opening " << fnCluster << "." << std::endl;
			printf("fopen failed, errno = %d\n", errno);
			exit (EXIT_FAILURE);
		}
		InDA[t_id] = fopen(fnDA.c_str(), "rb");
		if ((InDA[t_id]==NULL)){
			std::cerr << "Error opening " << fnDA << "." << std::endl;
			printf("fopen failed, errno = %d\n", errno);
			exit (EXIT_FAILURE);
		}
		#if EBWT
		InEBWT[t_id] = fopen(fnBWT.c_str(), "rb");
		if ((InEBWT[t_id]==NULL)){
			std::cerr << "Error opening " << fnBWT << "." << std::endl;
			printf("fopen failed, errno = %d\n", errno);
			exit (EXIT_FAILURE);
		}
		#endif
	}
	
    //Similarity arrays
    dataTypeSim **SimArray_ = (dataTypeSim**) malloc(numRead * sizeof(dataTypeSim*)); 
    for(dataTypeNSeq i=0; i<numRead; ++i){
      SimArray_[i]  = (dataTypeSim*) malloc(numRef * sizeof(dataTypeSim));
      for(dataTypeNSeq j=0; j<numRef; ++j)
        SimArray_[i][j] = 0;
    }
    
	#if OMP
		d_total = omp_get_wtime();
	#else
		time_start(&t_total, &c_total); //start time
	#endif
  
	cerr << "Computing similarity arrays SimArray_i[1,numRead]..." << endl;
    
    #if OMP
      d_refine = omp_get_wtime();
    #else
      time_start(&t_refine, &c_refine); //start time for clusterRefine
    #endif
  
	dataTypeNChar AnaCluster=0;
  
	//START PARALLEL
    #if OMP
    #pragma omp parallel default(shared) reduction(+:AnaCluster)
    #endif
    {
        int tid = 0;
		int numthreads = 1;
		#if OMP
        tid=omp_get_thread_num();//id_thread
        numthreads = omp_get_num_threads();
		#endif

        dataTypeNChar chunk=(nClusters/numthreads);
        dataTypeNChar startRead=tid*chunk;
        
        dataTypeNChar endRead=(tid+1)*chunk;
        if(tid==numthreads-1)
            endRead=nClusters;
		
		chunk= endRead-startRead;
        
        double start=omp_get_wtime();
		
		fseek(InFileCluster[tid], startRead*sizeof(ElementCluster), SEEK_SET);
		
		/**Start analysis**/
		#if EBWT
			AnaCluster = clusterAnalyze(chunk,InFileCluster[tid],InDA[tid], InEBWT[tid],numRead,numRef, SimArray_);
		#else
			AnaCluster = clusterAnalyze(chunk,InFileCluster[tid],InDA[tid], numRead,numRef, SimArray_);
		#endif
		/**End analysis**/
		
		#pragma omp critical
        {
            std::cerr << "TIME THREAD " << tid << " = " << omp_get_wtime()-start << "(in seconds)\n";
        }
    }//end-pragma
  
	#if SMALL
    cout << "***FINAL***\n";
    for (dataTypeNSeq i=0; i<numRead; i++)
    {
        for(dataTypeNSeq j=0; j<numRef; j++)
            cout << (int)SimArray_[i][j] << "\t";
        cout << "\n";
    }
    cout << "***********\n";
    #endif

    #if OMP
      fprintf(stderr,"TIME clusterAnalyze: %.6lf\n", omp_get_wtime()-d_refine);
    #else
      fprintf(stderr,"TIME clusterAnalyze: %.6lf\n", time_stop(t_refine, c_refine));
    #endif
	
    //Close .clrs and .da (and .ebwt #if EBWT)
	t_id=0;
    #if OMP
    for(;t_id<num_threads; t_id++)
    #endif
    {
        fclose(InFileCluster[t_id]);
        fclose(InDA[t_id]);
		#if EBWT
			fclose(InEBWT[t_id]);
		#endif
    }
    //Open Output file
	string fnF=fileFasta+".res";
	
    #if BIN
		FILE *fdResultPos, *fdResultBin;
		string sResultPos(fnF.c_str());
		string sResultBin(fnF.c_str());
		sResultPos.append(".pos");
		sResultBin.append(".bin");

		fdResultPos = fopen(sResultPos.c_str(), "wb");
		cerr << "Writing " << sResultPos << endl;
		fdResultBin = fopen(sResultBin.c_str(), "wb");
		cerr << "Writing " << sResultBin << endl;
    #else
		FILE *fdResultF;
		string sResultTxt(fnF.c_str());
		sResultTxt.append(".txt");
		fdResultF=fopen(sResultTxt.c_str(), "w");
		if(fdResultF==NULL) {
			std::cerr << "Error opening " << sResultTxt << "." << std::endl;
			printf("fopen failed, errno = %d\n", errno);
			exit (EXIT_FAILURE);
		}
		cerr << "Writing " << sResultTxt << endl;
    #endif
	
    #if OMP
      d_refine = omp_get_wtime();
    #else
      time_start(&t_refine, &c_refine); //start time for clusterChoose
    #endif
	
	
    /**Start writing**/
	#if BIN
		clusterChoose(fdResultBin, fdResultPos, norm, beta, numRead, numRef, SimArray_);
	#else
		clusterChoose(fdResultF, norm, beta, numRead, numRef, SimArray_);
	#endif
    /**End writing**/
    
    #if OMP
      fprintf(stdout,"Time: %.6lf\n", omp_get_wtime()-d_refine);
    #else
      fprintf(stderr,"TIME clusterChoose: %.6lf\n", time_stop(t_refine, c_refine));
    #endif
    
    #if BIN
      fclose(fdResultBin);
      fclose(fdResultPos);
    #else
      fclose(fdResultF);
    #endif


    cout << "Cluster analysis completed with beta=" << beta << "." << endl;
    cout << "Number of clusters: " << AnaCluster << "." << endl;

    #if OMP
      fprintf(stdout,"Time: %.6lf\n", omp_get_wtime()-d_total);
    #else
      fprintf(stdout,"Time: %.6lf\n", time_stop(t_total, c_total));
    #endif

    #if FAST
    for(dataTypeNSeq i=0; i<numRead; ++i)
      free(SimArray_[i]);
    free(SimArray_);
    #endif

	return 0;
}
