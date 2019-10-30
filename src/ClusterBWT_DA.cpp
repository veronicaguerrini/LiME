#include "Tools.h"

/*
Step 2. analyzes alpha-clusters in order to evaluate a degree of similarity between any read and any genome in the collection.
If EBWT=1, the similarity scores are computed by using both the document array and the ebwt. 
If EBWT=0, the similarity scores are computed by using only the document array.

As preprocessing, it need to have the two data structures fileFasta.da and fileFasta.ebwt computed.

Input: fileFasta, read length, minimum similarity value beta, output, threads
 
Output: A .txt/.bin file containing the similarity scores of those reads whose maximum score is greater than beta.
 
*/

#ifndef EBWT
	#define EBWT 1
#endif

#ifndef CHECK
	#define CHECK 0
#endif

#ifndef SMALL
	#define SMALL 0
#endif

typedef std::pair<vector<dataTypeNSeq>, vector<dataTypeNSeq>> PairNReadRef;

struct entryDAeBWT{
	
	dataTypedimAlpha bwt;
	dataTypeNSeq da;

	// partial order comparison
	bool operator< (const entryDAeBWT& a) const
	{
		if (bwt == a.bwt) return da < a.da;
		return bwt < a.bwt;
	}
    
    static bool comp_bwt (const entryDAeBWT& a, const char& s)
    {
        return a.bwt < s;
    }
	
	operator int() const	{
		return da;
	}
	
	operator uchar() const	{
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
	if (i.bwt != j.bwt)
		return (i.bwt < j.bwt);
	else
		return (i.da < j.da);
}
#if EBWT
PairNReadRef Analysis(std::vector<entryDAeBWT>::iterator &lowBounds, std::vector<entryDAeBWT>::iterator &upBounds, std::vector<entryDAeBWT>::iterator &pos_nReadInCluster, vector<dataTypeNSeq> &read_symb, vector<dataTypeNSeq> &ref_symb, dataTypeNSeq nRead)
#else
PairNReadRef Analysis(std::vector<dataTypeNSeq>::iterator &lowBounds, std::vector<dataTypeNSeq>::iterator &upBounds, std::vector<dataTypeNSeq>::iterator &pos_nReadInCluster, vector<dataTypeNSeq> &read_symb, vector<dataTypeNSeq> &ref_symb, dataTypeNSeq nRead)
#endif
{
    PairNReadRef Pairs;
    
    dataTypeNSeq da, da_next;
	//Reads
    while ( lowBounds < pos_nReadInCluster ){
		#if EBWT
			da = lowBounds->da;
			da_next = lowBounds->da;
		#else
			da = *lowBounds;
            da_next = *lowBounds;
		#endif
        read_symb[ da ]=0;
        (Pairs.first).push_back(da);
        while ( ( da == da_next ) && ( lowBounds < pos_nReadInCluster ) )
        {
            read_symb[da] ++;
            lowBounds++;
            da = da_next;
			#if EBWT
				da_next = lowBounds->da;
			#else
				da_next = *lowBounds;
			#endif
        }
    }
    //Refs      
    while ( lowBounds < upBounds ){
		#if EBWT
			da = lowBounds->da;
            da_next = lowBounds->da;
		#else
			da = *lowBounds;
            da_next = *lowBounds;
		#endif

        ref_symb[ da -nRead ]=0;
        (Pairs.second).push_back(da-nRead);
        while ( ( da == da_next ) && ( lowBounds < upBounds ))
        {
            ref_symb[ da-nRead ]++;
            lowBounds++;
            da = da_next;
			#if EBWT
				da_next = lowBounds->da;
			#else
				da_next = *lowBounds;
			#endif
        }
    }
    
    return Pairs;
}

void Updating(PairNReadRef &Pairs, vector<dataTypeNSeq> &read_symb, vector<dataTypeNSeq> &ref_symb, dataTypeSim **SimArray_)
{
	//counter (t cannot be greater than 255, since each read has length 101)
    dataTypeSim t=0;
	
	for (dataTypeNSeq i=0; i<Pairs.first.size();i++) //range over reads
    {
        for (dataTypeNSeq j=0; j<Pairs.second.size();j++)//range over references
        {
            //Take the minimum
            if (read_symb[Pairs.first[i]] > ref_symb[Pairs.second[j]] )
            {
                #if CHECK
                assert(ref_symb[Pairs.second[j]]<USim_MAX);
                #endif
                t= ref_symb[Pairs.second[j]];
            }
            else
            {
                #if CHECK
                assert(read_symb[Pairs.first[i]]<USim_MAX);
                #endif
                t= read_symb[Pairs.first[i]];
            }
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
		std::vector<entryDAeBWT>::iterator lowBounds, upBounds, pos_nReadInCluster;
	#else
		vector<dataTypeNSeq> elebuffer;
		std::vector<dataTypeNSeq>::iterator lowBounds, upBounds, pos_nReadInCluster;
	#endif
	
	elebuffer.resize(sizeMaxBuf);
    PairNReadRef  idInCluster;  
	vector<dataTypeNSeq> read_symb;
    vector<dataTypeNSeq> ref_symb;
    read_symb.resize(nRead);
    ref_symb.resize(nRef);
	
	#if EBWT
        uchar s='\0', s_next='\0';
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
				for(dataTypeNChar e=0;e<numcharBWT;e++)
				{
					elebuffer[e].bwt=BWTbuffer[e];
					elebuffer[e].da=DAbuffer[e];
				}
                std::sort (elebuffer.begin(), elebuffer.begin()+ numcharBWT, compareByBWT);//w.r.t. ebwt, and then da
			#else
				std::sort (elebuffer.begin(), elebuffer.begin()+ numcharDA);//w.r.t. da
			#endif
                
			//Set lowBounds and upBounds
            #if EBWT
                lowBounds=elebuffer.begin();
                upBounds=elebuffer.begin();
                s=lowBounds->bwt;
                s_next=s+1;
                upBounds= std::lower_bound (lowBounds, elebuffer.begin()+ numcharBWT, uchar(entryDAeBWT{s_next}), *entryDAeBWT::comp_bwt);
			#else
				lowBounds=elebuffer.begin();
				upBounds=elebuffer.begin()+ numcharDA;
			#endif
                
			//Find the position of the first genome
			pos_nReadInCluster=std::lower_bound (lowBounds, upBounds, nRead);
            
			#if EBWT
			if ( (pos_nReadInCluster < upBounds) &&  (pos_nReadInCluster > lowBounds) )//At least one read and one ref
			#endif
            {
                idInCluster=Analysis(lowBounds, upBounds, pos_nReadInCluster, read_symb, ref_symb, nRead);
				Updating(idInCluster, read_symb, ref_symb, SimArray_);                   
            }
			
			#if EBWT	
                while (upBounds<elebuffer.begin()+ numcharBWT)
                {//There are still other symbols
					s = upBounds->bwt;
                    s_next = s+1;
                    lowBounds=upBounds;
                    upBounds= std::lower_bound (lowBounds, elebuffer.begin()+ numcharBWT, uchar(entryDAeBWT{s_next}), *entryDAeBWT::comp_bwt);
                    pos_nReadInCluster=std::lower_bound (lowBounds, upBounds, nRead);
                        
                    //Analyze iff there is at least one read vs one ref
                    if ( (pos_nReadInCluster < upBounds) &&  (pos_nReadInCluster > lowBounds) )
                    {
                        idInCluster=Analysis(lowBounds, upBounds, pos_nReadInCluster, read_symb, ref_symb, nRead);
                        Updating(idInCluster, read_symb, ref_symb, SimArray_);   
                    }//end-if analyzed
                }//end-while symbol
			#endif
                
        }//end-for clusters
        counter+=numcharCluster;
        chk-=numcharCluster;
            
        if(chk>0 && BUFFERCLUSTER<chk)
            numcharCluster=fread(clusterbuffer,sizeof(ElementCluster),BUFFERCLUSTER,InFileCluster);
        else
            numcharCluster=fread(clusterbuffer,sizeof(ElementCluster),chk,InFileCluster);

    }//end-while
    idInCluster.first.shrink_to_fit();
    idInCluster.second.shrink_to_fit();
    ref_symb.shrink_to_fit();
    read_symb.shrink_to_fit();
    elebuffer.shrink_to_fit();

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
			for (std::vector<type_cluster>::iterator it = outputPairs.begin() ; it != outputPairs.end(); ++it)
			{
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
		else
		{
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
    
	#if DEBUG
		cout << "numRead: " << numRead << " numRef: " << numRef << " minLCP: " << minLCP << endl;
	#endif
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
		#if CHECK
			stringstream ss;
			ss << "Id thread: " << tid << " working on " << chunk << " elements.\n";
			ss << "startRead " << startRead << "\tendRead " << endRead << ".\n";
			cout << ss.str();
        #endif
		
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
      fprintf(stderr,"TIME clusterRefine: %.6lf\n", omp_get_wtime()-d_refine);
    #else
      fprintf(stderr,"TIME clusterRefine: %.6lf\n", time_stop(t_refine, c_refine));
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
	string fnF="Res_"+fileFasta.substr(0,fileFasta.find(".fasta"));
	#if EBWT
		fnF+="_EBWT1";
	#else
		fnF+="_EBWT0";
	#endif
	
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
