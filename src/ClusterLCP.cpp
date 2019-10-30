#include "Tools.h"

/*
Step 1. detects alpha-clusters, i.e., blocks containing symbols belonging both to reads and to genomes, and whose associated suffixes share a common context of minimum length alpha.

As preprocessing, it need to have the two data structures fileFasta.lcp and fileFasta.da computed.

Input: fileFasta, total number of reads, total number of genomes, alpha, threads.
		
Output: fileFasta.alpha.clrs containing a pair ElementCluster (pStart,len) for each alpha-cluster detected; one auxiliary file.

*/

void StartOrRemain(dataTypeNSeq nRead, dataTypeNSeq* bufferEle, dataTypeNChar indexbuffer, dataTypeNChar ind,bool &init, ElementCluster &cluster, bool &nR, bool &nG)
{	
	if (not init) //Start a new cluster
	{
		init=true;
		cluster.pStart = ind-1;
		if (bufferEle[indexbuffer]<nRead)
			nR=true;
		else
			nG=true;
	}
	if (not (nR && nG))//Check if a read and a genome is in
	{
		if (bufferEle[indexbuffer+1]<nRead)
			nR=true;
		else
			nG=true;
	}
}

int Close(ElementCluster &cluster, dataTypeNChar ind, dataTypeNChar &lengthClust,vector<ElementCluster> &vOutput)
{
	cluster.len = ind - cluster.pStart;
	if (lengthClust<cluster.len)
		lengthClust=cluster.len;
        
	vOutput.push_back(cluster);
  
	return 1;
}
                          


int main(int argc, char **argv) {

	#if OMP
		double d_total;
	#else
		time_t t_refine=0, t_total=0;
		clock_t c_refine=0, c_total=0;
	#endif

	if( argc != 6 )
	{
		std::cerr << "Error usage: " << argv[0] << " fileFasta numReads numGenomes alpha threads" << std::endl;
		exit(1);
	}
    
	string fileFasta=argv[1]; 
	dataTypelenSeq alpha;
	dataTypeNSeq numReads, numGenomes;
	
	sscanf(argv[2], "%u", &numReads);
	sscanf(argv[3], "%u", &numGenomes);
	sscanf(argv[4], "%u", &alpha);
	
    int num_threads=1;
	sscanf(argv[5], "%d", &num_threads);
	
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
	
	//Open files
	string fnLCP, fnDA, fileOutput;
	fnLCP = fileFasta+".lcp\0";
	fnDA =fileFasta+".da\0";
	std::stringstream ssout;
	ssout << fileFasta << "." << alpha << ".clrs\0";
	fileOutput=ssout.str();
	
	FILE *OutCluster=fopen(fileOutput.c_str(), "w");
	if(OutCluster==NULL) {
		cerr << "Error opening " << fileOutput << ".";
		exit(1);
	}
	
	FILE *InLCP[num_threads];
	FILE *InDA[num_threads];
	int tid=0;
	#if OMP
    for(;tid<num_threads; tid++)
    #endif
    {
        InLCP[tid] = fopen(fnLCP.c_str(), "rb");
		InDA[tid] = fopen(fnDA.c_str(), "rb");
        if(!tid){
          std::cout << "\n\t" << fnLCP;
          std::cout << "\n\t" << fnDA;
        }
    }
	
	if ((InLCP[0]==NULL)){
		std::cerr << "Error opening " << fnLCP << "." << std::endl;
		printf("fopen failed, errno = %d\n", errno);
		exit (EXIT_FAILURE);
	}
	if ((InDA[0]==NULL)){
		std::cerr << "Error opening " << fnDA << "." << std::endl;
		printf("fopen failed, errno = %d\n", errno);
		exit (EXIT_FAILURE);
	}
	//InLCP dimension
	fseek(InLCP[0], 0, SEEK_END);
	dataTypeNChar sizeInLCP=ftell(InLCP[0])/sizeof(dataTypeNSeq);
	
	//Clustering
	#if OMP
		d_total = omp_get_wtime();
	#else
		time_start(&t_total, &c_total); //start time
	#endif
	
	dataTypeNChar maxLen=0, nClusters=0;
	
	//START PARALLEL
    #if OMP
    #pragma omp parallel default(shared) reduction(max:maxLen) reduction(+:nClusters)
    #endif
    {
        int t_id = 0;
		int numthreads = 1;
    #if OMP
        t_id=omp_get_thread_num();//id_thread
        numthreads = omp_get_num_threads();
    #endif
        
        dataTypeNChar chunk=(sizeInLCP/numthreads);
		dataTypeNChar bufferSize=(BUFFERLCPSIZE/numthreads);
		assert(bufferSize>0);
        
        dataTypeNChar startRead=t_id*chunk-1;
        dataTypeNChar endRead=(t_id+1)*chunk;
		
		if(t_id==0)
			startRead=0;
		
        if(t_id==numthreads-1)
            endRead=sizeInLCP;
    #if OMP
		double start=omp_get_wtime();
    #endif
        
		fseek(InLCP[t_id], startRead*sizeof(dataTypelenSeq), SEEK_SET);
		fseek(InDA[t_id], startRead*sizeof(dataTypeNSeq), SEEK_SET);

		//Output file contains a collection of pairs (pStart, len) each one corresponding to one alpha-cluster
		ElementCluster cluster;
		cluster.pStart=0;
		cluster.len=0;
		vector<ElementCluster> vOutput;
        
		bool init=false; //init is true if a cluster is open
		bool nR=false;	//nR is true when at least a read is in
		bool nG=false;	//nG is true when at least a genome is in
		
		//To read LCP and DA files
		dataTypeNChar numcharLCP;
		dataTypelenSeq* bufferLCP = new dataTypelenSeq[bufferSize];

		dataTypeNChar numcharDA;
		dataTypeNSeq *bufferEle = new dataTypeNSeq[bufferSize+1];
		
		dataTypeNChar index=1;
		if(tid!=0)
			index=tid*chunk;
		

		numcharLCP=fread(&bufferLCP[0],sizeof(dataTypelenSeq),1,InLCP[t_id]);
		numcharDA=fread(&bufferEle[0],sizeof(dataTypeNSeq),1,InDA[t_id]);
		
		chunk=endRead-startRead-1;
		
		while ((chunk>0) && (bufferLCP[0]>=alpha))
		{
			numcharLCP=fread(&bufferLCP[0],sizeof(dataTypelenSeq),1,InLCP[t_id]);
			numcharDA=fread(&bufferEle[0],sizeof(dataTypeNSeq),1,InDA[t_id]);
			chunk--;
			index++;
		}
		
		if(bufferLCP[0]<alpha)
		{
			while(numcharLCP>0)
			{
				if(bufferSize>=chunk)
					bufferSize=chunk;
					
				numcharLCP = fread(bufferLCP,sizeof(dataTypelenSeq),bufferSize,InLCP[t_id]);
				numcharDA = fread(bufferEle+1,sizeof(dataTypeNSeq),bufferSize,InDA[t_id]);
                
				for(dataTypeNChar indexbuffer=0; indexbuffer<numcharLCP; indexbuffer++)
				{
					if(bufferLCP[indexbuffer]>=alpha) //Start or remain in a cluster
						StartOrRemain(numReads, bufferEle, indexbuffer, index, init, cluster, nR, nG);
					else    //End a cluster
					{
						if (init && nR && nG)
							nClusters+=Close(cluster, index, maxLen,vOutput);
						
						init=false, nR=false, nG=false;
						cluster.pStart=0, cluster.len=0;
					}
					index++;
				}//end-for
				
                #if OMP
					#pragma omp critical
				#endif
                {
                    for(dataTypeNChar i=0;i<vOutput.size();i++)
                        fwrite(&vOutput[i],sizeof(ElementCluster),1,OutCluster);
                }
                
				//Read the LCP and DA files
				bufferEle[0]=bufferEle[numcharDA];
				chunk-=numcharLCP;
                vOutput.clear();
			}//end-while
			
			//We need to close a possibly open cluster 
			if(init && (tid==num_threads-1) && nR && nG)//Last thread
				nClusters+=Close(cluster, index, maxLen,vOutput);
			else if (init && (t_id!=num_threads-1)) //Manage straddling clusters 
			{
				bool stayIn=true;
				do {
					numcharLCP=fread(&bufferLCP[0],sizeof(dataTypelenSeq),1,InLCP[t_id]);
					numcharDA=fread(&bufferEle[1],sizeof(dataTypeNSeq),1,InDA[t_id]);
					
					if(bufferLCP[0]>=alpha) //Remain in THE cluster
						StartOrRemain(numReads, bufferEle, 0, index, init, cluster, nR, nG);
					else//End THE cluster
					{
						stayIn=false;
						if (nR && nG)
							nClusters+=Close(cluster, index, maxLen,vOutput);
					}//end-else
					index++;
				} while(stayIn);
				
			}//end-else-if
            #if OMP
				#pragma omp critical
			#endif
            {
                for(dataTypeNChar i=0;i<vOutput.size();i++)
                    fwrite(&vOutput[i],sizeof(ElementCluster),1,OutCluster);
            }
		}//end-if
		#if OMP
			#pragma omp critical
			{
				std::cerr << "TIME THREAD " << t_id << " = " << omp_get_wtime()-start << "(in seconds)\n";
			}
		#endif
		
		delete[] bufferLCP;
		delete[] bufferEle;
		
	}//end-pragma
	tid=0;
    #if OMP
    for(;tid<num_threads; tid++)
    #endif
    {
        fclose(InLCP[tid]);
        fclose(InDA[tid]);
    }
	fclose(OutCluster);
	
	string fileaux=fileFasta.substr(0,fileFasta.find(".fasta"))+".out";
	
   //Write auxiliary file
	FILE * outAux = fopen(fileaux.c_str(), "w");
	if (outAux==NULL) {
    std::cerr << "Error opening " << fileaux << "." << std::endl;
    printf("fopen failed, errno = %d\n", errno);
    exit (EXIT_FAILURE);
	}
  
	fwrite(&numReads,sizeof(dataTypeNSeq),1,outAux);
	fwrite(&numGenomes,sizeof(dataTypeNSeq),1,outAux);
	fwrite(&alpha,sizeof(dataTypelenSeq),1,outAux);
	fwrite(&maxLen,sizeof(dataTypeNChar),1,outAux);
	fwrite(&nClusters,sizeof(dataTypeNChar),1,outAux);
                                  
	fclose(outAux);
                                  
	cout << "Clustering process with alpha=" << alpha << " completed.\nTotal number of clusters: " << nClusters << ".\nMaximum cluster size: " << maxLen << "." << endl;
	#if OMP
      fprintf(stdout,"Time: %.6lf\n", omp_get_wtime()-d_total);
    #else
      fprintf(stdout,"Time: %.6lf\n", time_stop(t_total, c_total));
    #endif

return 0;
}
