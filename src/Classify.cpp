#include "Tools.h"

/*
Step 3. performs the read assignment. 
Every read is either assigned to a particular taxon, or is reported as not classified.
The taxonomic rank to which the read has to assigned can be fixed, and in this case, the reads that can be assigned to more than one taxon are reported as ambiguous. 

Input: The number of .txt\.bin files, list of names of .txt\.bin files obtained from ClusterBWT, total number of reads, total number of genomes, output name file, taxonomy file, rank, threads

The database taxonomy information must be provided in a ';'-separated file mapping reference sequence IDs to their lineage, whose first line is:
Seq_ID;Kingdom_TaxID;Phylum_TaxID;Class_TaxID;Order_TaxID;Family_TaxID;Genus_TaxID;Species_TaxID
Ex.
NC_011750.1;2;1224;1236;91347;543;561;562	
	
Output: A .txt file containing classification assignments for each read. The assignment output has 4 columns.
  The first column specifies if the read is not classified (U), classified ambiguously (A), assigned to a taxon at the given rank (C) or at a higher rank (H).
  The second column is the read index.
  The third column is either NA (in case of U or A) or the taxID to which the read is assigned.
  The fourth column is the similarity score for the classification.
*/

#define MAXFILES 4 //number of .txt/.bin files
#define RANK 6

#ifndef HIGHER
	#define HIGHER 0
#endif

//To fix taxonomic rank
#if HIGHER
void FixRank(int ind,std::vector<dataTypeSet> &v_corRef, dataTypeNSeq **m_corRef, string fileTaxID){
#else
void FixRank(int ind,std::vector<dataTypeSet> &v_corRef, string fileTaxID){
#endif
	std::ifstream f_taxid;
	f_taxid.open(fileTaxID.c_str(), std::ifstream::in);
	if (!f_taxid.is_open()){
		std::cerr << "Error opening file " << fileTaxID << ".\n";
		exit(1);
	}
	int n=0;
    
	vector<string> lin;
	lin.resize(7);
	//read header
	getline(f_taxid,lin[0],';');
	getline(f_taxid,lin[1],';');
    getline(f_taxid,lin[2],';');
	getline(f_taxid,lin[3],';');
    getline(f_taxid,lin[4],';');
	getline(f_taxid,lin[5],';');
    getline(f_taxid,lin[6],'\n');
	cout << "Taxonomic rank: " << lin[ind] << endl;
  
	getline(f_taxid,lin[0],';');
	getline(f_taxid,lin[1],';');
    getline(f_taxid,lin[2],';');
	getline(f_taxid,lin[3],';');
    getline(f_taxid,lin[4],';');
	getline(f_taxid,lin[5],';');
    getline(f_taxid,lin[6],'\n');
	while (f_taxid.good()) {
    #if TaxLevel == 1
		if(lin[ind]!="")
			v_corRef.push_back(atoi(lin[ind].c_str()));
		#if HIGHER
		for(int i=ind-1; i>=0; i--){
			if(lin[i+1]!="")
				m_corRef[i][n]=atoi(lin[i+1].c_str());
		}
		#endif
    #else
		v_corRef.push_back(lin[0]);
    #endif
		getline(f_taxid,lin[0],';');
		getline(f_taxid,lin[1],';');
		getline(f_taxid,lin[2],';');
		getline(f_taxid,lin[3],';');
		getline(f_taxid,lin[4],';');
		getline(f_taxid,lin[5],';');
		getline(f_taxid,lin[6],'\n');
		n++;
	}
	
	return;
}

bool compareFirst (std::pair<float,uint> a,std::pair<float,uint> b) { return isgreater(b.first,a.first); }

/*****BIN*****/
void ReadMax(FILE *fIn, size_t pos, float *maxSim, std::vector<std::pair<float,dataTypeNSeq>>& vectorLine){
	
	dataTypeNSeq ReadIdRef=0, size=0;
	float Sim = 0.0;
	
	fseek(fIn, pos*sizeof(pair_sim), SEEK_SET);
	fread(maxSim, sizeof(float), 1, fIn);
	fread(&size, sizeof(dataTypeNSeq), 1, fIn);

	vectorLine.push_back(std::make_pair(*maxSim,size));
	for(uint i=0; i<size; i++){
      fread(&Sim, sizeof(float), 1, fIn);
      fread(&ReadIdRef, sizeof(uint), 1, fIn);
      vectorLine.push_back(std::make_pair(Sim,ReadIdRef));
    }
}
/***** NOT BIN*****/
void ReadMax(const string & s, float *maxSim){
  istringstream is( s );
  is >> (*maxSim);
}

/*****BIN*****/
void ReadSimilarity(std::vector<std::pair<float,dataTypeNSeq>> &vectorLine, std::set<dataTypeSet> & setTarg, std::vector<dataTypeSet> &v_corRef ){

	dataTypeNSeq ReadIdRef=0;
	float Sim = 0.0, maxSim = vectorLine[0].first;
	dataTypeNSeq j, size = vectorLine[0].second;

	for(j=1; j<=size; j++)
	{      
		Sim = vectorLine[j].first;
		ReadIdRef = vectorLine[j].second;
		if (maxSim - Sim < (static_cast<float>(ERROR)) ){
			setTarg.insert(v_corRef[ReadIdRef]);//store only genome identifiers for which the similarity value is close to the maximum value maxSim
		}
	}
}

/***** NOT BIN*****/
void ReadSimilarity(const string & s, std::set<dataTypeSet> & setTarg, std::vector<dataTypeSet> &v_corRef){

	float maxSim;
	dataTypeNSeq ReadIdRef;
	bool val;
	istringstream is( s );
	val= static_cast<bool> (is >> maxSim);

	while( val == 1 ) {
		val = static_cast<bool> (is >> ReadIdRef); //genome identifier index
		if ( val == 1 ) {
			float Sim;
			val = static_cast<bool> (is >> Sim);			
			if (maxSim - Sim < (static_cast<float>(ERROR)) ){
				setTarg.insert(v_corRef[ReadIdRef]);//store only genome identifiers for which the similarity value is close to the maximum value maxSim
			}	
		}
	}
}

#if BIN
float reExamination_2(std::vector<std::pair<float,dataTypeNSeq>> *vectorLine, uint &numFile, dataTypeNSeq &numTarg, std::set<dataTypeNSeq> *setGen){
#else
float reExamination_2(std::vector<string> &lineInput, uint &numFile, dataTypeNSeq &numTarg, std::set<dataTypeNSeq> *setGen){
#endif
	
	float highest=0.0;

  	float *SimArray_[MAXFILES+1];//matrix numTarg x numFile+1

	for(uint i=0; i<numFile+1; ++i){
		SimArray_[i]  = (float*) malloc(numTarg * sizeof(float));
		for(dataTypeNSeq j=0; j<numTarg; ++j)
			SimArray_[i][j] = 0.0;
	}

	float sim_value;
    dataTypeNSeq j;
	dataTypeNSeq IdRef;
	#if BIN
		for(uint i=0; i<numFile; ++i){
			dataTypeNSeq size = vectorLine[i][0].second;
			for(j=1; j<=size; j++){	
				sim_value = vectorLine[i][j].first;
				IdRef = vectorLine[i][j].second;
				assert(IdRef<numTarg);
				SimArray_[i][IdRef]+=sim_value;//store all the sim_values of genomes
			}
		}
	#else
		float max[MAXFILES];
		bool val;
		for(uint i=0; i<numFile; ++i){//Keep in memory all the sim values
			max[i]=0;
			istringstream is(lineInput[i]);
			val= static_cast<bool> (is >> max[i]);
			while( val == 1 ) {
				val =static_cast<bool> (is >> IdRef); //genome identifier index
				if ( val == 1 ) {
					val =static_cast<bool>(is >> sim_value);
					assert(IdRef<numTarg);
					SimArray_[i][IdRef]+=sim_value;//store all the sim_values of genomes
				}
			}
		}
	#endif
	
	//Take max(F,RC)
	int row=0;
	do {
		for(j=0; j<numTarg; ++j){
			if(SimArray_[row][j]<SimArray_[row+1][j])
				SimArray_[row][j]=SimArray_[row+1][j];
			if (numFile==2) 
				SimArray_[numFile][j]=SimArray_[row][j];
		}
		row+=2;
		if (numFile==2)
			row=3;
	} while(row<3);
	
	//Find maximum
	for(j=0; j<numTarg; ++j){
		if (numFile==4)
			SimArray_[numFile][j]=SimArray_[0][j]+SimArray_[2][j];
		
		if (highest<SimArray_[numFile][j])
			highest=SimArray_[numFile][j];
	}
	//Classify
	for(j=0; j<numTarg; ++j){
		if(highest-SimArray_[numFile][j]< (static_cast<float>(ERROR)) )
		  (*setGen).insert(j);
	}
		
	assert((*setGen).size()>0);
	
	for(uint i=0; i<numFile+1; ++i)
    free(SimArray_[i]);

	return highest;
}
#if BIN
float reExamination_1(std::set<dataTypeNSeq> &setGen, std::vector<std::pair<float,dataTypeNSeq>> *vectorLine, uint &numFile, dataTypeNSeq idSeqRead,std::vector<dataTypeSet> &v_corRef, ofstream &out, std::vector<typeOutput> &vOutput, dataTypeNSeq &numC, dataTypeNSeq numTarg){
#else
float reExamination_1(std::set<dataTypeNSeq> &setGen, std::vector<string> &lineInput, uint &numFile, dataTypeNSeq idSeqRead,std::vector<dataTypeSet> &v_corRef, ofstream &out, std::vector<typeOutput> &vOutput, dataTypeNSeq &numC, dataTypeNSeq numTarg){
#endif

	float maxSum=0.0;
	std::vector<occ> CompareMates;
	occ el;
	el.TaxID=0, el.maxRead=0, el.maxMate=0;
  
	float maxSim, Sim;
	dataTypeNSeq ReadIdRef;
	dataTypeNSeq index=0;	
	for(uint i=0;i<numFile;i++ )
    { 
		#if BIN
			maxSim = vectorLine[i][0].first;
			dataTypeNSeq j, size = vectorLine[i][0].second;
		#else
			istringstream is(lineInput[i]);
			is >> maxSim; //Read the maximum similarity value
		#endif

		#if BIN
		for(j=1; j<=size; j++){
			Sim = vectorLine[i][j].first;
			ReadIdRef = vectorLine[i][j].second;
		#else
		while ((is >> ReadIdRef) && (is >> Sim)){//Read line	
		#endif
			if (maxSim - Sim < (static_cast<float>(ERROR)) )
			{//If the similarity value is close to the maximum --> insert an entry in CompareMates
				el.TaxID=v_corRef[ReadIdRef];
				if(i<numFile/2 || numFile==2)
				{
					el.maxRead=Sim;
					el.maxMate=0;
				}
				else
				{
					el.maxRead=0;
					el.maxMate=Sim;
				}
				index=0;//index for running on the list of TaxID
				while ((index<CompareMates.size()) && (el.TaxID <CompareMates[index].TaxID))
					index++;

				//Possibly update the similarity value 
                if ( (index < CompareMates.size()) && (el.TaxID == CompareMates[index].TaxID) )
                {
                    if (i<numFile/2 || numFile==2){
						if (el.maxRead > CompareMates[index].maxRead)
							CompareMates[index].maxRead=el.maxRead;
                    }
                    else if (el.maxMate > CompareMates[index].maxMate){
                        CompareMates[index].maxMate=el.maxMate;
                    }
                }
                else{//TaxID not in the list --> insert it
                    if ( index < CompareMates.size()) //index is not pointing the end of CompareMates
                        CompareMates.insert(CompareMates.begin()+index,el);
                    else    //append the new entry
                        CompareMates.push_back(el);
                }
	        }
	    }//end-for
	}//end-for
	
	maxSum=0.0;
	dataTypeNSeq TaxIDMaxSum=0;
	float secondMaxSum=0;
  
	for(uint i=0; i<CompareMates.size(); i++){
		if ((CompareMates[i].maxRead>0 && (CompareMates[i].maxMate>0)) || (numFile==2)){//Take species appearing in both mates (for paired-end reads)
			if (maxSum <= CompareMates[i].maxRead + CompareMates[i].maxMate)
			{
				secondMaxSum = maxSum;
				maxSum = CompareMates[i].maxRead + CompareMates[i].maxMate;
				TaxIDMaxSum = CompareMates[i].TaxID;
			}
			else
			{
				if (secondMaxSum < CompareMates[i].maxRead + CompareMates[i].maxMate)
					secondMaxSum = CompareMates[i].maxRead + CompareMates[i].maxMate;
			}
		}//end-if
	}//end-for
  
	if ( maxSum > secondMaxSum + 0.001){
		#if OMP
			vOutput[idSeqRead].type='C';
			vOutput[idSeqRead].TaxID=TaxIDMaxSum+1;
			vOutput[idSeqRead].sim=maxSum;
		#else
			out << "C," << idSeqRead << "," << TaxIDMaxSum << "," << maxSum << "\n";
		#endif
		numC++;
	} 
	else {
		#if BIN
			maxSum=reExamination_2(vectorLine, numFile, numTarg, &setGen);
		#else
			maxSum=reExamination_2(lineInput,numFile,numTarg, &setGen);
		#endif
		std::set<dataTypeNSeq>::iterator it = setGen.begin();
		dataTypeSet taxon=v_corRef[*it];
		for(;it != setGen.end(); it++)
		{
			if(taxon!=v_corRef[*it])
				return maxSum;
		}
		
		#if OMP
			vOutput[idSeqRead].type='C';
			vOutput[idSeqRead].TaxID=taxon;
			vOutput[idSeqRead].sim=maxSum;
		#else
			out << "C," << idSeqRead << "," << taxon << "," << maxSum << "\n";
		#endif
		numC++;
		setGen.clear();
	}
	CompareMates.clear();
	return maxSum;
}

int AssignToHigher(int index, dataTypeNSeq idSeqRead, std::set<dataTypeNSeq> &setTarg, dataTypeNSeq **m_corRef, ofstream &out, std::vector<typeOutput> &vOutput, dataTypeNSeq &numH, float maxSum){
	//Assign to higher taxonomic ranks
	bool flagAssign;
	
	while(index>=0){
		flagAssign=true;
		std::set<dataTypeNSeq>::iterator it=setTarg.begin();
		dataTypeNSeq taxon=m_corRef[index][*it];
		for (std::set<dataTypeNSeq>::iterator it=setTarg.begin(); it!=setTarg.end(); it++)
		{
			if(taxon!=m_corRef[index][*it]){
				flagAssign=false;
				break;
			}
		}
		if(flagAssign && taxon!=0){//Assign
			#if OMP 
				vOutput[idSeqRead].type='H';
				vOutput[idSeqRead].TaxID=taxon;
				vOutput[idSeqRead].sim=maxSum;
			#else
				out << "H," << idSeqRead << "," << taxon << "," << maxSum << "\n";
			#endif
				numH++;
			return 0;
		}
		else
			index--;
	}
	#if OMP
		vOutput[idSeqRead].type='A';
	#else
		out << "A," << idSeqRead << "," << "NA,0\n";
	#endif
	
	return 1;
}

int main(int argc, char **argv) {

	#if OMP
		double d_total;
		d_total = omp_get_wtime();
	#else
		time_t t_total=0;
		clock_t c_total=0;
		time_start(&t_total, &c_total);
	#endif

	#if DEBUG
		cout<<"##\n";
		cout<<"DEBUG == "<<DEBUG<<endl;
		cout<<"##\n";
	#endif
  
	if (argc < 2){
		std::cerr << "Error usage " << argv[0] << " N fileInput1 fileInput2 ... fileInputN numReads numGenomes fileOutput fileTaxo taxRank numThreads" << std::endl;
		exit(1);
	}
  
	uint numFile;
	sscanf(argv[1], "%u", &numFile);
  
	if( argc != (int)numFile+8){
		std::cerr << "Error usage " << argv[0] << " N fileInput1 fileInput2 ... fileInputN numReads numGenomes fileOutput fileTaxo taxRank numThreads" << std::endl;
		exit(1);
	}
  
	if (numFile!=2 && numFile!=4 ){ 
		//MAXFILES==4
		std::cerr << "Error usage " << argv[0] << ": the allowed number of input files is 2 (single-end reads), or 4 (paired-end reads)" << std::endl;
		exit(1);
	}
  
	vector<string> fileInput;
	for(uint i=2;i<numFile+2;i++){
		fileInput.push_back(argv[i]);
	}
  
	dataTypeNSeq numReads; 
	sscanf(argv[numFile+2], "%u", &numReads);

	dataTypeNSeq numTarg; 
	sscanf(argv[numFile+3], "%u", &numTarg);
  
	string fileOutput;
	fileOutput = argv[numFile+4];
   
	string fileTaxID= argv[numFile+5];
    
    int taxRank;
	sscanf(argv[numFile+6], "%d", &taxRank);
	if(taxRank > RANK || taxRank < 1)
		std::cerr << "Error usage: taxRank 1=Phylum, 2=Class, 3=Order, 4=Family, 5=Genus, 6=Species." << std::endl;
	
	int numThreads;
	sscanf(argv[numFile+7], "%d", &numThreads);

	#if OMP
		omp_set_num_threads(numThreads);

		int threads=0;
		#pragma omp parallel
		{
			threads = omp_get_num_threads();
			if(omp_get_thread_num()==0)
				printf("Number of threads: %d\n", threads);
		}
		printf("Number of processors: %d\n", omp_get_num_procs());
	  #endif
	  
	//Make correspond genomes to TaxID at a given taxonomic rank
	cout << "Reading " << fileTaxID << endl;
  
	std::vector<dataTypeNSeq> v_corRef;
	
	#if HIGHER
		dataTypeNSeq **m_corRef = (dataTypeNSeq**) malloc(RANK * sizeof(dataTypeNSeq*)); 
		for(dataTypeNSeq i=0; i<RANK; ++i){
		  m_corRef[i]  = (dataTypeNSeq*) malloc(numTarg*sizeof(dataTypeNSeq));
		  for(dataTypeNSeq j=0; j<numTarg; ++j)
			m_corRef[i][j] = 0;
		}
        FixRank(taxRank, v_corRef, m_corRef, fileTaxID);
	#else
        FixRank(taxRank, v_corRef, fileTaxID);
	#endif
	
	if(numTarg!=v_corRef.size()){
		cerr << "Number of taxIDs lower than genome number: poor taxonomy information to classify." << endl;
		exit(1);
	}
  
	//Open .txt/.bin files
	std::cout << "Reading files:";

	#if BIN
		FILE* fdResult_Bin[numThreads][MAXFILES];
		FILE* fdResult_Pos[numThreads][MAXFILES];
	#else
		std::vector <ifstream*> fdResult_;
		std::vector<string> lineInput;
	#endif
  
	int tid=0;
	#if BIN
		#if OMP
		for(;tid<numThreads; tid++)
		#endif
		{
			for(uint i=0;i<numFile;i++ ){
				string sFileBin(fileInput[i]);
				string sFilePos(fileInput[i]);
				sFileBin.append(".bin");
				sFilePos.append(".pos");
				fdResult_Bin[tid][i] = fopen(sFileBin.c_str(), "rb");
				fdResult_Pos[tid][i] = fopen(sFilePos.c_str(), "rb");
				if(!tid){
					std::cout << "\n\t" << sFileBin;
					std::cout << "\n\t" << sFilePos;
				}
			}	
		}
	#else
		for(uint i=0;i<numFile;i++ ){
			string sFileTxt(fileInput[i]);
			sFileTxt.append(".txt");
			lineInput.push_back("");
			std::cout << "\n\t" << sFileTxt;
			ifstream* f = new ifstream;
			fdResult_.push_back(f);
			fdResult_[i]->open(sFileTxt.c_str());
			if (!fdResult_[i]->is_open()) {
				cerr << "Error opening " << sFileTxt << endl;
				exit(1);
			}
		}
	#endif
	std::cout << std::endl;
    
	//Output
	dataTypeNSeq numC=0, numNC=0, numA=0, numH=0;
	dataTypeNSeq idSeqRead=0;
	std::vector< std::pair<float,dataTypeNSeq> > vectorLine[MAXFILES];  
  
	std::ofstream out;
	out.open(fileOutput.c_str(), std::ofstream::out);
	if (!out.is_open()) cerr << "ERROR: File Output not Open" << endl;
	out << "C/U/A/H,IdSeqRead,TaxID,maxSim\n"; 
	
	#if OMP
		typeOutput aux;
		aux.type='U';
		aux.TaxID=0;
		aux.sim=0.0;
		vector<typeOutput> vOutput(numReads, aux); 
    #endif
	//START
	cerr << "Start comparing..." << endl;
  
	#if BIN == 0
		for(uint i=0;i<numFile;i++){
			getline(*fdResult_[i], lineInput[i]);
		}
	#endif


	#if OMP
		#pragma omp parallel for default(shared) firstprivate(vectorLine) reduction(+:numC) reduction(+:numH) reduction(+:numNC) reduction(+:numA)
	#endif
	for(idSeqRead=0; idSeqRead<numReads; idSeqRead++){
		int tid = 0;
		#if OMP
			tid = omp_get_thread_num();
		#endif

		size_t pos[MAXFILES];
		std::vector<float> maxSim_;//store the maximum similarity value

		for(uint i=0;i<numFile;i++){
			maxSim_.push_back(0);
			#if BIN
				fseek(fdResult_Pos[tid][i], (idSeqRead)*sizeof(size_t), SEEK_SET);
				fread(&pos[i], sizeof(size_t), 1, fdResult_Pos[tid][i]);
			#endif
		}
		
		std::vector< std::pair<float,uint> > vectorFileMax;  //vector storing in which file the max similarity value is
		bool notClass=true;
		
		for(uint i=0;i<numFile;i++){
			#if BIN
				if(pos[i])
					ReadMax(fdResult_Bin[tid][i],pos[i], &maxSim_[i],vectorLine[i]);
				else
					vectorLine[i].push_back(std::make_pair(0.0,0));
			#else
				ReadMax(lineInput[i], &maxSim_[i]);
			#endif
			if (maxSim_[i]) {
				vectorFileMax.push_back(std::make_pair(maxSim_[i],i));
				notClass=false;
			}
	
		}
      
		if (notClass){//Case read NOT CLASSIFIED
			#if OMP == 0
			  out << "U," << idSeqRead << ",NA,0\n";
			#endif
			numNC++;
		}
		else{
			std::set<dataTypeSet> setTax;
			uint sizeMax_multiset=0;
			
			std::sort (vectorFileMax.begin(), vectorFileMax.end(), compareFirst);	
			
			std::pair<float,uint> curr_pair;
			std::pair<float,uint> highest_pair =*vectorFileMax.rbegin();
			float highest = highest_pair.first;
			
			for (std::vector<std::pair<float,uint> >::iterator it=vectorFileMax.begin(); it<vectorFileMax.end(); it++){
				curr_pair=*it;
				if(highest-curr_pair.first < static_cast<float>(ERROR) )
					sizeMax_multiset++;
			}
			assert(sizeMax_multiset>0 && sizeMax_multiset<=numFile);

			for (uint i=0; i<sizeMax_multiset;i++){
				highest_pair=*(vectorFileMax.rbegin()+i);
				#if BIN
					ReadSimilarity(vectorLine[highest_pair.second], setTax, v_corRef);
				#else
					ReadSimilarity(lineInput[highest_pair.second], setTax, v_corRef);
				#endif
			}
			
			if (setTax.size()==1){//Assign the read to a UNIQUE TAXID --> CLASSIFIED
				std::set<dataTypeSet>::iterator it = setTax.begin();
				#if OMP 
					vOutput[idSeqRead].type='C';
					vOutput[idSeqRead].TaxID=*it;
					vOutput[idSeqRead].sim=highest;
				#else
					out << "C," << idSeqRead << "," << *it << "," << highest << "\n";
				#endif
				numC++;
			}
			else{//Ambiguous re-examination
				set<dataTypeNSeq> setGen;
				#if BIN
				  highest=reExamination_1(setGen, vectorLine, numFile, idSeqRead, v_corRef, out, vOutput, numC, numTarg);
				#else
				  highest=reExamination_1(setGen, lineInput, numFile, idSeqRead, v_corRef, out, vOutput, numC, numTarg);
				#endif
				if (setGen.size()>0)
				{	
					#if HIGHER
						numA+=AssignToHigher(taxRank-1, idSeqRead, setGen, m_corRef, out, vOutput, numH, highest); 
					#else
						numA++;
						#if OMP
							vOutput[idSeqRead].type='A';
						#else
							out << "A," << idSeqRead << "," << "NA,0\n";
						#endif
					#endif
				}
			}//end-else   
		}//end if-else

		vectorFileMax.shrink_to_fit();
		maxSim_.shrink_to_fit();

		for(uint i=0;i<numFile;i++){
			#if BIN
			  vectorLine[i].clear();
			#else
			  getline(*fdResult_[i], lineInput[i]);
			#endif
		}
	}//end-for (end reading files)

	tid=0;
	#if BIN
		#if OMP
		for(;tid<numThreads; tid++)
		#endif
		{
			for(uint i=0;i<numFile;i++ ){
				fclose(fdResult_Bin[tid][i]);
				fclose(fdResult_Pos[tid][i]);
			}
		}
	#else
		for(uint i=0;i<numFile;i++ ){
			assert(fdResult_[i]->eof());
			fdResult_[i]->close();
		}
	#endif

	#if OMP
    for(idSeqRead=0; idSeqRead<numReads; idSeqRead++){
		if(vOutput[idSeqRead].type=='U')
			out << "U," << idSeqRead << ",NA,0\n";
		else if (vOutput[idSeqRead].type=='A')
			out << "A," << idSeqRead << "," << "NA,0\n";
		else if (vOutput[idSeqRead].type=='C')
			out << "C," << idSeqRead << "," << vOutput[idSeqRead].TaxID << "," << vOutput[idSeqRead].sim << "\n";
		else
			out << "H," << idSeqRead << "," << vOutput[idSeqRead].TaxID << "," << vOutput[idSeqRead].sim << "\n";
    }
    vOutput.shrink_to_fit();
	#endif

	out.close();

	std::cout << "Classification process at level " << taxRank << " completed.\nNumber of successfully classified reads: " << numC << "/" << numReads<< ";" << std::endl;
	#if HIGHER
	std::cout << "\tClassified at higher taxonomic ranks: " << numH << "." << endl;
	#endif
	std::cout << "\tAmbiguously classified reads: " << numA << "." << std::endl;
	std::cout << "\tNot classified reads: " << numNC << "." << endl;
	#if OMP
		fprintf(stdout,"Time: %.6lf\n", omp_get_wtime()-d_total);
	#else
		fprintf(stdout,"Time: %.6lf\n", time_stop(t_total, c_total));
	#endif
	return 0;
}
