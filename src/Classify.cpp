#include "Tools.h"

/*
Step 3. performs the read assignment. 
Every read is either assigned to a particular taxon, or is reported as not classified.
The taxonomic rank to which the read is possibly assigned can be fixed, and in this case, the reads that can be assigned to more than one taxon are reported as ambiguous.

Input: The number of .txt\.bin files, list of names of .txt\.bin files obtained from ClusterBWT, total number of reads, total number of genomes, output name file, taxonomy file, rank, threads

The database taxonomy information must be provided in a ';'-separated file mapping reference sequence IDs to their lineage, whose first line is:
 Seq_ID;Kingdom_TaxID;Species_TaxID
 Ex.;Genus_TaxID:Family_TaxID;Order_TaxID;Class_TaxID;Phylum_TaxID
Ex.
NC_011750.1;562;561;543;91347;1236;1224
	
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
void FixRank(int ind,std::vector<dataTypeSet> &v_corRef, std::vector<string> &names, dataTypeNSeq **m_corRef, string fileTaxID){
#else
void FixRank(int ind,std::vector<dataTypeSet> &v_corRef, std::vector<string> &names, string fileTaxID){
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
	for(int s=0;s<RANK;s++)
		getline(f_taxid,lin[s],';');
    getline(f_taxid,lin[6],'\n');
	cout << "Taxonomic rank: " << lin[ind] << endl;
  
	for(int s=0;s<RANK;s++)
		getline(f_taxid,lin[s],';');
    getline(f_taxid,lin[6],'\n');
	
	if(ind>0){
		while (f_taxid.good()) {
			if(lin[ind]!="")
				v_corRef.push_back(atoi(lin[ind].c_str()));
			#if HIGHER
				for(int i=ind-1; i<RANK; i++){
					if(lin[i+1]!="")
						m_corRef[i][n]=atoi(lin[i+1].c_str());
				}
				n++;
			#endif
			for(int s=0;s<RANK;s++)
				getline(f_taxid,lin[s],';');
			getline(f_taxid,lin[6],'\n');
		}
	}
	else{//ind==0
		while (f_taxid.good()) {
			v_corRef.push_back(n);
			names.push_back(lin[0]);
			#if HIGHER
				for(int i=ind-1; i<RANK; i++){
					if(lin[i+1]!="")
						m_corRef[i][n]=atoi(lin[i+1].c_str());
				}
			#endif
			for(int s=0;s<RANK;s++)
				getline(f_taxid,lin[s],';');
			getline(f_taxid,lin[6],'\n');
			n++;
		}
	}
}

bool compareFirst (std::pair<float,uint> a,std::pair<float,uint> b) { return isgreater(b.first,a.first); }

#if BIN
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
    void ReadSimilarity(std::vector<std::pair<float,dataTypeNSeq>> &vectorLine, std::vector<dataTypeNSeq> & setTarg, std::vector<dataTypeSet> &v_corRef ){
        
        dataTypeNSeq ReadIdRef=0;
        float Sim = 0.0, maxSim = vectorLine[0].first;
        dataTypeNSeq j, size = vectorLine[0].second;
        
        for(j=1; j<=size; j++)
        {
            Sim = vectorLine[j].first;
            ReadIdRef = vectorLine[j].second;
            if (maxSim - Sim < (static_cast<float>(ERROR)) ){
                std::vector<dataTypeNSeq>::iterator it = find (setTarg.begin(), setTarg.end(), ReadIdRef);
                if (it == setTarg.end())
                    setTarg.push_back(ReadIdRef);//store only genome identifiers for which the similarity value is close to the maximum value maxSim
            }
        }
    }
#else
    void ReadMax(const string & s, float *maxSim){
      istringstream is( s );
      is >> (*maxSim);
    }

    void ReadSimilarity(const string & s, std::vector<dataTypeNSeq> & setTarg, std::vector<dataTypeSet> &v_corRef){

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
                    std::vector<dataTypeNSeq>::iterator it = find (setTarg.begin(), setTarg.end(), ReadIdRef);
                    if (it == setTarg.end())
                        setTarg.push_back(ReadIdRef);//store only genome identifiers for which the similarity value is close to the maximum value maxSim
                }
            }
        }
    }
#endif
    
#if OMP 
	void Assign(std::vector<typeOutput> &vOutput, set<dataTypeSet> &setTax, dataTypeNSeq idSeqRead, float h){//Assign the read to a UNIQUE TAXID --> CLASSIFIED
		std::set<dataTypeSet>::iterator it = setTax.begin();		
		vOutput[idSeqRead].type='C';
		vOutput[idSeqRead].TaxID=*it;
		vOutput[idSeqRead].sim=h;
	}
#else
	void Assign(std::ofstream &out, std::set<dataTypeSet> &setTax, dataTypeNSeq idSeqRead, float h){//Assign the read to a UNIQUE TAXID --> CLASSIFIED
		std::set<dataTypeSet>::iterator it = setTax.begin();
		out << "C," << idSeqRead << "," << *it << "," << h << "\n";
	}
#endif

#if BIN
float Exam_2(std::vector<std::pair<float,dataTypeNSeq>> *vectorLine, uint &numFile, dataTypeNSeq &numTarg, std::set<dataTypeNSeq> *setGen){
#else
float Exam_2(std::vector<string> &lineInput, uint &numFile, dataTypeNSeq &numTarg, std::set<dataTypeNSeq> *setGen){
#endif
	
	float h=0.0;
	vector<float> highest (2,0.0);

	float *SimArray_[MAXFILES];//matrix numTarg x numFile

	for(uint i=0; i<numFile; ++i){
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
	
	//Find max({1F+2RC}) and max({1RC+2RC})
	for(uint k=0; k<2; k++){
		for(j=0; j<numTarg; ++j){
			if (k==0 && numFile==4){
				SimArray_[0][j]+=SimArray_[3][j];
				SimArray_[1][j]+=SimArray_[2][j];
			}	
			if (highest[k]<SimArray_[k][j])
				highest[k]=SimArray_[k][j];
		}
	}
	
	if(highest[0] > highest[1])
	{
		h=highest[0];
		for(j=0; j<numTarg; ++j){
			if(h-SimArray_[0][j] < (static_cast<float>(ERROR)) )
			(*setGen).insert(j);
		}
	}
	else if(highest[0] < highest[1])
	{
		h=highest[1];
		for(j=0; j<numTarg; ++j){
			if(h-SimArray_[1][j] < (static_cast<float>(ERROR)) )
			(*setGen).insert(j);
		}
	}
	else //highest[0] == highest[1]
	{
		h=highest[0];
		for(j=0; j<numTarg; ++j){
			if( (h-SimArray_[0][j] < (static_cast<float>(ERROR))) || (h-SimArray_[1][j]< (static_cast<float>(ERROR))) )
				(*setGen).insert(j);
		}
	}
	assert((*setGen).size()>0);
	
	for(uint i=0; i<numFile; ++i)
		free(SimArray_[i]);

	return h;
}
#if OMP 
int AssignToHigher(int index, dataTypeNSeq idSeqRead, std::set<dataTypeNSeq> &setTarg, dataTypeNSeq **m_corRef, std::vector<typeOutput> &vOutput, dataTypeNSeq &numH, float maxSum){
#else
int AssignToHigher(int index, dataTypeNSeq idSeqRead, std::set<dataTypeNSeq> &setTarg, dataTypeNSeq **m_corRef, ofstream &out, dataTypeNSeq &numH, float maxSum){
#endif
	//Assign to higher taxonomic ranks
	bool flagAssign, search=true;
	int amb=1;
	
	while(search && index<RANK){
		flagAssign=true;
		std::set<dataTypeNSeq>::iterator it=setTarg.begin();
		dataTypeNSeq taxon=m_corRef[index][*it];
		for (std::set<dataTypeNSeq>::iterator it=setTarg.begin(); it!=setTarg.end() && flagAssign; it++)
		{
			if(taxon!=m_corRef[index][*it])
				flagAssign=false;
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
			amb=0;
			search=false;
		}
		else
			index++;
	}
	if(search){
		#if OMP
			vOutput[idSeqRead].type='A';
		#else
			out << "A," << idSeqRead << "," << "NA,0\n";
		#endif
	}	
	return amb;
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
	if(taxRank > RANK || taxRank < 0)
		std::cerr << "Error usage: taxRank 0=Genome, 1=Species, 2=Genus, 3=Family, 4=Order, 5=Class, 6=Phylum." << std::endl;
	
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
	std::vector<string> names;
	
	#if HIGHER
		dataTypeNSeq **m_corRef = (dataTypeNSeq**) malloc(RANK * sizeof(dataTypeNSeq*)); 
		for(dataTypeNSeq i=0; i<RANK; ++i){
		  m_corRef[i]  = (dataTypeNSeq*) malloc(numTarg*sizeof(dataTypeNSeq));
		  for(dataTypeNSeq j=0; j<numTarg; ++j)
			m_corRef[i][j] = 0;
		}
        FixRank(taxRank, v_corRef, names, m_corRef, fileTaxID);
	#else
        FixRank(taxRank, v_corRef, names, fileTaxID);
	#endif
	
	if(numTarg!=v_corRef.size()){
        cerr << "Number of taxIDs = " << v_corRef.size() << " lower than genome number: poor taxonomy information to classify." << endl;
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
		#pragma omp parallel for default(shared) firstprivate(vectorLine) reduction(+:numC) reduction(+:numH) reduction(+:numA) reduction(+:numNC)
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
			std::vector<dataTypeNSeq> setI;
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
					ReadSimilarity(vectorLine[highest_pair.second], setI, v_corRef);
				#else
					ReadSimilarity(lineInput[highest_pair.second], setI, v_corRef);
				#endif
			}
			assert(setI.size()>0);
			std::set<dataTypeSet> setTax;
			for(uint e = 0; e< setI.size(); e++)
				setTax.insert(v_corRef[setI[e]]);
			if (setTax.size()==1){
				#if OMP 
					Assign(vOutput,setTax,idSeqRead,highest);
				#else
					Assign(out, setTax,idSeqRead,highest);
				#endif
				numC++;
				//class_1++;
			}
			else{//Ambiguous exam
				setTax.clear();
				std::vector<float> CompareGMax[MAXFILES]; 
				for(uint e = 0; e< setI.size(); e++) {
					for(uint i=0;i<numFile;i++){
						//Search genome setI[e]
						#if BIN
							dataTypeNSeq size = vectorLine[i][0].second;
							bool fnd=false;
							if(size>0){
								for(dataTypeNSeq j=1; (j<=size && !fnd) ;j++){
									if(vectorLine[i][j].second==setI[e]){
										fnd=true;
										CompareGMax[i].push_back(vectorLine[i][j].first);
									}
								}
							}
							if(!fnd)
								CompareGMax[i].push_back(0.0);
						#else
							float sim=0.0;
							dataTypeNSeq ref;
							bool val,fnd=false;
							istringstream isLine ( lineInput[i] );
							val= static_cast<bool> (isLine >> sim);
							while( val == 1 && !fnd) {
								val = static_cast<bool> (isLine >> ref); //genome identifier index
								if ( val == 1 ) {
									val = static_cast<bool> (isLine >> sim);			
									if (ref==setI[e]){
										fnd=true;
										CompareGMax[i].push_back(sim);
									}	
								}
							}
							if(!fnd)
								CompareGMax[i].push_back(0.0);
						#endif
					}
				}
				vector<float> maxsum (2,0.0);
				for (uint k=0;k<2; k++){
					for (uint j=0;j<CompareGMax[0].size();j++){
						if(k==0 && numFile==4){//Sum mates and find max values
							CompareGMax[0][j]+=CompareGMax[3][j];
							CompareGMax[1][j]+=CompareGMax[2][j];
						}
						if(maxsum[k] < CompareGMax[k][j])
							maxsum[k]= CompareGMax[k][j];
					}
				}
				bool Exam=true;
				if (maxsum[0] > maxsum[1] + static_cast<float>(ERROR)){
					for (dataTypeNSeq j=0;j<CompareGMax[0].size();j++){
						if(maxsum[0]==CompareGMax[0][j]){
							setTax.insert(v_corRef[setI[j]]);
						}
					}
					if (setTax.size()==1){//Assign 
						#if OMP 
							Assign(vOutput,setTax,idSeqRead,maxsum[0]);
						#else
							Assign(out, setTax,idSeqRead,maxsum[0]);
						#endif
						numC++;
						Exam=false;
					}
				}
				else if(maxsum[1] > maxsum[0] + static_cast<float>(ERROR)){
					for (uint j=0;j<CompareGMax[1].size();j++){
						if(maxsum[1]==CompareGMax[1][j]){
							setTax.insert(v_corRef[setI[j]]);
						}
					}
					if (setTax.size()==1){//Assign 
						#if OMP 
							Assign(vOutput,setTax,idSeqRead,maxsum[1]);
						#else
							Assign(out, setTax,idSeqRead,maxsum[1]);
						#endif
						numC++;
						Exam=false;
					}
				}

				if(Exam)
				{
					setTax.clear();
					set<dataTypeNSeq> setGen;
					#if BIN
						highest=Exam_2(vectorLine, numFile, numTarg, &setGen);
					#else
						highest=Exam_2(lineInput,numFile,numTarg, &setGen);
					#endif
					
					bool is_class=true;
					std::set<dataTypeNSeq>::iterator it = setGen.begin();
					dataTypeSet taxon=v_corRef[*it];
                    setTax.insert(taxon);
					for(;it != setGen.end(); it++)
					{
						if(taxon!=v_corRef[*it])
							is_class=false;
					}
					if(is_class){
						#if OMP 
							Assign(vOutput,setTax,idSeqRead,highest);
						#else
							Assign(out, setTax,idSeqRead,highest);
						#endif
						numC++;
						//class_3++;
						setGen.clear();
					}
					else{	
						#if HIGHER
						    #if OMP
							numA+=AssignToHigher(taxRank-1, idSeqRead, setGen, m_corRef, vOutput, numH, highest);
						    #else
							numA+=AssignToHigher(taxRank-1, idSeqRead, setGen, m_corRef, out, numH, highest);
						    #endif
						#else
							numA++;
							#if OMP
								vOutput[idSeqRead].type='A';
							#else
								out << "A," << idSeqRead << "," << "NA,0\n";
							#endif
						#endif
					}					
				}//end-if Exam
				
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
