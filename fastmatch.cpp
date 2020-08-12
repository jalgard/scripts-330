#include "fastmatch.h"
#include "dnacommon.h"

void Preprocess_Collection(int M, int Q, std::unordered_map<int, std::string>& read_collection,
        std::unordered_map<std::string, std::vector<int> >& qgramm_lib)
{
    std::cout << "In Preprocess_Collection\n";
    for(auto& x : read_collection)
    {
        if(M * Q >= x.second.size()) continue;
        std::string Mr = "";
        for(size_t i = 0; i < x.second.size(); i += M)
            Mr += x.second[i];

        for(size_t i = 0; i < Mr.size()-Q; i++)
        {
            std::string qgramm = Mr.substr(i, Q);
            qgramm_lib[qgramm].push_back(x.first);	// each value of Q-gramm library
            qgramm_lib[qgramm].push_back(i * M);	// is a pair of int-s with indecies  i, i+1
							// where [i] is read ID and [i+1] is position
							// of Q-gramm in read's sequence
        }
    }
    std::cout << "Out Preprocess_Collection\n";
    // DEBUG mode lists Q-gramm lib
    // printing out Q-gramm and its parent ID - position pairs

    #ifdef DEBUG_FASTMATCH
    for(auto x : qgramm_lib)
    {
        std::cout << "Qgramm " << x.first;
	int o = 1;
        for(auto y : x.second) {
            std::cout << " " << y; if(++o%2) std::cout << ";";
        }
        std::cout << std::endl;
    }
    #endif

    // DEBUG COMPLETE
}

/*
 * MULTITHREADED IMPLEMENTATION OF UNGAPPED PATTERN PRELOCATION
 */

std::vector<int> Ungapped_Find_Pattern_For_One_Polyphase_Task(InnerPatternMatchTask task)
{
    std::vector<int> results;
    std::vector<std::string> QSP;
    for(int i = 0; i < task.PM.size(); i += task.Q)
    {
	std::string qgr = task.PM.substr(i, task.Q);
	if(qgr.size()==task.Q) QSP.push_back(qgr);
    }

    std::vector<std::vector<int> > hits;
    for(int j = 0; j < QSP.size(); j++)
    {
        auto match = task.qgramm_lib->find(QSP[j]);
	if(match != task.qgramm_lib->end())
	    hits.push_back(match->second);
	else break;
    }

    if(hits.size() != QSP.size()) return results;

    else for(int i = 0; i < hits.size()-1; i++) {
	std::vector<int> accepted_hits;
	for(int x1 = 0; x1 < hits[i+1].size(); x1 += 2) {
	    bool accepted = false;
	    for(int x2 = 0; x2 < hits[i].size(); x2 += 2) {
		if(hits[i+1][x1] == hits[i][x2] &&
		   hits[i+1][x1+1] - task.M*task.Q == hits[i][x2+1]) {
		    accepted = true; break;
		}
	    }
	    if(accepted) {
		accepted_hits.push_back(hits[i+1][x1]);
		accepted_hits.push_back(hits[i+1][x1+1]);
	    }
	}
	hits[i+1] = accepted_hits;
    }
    for(int i = 0; i < hits[hits.size()-1].size(); i+=2) {
	if(hits[hits.size()-1][i+1]-task.phase-static_cast<int>(hits.size()-1)*task.M*task.Q >= 0)
	{
	    results.push_back(hits[hits.size()-1][i]);
	    results.push_back(hits[hits.size()-1][i+1]-task.phase-static_cast<int>(hits.size()-1)*task.M*task.Q);
	}
    }
    return results;
}

std::vector<int> Ungapped_Find_Pattern(const std::string& pattern, int M, int Q,
        std::unordered_map<std::string, std::vector<int> >& qgramm_lib)
{
    std::vector<int> results;
    std::vector<std::future<std::vector<int> > > job_results;
    for(int phase = 0; phase < M; phase++)
    {
	std::string PM = "";
	for(int i = phase; i < pattern.size(); i += M)
	    PM += pattern[i];
	InnerPatternMatchTask task;
	task.Q=Q; task.M=M;
	task.qgramm_lib = &qgramm_lib;
	task.PM = PM;
	task.phase = phase;

	job_results.push_back(std::async(std::launch::async,
	    Ungapped_Find_Pattern_For_One_Polyphase_Task, task));
    }
    for(auto &x : job_results)
    {
	auto result = x.get();
	results.insert(results.begin(), result.begin(), result.end());
    }
    return results;
}

bool Ungapped_Match_Pattern(const std::string& pattern, int collection_key, int pattern_start_pos,
        std::unordered_map<int, std::string>& read_collection,
        std::unordered_map<std::string, std::vector<int> >& qgramm_lib, int mm_count)
{
    std::string putative_region = read_collection[collection_key].substr(pattern_start_pos);
    if(pattern.size() > putative_region.size()) return false;
    int mm_seen = 0;
    for(size_t pos = 0; pos < pattern.size(); pos++)
    {
    	if(pattern[pos] != putative_region[pos])
    	{
    	    mm_seen++;
    	    if(mm_seen > mm_count) return false;
    	}
    }
    return true;
}


// INTERFACES

void Process_Reads_FASTQ(const char* fastq, int M, int Q,
        std::unordered_map<int, std::string>& read_collection,
        std::unordered_map<int, std::string>& read_names,
        std::unordered_map<std::string, std::vector<int> >& qgramm_lib)
{
    std::ifstream ifastq(fastq); std::string buffer_line;
    int line_counter = 0; int read_counter = 0;
    while(std::getline(ifastq, buffer_line, '\n')) {
	if(line_counter % 4 == 0) {
	    read_counter++;
	    read_names[read_counter] = buffer_line;

	}
	if(line_counter % 4 == 1) {
	    read_collection[read_counter] = buffer_line;
        read_counter++;
        read_names[read_counter] = libdna::rcDNA(buffer_line);
	}
	line_counter++;
    }
    #ifdef DEBUG_FASTMATCH
	std::cout << "Reads read from file " << fastq << ": " << read_counter << std::endl;
    #endif

    Preprocess_Collection(M, Q, read_collection, qgramm_lib);
}

std::vector<std::vector<int> > Locate_Pattern_With_MM(const std::string& P, int M, int Q,
        std::unordered_map<int, std::string>& read_collection,
        std::unordered_map<std::string, std::vector<int> >& qgramm_lib, int MM)
{
    std::vector<std::vector<int> > results;

    std::vector<int> preliminary_location = Ungapped_Find_Pattern(P, M, Q,
        qgramm_lib);

    for(int i = 0; i < preliminary_location.size(); i+=2)
    {
    	if(Ungapped_Match_Pattern(P, preliminary_location[i],
    	    preliminary_location[i+1], read_collection, qgramm_lib, MM))
    	{
    	    std::vector<int> hit(2,0);
    	    hit[0] = preliminary_location[i];
    	    hit[1] = preliminary_location[i+1];
    	    results.push_back(hit);
    	}
    }
    return results;
}




#ifdef FASTMATCH


int main(int argc, char** argv)
{
    // TEST 1
    /*
    std::cout << "Unit test for T-less aligner" << std::endl;

    std::unordered_map<int, std::string> collection;
    collection[1] = "TESTSTRINGOFSOMETHINGGOFSOMETHING";
    collection[2] = "TESTSTRINGOFNOTHING";

    std::unordered_map<std::string, std::vector<int> > qgramm_lib;
    Preprocess_Collection(2, 3, collection, qgramm_lib);

    std::string P = "GOFSOMETHINGGOFSOMETHING";
    std::vector<int> prelim_results = Ungapped_Find_Pattern(P, 2, 3, qgramm_lib);

    for(auto x : prelim_results)
	std::cout << x << ",";
	std::cout << std::endl;

    */

    // TEST 2
    /*
    std::unordered_map<int, std::string> collection, patterns;

    std::ifstream ifastq("/Users/jalgard/Bioinf/testlib.fastq");
    int line_counter = 0; int read_counter = 0; std::string buffer_line;
    while(!ifastq.eof()) {
	std::getline(ifastq, buffer_line, '\n');
	if(line_counter % 4 == 1) {
	    read_counter++; collection[read_counter] = buffer_line;
	}
	line_counter++;
    }

    std::cout << "Total reads read: " << read_counter << std::endl;

    std::ifstream ipattern("/Users/jalgard/Bioinf/testpatterns.fastq");
    line_counter = 0; read_counter = 0;
    while(!ipattern.eof()) {
	std::getline(ipattern, buffer_line, '\n');
	if(line_counter % 4 == 1) {
	    read_counter++; patterns[read_counter] = buffer_line;
	}
	line_counter++;
    }

    std::cout << "Total patterns read: " << read_counter << std::endl;


    std::unordered_map<std::string, std::vector<int> > qgramm_lib;

    int Q = std::atoi(argv[1]); int M = std::atoi(argv[2]); int MM = std::atoi(argv[3]);

    Preprocess_Collection(Q, M, collection, qgramm_lib);


    for(auto x : patterns)
    {
	std::cout << "Pattern " << x.first << std::endl;
	std::vector<int> prelim_results = Ungapped_Find_Pattern(x.second, Q, M, qgramm_lib);
	std::cout << prelim_results.size() << " : ";
	std::vector<int> true_matched;
	for(int i = 0; i < prelim_results.size(); i+=2)
	{
	    if(Ungapped_Match_Pattern(x.second, prelim_results[i], prelim_results[i+1],
		collection,qgramm_lib, MM))
	    {
		true_matched.push_back(prelim_results[i]);
		true_matched.push_back(prelim_results[i+1]);
	    }

	}
	std::cout << true_matched.size();
	std::cout << std::endl;
    }
    */

    // TEST 3
    std::unordered_map<int, std::string> collection, names, patterns;
    std::unordered_map<std::string, std::vector<int> > qgramm_lib;

    std::ifstream ipattern("/Users/jalgard/Bioinf/testpatterns.fastq");
    int line_counter = 0; int read_counter = 0; std::string buffer_line;
    while(!ipattern.eof()) {
	std::getline(ipattern, buffer_line, '\n');
	if(line_counter % 4 == 1) {
	    read_counter++; patterns[read_counter] = buffer_line;
	}
	line_counter++;
    }
    std::cout << "Total patterns read: " << read_counter << std::endl;
    int Q = std::atoi(argv[1]); int M = std::atoi(argv[2]); int MM = std::atoi(argv[3]);
    Process_Reads_FASTQ("/Users/jalgard/Bioinf/testlib.fastq", M, Q, collection,
        names, qgramm_lib);
    for(auto x : patterns)
    {
	std::cout << "Pattern " << x.first << std::endl;
	auto results = Locate_Pattern_With_MM(x.second, M, Q, collection, qgramm_lib, MM);
	std::cout << "Matches: " << results.size() << std::endl;
    }
    return 0;
}

#endif
