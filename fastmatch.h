#ifndef TALIGNER_FASTMATCH
#define TALIGNER_FASTMATCH

#include<string>
#include<vector>
#include <unordered_map>
#include<future>

#include<iostream>
#include<fstream>


/*
 *  Preprocess read collection to Q-gramm library
 *  M - resize factor
 *  Q - qgramm length
 */
void Preprocess_Collection(int M, int Q, std::unordered_map<int, std::string>& read_collection,
        std::unordered_map<std::string, std::vector<int> >& qgramm_lib);

/*
 *  Locates putative places in reads collection, 
 *  where the pattern can be found (strongly matches each M-th nucleotide
 *  in all possible phase shift)
 */

//std::vector<int> Ungapped_Find_Pattern(const std::string& pattern, int M, int Q,
//        std::unordered_map<std::string, std::vector<int> >& qgramm_lib);

/*
 * MULTITHREADED IMPLEMENTATION
 */

/*
 * Wrapper for task parameters
 * Task is each polyphase decomposition
 */

struct InnerPatternMatchTask
{
    int phase;
    std::unordered_map<std::string, std::vector<int> >* qgramm_lib;
    std::string PM;
    int Q;
    int M;
};

/*
 * One task matching funtion. Task is one polyphase decomposition
 */

std::vector<int> Ungapped_Find_Pattern_For_One_Polyphase_Task(InnerPatternMatchTask task);

/*
 * MULTITHREADED matching function, threads number = M;
 */

std::vector<int> Ungapped_Find_Pattern(const std::string& pattern, int M, int Q,
        std::unordered_map<std::string, std::vector<int> >& qgramm_lib);

/*
 *  Verifies putative match for pattern with at most mm_count mismatches
 *
 */

bool Ungapped_Match_Pattern(const std::string& pattern, int collection_key, int pattern_start_pos,
        std::unordered_map<int, std::string>& read_collection,
        std::unordered_map<std::string, std::vector<int> >& qgramm_lib, int mm_count = 0);

/*
 * Interface function!
 * Returns the vector of vectors <read ID, position> of verified matches
 * of pattern P
 * with at most MM mismatches allowed
 *
 */

std::vector<std::vector<int> > Locate_Pattern_With_MM(const std::string& P, int M, int Q,
        std::unordered_map<int, std::string>& read_collection,
        std::unordered_map<std::string, std::vector<int> >& qgramm_lib, int MM = 0);

/*
 * Interface function!
 * Processes FASTQ file and returns reads collection, read ID 2 reads name table,
 * qgramm lib with M,Q
 *
 */

void Process_Reads_FASTQ(const char* fastq, int M, int Q,
        std::unordered_map<int, std::string>& read_collection,
        std::unordered_map<int, std::string>& read_names,
        std::unordered_map<std::string, std::vector<int> >& qgramm_lib);




#endif