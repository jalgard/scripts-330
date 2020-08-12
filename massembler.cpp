// libdna common tools
#include "dnacommon.h"

// libdna fastmatch
#include "fastmatch.h"
#include<algorithm>

using namespace std;
using namespace libdna;

// for min read length 100
const int M = 8;    // 11
const int Q = 8;    // 8


const int L = 70; // set L, K and Z according to the assembly protocol
const int K = 60;
const int Z = 100;

int Verify_Overlap(const string& Left,  int match_at_L,
                   const string& Right, int match_at_R, int Ol)
{
    int l_ov_s = match_at_L - match_at_R;
    if(l_ov_s < 0) return -1;                     // L str is included in R

    int max_overlength = Left.size() - l_ov_s;
    if(max_overlength >= Right.size()) return -2; // R str is included in L

    if(Left.substr(l_ov_s, max_overlength) == Right.substr(0, max_overlength))
        return max_overlength;
    return 0;
}

// will extend seed sequence 1 step forward
void Extend_Seed_At_Prime(string& seed, unordered_map<string, vector<int> >& qgramm_lib,
    unordered_map<int, string>& read_collection)
{
    // try to extend leftward ----> (seed)
    //                           -------- (read)
    if(seed.size() < L) return;

    auto seed_suffix = seed.substr(seed.size() - L, L);

    auto hits = Locate_Pattern_With_MM(seed_suffix, M, Q, read_collection, qgramm_lib, 0);

    map<string, int> kmer_counts;
    int index = 0;
    map<int, string> index2kmer;

    for(auto& x : hits)
    {
        if(read_collection[x[0]].size() > x[1] + L + K)
        {
            auto kmer = read_collection[x[0]].substr(x[1]+L, K);
            kmer_counts[kmer] += 1;
            index2kmer[index] = kmer;
        }
        index++;
    }
    index = 0;
    for(auto& x : hits)
    {
        x.push_back(kmer_counts[index2kmer[index]]);
        index++;
    }
    sort(hits.begin(), hits.end(), [](vector<int> a, vector<int> b)
			{ return a[1] > b[1]; });
    stable_sort(hits.begin(), hits.end(), [](vector<int> a, vector<int> b)
        	{ return a[2] > b[2]; });

    bool left_extended = false; int j = 0;

    while(!left_extended && j < hits.size())
    {
        auto rd = read_collection[hits[j][0]];

        int lol = Verify_Overlap(seed, seed.size()-L, rd, hits[j][1], L);
        if(lol > L && rd.size()-lol > K) {
            left_extended = true;
            seed = seed + rd.substr(lol, K);
        }
        j++;
    }
}

vector<vector<string> > Load_Seeds(const char* seedfile)
{
    vector<vector<string> > primers;
    map<string, string> fae;
    read_fasta_upper(seedfile, fae);
    for(auto& i : fae) {
        vector<string> pri(2,"");
        pri[0] = i.second; pri[1] = i.first;
        primers.push_back(pri);
    }
    return primers;
}

int Check_Circle(string seed, const string& prime)
{
    int msz = -1;
    while(1)
    {
        size_t sp = seed.rfind(prime);
        if(sp == string::npos || sp == 0) return msz;
        msz = sp;
        seed = seed.substr(0, msz);
    }

}

int main(int argc, char** argv)
{

    unordered_map<int, string> read_collection, read_names;
    unordered_map<string, vector<int> > qgramm_lib;

    auto primers = Load_Seeds(argv[2]);
    cout << primers.size() << " primers to extned\n";
    // Load and hash reads
    Process_Reads_FASTQ(argv[1], M, Q, read_collection, read_names, qgramm_lib);
    cout << "Reads loaded!\n";


    cout << "Reads from file " << argv[1] << " hashed\n";
    ofstream outlog("miniassmL.log");
    ofstream outfas("assemblyL.fa");
    int seed_count = 0;
    for(int j = 0; j < primers.size(); j++)
    {
        // extend each pattern

        string prime = primers[j][0];
        string seed = prime;
        string seed_name = primers[j][1];

        cout << "Extending seed " << seed_name << "\n";
        for(int k = 0; k < 35; k++)
        {
            Extend_Seed_At_Prime(seed, qgramm_lib, read_collection);
        }
        int circle_sz = Check_Circle(seed, prime);

        if(circle_sz > 1000)
        {
            outfas << ">" << seed_name << "\n" << seed.substr(0, circle_sz) << "\n";
            outlog << seed_name << "\tSz:" << circle_sz <<  "\n";
            cout << "Assembled! Size:" << circle_sz << "\n";
        }
        else
        {
            outlog << seed_name << "\tSz:" << circle_sz << "\tUnassembled(" << seed.size() << ")\n";
            cout << "Unassembled! Seed size: " << seed.size() << "\n";
            outlog << "\n" << seed << "\n...\n";
        }
    }
    return 0;
}
