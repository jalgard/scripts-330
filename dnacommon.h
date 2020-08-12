#ifndef DNACOMMON_H
#define DNACOMMON_H

#include <string>
#include <locale>
#include <vector>
#include <sstream>
// for work with FASTA/FASTQ
#include <map>
#include <fstream>


namespace libdna {

static char libdnaComplement(char a)
{
    switch(a) {
        case 'A' : return 'T';
        case 'T' : return 'A';
        case 'G' : return 'C';
        case 'C' : return 'G';
        case 'N' : return 'N';
        case 'a' : return 't';
        case 't' : return 'a';
        case 'g' : return 'c';
        case 'c' : return 'g';
        case 'n' : return 'n';
        case '-' : return '-';
        default: return 'N';
    }
}

static char toDNAchar(char a)
{
    switch(a) {
        case 'A' : return 'a';
        case 'T' : return 't';
        case 'G' : return 'G';
        case 'C' : return 'C';
        case 'N' : return 'N';
        case 'a' : return 'a';
        case 't' : return 't';
        case 'g' : return 'g';
        case 'c' : return 'c';
        case 'n' : return 'n';
        case '-' : return '-';
        default: return 'N';
    }
}

static std::string upDNA(const std::string& dna)
{
    std::string tmp=dna;
    for (size_t pos = 0; dna[pos] != '\0'; ++pos)
        tmp[pos] = std::toupper(tmp[pos]);
    return tmp;
}

static std::string lowDNA(const std::string& dna)
{
    std::string tmp=dna;
    for (size_t pos = 0; dna[pos] != '\0'; ++pos)
        tmp[pos] = std::tolower(tmp[pos]);
    return tmp;
}

static std::string comDNA(const std::string& dna)
{
    std::string tmp=dna;
    for(size_t pos = 0; dna[pos] != '\0'; ++pos)
	tmp[pos]=libdnaComplement(dna[pos]);
    return tmp;
}

static std::string revDNA(const std::string& dna)
{
    std::string tmp="";
    for(int pos=dna.length()-1; pos >= 0; --pos)
        tmp += dna[pos];
    return tmp;
}

static std::string rcDNA(const std::string& dna)
{
    std::string tmp = revDNA(dna);
    for(size_t pos = 0; tmp[pos] != '\0'; ++pos)
	tmp[pos]=libdnaComplement(tmp[pos]);
    return tmp;
}

static std::string toDNA(const std::string& dna)
{
    std::string tmp = dna;
    for(size_t pos = 0; tmp[pos] != '\0'; ++pos)
	   tmp[pos]=toDNAchar(tmp[pos]);
    return tmp;
}


static std::vector<std::string>& libdnaSplit(const std::string& s, char delim, std::vector<std::string>& elems)
{
    std::stringstream ss(s); std::string item;
    while (std::getline(ss, item, delim)) elems.push_back(item);
    return elems;
}

static std::vector<std::string> Tokenize(const std::string& s, char delim)
{
    std::vector<std::string> elems;
    libdnaSplit(s, delim, elems);
    return elems;
}

static std::string Token(const std::string& s, const std::string& start, const std::string& end)
{
    std::string token = "";
    std::size_t ss = s.find(start); std::size_t se = s.find(end, ss + start.size() + 1);
    if(ss != std::string::npos && se != std::string::npos && ss < se)
    {
	token = s.substr(ss + start.size(), se - ss - start.size());
    }
    return token;
}

static std::string MergeTokens(std::vector<std::string>& elems, size_t start, size_t end)
{
    std::string merged = "";
    if(end == 0) end = elems.size();
    for(size_t i = start; i < end; i++) merged+=elems[i];
    return merged;
}

static std::string i2s(const int& value)
{
    std::ostringstream tmp;
    tmp << value;
    return tmp.str();
}

static std::string d2s(const double& value, const int prec)
{
    std::ostringstream tmp;
    tmp.precision(prec); tmp  << value;
    return tmp.str();
}


// Generic FASTA reader
static void read_fasta(const char* file, std::map<std::string, std::string>& fae,
	bool name_space_split = false)
{
    std::ifstream ifile(file); std::string line, last;
	while(std::getline(ifile, line, '\n'))
	{
		if(line[0] == '>')
		{
			last = line.substr(1); if(name_space_split) last = Tokenize(last,' ')[0];
		}
		else fae[last]+=line;
	}
}

static void read_fasta_upper(const char* file, std::map<std::string, std::string>& fae,
	bool name_space_split = false)
{
    std::ifstream ifile(file); std::string line, last;
	while(std::getline(ifile, line, '\n'))
	{
		if(line[0] == '>')
		{
			last = line.substr(1); if(name_space_split) last = Tokenize(last,' ')[0];
		}
		else fae[last]+=upDNA(toDNA(line));
	}
}

// Generic FASTQ reader, by default loads reads in map [read name -> read seq]
// can save memory, discarding read's name

static void read_fastq(const char* file, std::map<std::string, std::string>& fqe,
	bool rename_reads_to_ids = false)
{
    std::ifstream ifile(file); std::string line1, line2, line3, line4;
	unsigned int id_count=0;
	while(  std::getline(ifile, line1, '\n') &&
			std::getline(ifile, line2, '\n') &&
			std::getline(ifile, line3, '\n') &&
			std::getline(ifile, line4, '\n'))
	{
		if(line1[0] == '@' && line2.size() > 1 && line2.size() == line4.size())
		{
			if(!rename_reads_to_ids) fqe[line1] = line2;
			else fqe[i2s(++id_count)] = line2;
		}
	}
}

static void read_fastq_quals(const char* file, std::map<std::string, std::string>& fqe,
    std::map<std::string, std::string>& fqq, bool rename_reads_to_ids = false)
{
    std::ifstream ifile(file); std::string line1, line2, line3, line4;
	unsigned int id_count=0; std::string rd_id;
	while(  std::getline(ifile, line1, '\n') &&
			std::getline(ifile, line2, '\n') &&
			std::getline(ifile, line3, '\n') &&
			std::getline(ifile, line4, '\n'))
	{
		if(line1[0] == '@' && line2.size() > 1 && line2.size() == line4.size())
		{
            rd_id = line1; if(rename_reads_to_ids) rd_id = i2s(++id_count);
			fqe[rd_id] = line2;
			fqq[rd_id] = line4;
		}
	}
}

} // libdna

#endif
