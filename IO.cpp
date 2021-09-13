/////////////////////////////////////////////////////////////////
// IO.cpp
/////////////////////////////////////////////////////////////////

#include "IO.hpp"

/////////////////////////////////////////////////////////////////
// IO::ReadMFA()
//
// Read MFA format file.
/////////////////////////////////////////////////////////////////

std::vector<Sequence> IO::ReadMFA(const std::string &filename)
{
    // attempt to open file
    
    std::ifstream infile(filename.c_str());
    if (infile.fail()) Error("Unable to open MFA input file: %s", filename.c_str());
    
    std::vector<Sequence> ret;
    std::string s;

    // read lines from file
    
    while (getline(infile, s))
    {
        if (s.length() == 0) continue;

        // check for header line
        
        if (s[0] == '>')
        { 
            ret.push_back(Sequence(s.substr(1), std::string("@")));
        }

        // process non-header lines
        
        else
        {
            if (ret.size() == 0) Error("Expected header line in MFA file: %s", filename.c_str());
            for (size_t i = 0; i < s.length(); i++)
            {
                // skip whitespace
                
                if (isspace(s[i])) continue;
                if (!isalpha(s[i]) && s[i] != '-' && s[i] != '.' && s[i] != '(' && s[i] != ')')
                    Error("Unexpected character '%c' in file \"%s\"", s[i], filename.c_str());
                
                // save characters
                
                ret.back().data.push_back(s[i]);
            }
        }
    }
    
    return ret;  
}

/////////////////////////////////////////////////////////////////
// IO::ReadUnalignedSequences()
//
// Read unaligned sequences in MFA format.
/////////////////////////////////////////////////////////////////

std::vector<Sequence> IO::ReadUnalignedSequences(const std::string &filename)
{
    std::vector<Sequence> ret = StripGaps(ReadMFA(filename));
    
    // check that all sequences have only letters
    
    for (size_t i = 0; i < ret.size(); i++)
    {
        for (size_t j = 1; j < ret[i].data.length(); j++)
        {
            if (!isalpha(ret[i].data[j]))
                Error("Unaligned sequences should contain only letters: %s.", ret[i].data.substr(1).c_str());
        }
    }
    
    return ret;
}

/////////////////////////////////////////////////////////////////
// IO::ReadAlignmentWithConsensus()
//
// Read aligned sequences in MFA format along with consensus
// structure.
/////////////////////////////////////////////////////////////////

std::vector<Sequence> IO::ReadAlignmentWithConsensus(const std::string &filename, std::string &consensus)
{
    std::vector<Sequence> ret = ReadMFA(filename);
    Assert(ret.size() > 0, "Input file does not contain consensus sequence.");

    // check that all sequences except last have only letters and gap characters;
    // set '-' as the gap character and '.' as the unpaired character
    
    for (size_t i = 0; i < ret.size() - 1; i++)
    {
        for (size_t j = 1; j < ret[i].data.length(); j++)
        {
            if (ret[i].data[j] == '.') ret[i].data[j] = '-';
            if (ret[i].data[j] == '(' || ret[i].data[j] == ')')
                Error("Consensus structure should be last in the input alignment file.");
        }
    }
    
    // check that the last sequence is a consensus structure

    for (size_t j = 1; j < ret.back().data.length(); j++)
    {
        if (isalpha(ret.back().data[j])) Error("Last sequence of input alignment file should be a consensus structure.");

        // convert -'s to .'s
        
        if (ret.back().data[j] == '-') ret.back().data[j] = '.';
    }
    
    consensus = ret.back().data;
    ret.pop_back();
    
    return ret;  
}

/////////////////////////////////////////////////////////////////
// IO::StripGaps()
//
// Project alignment to unaligned sequences.
/////////////////////////////////////////////////////////////////

std::vector<Sequence> IO::StripGaps(const std::vector<Sequence> &alignment)
{
    std::vector<Sequence> ret(alignment.size());
    
    for (size_t i = 0; i < alignment.size(); i++)
    {
        ret[i].name = alignment[i].name;
        ret[i].data = RemoveGaps(alignment[i].data);
    }
    
    return ret;  
}

/////////////////////////////////////////////////////////////////
// IO::WriteToFASTA()
//
// Write sequence to FASTA file
/////////////////////////////////////////////////////////////////

void IO::WriteToFASTA(const std::string &filename,
                      const std::string &name,
                      const std::string &sequence)
{
    std::ofstream outfile(filename.c_str());
    if (outfile.fail()) Error("Unable to open temporary file for writing.");
    outfile << ">" << name << std::endl
            << sequence.substr(1) << std::endl;
    outfile.close();
}
