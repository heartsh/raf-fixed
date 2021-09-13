/////////////////////////////////////////////////////////////////
// IO.hpp
//
// Routines for input/output.
/////////////////////////////////////////////////////////////////

#ifndef IO_HPP
#define IO_HPP

#include <vector>
#include <string>
#include "Utilities.hpp"
#include "SparseMatrix.hpp"

/////////////////////////////////////////////////////////////////
// struct Sequence
/////////////////////////////////////////////////////////////////

struct Sequence {
    std::string name;
    std::string data;
 
    Sequence() : name(""), data("") {}
    Sequence(const Sequence &rhs) : name(rhs.name), data(rhs.data) {}
    Sequence(const std::string &name, const std::string &data) : name(name), data(data) {}
    Sequence &operator=(const Sequence &rhs)
    {
        if (this != &rhs)
        {
            name = rhs.name;
            data = rhs.data;
        }
        return *this;
    }        
};

/////////////////////////////////////////////////////////////////
// class IO
/////////////////////////////////////////////////////////////////

class IO
{
public:
    static std::vector<Sequence> ReadMFA(const std::string &filename);
    static std::vector<Sequence> ReadUnalignedSequences(const std::string &filename);
    static std::vector<Sequence> ReadAlignmentWithConsensus(const std::string &filename, std::string &consensus);
    static std::vector<Sequence> StripGaps(const std::vector<Sequence> &alignment);
    
    static void WriteToFASTA(const std::string &filename,
                             const std::string &name,
                             const std::string &sequence);
    
    template<class T>
    static SparseMatrix<T> ReadBasePairingPosteriors(const std::string &filename, 
                                                     const T base_pairing_posterior_cutoff,
                                                     const int min_hairpin_length);
    
    template<class T>
    static SparseMatrix<T> ReadViennaRNABasePairingPosteriors(const std::string &filename, 
                                                              const T base_pairing_posterior_cutoff,
                                                              const int min_hairpin_length);
    
    template<class T>
    static SparseMatrix<T> ReadAlignmentPosteriors(const std::string &filename, 
                                                   const T alignment_posterior_cutoff);
    
    template<class T>
    static SparseMatrix<T> ComputeBasePairingPosteriorsCONTRAfold(const std::string &contrafold_dir,
                                                                  const std::string &tmp_dir_name,
                                                                  const int index,
                                                                  const bool toggle_noncomplementary,
                                                                  const bool toggle_verbose,
                                                                  const T base_pairing_posterior_cutoff,
                                                                  const int min_hairpin_length);
    
    template<class T>
    static SparseMatrix<T> ComputeBasePairingPosteriorsViennaRNA(const std::string &vienna_rna_dir,
                                                                 const std::string &tmp_dir_name,
                                                                 const int index,
                                                                 const bool toggle_verbose,
                                                                 const T base_pairing_posterior_cutoff,
                                                                 const int min_hairpin_length);
    
    template<class T>
    static SparseMatrix<T> ComputeAlignmentPosteriorsCONTRAlign(const std::string &probcons_rna_dir,
                                                                const std::string &tmp_dir_name,
                                                                const int index1,
                                                                const int index2,
                                                                const bool toggle_verbose,
                                                                const T alignment_posterior_cutoff);
};

#include "IO.ipp"

#endif
