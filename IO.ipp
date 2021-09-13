/////////////////////////////////////////////////////////////////
// IO.ipp
/////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////
// IO::ReadBasePairingPosteriors()
//
// Read base-pairing probabilities.
/////////////////////////////////////////////////////////////////

template<class T>
SparseMatrix<T> IO::ReadBasePairingPosteriors(const std::string &filename, 
                                              const T base_pairing_posterior_cutoff,
                                              const int min_hairpin_length)
{
    Assert(base_pairing_posterior_cutoff >= 0, "Base-pairing posterior cutoff should be nonnegative: %lf", double(base_pairing_posterior_cutoff));
    
    std::ifstream infile(filename.c_str());
    if (infile.fail()) Error("Unable to load base-pairing probability file: %s", filename.c_str());
    
    std::map<std::pair<int,int>, T> elems;
    std::string token;
    int row = 0;
    
    // read input file
    
    while (infile >> token)
    {
        std::string::size_type colon_pos = token.find(':');
        if (colon_pos == std::string::npos)
        {
            
            // read row
            
            int index = 0;
            if (!ConvertToNumber(token, index)) Error("Could not read row: %s", filename.c_str());
            if (index <= 0) Error("Row numbers must be positive: %s", filename.c_str());
            if (index != row+1) Error("Rows of base-pairing probability matrix must occur in increasing order: %s", filename.c_str());
            row = index;
            
            // read sequence letter
            
            if (!(infile >> token)) Error("Expected sequence letter after row number: %s", filename.c_str());
            if (token.length() != 1) Error("Expected sequence letter after row number: %s", filename.c_str());      
            if (!isalpha(token[0])) Error("Unknown character '%c' in base pair probability file: %s", token[0], filename.c_str());
            
        }
        else
        {
            
            // read column:value
            
            int col;
            T value;
            if (!ConvertToNumber(token.substr(0, colon_pos), col)) Error("Could not read column number: %s", filename.c_str());
            if (col <= 0) Error("Column index must be positive: %s", filename.c_str());
            if (!ConvertToNumber(token.substr(colon_pos+1), value)) Error("Could not read matrix entry value: %s", filename.c_str());
            if (value < T(-1e-6)) Error("Alignment match probability matrix should not contain negative values: %s", filename.c_str());
            if (value > T(1.0+1e-6)) Error("Alignment match probability matrix should not contain values over 1: %s", filename.c_str());
            if (value >= base_pairing_posterior_cutoff)
            {
                if (row + min_hairpin_length <= col)
                    elems[std::make_pair(row,col)] = value;
            }
        }
    }
    infile.close();
    
    return SparseMatrix<T>(elems, row+1, row+1, 0);
}

/////////////////////////////////////////////////////////////////
// IO::ReadViennaRNABasePairingPosteriors()
//
// Read ViennaRNA base-pairing probabilities.
/////////////////////////////////////////////////////////////////

template<class T>
SparseMatrix<T> IO::ReadViennaRNABasePairingPosteriors(const std::string &filename, 
                                                       const T base_pairing_posterior_cutoff,
                                                       const int min_hairpin_length)
{
    Assert(base_pairing_posterior_cutoff >= 0, "Base-pairing posterior cutoff should be nonnegative: %lf", double(base_pairing_posterior_cutoff));
    
    std::ifstream infile(filename.c_str());
    if (infile.fail()) Error("Unable to load ViennaRNA base-pairing probability file: %s", filename.c_str());
    
    std::map<std::pair<int,int>, T> elems;
    std::string s;
    int row, col, L = 0;
    T value;

    // read input file

    while (std::getline(infile, s))
    {
        std::istringstream iss(s);
        if (iss >> s)
        {
            if (s == "/sequence")
            {
                if (L != 0) Error("Too many sequences in file.");
                while (std::getline(infile, s))
                {
                    if (s.length() == 0) continue;
                    if (s[0] == ')') break;
                    L += int(s.length()) - 1;
                }
            }
            else if (ConvertToNumber(s, row))
            {
                if (iss >> col >> value >> s)
                {
                    if (s == "ubox")
                    {
                        if (row <= 0) Error("Row numbers must be positive: %s", filename.c_str());
                        if (row > L) Error("Row numbers must be at most the length of the sequence: %s", filename.c_str());
                        if (col <= row) Error("Column numbers must be greater than row numbers: %s", filename.c_str());
                        if (col > L) Error("Column numbers must be at most the length of the sequence: %s", filename.c_str());
                        if (value < T(-1e-6)) Error("Alignment match probability matrix should not contain negative values: %s", filename.c_str());
                        if (value > T(1.0+1e-6)) Error("Alignment match probability matrix should not contain values over 1: %s", filename.c_str());

                        if (value >= base_pairing_posterior_cutoff)
                        {
                            if (row + min_hairpin_length <= col)
                                elems[std::make_pair(row,col)] = value;
                        }
                    }
                }                
            }
        }
    }
    
    infile.close();
    
    return SparseMatrix<T>(elems, L+1, L+1, 0);
}

/////////////////////////////////////////////////////////////////
// IO::ReadAlignmentPosteriors()
//
// Read alignment match probabilities.
/////////////////////////////////////////////////////////////////

template<class T>
SparseMatrix<T> IO::ReadAlignmentPosteriors(const std::string &filename, 
                                            const T alignment_posterior_cutoff)
{
    Assert(alignment_posterior_cutoff >= 0, "Alignment posterior cutoff should be nonnegative: %lf", double(alignment_posterior_cutoff));
    
    std::ifstream infile(filename.c_str());
    if (infile.fail()) Error("Unable to load alignment match probability file: %s", filename.c_str());

    // get sequence lengths
    
    std::vector<int> lengths;
    std::string s;

    // parse MFA header
    
    while (getline(infile, s))
    {
        if (s.length() == 0) continue;
        if (s == "#") break;

        if (s[0] == '>')
            lengths.push_back(0);
        else
        {
            if (lengths.size() == 0) Error("Expected MFA header line in file: %s", filename.c_str());
            for (size_t i = 0; i < s.length(); i++)
                if (isalpha(s[i])) ++(lengths.back());
        }
    }

    if (lengths.size() != 2) Error("Expected to read pairwise MFA first in file: %s", filename.c_str());
    if (s != "#") Error("Expected MFA to be followed by separator character '#'.");

    // retain only sequence lengths

    const int LA = lengths[0];
    const int LB = lengths[1];

    /* TODO: get rid of probcons
    int LA, LB;
    infile >> LA >> LB;
    */

    if (LA <= 0) Error("Length of first sequence in alignment match probability file must be positive: %s", filename.c_str());
    if (LB <= 0) Error("Length of second sequence in alignment match probability file must be positive: %s", filename.c_str());
    
    std::map<std::pair<int,int>, T> elems;
    std::string token;
    int row = 0;
    
    // read input file
    
    while (infile >> token)
    {
        std::string::size_type colon_pos = token.find(':');
        if (colon_pos == std::string::npos)
        {
            // read row
            
            int index = 0;
            if (!ConvertToNumber(token, index)) Error("Could not read row: %s", filename.c_str());
            if (index <= 0) Error("Row numbers must be positive: %s", filename.c_str());
            if (index <= row) Error("Rows of alignment match probability matrix must occur in increasing order: %s", filename.c_str());
            row = index;
            if (row > LA) Error("Row index out-of-range: %s", filename.c_str());
        }
        else
        {
            // read column:value
            
            int col;
            T value;
            if (!ConvertToNumber(token.substr(0, colon_pos), col)) Error("Could not read column number: %s", filename.c_str());
            if (col <= 0 || col > LB) Error("Column index out-of-range: %s", filename.c_str());
            if (!ConvertToNumber(token.substr(colon_pos+1), value)) Error("Could not read matrix entry value: %s", filename.c_str());
            if (value < T(-1e-6)) Error("Alignment match probability matrix should not contain negative values: %s", filename.c_str());
            if (value > T(1.0+1e-6)) Error("Alignment match probability matrix should not contain values over 1: %s", filename.c_str());
            if (value >= alignment_posterior_cutoff)
                elems[std::make_pair(row,col)] = value;
        }
    }
    infile.close();
    
    return SparseMatrix<T>(elems, LA+1, LB+1, 0);
}

/////////////////////////////////////////////////////////////////
// IO::ComputeBasePairingPosteriorsCONTRAfold()
//
// Compute base-pairing probabilities.
/////////////////////////////////////////////////////////////////

template<class T>
SparseMatrix<T> IO::ComputeBasePairingPosteriorsCONTRAfold(const std::string &contrafold_dir,
                                                           const std::string &tmp_dir_name,
                                                           const int index,
                                                           const bool toggle_noncomplementary,
                                                           const bool toggle_verbose,
                                                           const T base_pairing_posterior_cutoff,
                                                           const int min_hairpin_length)
{
    const double starting_time = GetSystemTime();
    const std::string sequence_filename = SPrintF("%s/seq%d.fasta", tmp_dir_name.c_str(), index);
    const std::string posteriors_filename = SPrintF("%s/seq%d.posteriors", tmp_dir_name.c_str(), index);  
    const std::string command = SPrintF("%s/contrafold predict %s %s --posteriors %lf %s",
                                        contrafold_dir.c_str(),
                                        (toggle_noncomplementary ? "--noncomplementary" : ""), 
                                        sequence_filename.c_str(), 
                                        double(base_pairing_posterior_cutoff),
                                        posteriors_filename.c_str());
    
    if (toggle_verbose) WriteProgressMessage(SPrintF("Running: %s", command.c_str()));
    if (system(command.c_str())) Error("Error running CONTRAfold.");
    const SparseMatrix<T> res = ReadBasePairingPosteriors<T>(posteriors_filename, base_pairing_posterior_cutoff, min_hairpin_length);
    if (toggle_verbose)
    {
        WriteProgressMessage("");
        std::cerr << SPrintF("Folded %s (%d): %lf seconds", 
                             sequence_filename.c_str(), res.GetNumRows()-1, 
                             GetSystemTime() - starting_time) << std::endl;
    }
    return res;
}

/////////////////////////////////////////////////////////////////
// IO::ComputeBasePairingPosteriorsViennaRNA()
//
// Compute base-pairing probabilities.
/////////////////////////////////////////////////////////////////

template<class T>
SparseMatrix<T> IO::ComputeBasePairingPosteriorsViennaRNA(const std::string &vienna_rna_dir,
                                                          const std::string &tmp_dir_name,
                                                          const int index,
                                                          const bool toggle_verbose,
                                                          const T base_pairing_posterior_cutoff,
                                                          const int min_hairpin_length)
{
    const double starting_time = GetSystemTime();
    const std::string sequence_filename = SPrintF("%s/seq%d.fasta", tmp_dir_name.c_str(), index);
    const std::string posteriors_filename = SPrintF("%s/seq%d.posteriors", tmp_dir_name.c_str(), index);  
    const std::string intermediate_filename = SPrintF("%s/seq%d_dp.ps", tmp_dir_name.c_str(), index);
    const std::string relative_sequence_filename = (sequence_filename[0] == '/' ? sequence_filename : "../" + sequence_filename);
    const std::string command = SPrintF("cd %s; %s/RNAfold -p < %s > /dev/null; cd ..",
                                        tmp_dir_name.c_str(), vienna_rna_dir.c_str(), relative_sequence_filename.c_str());
    if (toggle_verbose) WriteProgressMessage(SPrintF("Running: %s", command.c_str()));
    if (system(command.c_str())) Error("Error running ViennaRNA.");
    const SparseMatrix<T> res = ReadViennaRNABasePairingPosteriors<T>(intermediate_filename.c_str(), base_pairing_posterior_cutoff, min_hairpin_length);

    if (toggle_verbose)
    {
        WriteProgressMessage("");
        std::cerr << SPrintF("Folded %s (%d): %lf seconds", 
                             sequence_filename.c_str(), res.GetNumRows()-1, 
                             GetSystemTime() - starting_time) << std::endl;
    }
    return res;
}

/////////////////////////////////////////////////////////////////
// IO::ComputeAlignmentPosteriorsCONTRAlign()
//
// Compute alignment match probabilities.
/////////////////////////////////////////////////////////////////

template<class T>
SparseMatrix<T> IO::ComputeAlignmentPosteriorsCONTRAlign(const std::string &contralign_dir,
                                                         const std::string &tmp_dir_name,
                                                         const int index1,
                                                         const int index2,
                                                         const bool toggle_verbose,
                                                         const T alignment_posterior_cutoff)
{
    const double starting_time = GetSystemTime();
    const std::string sequence1_filename = SPrintF("%s/seq%d.fasta", tmp_dir_name.c_str(), index1);
    const std::string sequence2_filename = SPrintF("%s/seq%d.fasta", tmp_dir_name.c_str(), index2);
    const std::string posteriors_filename = SPrintF("%s/seq%d-%d.posteriors", tmp_dir_name.c_str(), index1, index2);

    const std::string command = SPrintF("%s/contralign predict %s %s --conflate --posteriors %lf %s 2> /dev/null",
                                        contralign_dir.c_str(),
                                        sequence1_filename.c_str(),
                                        sequence2_filename.c_str(),
                                        double(alignment_posterior_cutoff),
                                        posteriors_filename.c_str());
    
    if (toggle_verbose) WriteProgressMessage(SPrintF("Running: %s", command.c_str()));
    if (system(command.c_str())) Error("Error running CONTRALIGN.");
    const SparseMatrix<T> res = ReadAlignmentPosteriors<T>(posteriors_filename, alignment_posterior_cutoff);
    if (toggle_verbose)
    {
        WriteProgressMessage("");
        std::cerr << SPrintF("Aligned %s (%d) with %s (%d): %lf seconds", 
                             sequence1_filename.c_str(), res.GetNumRows()-1, 
                             sequence2_filename.c_str(), res.GetNumCols()-1,
                             GetSystemTime() - starting_time) << std::endl;
    }
    return res;
}


/*
/////////////////////////////////////////////////////////////////
// CompleteStructure()
//
// Complete partially filled secondary structure.
/////////////////////////////////////////////////////////////////

std::string CompleteStructure (int index,
			       const std::string &sequence,
			       const std::string &consensus,
			       const bool toggle_complementary){

  if (!toggle_complete) return consensus;
  
  double starting_time = GetSystemTime();
  std::string input_filename = SPrintF ("%s/seq%d.consensus.bpseq", tmp_dir_name.c_str(), index);
  std::ofstream outfile (input_filename.c_str());
  if (outfile.fail()) Error("Unable to open temporary file for writing.");

  Assert(sequence.length() == consensus.length(), "Size mismatch.");

  // build secondary structure mapping

  std::vector<std::pair<char,int> > mapping (1, std::make_pair('@', -1));
  std::vector<int> stack;
  
  for (size_t i = 1; i < consensus.length(); i++){
    if (isalpha(sequence[i])){
      mapping.push_back (std::make_pair(sequence[i],-1));
      if (consensus[i] == ')' && stack.back()){
	mapping[stack.back()].second = int(mapping.size()) - 1;
	mapping[int(mapping.size()) - 1].second = stack.back();
      }
    }
    if (consensus[i] == '('){
      stack.push_back (isalpha(sequence[i]) ? int(mapping.size()) - 1 : 0);
    } else if (consensus[i] == ')'){
      stack.pop_back();
    }
  }
  
  for (size_t i = 1; i < mapping.size(); i++){
    outfile << i << " " << mapping[i].first << " " << mapping[i].second << std::endl;
  }
  outfile.close();

  // run CONTRAfold

  std::string output_filename = SPrintF ("%s/seq%d.completed.bpseq", tmp_dir_name.c_str(), index);

  std::string command = SPrintF ("~/rna_alignment/programs/raf/contrafold predict %s %s --constraints --bpseq %s", 
				 (toggle_complementary ? "" : "--noncomplementary"), 
				 input_filename.c_str(), 
				 output_filename.c_str());
  
  if (toggle_verbose) WriteProgressMessage (SPrintF ("Running: %s", command.c_str()));
  if (system (command.c_str())) Error("Error running CONTRAfold.");

  // parse output
  
  std::ifstream infile (output_filename.c_str());
  if (infile.fail()) Error("Unable to open temporary file for reading.");

  for (size_t i = 1; i < mapping.size(); i++){
    int left, right;
    char ch;
    if (!(infile >> left >> ch >> right)) Error("Temporary file has incorrect format.");
    Assert(left == int(i), "Temporary file has incorrect format.");
    mapping[left].second = right;
  }

  // generate projected secondary structure

  std::string new_struct = "@";
  int pos = 0;
  
  for (size_t i = 1; i < consensus.length(); i++){
    if (isalpha(sequence[i])){
      pos++;
      if (mapping[pos].second == 0){
	new_struct.push_back ('.');
      } else if (mapping[pos].second > pos){
	new_struct.push_back ('(');
      } else if (mapping[pos].second < pos){
	new_struct.push_back (')');
      } else {
	Assert(false, "Should not get here.");
      }
    } else {
      new_struct.push_back ('.');
    }
  }

  if (toggle_verbose) {
    WriteProgressMessage ("");
    std::cerr << SPrintF ("Completed structure %d (%d): %lf seconds", index, int(mapping.size()) - 1, GetSystemTime() - starting_time) << std::endl;
  }

  return new_struct;
}
*/
