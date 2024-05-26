class Nussinov():
    """The algorithm of Nussinov is an RNA secondary structure folding algorithm. It was developed by Ruth Nussinov et al.
    and was published in 1978:
            Nussinov, Ruth, et al. "Algorithms for loop matchings."
            SIAM Journal on Applied mathematics 35.1 (1978): 68-82.
            http://rci.rutgers.edu/~piecze/GriggsNussinovKleitmanPieczenik.pdf
    """
    def __init__(self, rnaSequence):
        """rnaSequence: The RNA sequence for which the folding should be computed."""
        self.sequence = rnaSequence
        self.pairedBases = {}
        self.computationMatrix = []

    def computeMatrix(self):
        """This function computes the matrix which the Nussinov-algorithm is based on."""
        seq_len = len(self.sequence)
        self.computationMatrix = [[0] * seq_len for _ in range(seq_len)]

        for length in range(1, seq_len):  # length is the span of the subsequence
            for i in range(seq_len - length):  # i is the start index
                j = i + length  # j is the end index
                self.computeMatrixCell(i, j)
          
    def computeMatrixCell(self, i, j):
        """This function computes the value for every cell of the matrix for the Nussinov-algorithm.
            i:  First index of cell of the Nussinov-matrix
            j:  Second index of cell of the Nussinov-matrix
        """
        if i >= j:
            return
        max_value = self.computationMatrix[i][j-1]
        for k in range(i, j):
            if self.complementary(self.sequence[k], self.sequence[j]):
                pairing_value = (self.computationMatrix[i][k-1] if k > i else 0) + \
                               (self.computationMatrix[k+1][j-1] if k+1 <= j-1 else 0) + 1
                max_value = max(max_value, pairing_value)

        self.computationMatrix[i][j] = max_value
       
    def complementary(self, characterA, characterB):
        """Returns True if two RNA nucleotides are complementary, False otherwise.
        Nucleotides are complementary if they are "A" and "U" or "C" and "G".
            characterA: First nucleotide
            characterB: Second nucleotide"""
        complementary_pairs = {('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C')}
        return (characterA, characterB) in complementary_pairs

    def traceback(self, i, j):
        """Computes the traceback for the Nussinov algorithm.
            i: First index of cell of the Nussinov matrix
            j: Second index of cell of the Nussinov matrix
        """
        if j <= i:
            return
        if self.computationMatrix[i][j] == self.computationMatrix[i][j-1]:
            self.traceback(i, j-1)
        else:
            for k in range(i, j):
                if self.complementary(self.sequence[k], self.sequence[j]):
                    if self.computationMatrix[i][j] == (self.computationMatrix[i][k-1] if k > i else 0) + \
                                                       (self.computationMatrix[k+1][j-1] if k+1 <= j-1 else 0) + 1:
                        self.pairedBases[k] = j
                        self.traceback(i, k-1)
                        self.traceback(k+1, j-1)
                        return

    def execute(self):
        """To compute the Nussinov algorithm, execute this method. It returns a dictionary with the paired bases."""
        user_input = input("Enter RNA sequence or press 1 to use the default sequence: ")
        if user_input == '1':
            self.sequence = "GGCAGUACCAAGUCGCGAAAGCGAUGGCCUUGCAAAGGGUAUGGUAAUAAGCUGCC"
        else:
            self.sequence = user_input
        
        self.computeMatrix()
        self.traceback(0, len(self.sequence) - 1)
        print(self.pairedBases)
        print(len(self.pairedBases))
        return self.pairedBases


# Initialize the Nussinov algorithm without a predefined RNA sequence
nussinov_instance = Nussinov("")

# Execute the algorithm to compute the RNA folding.
paired_bases = nussinov_instance.execute()
