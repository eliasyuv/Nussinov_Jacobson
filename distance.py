import RNA

# biological result (NMR)
sequence1 = "GGCAGUACCAAGUCGCGAAAGCGAUGGCCUUGCAAAGGGUAUGGUAAUAAGCUGCC"
structure1 = "(((((((((..(((((....))))).(((((....)))))..))).....))))))"

# nussinov
sequence2 = "GGCAGUACCAAGUCGCGAAAGCGAUGGCCUUGCAAAGGGUAUGGUAAUAAGCUGCC"
structure2 ="(((.(()((..)((()))..)))(((())((())))(..(((..())))..)).))"


distance = RNA.bp_distance(structure1, structure2)

print(f"Base pair distance: {distance}")
