import pandas as pd

def find_palindromic_sequences(sequence, min_length=3, max_length=None):
    """
    Find palindromic sequences within the given sequence.

    This function searches for palindromic sequences within the provided sequence.
    A palindromic sequence is a sequence that reads the same forward and backward.

    Parameters:
        sequence (str): The input sequence in which to search for palindromic sequences.
        min_length (int): The minimum length of palindromic sequences to search for (default is 3).
        max_length (int): The maximum length of palindromic sequences to search for.
            If None, the maximum length is set to half the length of the input sequence.

    Returns:
        list: A list of tuples containing palindromic sequences found in the input sequence.
            Each tuple contains the start position, palindromic fragment, end position, and palindromic fragment.
    """
    
    if max_length is None:
        max_length = len(sequence) // 2
    
    palindromic_sequences = []
    
    for length in range(min_length, max_length + 1):
        for i in range(len(sequence) - length + 1):
            fragment = sequence[i:i+length]
            if fragment == fragment[::-1]: 
                for j in range(i + length, len(sequence) - length + 1):
                    if sequence[j:j+length] == fragment:
                        palindromic_sequences.append((i, fragment, j, fragment))
    
    return palindromic_sequences

def find_direct_repeats(sequence, min_length=3, max_length=None):
    """
    Find direct repeat sequences within the given sequence.

    This function searches for direct repeat sequences within the provided sequence.
    A direct repeat sequence is a sequence that is repeated identically in the same orientation.

    Parameters:
        sequence (str): The input sequence in which to search for direct repeat sequences.
        min_length (int): The minimum length of direct repeat sequences to search for (default is 3).
        max_length (int): The maximum length of direct repeat sequences to search for.
            If None, the maximum length is set to half the length of the input sequence.

    Returns:
        list: A list of tuples containing direct repeat sequences found in the input sequence.
            Each tuple contains the start position, repeat fragment, end position, and repeat fragment.
    """
    
    if max_length is None:
        max_length = len(sequence) // 2
    
    direct_repeats = []
    
    for length in range(min_length, max_length + 1):
        for i in range(len(sequence) - length + 1):
            fragment = sequence[i:i+length]
            for j in range(i + length, len(sequence) - length + 1):
                if sequence[j:j+length] == fragment:
                    direct_repeats.append((i, fragment, j, fragment))
    
    return direct_repeats

def find_inverted_repeats(sequence, min_length=3, max_length=None):
    """
    Find inverted repeat sequences within the given sequence.

    This function searches for inverted repeat sequences within the provided sequence.
    An inverted repeat sequence is a sequence that forms a hairpin or stem-loop structure when folded.

    Parameters:
        sequence (str): The input sequence in which to search for inverted repeat sequences.
        min_length (int): The minimum length of inverted repeat sequences to search for (default is 3).
        max_length (int): The maximum length of inverted repeat sequences to search for.
            If None, the maximum length is set to half the length of the input sequence.

    Returns:
        list: A list of tuples containing inverted repeat sequences found in the input sequence.
            Each tuple contains the start position, repeat fragment, end position, and inverted repeat fragment.
    """
    
    if max_length is None:
        max_length = len(sequence) // 2
    
    inverted_repeats = []
    
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    
    for length in range(min_length, max_length + 1):
        for i in range(len(sequence) - length + 1):
            fragment = sequence[i:i+length]
            inverted_fragment = "".join(complement[b] for b in fragment[::-1])
            for j in range(i + length, len(sequence) - length + 1):
                if sequence[j:j+length] == inverted_fragment:
                    inverted_repeats.append((i, fragment, j, inverted_fragment))
    
    return inverted_repeats

def score_sequence(fragment, sequence):
    """
    Score the given fragment within the provided sequence.

    This function calculates a score for the given fragment based on its length, frequency, and distance within the sequence.

    Parameters:
        fragment (str): The fragment to be scored.
        sequence (str): The sequence in which the fragment occurs.

    Returns:
        int: The score calculated for the fragment within the sequence.
    """
    
    # Score based on the length of the fragment
    length_score = len(fragment)
    
    # Score based on the frequency of the fragment
    frequency_score = sequence.count(fragment)
    
    # Score based on the distance between the fragments
    positions = [i for i in range(len(sequence)) if sequence.startswith(fragment, i)]
    if len(positions) > 1:
        distance_score = min(positions[i+1] - positions[i] for i in range(len(positions) - 1))
    else:
        distance_score = len(sequence)
    
    # Total score
    total_score = length_score + frequency_score - distance_score
    return total_score

def find_and_score_sequences(sequence, meme_output, min_length=4, max_length=None):
    """
    Find and score different types of sequences within the given sequence.

    This function finds and scores palindromic sequences, direct repeat sequences,
    and inverted repeat sequences within the provided sequence. It then saves the results to an Excel file.

    Parameters:
        sequence (str): The input sequence in which to search for sequences.
        meme_output (str): The path to the output directory where the Excel file will be saved.
        min_length (int): The minimum length of sequences to search for (default is 4).
        max_length (int): The maximum length of sequences to search for.
            If None, the maximum length is set to half the length of the input sequence.

    Returns:
        None
    """
    
    palindromic_sequences = find_palindromic_sequences(sequence, min_length, max_length)
    direct_repeats = find_direct_repeats(sequence, min_length, max_length)
    inverted_repeats = find_inverted_repeats(sequence, min_length, max_length)
    
    scored_palindromic_sequences = [("palindromic", pos1, frag1, pos2, frag2, score_sequence(frag1, sequence)) for pos1, frag1, pos2, frag2 in palindromic_sequences]
    scored_direct_repeats = [("direct_repeat", pos1, frag1, pos2, frag2, score_sequence(frag1, sequence)) for pos1, frag1, pos2, frag2 in direct_repeats]
    scored_inverted_repeats = [("inverted_repeat", pos1, frag1, pos2, frag2, score_sequence(frag1, sequence)) for pos1, frag1, pos2, frag2 in inverted_repeats]
    
    results = {
        "palindromic_sequences": scored_palindromic_sequences,
        "direct_repeats": scored_direct_repeats,
        "inverted_repeats": scored_inverted_repeats
    }
    
    # Make a excel file with the results
    df = pd.DataFrame()
    for key, value in results.items():
        if value:
            df = pd.concat([df, pd.DataFrame(value, columns=["Type", "Position1", "Sequence1", "Position2", "Sequence2", "Score"])])
    
    df.to_excel(meme_output+"/patterns.xlsx", index=False)
    
    return None

# if __name__ == "__main__":
#     sequence = "AGTACGTAGCCGAGCATGATCATCGAGCCGACGATCGA"
#     patterns = find_and_score_sequences(sequence, meme_output="/home/maria/TFM/Output_Ivan_Examples")
#     print(patterns)
