from Bio import Align
from calculate_identity import calculate_identity

def filter_sequences(sequences, configParams):
    """
    Filter sequences based on pairwise sequence identity.

    This function performs pairwise sequence alignment between the input sequences
    and filters out those that have a sequence identity greater than or equal to
    the specified maximum identity threshold.

    Parameters:
        sequences (list): A list of dictionaries containing sequence information.
            Each dictionary should contain a key "seq" with a BioPython SeqRecord object
            representing the sequence.
        configParams (Variables): An instance of the Variables class containing configuration parameters.
            The configuration parameters should include the maximum identity threshold.

    Returns:
        tuple: A tuple containing two elements:
            - A list of dictionaries containing the filtered sequences that meet the identity threshold.
            - A dictionary containing pairwise sequence identities for the retained sequences.
              The keys are sequence IDs, and the values are lists of tuples.
              Each tuple contains the ID of the aligned sequence, the identity score, and the aligned sequence itself.

    Note:
        The filtering process is performed by pairwise sequence alignment using the Needleman-Wunsch algorithm.
        Sequences with an identity greater than or equal to the maximum identity threshold are filtered out.
        The output retains only sequences that meet the identity threshold, along with their pairwise sequence identities.
    """
    
    max_identity = configParams["max_identity"]
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    idx_to_remove = set()
    identities_dict, identities, sequences_to_keep = {}, [], []
    
    sequences_list = [d for subseq in sequences for d in subseq]

    for idx, sequence in enumerate(sequences_list):
        if idx not in idx_to_remove:
            seq1 = str(sequence["seq"].seq)
            idx2 = idx
            for _,sequence2 in enumerate(sequences_list[(idx + 1):]):
                idx2 = idx2 + 1
                seq2 = str(sequence2["seq"].seq)
                if idx2 not in idx_to_remove:
                    if seq1 != seq2:
                        if len(seq1) == 0:
                            idx_to_remove.add(idx)
                            continue
                        elif len(seq2) == 0:
                            idx_to_remove.add(idx2)
                            continue
                        else:
                            alignments = aligner.align(seq1, seq2)
                            identity = calculate_identity(alignments)

                            if identity >= max_identity:
                                idx_to_remove.add(idx2)
                            else:
                                identities.append((sequence["seq"].id, sequence2["seq"].id, identity))
            identities_dict[sequence["seq"].id] = identities
        else:
            continue

    for idx, seq in enumerate(sequences_list):
        if idx not in idx_to_remove:
            sequences_to_keep.append(seq)
    
    return sequences_to_keep, identities_dict
