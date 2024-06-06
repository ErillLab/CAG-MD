
def create_fasta(sequences):
    """
    Function that receives a dictionary with the sequences of the intergenic regions and creates a FASTA file with the sequences.
    
    Parameters:
        sequences (dict): a dictionary containing the sequences of the intergenic regions.
        
    Returns:
        records (list): a list of SeqRecord objects containing the sequences.
    """
    
    records = []
    for id, seqs in sequences.items():
        for seq in seqs[0]:
            seq["seq"].id = id+"$"+seq["seq"].id
            records.append(seq["seq"])

    return records