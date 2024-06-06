from Bio.SeqRecord import SeqRecord
from get_genome import fetch_genome_sequence
from logger import logger

def fetch_sequences(intergenic_positions, configParams):
    """
    Extracts the DNA sequence of the intergenic regions in the operon.

    Parameters:
        intergenic_positions (dict): a dictionary containing information about the intergenic positions.
        configParams (Variables): an instance of the Variables class containing config parameters.

    Returns:
        intergenic_sequences (list): a list of SeqRecord objects containing the DNA sequences of the intergenic regions.
    """
    intergenic_sequences = []
    
    gene_record = intergenic_positions[0]["gene_record"]
    promoter_info = intergenic_positions[0]["promoter_info"]
    
    genome = fetch_genome_sequence(gene_record, configParams)
    i = 1
    for neighboring_gene, intergenic_region in promoter_info:
        start = intergenic_region["start"]
        stop = intergenic_region["stop"]
        
        intergenic_sequence = genome[start:stop]
        
        id_name = f"seq{i}:{gene_record["acc"]}"
        
        intergenic_seqRecord = SeqRecord(intergenic_sequence.seq, id=id_name, description=f"Region between {start} and {stop}.")
        intergenic_sequences.append({"start": start, "stop": stop, "seq": intergenic_seqRecord, "length": intergenic_region["length"]})
    
        i += 1
        
    return intergenic_sequences