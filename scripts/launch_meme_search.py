import subprocess
import tempfile
from Bio import SeqIO
import os
from create_fasta import create_fasta

def launch_meme_search(sequences_dict, configParams, output, meme_output = False):
    """
    Launch MEME motif discovery search.

    This function generates a FASTA file from the provided sequences, launches a MEME search
    using the specified parameters, and saves the results to the specified output directory.

    Parameters:
        sequences_dict (dict): a dictionary containing the sequences to search for motifs.
            The keys represent sequence identifiers, and the values are BioPython SeqRecord objects.
        configParams (Variables): an instance of the Variables class containing MEME search parameters.
            The MEME search parameters include mod, nmotifs, minw, maxw, revcomp, and pal.
        output (str): the path to the output directory where MEME results will be saved.
        meme_output (bool or str): if False, the FASTA file with the sequences is created.
            If not False, it indicates the output folder for the new MEME execution.

    Returns:
        list: a list of SeqRecord objects containing the input sequences in FASTA format.

    Note:
        The function generates a temporary FASTA file from the input sequences, launches the MEME search
        using the specified parameters, and saves the results to the specified output directory.
        If an error occurs during the process, appropriate error messages are printed.
    """

    fasta_records = None
    try:
        output_file_path = os.path.join(output,"meme_sequences_filtered.fasta")
        
        if not meme_output: 
            # Generate FASTA content      
            fasta_records = create_fasta(sequences_dict)
            
            with open(output_file_path, mode="w") as output_fasta:
                SeqIO.write(fasta_records, output_fasta, "fasta")
        else:
            output = meme_output
        
        # Get MEME parameters
        mod = str(configParams.meme_mod)
        nmotifs = str(configParams.meme_nmotifs)
        minw = str(configParams.meme_minw)
        maxw = str(configParams.meme_maxw)
        options = ["-dna", "-mod", mod, "-nmotifs", nmotifs, "-minw", minw, "-maxw", maxw]
        
        if configParams.meme_revcomp:
            options.append("-revcomp")
        
        if configParams.meme_pal:
            options.append("-pal")
        
        # Launch MEME search
        command = ["meme", output_file_path, "-oc", output] + options
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        stdout, stderr = process.communicate()
        
        # Check for errors
        if stderr:
            print("Error while running MEME:", stderr)
    except Exception as e:
        print("An error occurred:", e)
    
    return fasta_records