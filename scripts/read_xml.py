from Bio import motifs

def read_xml(path):
    """
    Reads a MEME XML file and parses it to extract motifs.

    Parameters:
    path (str): the file path to the MEME XML file.

    Returns:
    meme_motifs (Bio.motifs.Motif): a Motif object containing the parsed motifs from the MEME XML file. Returns None if an error occurs during reading.
    """
    
    try:
        with open(path) as handle:
            meme_motifs = motifs.parse(handle, "meme")
    except Exception as e:
        print("Error reading MEME XML file:", e)
        return None
    
    return meme_motifs