from read_and_validate_config import read_config
from extract_params_config import extract_params
from search_blast import search_blast
from get_gene_records import get_gene_records
from fetch_sequences import fetch_sequences
from filter_sequences import filter_sequences
from launch_meme_search import launch_meme_search
from read_xml_meme import read_xml_meme, total_coverage
from find_pattern import find_and_score_sequences
from Bio import Entrez
from divergence import calculate_motif_divergence
from variables import *
from search_promoters import search_promoters
from logger import logger
import os
import json
import pickle

def main():
    root = "/home/maria/TFM"
    pathConfig = root+"/config/config.json"
    
    config = read_config(pathConfig)

    configParams = extract_params(config)
    
    Entrez.email = configParams.EMAIL
    Entrez.api_key = configParams.API_KEY
    
    meme_output = root+"/"+config["output_parameters"]["folder_name"]
    
    if not os.path.isdir(meme_output):
        os.mkdir(meme_output)
    
    # Commented to avoid running the BLAST search bc it takes a lot of time
    # hits = search_blast(params=configParams)
    # Read the JSON file named hits.json
    with open(root+"/hits.json", "r") as file:
        hits = json.load(file)
    
    id_sequences, total_sequences, identities_id_dict = {}, {}, {}
    all_sequences = []
    id_sequences_pickle_path = meme_output + "/id_sequences.pickle"
    usar_cache = True
    # Provisional code to avoid this part of the code
    if not os.path.exists(id_sequences_pickle_path) or not usar_cache:
        for id,hits in hits.items():
            hit_sequences, sequences = [], []
            for hit in hits:
                records = get_gene_records(hit, configParams)
                intergenic_positions = search_promoters(hit = records, params = configParams)
                sequences.append(fetch_sequences(intergenic_positions, configParams))
            sequences_filtered, identities_dict = filter_sequences(sequences, configParams)
            identities_id_dict[id] = identities_dict
            hit_sequences.append(sequences_filtered)
            id_sequences[id] = hit_sequences
            total_sequences[id] = len(sequences_filtered)
            all_sequences.append(sequences_filtered)
        
        to_pickle(id_sequences, id_sequences_pickle_path)
    
    else:
        id_sequences = from_pickle(id_sequences_pickle_path)
        for id, sequences in id_sequences.items():
            total_sequences[id] = len(sequences[0])
            
    logger.info("- Launching MEME search...")
    fasta_records = launch_meme_search(id_sequences, configParams, meme_output)
    logger.info("\t- MEME search finished.")
    logger.info("- Reading MEME output and creating output files...")
    motifs_dict = read_xml_meme(meme_output, total_sequences)
    motifs_coverage = total_coverage(motifs_dict, total_sequences, meme_output)
    motif_divergence = calculate_motif_divergence(motifs_dict, id_sequences, meme_output)
    
    for motif, sequences in motifs_dict.items():
        sequence = str(sequences["motif"])
        find_and_score_sequences(str(sequence), meme_output)
    logger.info("\t- Output files created.")
    
    pass

def to_pickle(obj, path):
    with open(path, 'wb') as file:
        pickle.dump(obj, file)

def from_pickle(path):
    with open(path, 'rb') as file:
        return pickle.load(file)

if __name__ == "__main__":
    main()
