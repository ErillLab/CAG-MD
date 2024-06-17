import os
import argparse
from tqdm import tqdm
from Bio import Entrez
from variables import *
from logger import logger
from search_blast import search_blast
from final_summary import final_summary
from fetch_sequences import fetch_sequences
from search_promoters import search_promoters
from get_gene_records import get_gene_records
from filter_sequences import filter_sequences
from read_and_validate_config import read_config
from extract_params_config import extract_params
from launch_meme_search import launch_meme_search
from divergence import calculate_motif_divergence
from coverages import query_coverage, total_coverage


def main(root):
    pathConfig = os.path.join(root,"config/config.json")
    
    config = read_config(pathConfig)

    configParams = extract_params(config)
    
    Entrez.email = configParams.EMAIL
    Entrez.api_key = configParams.API_KEY
    run_meme_only = configParams.run_only_meme
    
    output = os.path.join(root, config["output_parameters"]["folder_name"])
    
    if not os.path.isdir(output):
        os.mkdir(output)
    id_sequences, total_sequences, identities_id_dict = {}, {}, {}
    
    if run_meme_only == False:
        meme_output = False
        hits = search_blast(params=configParams)
        
        for id, hits in hits.items():
            hit_sequences, sequences = [], []
            print(f"\tObtaining the gene records for the {id} hits. See the log file for more information.")
            
            for hit in tqdm(hits):
                records = get_gene_records(hit, configParams)
                if records is None:
                    logger.info(f"Skipping hit {hit} due to failed record retrieval.")
                    continue
                intergenic_positions = search_promoters(hit=records, configParams=configParams)
                sequences.append(fetch_sequences(intergenic_positions, configParams))
                
            print("\t\tFiltering the sequences...")
            logger.info("\t\tFiltering the sequences...")
            sequences_filtered, identities_dict = filter_sequences(sequences, configParams)
            
            identities_id_dict[id] = identities_dict
            hit_sequences.append(sequences_filtered)
            id_sequences[id] = hit_sequences
            total_sequences[id] = len(sequences_filtered)
    else:
        meme_output = os.path.join(output, config["output_parameters"]["folder_rerun_meme"])
    
    logger.info("Running MEME...")
    print("Running MEME...")
    fasta_records = launch_meme_search(id_sequences, configParams, output, meme_output)
    logger.info("\tMEME search finished.")
    print("\tMEME search finished.")
    
    logger.info("Reading MEME output and creating output files...")
    print("Reading MEME output and creating output files...")
    motifs_dict = query_coverage(output, total_sequences)
    motifs_coverage = total_coverage(motifs_dict, total_sequences, output)
    motif_divergence = calculate_motif_divergence(motifs_dict, id_sequences, output)
    final_summary(output)
    logger.info("Output files created.")
    print("Output files created.")
    
    return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the main pipeline with the specified root directory.")
    parser.add_argument("root", type=str, help="The root directory for the pipeline.")
    args = parser.parse_args()
    
    main(args.root)