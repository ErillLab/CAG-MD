import pandas as pd
import numpy as np
from Bio import Align
from calculate_identity import calculate_identity

def calculate_motif_divergence(motifs_dict, id_sequences, path_output):
    """
    Calculate motif divergence statistics based on pairwise sequence alignment.

    This function calculates the motif divergence statistics for each motif based on the pairwise alignment
    of sequences associated with each motif instance. It calculates the mean and standard deviation of distances
    between sequences, as well as unique sequence lengths.

    Parameters:
        motifs_dict (dict): a dictionary containing motif data, including instances and sequence names.
        id_sequences (dict): a dictionary containing sequences associated with each motif instance.
        path_output (str): the path to the output directory where the CSV file will be saved.

    Returns:
        dict: a dictionary containing motif divergence statistics for each motif.
            Each motif key is associated with a dictionary containing statistics for each query sequence.
    """
    
    motif_divergence_stats = {}
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    
    # Get all query ids
    all_query_ids = {id_query_sequence for motif_data in motifs_dict.values() for id_query_sequence in motif_data["instances"].keys()}
    
    # Create a dictionary with the sequences for each motif
    sequences_dict = {}
    for motif_id, motif_data in motifs_dict.items():
        prov_dict = {}
        motif = motif_data["motif"]
        
        for id_query, instances in motif_data["instances"].items():
            if not instances:
                continue
            id_promotores = id_sequences[id_query][0]
            for instance in instances:
                id_query_sequence = instance.sequence_name.split("$")[0]
                for seq in id_promotores:
                    if id_query_sequence in seq["seq"].id:
                        seq_name = instance.sequence_name
                    else:
                        seq_name = instance.sequence_name.split("$")[-1]
                    if seq["seq"].id == seq_name:
                        if id_query_sequence not in prov_dict:
                            prov_dict[id_query_sequence] = [(seq["seq"], len(seq["seq"]))]
                        else:
                            prov_dict[id_query_sequence].append((seq["seq"], len(seq["seq"])))
        sequences_dict[str(motif)] = prov_dict
    
    # Calculate the mean and std distance between sequences for each motif
    for motif, sequences in sequences_dict.items():
        motif_stats = {}
        for id_query_sequence in all_query_ids:
            if id_query_sequence in sequences:
                seq_list = sequences[id_query_sequence]
                if len(seq_list) < 2:
                    motif_stats[id_query_sequence] = {
                        "mean_distance": float("nan"),
                        "std_distance": float("nan"),
                        "sequence_lengths": [seq[1] for seq in seq_list]
                    }
                else:
                    distances = []
                    for i in range(len(seq_list)):
                        for j in range(i + 1, len(seq_list)):
                            alignment = aligner.align(seq_list[i][0].seq, seq_list[j][0].seq)
                            distance = 1 - calculate_identity(alignment)
                            distances.append(distance)
                    
                    if distances:
                        mean_distance = np.mean(distances)
                        std_distance = np.std(distances)
                        motif_stats[id_query_sequence] = {
                            "mean_distance": mean_distance,
                            "std_distance": std_distance,
                            "sequence_lengths": [seq[1] for seq in seq_list]
                        }
                    else:
                        motif_stats[id_query_sequence] = {
                            "mean_distance": float("nan"),
                            "std_distance": float("nan"),
                            "sequence_lengths": [seq[1] for seq in seq_list]
                        }
            else:
                motif_stats[id_query_sequence] = {
                    "mean_distance": float("nan"),
                    "std_distance": float("nan"),
                    "sequence_lengths": []
                }

        motif_divergence_stats[motif] = motif_stats
    
    save_results_to_csv(motif_divergence_stats, path_output + "/motif_divergence_stats.csv")
    
    return motif_divergence_stats

def save_results_to_csv(motif_divergence_stats, output_filename):
    """
    Save motif divergence statistics to a CSV file.

    This function saves the motif divergence statistics, including mean distance, standard deviation of distance,
    and unique sequence lengths, to a CSV file. The statistics are organized by motif and query sequence.

    Parameters:
        motif_divergence_stats (dict): A dictionary containing motif divergence statistics.
            Each motif key is associated with a dictionary containing statistics for each query sequence.
        output_filename (str): The filename (including path) for the output CSV file.

    Returns:
        None
    """
    
    rows = []
    for motif, stats in motif_divergence_stats.items():
        for query, data in stats.items():
            # Save unique sequence lengths as a single cell
            unique_lengths = ",".join(map(str, set(data.get("sequence_lengths", []))))
            rows.append({
                "Motif": motif,
                "Query": query,
                "Mean distance": data["mean_distance"],
                "Std distance": data["std_distance"],
                "Sequence Lengths": unique_lengths
            })
    
    df = pd.DataFrame(rows)
    df.to_csv(output_filename, index=False)

