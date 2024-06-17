from read_xml import read_xml
import csv
import os

class MyInstance:
    def __init__(self, sequence_name, start, end, pvalue, strand, data):
        self.sequence_name = sequence_name
        self.start = start
        self.end = end
        self.pvalue = pvalue
        self.strand = strand
        self.sequence = data

def query_coverage(original_path, total_instances):
    """
    Calculate motif coverage for each query sequence and generate summary CSV files.

    This function reads motifs from a MEME XML file, calculates the coverage of motifs across query sequences,
    and generates summary and instance CSV files. The coverage is defined as the ratio of sequences containing
    at least one instance of the motif to the total number of sequences.

    Parameters:
        original_path (str): the path to the directory containing the MEME XML file.
        total_instances (dict): a dictionary where keys are query sequence IDs and values are the total number
                                of sequences for each query.

    Returns:
        dict: a dictionary containing motif data, including instances and coverage information.
    """

    path = os.path.join(original_path, "meme.xml")

    meme_motifs = read_xml(path)

    motifs_dict = {}
    
    for motif in meme_motifs:
        myinstances = {id_name_query: [] for id_name_query in total_instances.keys()}
        seen_sequences = {id_name_query: set() for id_name_query in total_instances.keys()}
        
        for instance in motif.instances:
            sequence_name = instance.sequence_name
            start = instance.start
            end = start + instance.length
            strand = instance.strand
            pvalue = instance.pvalue
            sequence = instance._data
            
            id_name_query = sequence_name.split("$")[0]
            
            # Ensure only one instance per sequence
            if sequence_name not in seen_sequences[id_name_query]:
                seen_sequences[id_name_query].add(sequence_name)
                myinstances[id_name_query].append(MyInstance(sequence_name, start, end, pvalue, strand, sequence))
            
        motifs_dict[motif.alt_id] = {
            "motif": motif.consensus, 
            "instances": myinstances,
            "evalue": motif.evalue, 
            "length": motif.length
        }

    # Prepare data for the summary CSV
    summary_rows = []
    instance_rows = []

    for motif, data in motifs_dict.items():
        coverages = {}
        for id_name, instances in data["instances"].items():
            total = total_instances[id_name]
            coverage = len(instances) / total
            coverages[id_name] = {
                "coverages": coverage,
                "num_instances": len(instances),
                "total_instances": total
            }
            summary_row = {
                "Query": id_name,
                "Motif": str(data["motif"]),
                "Coverage": coverage,
                "Number of motif instances": len(instances),
                "Total query sequences": total,
                "Motif length": data["length"]
            }
            summary_rows.append(summary_row)

            for instance in instances:
                instance_row = {
                    "Sequence name": instance.sequence_name,
                    "Start": instance.start,
                    "End": instance.end,
                    "P-value": instance.pvalue,
                    "Strand": instance.strand,
                    "Motif instance": instance.sequence
                }
                instance_rows.append(instance_row)
        
        data["coverage_info"] = coverages
    
    # Write the summary CSV
    summary_csv_path = os.path.join(original_path, "motifs_summary.csv")
    with open(summary_csv_path, mode="w", newline="") as summary_csv_file:
        fieldnames = ["Query", "Motif", "Coverage", "Number of motif instances", "Total query sequences", "Motif length"]
        writer = csv.DictWriter(summary_csv_file, fieldnames=fieldnames)
        
        writer.writeheader()
        writer.writerows(summary_rows)
    
    # Write the instances CSV
    instances_csv_path = os.path.join(original_path, "motif_instances.csv")
    with open(instances_csv_path, mode="w", newline="") as instances_csv_file:
        fieldnames = ["Sequence name", "Start", "End", "P-value", "Strand", "Motif instance"]
        writer = csv.DictWriter(instances_csv_file, fieldnames=fieldnames)
        
        writer.writeheader()
        writer.writerows(instance_rows)

    return motifs_dict

def total_coverage(total_motifs, total_instances, output_path):
    """
    Calculate and summarize total motif coverage across all query sequences.

    This function calculates the total coverage of motifs across all query sequences, taking into account
    the ratio of instances in each query. It generates a CSV file with the summarized coverage information.

    Parameters:
        total_motifs (dict): a dictionary containing motif data and coverage information for each query sequence.
        total_instances (dict): a dictionary where keys are query sequence IDs and values are the total number
                                of sequences for each query.
        output_path (str): the path to the directory where the total coverage CSV file will be saved.

    Returns:
        dict: a dictionary containing the summed coverage for each motif.
    """
    
    total_sequences = sum(total_instances.values())
    
    motifs_coverage_sums = {}
    rows = []
    
    # First pass to calculate the total coverage for each motif
    for memeid, data in total_motifs.items():
        motif = str(data["motif"])
        if motif not in motifs_coverage_sums:
            motifs_coverage_sums[motif] = 0
            
        for id_name, coverage_info in data["coverage_info"].items():
            coverage = coverage_info["coverages"]
            id_ratio = coverage_info["total_instances"] / total_sequences
            id_coverage = coverage * id_ratio
            
            motifs_coverage_sums[motif] += id_coverage
    
    # Second pass to create the rows for the CSV
    for memeid, data in total_motifs.items():
        motif = str(data["motif"])
        for id_name, coverage_info in data["coverage_info"].items():
            coverage = coverage_info["coverages"]
            id_ratio = coverage_info["total_instances"] / total_sequences
            id_coverage = coverage * id_ratio
            
            row = {
                "Query": id_name,
                "Motif": motif,
                "Total coverage (sum of querys)": motifs_coverage_sums[motif],
                "Motif coverage": id_coverage,
                "Query coverage": coverage,
                "Total motif instances": coverage_info["total_instances"],
                "Total MEME sequences": total_sequences,
                "Query ratio": id_ratio,
            }
            rows.append(row)
    
    total_coverage_path = os.path.join(output_path, "total_coverage.csv")
    with open(total_coverage_path, mode="w", newline="") as coverage_file:
        fieldnames = ["Query", "Motif", "Total coverage (sum of querys)",
                      "Motif coverage", "Query coverage", "Total motif instances",
                      "Total MEME sequences", "Query ratio"]
        writer = csv.DictWriter(coverage_file, fieldnames=fieldnames)
        
        writer.writeheader()
        writer.writerows(rows)
    
    return motifs_coverage_sums