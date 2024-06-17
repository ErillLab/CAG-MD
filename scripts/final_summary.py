import pandas as pd
import csv
import os

def final_summary(path):
    """
    Combine motif summary, total coverage, and divergence statistics into a single CSV file.

    This function reads the motif summary, total coverage, and motif divergence statistics from their respective
    CSV files and combines them into a single CSV file with relevant columns.

    Parameters:
        path (str): the path with the output files.

    Returns:
        None
    """
    summary_path = os.path.join(path, "motifs_summary.csv")
    total_coverage_path = os.path.join(path, "total_coverage.csv")
    divergence_stats_path = os.path.join(path, "motif_divergence_stats.csv")
    
    # Read the CSV files
    summary_df = pd.read_csv(summary_path)
    total_coverage_df = pd.read_csv(total_coverage_path)
    divergence_stats_df = pd.read_csv(divergence_stats_path)
    
    # Select relevant columns from each DataFrame
    summary_df = summary_df[["Query", "Motif", "Coverage"]]
    total_coverage_df = total_coverage_df[["Query", "Motif", "Total coverage (sum of querys)"]]
    divergence_stats_df = divergence_stats_df[["Motif", "Query", "Mean distance"]]
    
    # Merge the DataFrames on "Query" and "Motif"
    combined_df = summary_df.merge(total_coverage_df, on=["Query", "Motif"], how="left")
    combined_df = combined_df.merge(divergence_stats_df, on=["Query", "Motif"], how="left")
    
    # Write the combined DataFrame to a new CSV file
    combined_csv_path = os.path.join(path, "final_summary.csv")
    combined_df.to_csv(combined_csv_path, index=False)