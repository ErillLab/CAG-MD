# CAG-MD
Comparative genomics-assisted motif discovery

# Configuration file guidelines

This document provides detailed information on the configuration parameters used in the project. The configuration file follows a schema defined by Cerberus, ensuring that all inputs meet the specified criteria. Below are the parameters and their validation rules.

## 1. Entrez Parameters

The `entrez_parameters` section configures parameters related to Entrez queries.

- **request_limit** (integer, required): the maximum number of requests allowed. Must be at least 1.
- **sleep_time** (number, required): the time (in seconds) to wait between requests. Must be a non-negative number.
- **email** (string, required): a valid email address used for Entrez queries. Must match the regex pattern for email addresses.
- **api_key** (string, optional): the API key for Entrez. Can be `null`.

## 2. BLAST Parameters

The `blast_parameters` section configures parameters related to BLAST searches.

- **db** (string, optional): the database to search against. Can be `null`.
- **tax_IDs** (list, optional): a list of taxonomy IDs to include in the search. Can be `null`.
- **e_value** (number, optional): the E-value threshold for the BLAST search. Can be `null`.
- **query_coverage** (number, optional): the query coverage percentage. Can be `null`.
- **max_hits** (integer, optional): the maximum number of hits to return. Must be at least 1 if specified. Can be `null`.

## 3. Sequences Parameters

The `sequences_parameters` section configures parameters related to sequence processing.

- **max_intergenic_size** (number, required): the maximum intergenic size. Cannot be `null`.
- **min_intergenic_size** (number, required): the minimum intergenic size. Cannot be `null`.
- **upstream_size_region** (number, required): the size of the upstream region. Cannot be `null`.
- **downstream_size_region** (number, required): the size of the downstream region. Cannot be `null`.
- **max_sequence_length** (integer, optional): the maximum sequence length. Can be `null`.

## 4. Maximum Identity

- **max_identity** (float, required): the maximum identity threshold. Must be between 0.0 and 1.0 inclusive. Cannot be `null`.

## 5. MEME Parameters

The `meme_parameters` section configures parameters for MEME motif discovery.

- **mod** (string, required): the model to use. Cannot be `null`.
- **nmotifs** (integer, required): the number of motifs to find. Cannot be `null`.
- **minw** (integer, required): the minimum motif width. Must be at least 0. Cannot be `null`.
- **maxw** (integer, required): the maximum motif width. Cannot be `null`.
- **revcomp** (boolean, required): whether to search both strands. Cannot be `null`.
- **pal** (boolean, required): whether to search for palindromic motifs. Cannot be `null`.

## 6. Run Only MEME

- **run_only_meme** (boolean, required): flag to indicate if only MEME should be run. Cannot be `null`.

## 7. Output Parameters

The `output_parameters` section configures parameters related to output directories.

- **folder_name** (string, required): the name of the output folder. Cannot be `null`.
- **folder_rerun_meme** (string, optional): the folder name for rerunning MEME. Can be `null`.

## 8. Input Records

- **input_records** (list, required): a list of input records. Must contain at least one record. Cannot be `null`.

# Execution Instructions

## Step-by-Step Guide

### 1. Open the Terminal

Ensure you have access to the command line interface on your operating system. This could be Terminal on macOS or Linux, or Command Prompt/PowerShell on Windows.

### 2. Navigate to the Project Directory

Change your current directory to the root directory of the project where the main script is located. You can use the `cd` command followed by the path to the project directory. For example:

```sh
cd path/to/your/project
```

### 3. Execute the Main Script

To run the main pipeline, use the Python interpreter followed by the script name and the required root directory argument. The root directory is the main directory where the pipeline will operate. The command structure is as follows:

```sh
python main.py <root_directory>
```

#### Example

If your root directory is located at `/data/project_root`, the command would be:

```sh
python main.py /data/project_root
```
