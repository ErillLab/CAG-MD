from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from tqdm import tqdm
from logger import logger
import time


def search_blast(params, db="nr", max_hits=50, e_value=10E-10, taxIDs=None, min_coverage=None):
    '''
    Perform a BLAST search for a set of protein records.

    This function conducts a BLAST search for a set of protein records using the specified parameters.
    It retrieves protein sequences from NCBI using Entrez, performs BLAST searches via NCBIWWW,
    and filters the search results based on specified criteria.

    Parameters:
        params (Variables): an instance of the Variables class containing configuration parameters.
        db (str, optional): the name of the BLAST database to be queried. Defaults to 'nr'.
        max_hits (int, optional): the maximum number of hits to return. Defaults to 50.
        e_value (float, optional): the threshold for the E-value of the hits. Defaults to 10E-10.
        taxIDs (list of strings, optional): the taxonomic IDs limits to the BLAST search. Defaults to None.
        min_coverage (float, optional): the minimum coverage for the hits. Defaults to None.

    Returns:
        hits_dict (dict): A dictionary containing the input record ID as the key and a list of hit accession numbers as the value.
    
    Raises:
        Exception: if no protein ID is provided for the BLAST search.
        Exception: if the BLAST search fails after the maximum number of attempts (REQUEST_LIMIT).

    Note:
        This function logs various steps of the BLAST search process using the provided logger.
    '''
    
    REQUEST_LIMIT = params.REQUEST_LIMIT
    SLEEP_TIME = params.SLEEP_TIME
    
    input = params.input_records
    
    db = params.db
    taxIDs = params.tax_IDs
    min_coverage = params.query_coverage
    max_hits = params.max_hits
    e_value = params.e_value
    
    logger.info("Starting BLAST search:")
    
    if len(input) < 1:
        raise Exception("\tAt least one protein ID is needed to perform BLAST. Try again.")
    
    hits_dict = {}
    blast_records = None
    
    for id in tqdm(input):
        logger.info("\tGetting protein record")
        
        for i in range(REQUEST_LIMIT):
            
            try:
                handle = Entrez.efetch("protein", id=id, rettype="fasta", retmode="text")
                time.sleep(SLEEP_TIME)
                break
                
            except:
                logger.warning(f"\t\tNCBI exception raised on attempt {str(i)}. Reattempting.")
                
                if i == (REQUEST_LIMIT - 1):
                    logger.info(f"\t\tCould not download record after {str(REQUEST_LIMIT)} attempts.")
                    raise Exception(f"\t\tCould not download record after {str(REQUEST_LIMIT)} attempts.")
        
        logger.info("\tGetting protein sequence")
        input_seq = (SeqIO.read(handle, "fasta")).seq
        
        logger.info(f"\tPerforming BLAST search: {str(id)}")
        
        if taxIDs:
            taxon_query = "txid" + ",".join(map(str, taxIDs)) + "[ORGN]"
            logger.info(f"\t\tPerforming BLAST search on {",".join(map(str, taxIDs))}")
            for i in range(REQUEST_LIMIT):
                try:
                    result_handle = NCBIWWW.qblast(program="blastp",
                                                   database=db, sequence=input_seq,
                                                   entrez_query=taxon_query, expect=e_value,
                                                   hitlist_size=max_hits)
                    logger.info("\tGetting records.")
                    blast_records = list(NCBIXML.parse(result_handle))
                    time.sleep(SLEEP_TIME)
                    break
                    
                except:
                    logger.warning(f"\t\tNCBI exception raised on attempt {str(i)}. Reattempting.")

                    if i == (REQUEST_LIMIT - 1):
                        logger.error(f"\t\tCould not download record after {str(REQUEST_LIMIT)} attempts.")
        else:
            for i in range(REQUEST_LIMIT):
            
                try:
                    result_handle = NCBIWWW.qblast("blastp", db, input_seq, expect=e_value, hitlist_size=max_hits)
                    logger.info("\tGetting records.")
                    blast_records = list(NCBIXML.parse(result_handle))
                    time.sleep(SLEEP_TIME)
                    break
                    
                except:
                    logger.warning(f"\t\tNCBI exception raised on attempt {str(i)}. Reattempting.")
                
                if i == (REQUEST_LIMIT - 1):
                    logger.error(f"\t\tCould not download record after {str(REQUEST_LIMIT)} attempts.")
        
        hits_dict[id] = []
        
        for record in blast_records[0].alignments:
            
            current_hit = record.hit_id.split('|')[1]
            logger.info(f"\t\tAnalyzing hit {str(current_hit)}")
            
            for hit in record.hsps:
                if min_coverage:
                    cov = (hit.query_end - hit.query_start + 1) / len(input_seq)
                    if cov >= min_coverage:
                        if current_hit not in hits_dict[id]:
                            logger.info(f"\t\tAdding hit (coverage = {str(round(cov*100, 2))}): {str(current_hit)}.")
                            hits_dict[id].append(current_hit)
                    else:
                        logger.info(f"\t\tThe hit {str(current_hit)} has a coverage below the established threshold (coverage = {str(round(cov*100, 2))}).")
                else:
                    if current_hit not in hits_dict[id]:
                        logger.info(f"There is no minimum coverage. Adding hit: {str(current_hit)}")
                        hits_dict[id].append(current_hit)
                
    logger.info(f"\tReturning hits for {str(len(hits_dict))} records.")
    return hits_dict
