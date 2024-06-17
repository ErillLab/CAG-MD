from logger import logger
from get_genome import fetch_genome_annotations

def search_promoters(hit, configParams):
    """
    Search for promoters in the intergenic regions of a gene hit.

    This function identifies potential promoters in the intergenic regions of a gene hit based on the provided
    genome annotations and configuration parameters.

    Parameters:
        hit (dict): a dictionary containing information about the gene hit, including start and stop positions,
            strand, and other relevant details.
        configParams (Variables): an instance of the Variables class containing configuration parameters,
            such as the minimum and maximum intergenic distances and region sizes.

    Returns:
        list: a list containing dictionaries representing potential promoters found in the intergenic regions
            of the gene hit. Each dictionary contains information about the gene record and the corresponding
            intergenic regions where potential promoters are located.
    """
    
    
    promoters = []
    
    min_intergenic_distance = configParams.min_intergenic_size
    max_intergenic_distance = configParams.max_intergenic_size
    max_sequence_length = configParams.max_sequence_length

    if hit["acc"] is None:
        logger.info("No accession number provided in the gene hit.")
        return promoters
    
    logger.info(f"-Searching for possible promoters in the genome of {hit['acc']}")
    
    annotation, original_start, original_end = fetch_genome_annotations(hit, configParams)
    
    gene_length = abs(int(hit["stop"])-int(hit["start"]))+1
    gene_strand = hit["strand"]

    # Assign the correct gene start
    if gene_strand == -1:
        gene_start = configParams.downstream_size_region
        if int(hit["start"]) < gene_start:
            gene_start = int(hit["start"])-1
    else:
        gene_start = configParams.upstream_size_region
        if int(hit["start"]) < gene_start:
            gene_start = int(hit["start"])-1
    
    gene_end =  gene_start + gene_length
        
    intergenic_regions_list = []
    
    # If strand == 1, then we reversed the features
    if gene_strand == -1:
        annotations = annotation.features
    else:
        annotations = list(reversed(annotation.features))

    annotations_genes = [feature for feature in annotations if feature.type == "gene"]
    i, distance = 0, 0
    
    if gene_strand == -1:
        while i < len(annotations_genes) and annotations_genes[i].location.start <= gene_start:
            i += 1
    else:
        while i < len(annotations_genes) and annotations_genes[i].location.start >= gene_start:
            i += 1
    
    for feature in annotations_genes[i:]:
        next_start = int(feature.location.start)
        next_end = int(feature.location.end)
        
        if gene_strand == -1:
            distance = next_start - gene_end
        else:
            distance = gene_start - next_end
        
        if distance >= max_intergenic_distance or feature.strand != gene_strand:
            break
        
        # Save the start and end of the intergenic region
        if gene_strand == -1:
            start_intergenic_distance = gene_end
            end_intergenic_distance = next_start
        else: 
            start_intergenic_distance = next_end
            end_intergenic_distance = gene_start
        
        if distance >= min_intergenic_distance:
            neighboring_gene_info = {
                "start": next_start,
                "stop": next_end
            }
            intergenic_regions = {
                "start": start_intergenic_distance,
                "stop": end_intergenic_distance,
                "length": distance,
                "strand": gene_strand
            }

            intergenic_regions_list.append((neighboring_gene_info, intergenic_regions))

        # Actualize the gene_start and gene_end
        if gene_strand == -1:
            gene_start = next_start
            gene_end = next_end
        else:
            gene_start = next_start
            gene_end = next_end
    
    if min_intergenic_distance < distance < max_sequence_length:
        max_sequence_length = distance

    # Add putative promoter region after the last gene
    if gene_strand == -1:
        start_intergenic_distance = gene_end
        end_intergenic_distance = gene_end + max_sequence_length
    else:
        end_intergenic_distance = gene_start
        start_intergenic_distance = max(gene_start - max_sequence_length,0)
    
    intergenic_regions = {
        "start": start_intergenic_distance,
        "stop": end_intergenic_distance,
        "length": abs(end_intergenic_distance-start_intergenic_distance),
        "strand": gene_strand
    }
    
    neighboring_gene_info = {
        "start": gene_start,
        "stop": gene_end
    }
    
    intergenic_regions_list.append((neighboring_gene_info, intergenic_regions))
    
    promoter = {
        "gene_record": hit,
        "promoter_info": intergenic_regions_list
    }
    promoters.append(promoter)
    
    logger.info(f"\t-Putative promoters found for {hit['acc']}: {len(promoters)}")
    
    return promoters
