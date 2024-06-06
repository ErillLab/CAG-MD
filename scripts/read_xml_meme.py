"""
Reads the MEME XML file and returns the motifs in a dictionary format.
"""

from Bio import motifs
from xlwt import Workbook
import pickle
from divergence import calculate_motif_divergence

class MyInstance():
    def __init__(self, sequence_name, start, end, pvalue, strand, data):
        self.sequence_name = sequence_name
        self.start = start
        self.end = end
        self.pvalue = pvalue
        self.strand = strand
        self.sequence = data

def read_xml_meme(original_path, total_instances):

    path = original_path + "/meme.xml"
    wb = Workbook()
    
    try:
        with open(path) as handle:
            meme_motifs = motifs.parse(handle, "meme")
    except Exception as e:
            # Print the error message
            print("Error reading MEME XML file:", e)
            return None
    
    id_sheets = {}
    for name_id, total in total_instances.items():
        sheet = wb.add_sheet(name_id) 
        
        sheet.write(0, 0, "Motif")
        sheet.write(0, 1, "E-value")
        sheet.write(0, 2, "Length")
        sheet.write(0, 3, "Number of instances")
        sheet.write(0, 4, "Total instances")
        sheet.write(0, 5, "Coverage")

        sheet_instances = wb.add_sheet(name_id + "_instances")

        sheet_instances.write(0, 0, "Sequence name")
        sheet_instances.write(0, 1, "Start")
        sheet_instances.write(0, 2, "End")
        sheet_instances.write(0, 3, "P-value")
        sheet_instances.write(0, 4, "Strand")
        sheet_instances.write(0, 5, "Sequence")
        id_sheets[name_id] = {"sheet": sheet, "sheet_instances": sheet_instances, "row": 1}
        

    motifs_dict = {}
    
    for motif in meme_motifs:
        myinstances = {id_name_query:[] for id_name_query in total_instances.keys()}
        
        for instance in motif.instances:
            sequence_name = instance.sequence_name
            start = instance.start
            end = start + instance.length
            strand = instance.strand
            pvalue = instance.pvalue
            sequence = instance._data
            
            id_name_query = sequence_name.split("$")[0]
            myinstances[id_name_query].append(MyInstance(sequence_name, start, end, pvalue, strand, sequence))
            
        motifs_dict[motif.alt_id] = {"motif": motif.consensus, "instances": myinstances,
                                    "evalue": motif.evalue, "length": motif.length}
    
    c = 1
    for motif, data in motifs_dict.items():
        coverages = {}
        for id_name, instances in data["instances"].items():
            sheet_instances = id_sheets[id_name]["sheet_instances"]
            idx = id_sheets[id_name]["row"]
            for instance in instances:
                sheet_instances.write(idx, 0, instance.sequence_name)
                sheet_instances.write(idx, 1, int(instance.start))
                sheet_instances.write(idx, 2, int(instance.end))
                sheet_instances.write(idx, 3, int(instance.pvalue))
                sheet_instances.write(idx, 4, instance.strand)
                sheet_instances.write(idx, 5, instance.sequence)
                idx += 1
            id_sheets[id_name]["row"] = idx
            total = total_instances[id_name]
            coverages[id_name] = {"coverages": len(instances)/total * 100, "num_instances": len(instances), "total_instances": total}
        
        data["coverage_info"] = coverages
        for id, sheet_dict in id_sheets.items():
            sheet = sheet_dict["sheet"]
            sheet.write(c, 0, str(data["motif"]))
            sheet.write(c, 1, float(data["evalue"]))
            sheet.write(c, 2, float(data["length"]))

        for id_name, data_cov in coverages.items():
            sheet = id_sheets[id_name]["sheet"]
            sheet.write(c, 3, float(data_cov["num_instances"]))
            sheet.write(c, 4, float(data_cov["total_instances"]))
            sheet.write(c, 5, float(data_cov["coverages"]))
        c += 1
    wb.save(original_path + "/motifs_summary.xls")
    return motifs_dict
    

def total_coverage(total_motifs, total_instances, output_path):
    wb = Workbook()
    
    sheet = wb.add_sheet("Total coverage")
    
    #sheet.write(0, 0, "Query")
    sheet.write(1, 0, "Motif")
    #sheet.write(0, 2, "Coverage")
    id_cols = {}
    for col, id in enumerate(total_instances.keys()):
        sheet.write(1, col+1, id)
        id_cols[id] = col+1
    
    total_sequences = 0
    #for motif_dict in total_motifs.values():
    #    total_sequences += motif_dict["total_instances"]
    for id, nsequences in total_instances.items():
        total_sequences += nsequences
    
    motifs_coverage_sums = {}
    fila = 2
    for memeid, data in total_motifs.items():     
        sheet.write(fila, 0, str(data["motif"])) 
        for id_name, coverage_info in data["coverage_info"].items():
            coverage = coverage_info["coverages"]
            id_ratio = coverage_info["total_instances"] / total_sequences
            id_coverage = coverage * id_ratio
            sheet.write(fila, id_cols[id_name], id_coverage)
            
            if data["motif"] in motifs_coverage_sums:
                motifs_coverage_sums[data["motif"]] += id_coverage
            else:
                motifs_coverage_sums[data["motif"]] = id_coverage
        fila += 1
        
    c = 2
    nids = len(total_instances)
    sum_col = nids + 1
    sheet.write(1, sum_col, "Total coverage")
    for _, coverage_sum in motifs_coverage_sums.items():
        sheet.write(c, sum_col, float(coverage_sum))
        c += 1
        
    wb.save(f"{output_path}/total_coverage.xls")
    return motifs_coverage_sums

def from_pickle(path):
    with open(path, 'rb') as file:
        return pickle.load(file)

if __name__ == "__main__":
    paths = ['/home/maria/TFM/Output_Ivan_Examples']
    
    total_instances = {"QJF44568.1": 55, "WP_003238209.1": 38, "NP_388444.1": 65}
    #total_motifs = []
    total_motifs = read_xml_meme(paths[0], total_instances)
    #total_motifs.append(motif)
    
    # Read pickle file in path with the name id_sequences.pickle
    id_sequences = from_pickle(paths[0] + "/id_sequences.pickle")
    
    divergence = calculate_motif_divergence(total_motifs, id_sequences, paths[0])
    
    #motifs_coverage = total_coverage(total_motifs, total_instances, paths[0])
    