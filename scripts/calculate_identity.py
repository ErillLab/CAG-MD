
def calculate_identity(alignment):
    
    symbol_line = format(alignment[0]).split("\n")[1]
    
    identical_count = symbol_line.count("|")
    
    identity = identical_count/len(symbol_line)
    
    return identity