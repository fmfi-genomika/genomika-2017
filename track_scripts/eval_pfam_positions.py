import sys
import logging
from pprint import pprint, pformat
import copy

def get_gff3_data(gff3_filename, datatype):
    PROTEIN_NAME_COL = 0
    DATATYPE_COL = 1
    BEGIN = 3
    END = 4
    INFO = 8

    data = []

    f = open(gff3_filename, "r")
    for line in f:
        if line == "##FASTA\n":
            logging.debug("rest of gff3 file is just fasta")
            break
        if line[0] == "#":
            # comment
            continue
        line = line.strip()
        #logging.debug(pformat(line))
        row = line.split("\t")

        if row[DATATYPE_COL] != datatype: continue

        #logging.debug(pformat(row))
        try:
            dato = {
                "protein_name": row[PROTEIN_NAME_COL],
                "begin": int(row[BEGIN]),
                "end": int(row[END]),
                "info": row[INFO]
            }
        except Exception as e:
            logging.warning("failed to read dato: {}".format(e))
            continue

        # extract pfam dat from info
        if datatype == "Pfam":
            res = [s[5:] for s in dato["info"].split(";") if s[0:5] == "Name="]
            if (len(res) == 1):
                dato["pfam_id"] = res[0]
            else:
                logging.debug(pformat(dato))
                logging.warning("dato without Pfam ID in the info!")
                continue
        
        #logging.debug(pformat(dato))
        data.append(dato)
    f.close()
    return data

def get_protein_data(gene_positions_filename, exclude_non3=True):
    PROTEIN_NAME_COL = 1
    CHROMOSOME_NAME = 2
    STRAND = 3
    EXON_STARTS = 9
    EXON_ENDS = 10
    f = open(gene_positions_filename, "r")
    data = []
    nondivisible_proteins = 0
    for line in f:
        raw_dato = line.split("\t")
        dato = {
            "protein_name": raw_dato[PROTEIN_NAME_COL],
            "chr": raw_dato[CHROMOSOME_NAME],
            "strand": raw_dato[STRAND],
            "exons": []
        }
        for b, e in zip(raw_dato[EXON_STARTS].split(",")[:-1], raw_dato[EXON_ENDS].split(",")[:-1]):
            dato["exons"].append((int(b), int(e)))
        #logging.debug(pformat(dato))

        if raw_dato[STRAND] not in ("+", "-"):
            logging.warning("{} strand is corrupted: {}, excluding the dato".format(dato["protein_name"], raw_dato[STRAND]))
            continue

        # checking for divisibilty by 3
        total = sum([e-b for b,e in dato["exons"]])
        if total % 3 != 0:
            if nondivisible_proteins == 5:
                logging.warning("too many proteins with wrong amount of bases, further warnings are not shown")
            elif nondivisible_proteins > 5:
                pass
            else:
                logging.warning("total amount of bases in exons of gene {} is not divisible by 3: {}".format(dato["protein_name"], total))
                logging.warning(pformat(dato))
            nondivisible_proteins += 1
        if total % 3 == 0 or not exclude_non3:
            data.append(dato)

    f.close()
    if nondivisible_proteins > 0:
        logging.warning("AMOUNT OF 3-NONDIVISIBLE PROTEINS IS {}/{}".format(nondivisible_proteins, len(data)))
        if exclude_non3:
            logging.warning("NONDIVISIBLE PROTEINS ARE EXCLUDED (you can change this by setting parameter 'exclude_non3' to False)")
    return data

def get_coordinates_of_protein_part(segment_begin, segment_end, _exons, strand):
    exons = copy.deepcopy(_exons)
    if any([b < 0 or e < 0 for b,e in exons]):
        logging.warning("exons are negative?!")
        logging.warning(pformat(exons))

    #logging.debug("segment_begin: {}, segment_end: {}".format(segment_begin, segment_end))
    def transform_to_closed_intervals(exons):
        # because:
        # https://genome.ucsc.edu/FAQ/FAQtracks.html#tracks1 
        # https://www.biostars.org/p/84686/
        ans = []
        for b, e in exons:
            ans.append((b+1, e))
        return ans
    def transform_to_halfopen_intervals(exons):
        ans = []
        for b, e in exons:
            ans.append((b-1, e))
        return ans

    exons = transform_to_closed_intervals(exons)    
    
    #logging.debug("exons: " + pformat(exons))
    
    if strand == "-":
        #logging.debug("strand is negative")
        exons = [(-e, -b) for b,e in exons][::-1]
        #logging.debug("exons: " + pformat(exons))

    exon_status = []
    total_length = 0
    for b, e in exons:
        exon_length = e - b + 1
        #logging.debug("exon length {} (codons: {})".format(exon_length, exon_length / 3))
        # b - begin; e - end
        # c* - codon number (starting with 1)
        # r* - length of first/last codon in this exon  
        t = {"cb": 0, "ce": 0, "rb": 0, "re": 0}
        t["cb"] = 1 + total_length // 3
        t["ce"] = (total_length + exon_length) // 3 + (1 if (total_length + exon_length) % 3 != 0 else 0)
        t["rb"] = 3 - (total_length % 3)
        t["re"] = (total_length + exon_length) % 3
        if t["re"] == 0: t["re"] = 3
        #logging.debug(pformat(t))
        total_length += exon_length
        exon_status.append(t) 
    if total_length % 3 != 0:
        logging.warning("total_length {} is not divisible by 3".format(total_length))
    #logging.debug(pformat(exon_status))

    start_coordinates = None
    start_exon = None 
    end_coordinates = None
    end_exon = None 
    breakpoints = []
    for i in range(len(exon_status)):
        if exon_status[i]["ce"] < segment_begin: continue
        if exon_status[i]["cb"] > segment_end: break

        if (exon_status[i]["cb"] < segment_begin or (exon_status[i]["cb"] == segment_begin and exon_status[i]["rb"] == 3)) and segment_begin <= exon_status[i]["ce"]:
            #we have beginning
            start_exon = i
            start_coordinates = (segment_begin - exon_status[i]["cb"] - 1) * 3 + exon_status[i]["rb"]  +  exons[i][0]

        if exon_status[i]["cb"] <= segment_end and (segment_end < exon_status[i]["ce"] or (segment_end == exon_status[i]["ce"] and exon_status[i]["re"] == 3)):
            ## we have end
            end_exon = i
            end_coordinates = - (exon_status[i]["ce"] - segment_end - 1) * 3 - exon_status[i]["re"] + exons[i][1]

    if start_exon is None or end_exon is None:
        logging.warning("Ends of segment not found: {} {}".format(start_exon, end_exon))
        assert(false)

    # insert breakpoints
    for i in range(start_exon, end_exon):
        breakpoints.append(exons[i][1])
        breakpoints.append(exons[i+1][0])
    breakpoints = [start_coordinates] + breakpoints + [end_coordinates]
    breakpoints = [(breakpoints[2 * i], breakpoints[2 * i + 1]) for i in range(len(breakpoints)//2)]
    #logging.debug("breakpoints: {}".format(pformat(breakpoints)))
    
    if strand == "-":
        breakpoints = [(-b,-e) for e,b in breakpoints][::-1]
    #logging.debug(breakpoints)

    breakpoints = transform_to_halfopen_intervals(breakpoints)
    
    if any([x < 0 or y < 0 for x, y in breakpoints]):
        logging.warning("Some breakpoints are negative!!")
        logging.warning(pformat(breakpoints))
        logging.debug("{} {}".format(segment_begin, segment_end))
        logging.debug("start_exon: {}, start_coordinates: {}, end_exon: {}, end_coordinates: {}".format(start_exon, start_coordinates, end_exon, end_coordinates))
        logging.debug(exons)
        logging.debug(exon_status)

    return breakpoints


def get_positions(gff3_filename, gene_positions_filename, datatype):
    gff3_data = get_gff3_data(gff3_filename, datatype)
    protein_data = get_protein_data(gene_positions_filename)

    dict_protein_data = {dato["protein_name"]:dato for dato in protein_data}

    for d in gff3_data:
        prot_name = d["protein_name"]
        #logging.debug(prot_name)
        if not prot_name in dict_protein_data:
            logging.warning("no such protein in the file: {}".format(prot_name))
            continue
        resulting_positions = get_coordinates_of_protein_part(d["begin"], d["end"], 
                    dict_protein_data[prot_name]["exons"], dict_protein_data[prot_name]["strand"])
        d["coordinates"] = resulting_positions
        d["chromosome"] = dict_protein_data[prot_name]["chr"]
        d["strand"] = dict_protein_data[prot_name]["strand"]
        #logging.debug(pformat(d))
    return gff3_data

def print2bed(data):
    DELIMITER = "\t"
    for d in data:
        line = []
        # chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
        line.append(d["chromosome"])
        # chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
        line.append(d["coordinates"][0][0])
        # chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
        line.append(d["coordinates"][-1][1])
        # The 9 additional optional BED fields are:
        # name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
        line.append(d["pfam_id"])
        # score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray:
        line.append("0")
        # strand - Defines the strand. Either "." (=no strand) or "+" or "-".
        line.append(d["strand"])
        # thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, thickStart and thickEnd are usually set to the chromStart position.
        line.append(d["coordinates"][0][0])
        # thickEnd - The ending position at which the feature is drawn thickly (for example the stop codon in gene displays).
        line.append(d["coordinates"][-1][1])
        # itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
        line.append(0)
        # blockCount - The number of blocks (exons) in the BED line.
        line.append(len(d["coordinates"]))
        # blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
        line.append(",".join([str(e-b) for i, (b,e)  in enumerate(d["coordinates"])]))
        # blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.
        line.append(",".join([str(b - d["coordinates"][0][0]) for i, (b,e) in enumerate(d["coordinates"])]))

        print(DELIMITER.join([str(x) for x in line]))

def main():
    logging.basicConfig(level=logging.DEBUG)
    if (len(sys.argv) != 3):
        logging.error("Wrong argument count")
        exit(1)
    gff3_filename, gene_positions_filename = sys.argv[1:3]

    positions = get_positions(gff3_filename, gene_positions_filename, "Pfam")
    print2bed(positions)

if __name__ == "__main__":
    main()