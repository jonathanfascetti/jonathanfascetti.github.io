"""
Jonathan Fascetti wrote original
C. Matthew Sundling made edits
OChem Spring 2020 @ BCDS
Project 2 Solution
"""

from fascetti_sundling_collaboration_read_and_comp_fasta import *
import sys
import logging  # useful for print statements throughout your code

# Set logging level to:
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def combine_full_sequence_string(seq_contents):
    full_seq = "".join(seq_contents)  # join all seq pieces
    full_seq = full_seq.replace("\n", "")  # remove endlines
    return full_seq


def find_start(full_seq, full_list=True):  # provided default value
    if full_list == False:
        logging.debug("sequence from first ATG: (since full_list == False)")
        start_loc_in_line = full_seq.find("ATG")
        logging.debug("start_loc_in_line --> " + str(start_loc_in_line))
        cut_seq = full_seq.replace(full_seq[:start_loc_in_line], "")
        logging.debug("cut_seq --> " + cut_seq)
        return cut_seq
    else:
        logging.debug("whole sequence: (since full_list == True)")
        cut_seq = full_seq
        logging.debug("cut_seq --> " + cut_seq)
        return cut_seq


def split_codons(cut_seq):
    amino_acid_codons = []
    amino_acid_codons = list(
        cut_seq[0 + 3 * i : 3 + 3 * i] for i in range(0, len(cut_seq) // 3)
    )
    logging.debug("amino_acid_codons --> " + str(amino_acid_codons))
    return amino_acid_codons


def triplet_to_aa(amino_acid_codons, full_list=True):
    genetic_code = {}  # list to hold aa

    with open(
        "codon_table.txt", "r"
    ) as aa_file:  # ref list, make sure in same directroy
        lines = aa_file.readlines()
        for l in lines:
            index = l.find("#")
            codon = l[:index].split(":")
            if len(codon) == 2:  # found a codon
                genetic_code[codon[0].strip()] = codon[1].strip()

    position = 0
    if full_list == False:  # find location of first "ATG"
        position = amino_acid_codons.index("ATG")

    aa_seq = []
    for codon in amino_acid_codons[position:]:
        aa_seq.append(genetic_code[codon])
    return aa_seq


def find_stop(aa_list):
    indices = [str(i) for i, x in enumerate(aa_list) if x == "*"]

    if indices == []:
        return "Stop_location(s)= No stop codons found"

    stop_info = "Stop_location(s)= " + indices[0]  # there is at least one
    for i in indices[1:]:
        stop_info += ", " + i
    return stop_info  # final line to go into new fasta


def find_db(aa_list):
    base_list = ["R", "K", "H"]

    db_loc = []
    for i in range(len(aa_list) - 1):
        if aa_list[i] in base_list and aa_list[i + 1] in base_list:
            db_loc.append(i)

    if db_loc == []:
        return "DB_locations(s)= No double bases found"

    db_info = (
        "DB_locations= "
        + aa_list[db_loc[0]]
        + aa_list[db_loc[0] + 1]
        + ": "
        + str(db_loc[0])
    )  # there is at least one

    for i in db_loc[1:]:
        db_info += ", " + aa_list[i] + aa_list[i + 1] + ": " + str(i)
    return db_info  # final line to go into new fasta


def create_peptide_list(aa_list, line_length=70):
    logging.debug("peptide_list (individual) --> " + str(aa_list))

    joined = "".join(aa_list)
    total_length = len(joined)
    logging.debug("total_length --> " + str(total_length))

    aa_lines = []
    for i in range(0, total_length, line_length):
        aa_lines.append(joined[i : i + line_length])

    logging.debug("peptide_list (final) --> " + str(aa_lines))
    return aa_lines, total_length


def write_new_file(
    comment_line, stop_info, db_info, peptide_list, total_length
):
    comment_line = comment_line.strip()
    comment_line_ext = (
        comment_line
        + " | Length= "
        + str(total_length)
        + ", "
        + stop_info[:]
        + ", "
        + db_info
    )

    new_file = open(
        "fascetti_sundling_collaboration_project2.output.fasta", "w"
    )
    new_file.write(comment_line_ext)
    q = 0
    for lines in peptide_list:
        new_file.write(
            "\n" + peptide_list[q]
        )  # start printing out peptide_list line by line
        q = q + 1
    new_file.close()


comment_line = ret_comment(
    "pa1"
)  # put in fasta file name -.fasta in same directory
title = ret_title(comment_line)
seq_contents = ret_seq(
    "pa1"
)  # put in fasta file name -.fasta in same directory
full_seq = combine_full_sequence_string(seq_contents)

codon_list = split_codons(find_start(full_seq))
aa_list = triplet_to_aa(codon_list)
stop_info = find_stop(aa_list)
db_info = find_db(aa_list)

peptide_list, total_length = create_peptide_list(
    aa_list, len(seq_contents[0]) - 1
)
write_new_file(comment_line, stop_info, db_info, peptide_list, total_length)
