"""
Author  : D. R. R. Wijewardene      S14245
Date    : 12/01/2022
Project2 - Sequence Aligner - IMPLEMENTATION (MAIN PROGRAM)
"""

# implementation
import project_class

if __name__ == '__main__':
    sequences_file = input(
        "Insert a FASTA file with the desired sequences \n\t(for demonstration insert \"sequences_file.FASTA\") : ")
    seqs = project_class.SeqFile(sequences_file)

    # 1 - sequence type
    check_type_id = input("\nInsert sequence ID to check sequence type \n\t(for demonstration insert \"seq2_P\") : ")
    for i in range(len(seqs.records)):
        if seqs.records[i].id == check_type_id:
            print("Sequence Type of %s : %s\n" % (seqs.records[i].id, seqs.get_seq_type(seqs.records[i].seq)))

    # 5 - find a sequence from GenBank when its accession number is given
    acc = input(
        "Insert an accession number to get the sequence from GenBank \n\t(for demonstration insert \"AB044088.1\") : ")
    seqs.get_from_genbank(acc)
    print("GenBank record of %s was obtained and saved !" % acc)

    # 6 - calculate the melting temperature of a given nucleotide sequence
    mt_id = input(
        "\nInsert sequence ID of nucleotide sequence to calculate melting temperature \n\t(for demonstration insert \"seq7_N\") : ")
    for i in range(len(seqs.records)):
        if seqs.records[i].id == mt_id:
            seqs.melting_temp(seqs.records[i].seq)

    # 7 - calculate the isoelectric point of a given amino acid sequence
    ip_id = input(
        "\nInsert sequence ID of amino acid sequence to calculate isoelectric point \n\t(for demonstration insert \"seq3_P\") : ")
    for i in range(len(seqs.records)):
        if seqs.records[i].id == ip_id:
            seqs.isoelectric(seqs.records[i].seq)

    # 8 - calculate the aromaticity of a given amino acid sequence
    aro = input(
        "\nInsert sequence ID of amino acid sequence to calculate aromaticity \n\t(for demonstration insert \"seq4_P\") : ")
    for i in range(len(seqs.records)):
        if seqs.records[i].id == aro:
            seqs.aromaticity(seqs.records[i].seq)

    # 2 - global pairwise similarity score between all sequence pairs in the input file
    seqs.get_gps_score()
    print(
        "\nGlobal pairwise similarity scores between all sequence pairs in the input file were obtained and saved as global_pss.txt !\n")

    # 3 - scatter plot to represent pairwise similarity for two given sequences
    print("Insert sequence IDs of two sequences to obtain pairwise similarity scatter plot :")
    nuc1 = input("\tSequence 1 (for demonstration insert \"seq5_N\") : ")
    nuc2 = input("\tSequence 2 (for demonstration insert \"seq6_N\") : ")
    nuc1_record = ""
    nuc2_record = ""
    for i in range(len(seqs.records)):
        if seqs.records[i].id == nuc1:
            nuc1_record = seqs.records[i]
        elif seqs.records[i].id == nuc2:
            nuc2_record = seqs.records[i]
    seqs.scatter_plot(nuc1_record, nuc2_record)
    print("Scatter plot of pairwise similarity between %s and %s was obtained and saved !\n" % (nuc1_record.id, nuc2_record.id))

    # 4 - run BLAST searches on several sequence inputs in a single FASTA file
    print("BLAST search running...")
    seqs.seq_blast()
    print("BLAST searches on all sequence inputs in the input file were obtained and saved as blast_results.xml !")
