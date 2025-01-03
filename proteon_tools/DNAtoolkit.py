nucleotides = {
    "DNA": ["A", "T", "C", "G"],
    "RNA": ["A", "U", "C", "G"]}


def rndDNAstr(length):
    """random DNA string generator"""
    import random
    return ''.join([random.choice(nucleotides["DNA"])
                    for nuc in range(length)])

def rndRNAstr(length):
    """random RNA string generator"""
    import random
    return ''.join([random.choice(nucleotides["RNA"])
                    for nuc in range(length)])




def validateDNAseq (seq):
    """check the sequence to make sure it is a DNA string"""
    tmpseq = seq.upper()
    for nuc in tmpseq:
        if nuc not in nucleotides["DNA"]:
            return False
    return tmpseq

def validateRNAseq (seq):
    """check the sequence to make sure it is a rnA string"""
    tmpseq = seq.upper()
    for nuc in tmpseq:
        if nuc not in nucleotides["RNA"]:
            return False
    return tmpseq




def countnucfrq (seq):
    """counting nucleotide frequency"""
    tmpFreqDict={'A':0, 'T':0, 'C':0, 'G':0, 'U':0 }
    for nuc in seq:
        tmpFreqDict[nuc]+=1
    return tmpFreqDict

        


def transcription(seq):
    """replacing thymine with uracil""" 
    return seq.replace('T','U')

DNA_reverseComplement={'A':'T','T':'A','G':'C','C':'G'}

def reverse_complement(seq):
    """swapping A with T and G with C. reversing newly generated string"""
    return ''.join([DNA_reverseComplement[nuc] for nuc in seq])[::-1]




def GC_content(seq):
    """GC content in a DNA/RNA sequence"""
    return round((seq.count('C')+seq.count('G'))/len(seq)*100,2)

def GC_content_subseq(seq, k=20):
    """GC_content in a DNA/RNA subsequence length k. k=20 by defult"""
    res=[]
    for i in range(0,len(seq)-k+1,k):
        subseq=seq[i:i+k]
        res.append(GC_content(subseq))
    return res




DNA_codons = {
    # 'M' - START, '_' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}



def translate_DNAseq (seq, init_pos=0):
    """translating a DNA sequence into an aminoacid sequence"""
    return [DNA_codons[seq[pos:pos+3]] for pos in range(init_pos,len(seq)-2,3)]

def DNA_translation (seq, init_pos=0):
    """translating a DNA sequence into an aminoacid sequence"""
    aminoacids=[DNA_codons[seq[pos:pos+3]] for pos in range(init_pos,len(seq)-2,3)]
    polypeptide=[]
    for i in aminoacids:
        if i=='_':
            break
        else:
            polypeptide.append(i)
    return polypeptide



RNA_codons = {
    # 'M' - START, '_' - STOP
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UGU": "C", "UGC": "C",
    "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "UUU": "F", "UUC": "F",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAU": "H", "CAC": "H",
    "AUA": "I", "AUU": "I", "AUC": "I",
    "AAA": "K", "AAG": "K",
    "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUG": "M",
    "AAU": "N", "AAC": "N",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UGG": "W",
    "UAU": "Y", "UAC": "Y",
    "UAA": "_", "UAG": "_", "UGA": "_"
}


def translate_RNAseq (seq, init_pos=0):
    """translating a RNA sequence into an aminoacid sequence"""
    return [RNA_codons[seq[pos:pos+3]] for pos in range(init_pos,len(seq)-2,3)]


def RNA_translation (seq, init_pos=0):
    """translating a RNA sequence into an aminoacid sequence"""
    aminoacids=[RNA_codons[seq[pos:pos+3]] for pos in range(init_pos,len(seq)-2,3)]
    polypeptide=[]
    for i in aminoacids:
        if i=='_':
            break
        else:
            polypeptide.append(i)
    return polypeptide


def codon_usage(seq, aminoacid):
    """Provides the frequency of each codon encoding a given aminoacid in a DNA sequence"""
    tmplist=[]
    for i in range(0,len(seq)-2,3):
        if DNA_codons[seq[i:i+3]]==aminoacid:
            tmplist.append(seq[i:i+3])

    import collections

    freqdict=dict(collections.Counter(tmplist))
    totalwight=sum(freqdict.values())
    for codon in freqdict:
        freqdict[codon]=round(freqdict[codon]/totalwight,2)
    return freqdict

