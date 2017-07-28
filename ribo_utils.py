import os
import operator
import itertools
import gzip
import numpy as np
from scipy import stats
import cPickle as pickle
import math
import multiprocessing
from collections import defaultdict
'''
Colorblind safe colors from Bang Wong, Nature Methods 8. 441 (2011)
'''
black = (0,0,0)
orange = (230/255.0,159/255.0,0)
skyBlue = (86/255.0,180/255.0,233/255.0)
bluishGreen = (0,158/255.0,115/255.0)
yellow = (240/255.0,228/255.0,66/255.0)
blue = (0,114/255.0,178/255.0)
vermillion = (213/255.0,94/255.0,0)
reddishPurple = (204/255.0,121/255.0,167/255.0)
colors = [black, orange, skyBlue, bluishGreen, vermillion, blue, reddishPurple, yellow]
rainbow = [black, vermillion, orange, bluishGreen, blue, reddishPurple, 'violet']
markers = ['.', 'o', 'v', 's', '^', 'p', 'x', '+']
line_styles = ['solid', 'dashed', 'dotted']

bokeh_black = (0,0,0)
bokeh_orange = (230,159,0)
bokeh_skyBlue = (86,180,233)
bokeh_bluishGreen = (0,158,115)
bokeh_yellow = (240,228,66)
bokeh_blue = (0,114,178)
bokeh_vermillion = (213,94,0)
bokeh_reddishPurple = (204,121,167)

###############################
#Parralellization code from
# http://stackoverflow.com/questions/3288595/multiprocessing-using-pool-map-on-a-function-defined-in-a-class
###############################
def spawn(f):
    def fun(q_in,q_out):
        while True:
            i,x = q_in.get()
            if i == None:
                break
            q_out.put((i,f(x)))
    return fun

def parmap(f, X, nprocs = multiprocessing.cpu_count()):
    q_in   = multiprocessing.Queue(1)
    q_out  = multiprocessing.Queue()
    proc = [multiprocessing.Process(target=spawn(f),args=(q_in,q_out)) for _ in range(nprocs)]
    for p in proc:
        p.daemon = True
        p.start()
    sent = [q_in.put((i,x)) for i,x in enumerate(X)]
    [q_in.put((None,None)) for _ in range(nprocs)]
    res = [q_out.get() for _ in range(len(sent))]
    [p.join() for p in proc]
    return [x for i,x in sorted(res)]

##########
#FILE HANDLING
##########
def unPickle(fileName):
    #returns the pickled object stored in a pickle file
    f = open(fileName, 'r')
    o = pickle.load(f)
    f.close()
    return o

def makePickle(o, fileName, protocol=pickle.HIGHEST_PROTOCOL):
    f = open(fileName, 'w')
    pickle.dump(o, f, protocol=protocol)
    f.close()

def make_dir(dirname):
    """
    Makes the directory; doesn't throw an error if it exists.
    """
    if not os.path.exists(dirname):
        try:
            os.makedirs(dirname)
        except:
            print 'The directory was made by another thread extremely recently.'

def file_exists(fname):
    """
    makes sure a given file exists
    """
    if not os.path.exists(fname):
        return False
    fstats = os.stat(fname)
    if not fstats[6]:
        return False
    if not os.access(fname, os.R_OK):
        raise ValueError('Input File %s cannot be read' % fname)
    return True

def tsv_to_dict(filename, header = True, delimiter = '\t', key_column = 0, convert_to_float = False):
    """
    Will return a dict index first by the row labels, then by the column headers
    """
    return_dict = {}
    f = open(filename)
    lines  = f.readlines()
    headers = lines[0].strip('\n').split(delimiter)
    for line in lines[1:]:
        ll= line.strip('\n').split(delimiter)
        return_dict[ll[key_column]] = {}
        for i in range(0, len(ll)):
            if not i == key_column:
                if not convert_to_float:
                    return_dict[ll[key_column]][headers[i]]=ll[i]
                elif is_float(ll[i]):
                    return_dict[ll[key_column]][headers[i]]=float(ll[i])
    f.close()
    return return_dict

##########
#MATH
##########
def divideWithError(num, stdDevNum, denom, stdDevDenom):
    '''
    divides the two values with provided standard deviations, and returns the mean and error of the ratio using standard error propogation
    '''
    num = float(num)
    denom = float(denom)
    stdDevNum = float(stdDevNum)
    stdDevDenom = float(stdDevDenom)

    ratio = num/denom
    ratioError = ratio*math.sqrt((stdDevNum/num)**2+(stdDevDenom/denom)**2)

    return ratio, ratioError

def subtractWithError(num, stdDevNum, denom, stdDevDenom):
    '''
    divides the two values with provided standard deviations, and returns the mean and error of the ratio using standard error propogation
    '''
    num = float(num)
    denom = float(denom)
    stdDevNum = float(stdDevNum)
    stdDevDenom = float(stdDevDenom)

    ratio = num/denom
    ratioError = ratio*math.sqrt((stdDevNum/num)**2+(stdDevDenom/denom)**2)

    return ratio, ratioError

def next_square_number(number):
    return int(math.ceil(math.sqrt(number)))**2

def computePfromMeanAndStDevZscore(mean, standard_deviation, testValue):
    #computes probability that test value came from a gaussian with the given mean and standard deviation
    try:
        z = (float(testValue)-mean)/standard_deviation
        p = stats.norm.sf(z)
        return p, z
    except ZeroDivisionError:
        return 0.5, 0

def ranges_overlap(min1, max1, min2, max2):
    """

    :param min1:
    :param max1:
    :param min2:
    :param max2:
    :return: return True if the 2 ranges overlap (edge inclusive), else False
    """

    if min1 <= max2 and min2 <= max1:
        return True
    return False

def number_passing_cutoff(numbers, cutoff):
    i = 0
    for number in numbers:
        if number >= cutoff:
            i += 1
    return i

def significantly_enriched(xs, zthresh=2., scale='linear'):
    assert scale in ['linear', 'log']
    if scale =='log':
        xs = np.log2(xs)
    xs = stats.zscore(xs)
    return [x > zthresh for x in xs]

def filter_x_y_pairs(x, y, filter_list = [float('inf'), -1*float('inf')]):
    '''
    takes 2 paired arrays, and returns matched copies of them with any positions with values in
    filter_list removed from both arrays, to keep them synced.
    alos removes NaN (defined by testing if the entry equals itself, which fails for NaN)
    :param filter_list: list of values to remove
    :return:
    '''
    filtered_x, filtered_y = [], []
    assert len(x) == len(y)
    for i in range(len(x)):
        if x[i] not in filter_list and y[i] not in filter_list and x[i]==x[i] and y[i]==y[i]:
            filtered_x.append(x[i])
            filtered_y.append(y[i])
    return np.array(filtered_x), np.array(filtered_y)

def benjaminiHochbergCorrection(pValDict):
    '''
    takes a dictionary mapping key to p value
    returns a dictionary of Benjamini-Hochberg corrected Q values

    Q = p * n / k, where n is the # of observations, and k is the rank of the particular p-value among all p-values
    '''
    qValues = {}
    sorted_p = sorted(pValDict.iteritems(), key=operator.itemgetter(1))
    n = len(sorted_p)
    for i in range(n):
        k = i+1
        q = sorted_p[i][1] * n / k
        qValues[sorted_p[i][0]] = q
    return qValues

def bonferroniCorrection(pValDict):
    '''
    takes a dictionary mapping key to p value
    returns a dictionary of Bonferroni corrected Q values

    Q = p * n, where n is the # of observations
    '''
    qValues = {}
    sorted_p = sorted(pValDict.iteritems(), key=operator.itemgetter(1))
    n = len(sorted_p)
    for i in range(n):
        q = sorted_p[i][1] * n
        qValues[sorted_p[i][0]] = q
    return qValues

##################
#SEQUENCE HANDLING
##################

GENETIC_CODE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

def rna(dna_seq):
    return dna_seq.replace('T','U').replace('t', 'u')

def get_barcode(line):
    """
    - Extracts the barcode from the first line of a fastq quartet
        - Assumes the first line is of the form:
            @D5FF8JN1:4:1101:1220:2099#ACTTGA/1
    """
    return line.split('#')[-1].split('/')[0]

def convertFastaToDict(fastaFile):
    '''
    converts a fasta file to a dict of {sequenceName:sequence}
    can take extra files in * args
    '''
    if isinstance(fastaFile, list):
        files = fastaFile
    else:
        files = [fastaFile]
    currentName = None
    currentSequence = None
    seqDict = {}
    for currentFile in files:
        if currentFile.endswith('.gz'):
            f = gzip.open(currentFile)
        else:
            f = open(currentFile)
        for line in f:
            if not line.strip() == '' and not line.startswith('#'):  # ignore empty lines and commented out lines
                if line.startswith('>'):  # > marks the start of a new sequence
                    if not currentName == None:  # after we've reached the firtst > line, we know what the sequence corresponds to
                        seqDict[currentName] = currentSequence
                    currentName = line.strip()[1:].split()[
                        0]  # i've noticed the gencode names have extraneous numbering after some whitespace. This doens't match the GTF files, so I'm removing it.
                    currentSequence = ''
                else:
                    currentSequence += line.strip()
        f.close()
    seqDict[currentName] = currentSequence
    return seqDict

def hamming_N(str1, str2):
    if not len(str1) == len(str2):
        raise(ValueError, 'lengths don\'t match')
    str1 = str1.upper()
    str2 = str2.upper()
    str1 = str1.replace('N', '#')
    return sum(itertools.imap(operator.ne, str1, str2))

# from http://code.activestate.com/recipes/499304-hamming-distance/
def hamming_distance(str1, str2):
    assert len(str1) == len(str2)
    ne = operator.ne
    return sum(itertools.imap(ne, str1, str2))

def reverse_complement(seq, isRNA = False):
    seq = seq.upper()
    compDict = {'A':'T', 'T':'A', 'U':'A', 'C':'G', 'G':'C', 'N':'N', '-':'-', '.':'.', '*':'*'}
    revComp = ''.join([compDict[c] for c in seq[::-1]])
    if isRNA:
        return revComp.replace('T', 'U')
    return revComp

def close_float_value(a, b, max_percent=1.0):
    if a == 0 and b == 0:
        return True
    if not (a > 0 and b > 0):
        return False
    ratio = float(max(a, b)) / float(min(a, b))
    percent_increase = (ratio - 1.0) * 100.0
    return percent_increase < max_percent

def getAllMismatchedSeqs(kmer, mismatchPositions):
    nucs = ['A', 'C', 'G', 'T']
    #generate tuples of allowed nucs at each mismatch position using a recursive algorithm
    allowedNucs = {}
    mismatchPositions = np.array(mismatchPositions)
    assert len(set(mismatchPositions)) == len(mismatchPositions)
    if len(mismatchPositions) == 0:
        yield kmer
    else:
        mismatchNucs = [] + nucs
        #print kmer
        #print mismatchPositions
        #print mismatchPositions[0]
        #print kmer[mismatchPositions[0]]
        mismatchNucs.remove(kmer[mismatchPositions[0]])
        downstreamMismatchSeqs = getAllMismatchedSeqs(kmer[mismatchPositions[0]+1:], mismatchPositions[1:]-(mismatchPositions[0]+1))
        for mismatchNuc in mismatchNucs:
            for downstreamMismatchSeq in downstreamMismatchSeqs:
                returnSeq = kmer[:mismatchPositions[0]] + mismatchNuc +downstreamMismatchSeq
                assert len(returnSeq) == len(kmer)
                yield returnSeq

def getPaddedMismatchedAdjacentKmers(kmerSequence, padding, numMismatches):
    '''
    Yield all sequences of length (len(kmerSequence)+padding )that contain the given kmer, with exactly the given number of mismatches.
    The order yielded is as follows:
        First mismatches are allowed at position 0 to (numMismatches-1)
            For each register:
                Iterate through all possible nucs at mismatch position in alphabetical order
                    Iterate through each nucleotide in padding positions in alphabetical order.
                Shift to next register
            Move most 3' mismatch position down by one, but not past the end of the kmerSequence if end of KmerSequence
            is reached, shift secondmost 3' mismatch 1 nt 3', and reset most 3' mismatch to 1nt 3' of that one
    '''

    # for troubleshooting, want to check that no repeats are generated, so will assert that size of this list and set
    # must be the same
    kmer_set = set()
    kmer_list =[]
    upper_to_combined = {}
    nucs = 'ACGT'
    #initialize mismatchPositions
    #print numMismatches
    if numMismatches == 0:
        for mismatchedKmer in [kmerSequence]:
            for shift in range(padding+1):
                #generate all possible mismatches to the kmer
                for leftPaddingSeq in [''.join(i) for i in itertools.product(nucs, repeat = shift)]:
                    for rightPaddingSeq in [''.join(i) for i in itertools.product(nucs, repeat = padding-shift)]:
                        paddedSeq = leftPaddingSeq+mismatchedKmer+rightPaddingSeq
                        if paddedSeq not in kmer_set:
                            kmer_list.append(paddedSeq)
                        kmer_set.add(paddedSeq)
    else:
        mismatchPositionsList = itertools.combinations(range(len(kmerSequence)), numMismatches)
        for mismatchPositions in mismatchPositionsList:
            #print mismatchPositions
            for mismatchedKmer in getAllMismatchedSeqs(kmerSequence, mismatchPositions):
                for shift in range(padding+1):
                    #generate all possible mismatches to the kmer
                    for leftPaddingSeq in [''.join(i) for i in itertools.product(nucs, repeat = shift)]:
                        for rightPaddingSeq in [''.join(i) for i in itertools.product(nucs, repeat = padding-shift)]:
                            paddedSeq = leftPaddingSeq+mismatchedKmer+rightPaddingSeq
                            paddedUpper = paddedSeq.upper()
                            if paddedUpper not in kmer_set:
                                kmer_list.append(paddedUpper)
                            kmer_set.add(paddedUpper)

    #print kmer_list
    #print kmer_set
    #print len(kmer_list), len(kmer_set)
    #assert len(kmer_list) == len(kmer_set)

    return kmer_list


##################
#GTF file parsing and handling
##################

class genome_sequence():
    def __init__(self, fasta_file, *args):
        self.genome_sequence = convertFastaToDict(fasta_file, *args)

    def get_sequence(self, chromosome, start, end, strand):
        """
        returns a string of genome sequence at the given start and end positions, inclusive
        reverse-complemented for minus strand

        start and end are 1-indexed (first base-pair of genome is 1)
        """
        assert end >= start
        sequence = self.genome_sequence[chromosome][start - 1: end]
        if strand == '-':
            return reverse_complement(sequence)
        else:
            return sequence

class gtf_data():
    def __init__(self, gtf_file):
        self.gtf_entries = []
        self.transcript_to_entries = defaultdict(set)
        self.gene_to_entries = defaultdict(set)
        self.genes_to_tx = defaultdict(set)
        self.chr_to_entry = defaultdict(lambda : defaultdict(set))
        self.feature_type_summary = defaultdict(int)
        self.transcript_type_summary = defaultdict(int)
        self.add_gtf_data(gtf_file)

    def add_gtf_data(self, gtf_file):
        if gtf_file.endswith('.gz'):
            gtf = gzip.open(gtf_file)
        else:
            gtf = open(gtf_file)
        for line in gtf:
            if not line.startswith('#'):
                new_entry = gtf_entry(line)
                self.gtf_entries.append(new_entry)
                self.feature_type_summary[new_entry.get_value('type')] += 1
                self.transcript_type_summary[new_entry.get_value('transcript_type')] += 1
                gene_id = new_entry.get_value('gene_id')
                transcript_id = new_entry.get_value('transcript_id')
                strand = new_entry.get_value('strand')
                chromosome = new_entry.get_value('chr')
                self.chr_to_entry[strand][chromosome].add(new_entry)
                if gene_id != None:
                    self.gene_to_entries[gene_id].add(new_entry)
                    if transcript_id != None:
                        self.transcript_to_entries[transcript_id].add(new_entry)
                        self.genes_to_tx[gene_id].add(transcript_id)
        gtf.close()

    def print_transcript_multiplicity(self, gene_type=None):
        self.tx_counts_histogram = defaultdict(int)
        for gene_id in self.genes_to_tx:
            if gene_type == None or gene_type == sorted(self.gene_to_entries[gene_id])[0].get_value('gene_type'):
                self.tx_counts_histogram[len(self.genes_to_tx[gene_id])] += 1
        for count in sorted(self.tx_counts_histogram.keys()):
            print count, self.tx_counts_histogram[count]

    def spliced_length(self, transcript_id, exon_type='exon'):
        """
        exon_type can be CDS or exon.
        CDS wil start and end at CDS boundaries, so that's convenient
        Returns lenth of transcript or cds
        """
        ordered_exon_entries = self.sorted_exons(transcript_id, exon_type=exon_type)
        if len(ordered_exon_entries) == 0:
            return 0
        transcript_length = sum([exon_entry.length() for exon_entry in ordered_exon_entries])
        return transcript_length

    def sorted_exons(self, transcript_id, exon_type='exon'):
        """
        exon_type can be : CDS, exon, UTR, stop_codon, start_codon, or a list containing a combination therof.
        CDS wil start and end at CDS boundaries, but excludes the stop codon. So need to pass ['CDS','stop_codon'] to get the full coding sequence
        -Be careful not to mix annotation types that may overlap, for example exon, with any other, as you will get the wrong sequence, with duplicates.
        Returns exons in annotated order, based on start position of each exon
        Ordering is relative to the sense strand, so the first exon in the list will be the 5'-most exon in the transcript.
        However, the "end" of the exon boundary is always larger than the 'start'
        """
        transcript_entries = self.transcript_to_entries[transcript_id]
        ordered_exon_entries = sorted([entry for entry in transcript_entries if entry.is_type(exon_type)],
                                      key=lambda x: int(x.get_value('start')))
        # if this transcript is on the minus strand, the exon order needs to be flipped
        if len(ordered_exon_entries) > 0 and ordered_exon_entries[0].get_value('strand') == '-':
            ordered_exon_entries = ordered_exon_entries[::-1]
        return ordered_exon_entries

    def transcript_sequence(self, genome_sequence, transcript_id, exon_type='exon'):
        """
        exon_type can be CDS or exon.
        CDS will start and end at CDS boundaries, so that's convenient
        Returns sequence of transcript or cds
        """
        ordered_exon_entries = self.sorted_exons(transcript_id, exon_type=exon_type)
        transcript_sequence = ''.join([exon_entry.sequence(genome_sequence) for exon_entry in ordered_exon_entries])
        return transcript_sequence

    def tx_with_longest_CDS(self, gene_id, starting_subset=None):
        """
        starting_subset can be a list of transcript ids. If it is given, then only thos etranscripts will be considered
        """
        if starting_subset == None:
            transcripts = self.genes_to_tx[gene_id]
        else:
            transcripts = starting_subset
        if len(transcripts) == 1:
            return [sorted(transcripts)[0]]
        else:
            sorted_transcripts = sorted(transcripts,
                                        key=lambda x: int(self.spliced_length(x, exon_type=['CDS', 'stop_codon'])),
                                        reverse=True)
            longest_CDS_length = self.spliced_length(sorted_transcripts[0], exon_type=['CDS', 'stop_codon'])
            return [x for x in sorted_transcripts if
                    self.spliced_length(x, exon_type=['CDS', 'stop_codon']) == longest_CDS_length]

    def longest_tx(self, gene_id, starting_subset=None):
        if starting_subset == None:
            transcripts = self.genes_to_tx[gene_id]
        else:
            transcripts = starting_subset
        if len(transcripts) == 1:
            return [sorted(transcripts)[0]]
        else:
            sorted_transcripts = sorted(transcripts, key=lambda x: int(self.spliced_length(x, exon_type='exon')),
                                        reverse=True)
            longest_CDS_length = self.spliced_length(sorted_transcripts[0], exon_type='exon')
            return [x for x in sorted_transcripts if self.spliced_length(x, exon_type='exon') == longest_CDS_length]

    def pick_all_longest_CDS_transcripts(self):
        # picks trnascripts with longest CDS
        # If tied picks longest TX
        # Otherwise, pick the first one randomly and make note
        genes_with_ties = []
        chosen_tx = []
        for gene_id in self.genes_to_tx:
            tx_with_longest_CDS = self.tx_with_longest_CDS(gene_id)
            assert len(tx_with_longest_CDS) > 0
            if len(tx_with_longest_CDS) == 1:
                chosen_tx.append(tx_with_longest_CDS[0])
            else:
                tx_with_longest_tx = self.longest_tx(gene_id, starting_subset=tx_with_longest_CDS)
                assert len(tx_with_longest_tx) > 0
                if len(tx_with_longest_tx) == 1:
                    chosen_tx.append(tx_with_longest_tx[0])
                else:
                    genes_with_ties.append(gene_id)
                    chosen_tx.append(tx_with_longest_tx[0])
        print 'genes with ties for longest CDS and tx:', len(genes_with_ties)
        assert len(chosen_tx) == len(set(chosen_tx))
        return chosen_tx

    def filter_transcripts_by_value(self, key, allowed_values, starting_subset=None):
        # returns all entries for which the given key matches one of the allowed values
        chosen_tx = []
        if starting_subset == None:
            starting_subset = self.transcript_to_entries.keys()
        for transcript_id in starting_subset:
            if sorted(self.transcript_to_entries[transcript_id])[0].get_value(key) in allowed_values:
                chosen_tx.append(transcript_id)
        assert len(chosen_tx) == len(set(chosen_tx))
        return chosen_tx

    def write_transcript_entries_to_file(self, out_file, transcript_ids=None):
        if transcript_ids == None:
            transcript_ids = self.transcript_to_entries.keys()
        if out_file.endswith('.gz'):
            out_gtf = gzip.open(out_file, 'w')
        else:
            out_gtf = open(out_file, 'w')
        transcript_entries = []
        for transcript_id in transcript_ids:
            for transcript_entry in self.transcript_to_entries[transcript_id]:
                transcript_entries.append(transcript_entry)
        for transcript_entry in sorted(transcript_entries,
                                       key=lambda x: (x.get_value('chr'), int(x.get_value('start')))):
            out_gtf.write(transcript_entry.gtf_file_line)
        out_gtf.close()

    def find_smallest_annotation_at_position(self, chr, strand, position, type_restrictions=None):
        '''
        Finds the smallest (smallest end-start) entry at a given position
        :param chr: 
        :param position: 
        :return: 
        '''
        if type_restrictions == None:
            entries = self.entries_by_position[chr][strand][position]
        else:
            entries = [entry for entry in self.entries_by_position[chr][strand][position]
                       if entry.get_value('type') in type_restrictions]
        sorted_by_length = sorted(entries, key=lambda x: (int(x.get_value('end')) - int(x.get_value('start'))))
        return sorted_by_length[0]

class gtf_entry():
    # note: levels 1 and 2 are verified and manually annotated, repectively, 3 are automatically annotated
    def __init__(self, gtf_file_line):
        self.gtf_file_line = gtf_file_line
        self.fields = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'frame', 'additional']
        self.additional_mandatory_keys = ['gene_id', 'transcript_id', 'gene_type', 'gene_status', 'gene_name',
                                          'transcript_type', 'transcript_status', 'transcript_name', 'exon_number',
                                          'exon_id', 'level']
        ll = gtf_file_line.rstrip('\n').split('\t')
        self.primary_data = dict(zip(self.fields, ll))
        additional_pairs = self.primary_data['additional'].split('; ')
        self.secondary_data = dict([pair.split(' ') for pair in additional_pairs])
        for key in self.secondary_data:
            self.secondary_data[key] = self.secondary_data[key].strip('"')

    def __repr__(self):
        return str(self.primary_data)

    def is_type(self, entry_type):
        """
        Check if this entry is of the given primary type (third column of gtf file), or in the given list of types
        """

        if isinstance(entry_type, str):
            return self.primary_data['type'] == entry_type
        elif isinstance(entry_type, list):
            return self.primary_data['type'] in entry_type
        else:
            raise Exception("entry_type should be a string or list of strings, recieved type %s" % type(entry_type))

    def get_value(self, key):
        assert not (key in self.primary_data and key in self.secondary_data)
        if key in self.primary_data:
            return self.primary_data[key]
        elif key in self.secondary_data:
            return self.secondary_data[key]
        else:
            return None

    def length(self):
        return (int(self.get_value('end')) - int(self.get_value('start'))) + 1

    def sequence(self, genome_sequence):
        """
        return the sense strand sequence of this element
        This accounts for the strand information, so minus strand elements will be reverse complemented
        """
        return genome_sequence.get_sequence(self.get_value('chr'), int(self.get_value('start')),
                                            int(self.get_value('end')), self.get_value('strand'))