import os
import operator
import itertools
import gzip
import numpy as np
from scipy import stats
import cPickle as pickle
import math
import multiprocessing
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

#parralellization code from
# http://stackoverflow.com/questions/3288595/multiprocessing-using-pool-map-on-a-function-defined-in-a-class
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
    '''
    currentName = None
    currentSequence = None
    seqDict = {}
    f = open(fastaFile)
    for line in f:
        if not line.strip() == '' and not line.startswith('#'):#ignore empty lines and commented out lines
            if line.startswith('>'):#> marks the start of a new sequence
                if not currentName == None: #after we've reached the firtst > line, we know what the sequence corresponds to
                    seqDict[currentName] = currentSequence.upper()
                currentName = line.strip()[1:]
                currentSequence = ''
            else:
                currentSequence += line.strip()
    f.close()
    seqDict[currentName] = currentSequence.upper()
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

def revComp(seq, isRNA = False):
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






