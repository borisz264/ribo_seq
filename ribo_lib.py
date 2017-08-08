from collections import defaultdict
import re

import pysam
import ribo_utils
import numpy as np
import sys
import math

def assign_tx_reads(experiment_settings, lib_settings):
    lib_settings.write_to_log('counting reads or loading counts')
    if experiment_settings.get_property('force_recount') or not lib_settings.sequence_counts_exist():
        print "counting BAM reads"
        transcripts = {}
        tx_annotations = ribo_utils.tsv_to_dict(experiment_settings.get_property('canonical_tx_features'))
        tx_seqs = ribo_utils.convertFastaToDict(experiment_settings.get_property('canonical_tx_seqs'))
        samfile = pysam.AlignmentFile(lib_settings.get_transcript_mapped_reads(), "rb")
        for tx_id in [t for t in samfile.references if t in tx_annotations and t in tx_seqs]:
            transcripts[tx_id] = transcript(tx_id, tx_annotations[tx_id], tx_seqs[tx_id], samfile)
        samfile.close()
        ribo_utils.makePickle(transcripts, lib_settings.get_transcript_counts())
    lib_settings.write_to_log('done counting reads or loading counts')


class ribo_lib:
    def __init__(self, experiment_settings, lib_settings):
        """
        Constructor for Library class
        """
        self.experiment_settings = experiment_settings
        self.lib_settings = lib_settings
        self.get_property = self.experiment_settings.get_property
        self.get_rdir = experiment_settings.get_rdir
        self.get_wdir = experiment_settings.get_wdir

        print "unpickling %s counts" % lib_settings.sample_name
        self.transcripts = ribo_utils.unPickle(self.lib_settings.get_transcript_counts())
        #this summing has to happen after library initialization and unpickling
        self.total_mapped_fragments = sum([transcript.fragment_count for transcript in self.transcripts.values()])
        self.total_mapped_counts = {}  # maps a set of parameters to the total number of read counts under those parameters

    def __repr__(self):
        return "ribo-seq lib for %s" % (self.lib_settings.sample_name)

    def name_sorted_counts(self):
        #returns counts for each sequence in pool, sorted by their sequence names, alphabetically
        return np.array([self.transcripts[sequence_name].fragment_count for sequence_name in
                         sorted(self.transcripts)])

    def name_sorted_rpms(self):
        # returns reads per million for each sequence in pool, sorted by their sequence names, alphabetically
        return (10**6)*self.name_sorted_counts()/float(self.total_mapped_fragments)

    def sorted_names(self):
        return sorted(self.transcripts.keys())

    def get_transcript(self, sequence_name):
        if sequence_name in self.transcripts:
            return self.transcripts[sequence_name]
        else:
            return None

    def get_sample_name(self):
        return self.lib_settings.sample_name

    def get_counts(self, sequence_name):
        return self.transcripts[sequence_name].fragment_count

    def get_mappings_with_minimum_reads(self, minimum_reads, names_only = False):
        passing_mappings = set()
        for mapping in self.transcripts.values():
            if mapping.get_number_rt_stops() >= minimum_reads:
                passing_mappings.add(mapping)

        if names_only:
            return set([passing_mapping.sequence_name for passing_mapping in passing_mappings])
        else:
            return passing_mappings

    def get_all_fragment_length_counts(self):
        all_length_counts = defaultdict(int)
        for pool_sequence_mapping in self.transcripts.values():
            for fragment_length in pool_sequence_mapping.length_dist:
                all_length_counts[fragment_length] += pool_sequence_mapping.length_dist[fragment_length]
        return all_length_counts

    def get_fragment_count_distribution(self):
        return [self.transcripts[sequence].fragment_count for sequence in self.transcripts]

    def get_total_tx_reads(self, start_offset, stop_offset, read_end='5p', read_lengths='all'):
        parameter_string = 'tx_%d_%d_%s_%s' %(start_offset, stop_offset, read_end, str(read_lengths))
        if parameter_string not in self.total_mapped_counts:
            self.total_mapped_counts[parameter_string] = sum([tx.get_tx_read_count(start_offset, stop_offset, read_end=read_end,
                                                                                   read_lengths=read_lengths)
                                                              for tx in self.transcripts.values()])
        return  self.total_mapped_counts[parameter_string]


    def get_total_cds_reads(self, start_offset, stop_offset, read_end='5p', read_lengths='all'):
        parameter_string = 'cds_%d_%d_%s_%s' % (start_offset, stop_offset, read_end, str(read_lengths))
        if parameter_string not in self.total_mapped_counts:
            self.total_mapped_counts[parameter_string] = sum(
                [tx.get_cds_read_count(start_offset, stop_offset, read_end=read_end,
                                      read_lengths=read_lengths)
                 for tx in self.transcripts.values()])
        return self.total_mapped_counts[parameter_string]

    def get_rpm(self, sequence_name, start_offset, stop_offset, read_end='5p', read_lengths='all'):
        tx_counts = self.transcripts[sequence_name].get_tx_read_count(start_offset, stop_offset, read_end=read_end,
                                                                      read_lengths=read_lengths)
        total_counts = self.get_total_tx_reads(start_offset, stop_offset, read_end=read_end, read_lengths=read_lengths)

        return (10 ** 6) * float(tx_counts) / float(total_counts)


    def get_rpkm(self, sequence_name, start_offset, stop_offset, read_end='5p', read_lengths='all'):
        return self.get_rpm(sequence_name, start_offset, stop_offset, read_end=read_end, read_lengths=read_lengths) /\
               self.transcripts[sequence_name].tx_length


    def get_cds_rpm(self, sequence_name, start_offset, stop_offset, read_end='5p', read_lengths='all'):
        cds_counts = self.transcripts[sequence_name].get_cds_read_count(start_offset, stop_offset, read_end=read_end,
                                                                      read_lengths=read_lengths)

        total_counts = self.get_total_cds_reads(start_offset, stop_offset, read_end=read_end, read_lengths=read_lengths)
        return (10 ** 6) * float(cds_counts) / float(total_counts)


    def get_cds_rpkm(self, sequence_name, start_offset, stop_offset, read_end='5p', read_lengths='all'):
        return self.get_cds_rpm(sequence_name, start_offset, stop_offset, read_end=read_end, read_lengths=read_lengths) / \
               (self.transcripts[sequence_name].cds_length-start_offset+stop_offset)

class transcript:
    """
    Represents a single transcript from the genome
    Stores
        The transcript sequence
        The Trimmed sequence used for mapping
        The positions of all reads mapping to this sequence
        Total number of reads mapping to this sequence
        Fraction of library reads mapping here
        Enrichment relative to input library

    """
    def __init__(self, tx_id, tx_info, tx_seq, sam_file):
        self.strand = tx_info['strand']
        self.tx_length = int(tx_info['tx_length'])
        self.exon_starts = [int(start) for start in tx_info['exon_starts'].split(',')]
        self.exon_ends = [int(end) for end in tx_info['tx_exon_ends'].split(',')]
        self.cds_start = int(tx_info['cdsStart'])
        self.cds_end = int(tx_info['cdsEnd'])
        self.cds_length = self.cds_end-self.cds_start
        self.trailer_length = self.tx_length-self.cds_end
        self.leader_length =self.cds_start
        assert self.trailer_length + self.cds_length + self.leader_length == self.tx_length
        self.gene_id = tx_info['geneID']
        self.notes = [note for note in tx_info['notes'].split(',')]

        self.sequence_name = tx_id
        self.full_sequence = tx_seq
        self.fragment_5p_ends_at_position = defaultdict(int) #will map position to # of reads there
        self.fragment_3p_ends_at_position = defaultdict(int) #will map position to # of reads there
        self.fragment_5p_lengths_at_position = defaultdict(dict)  # will map position to a list of fragment lengths with 5' ends at that position
        self.fragment_3p_lengths_at_position = defaultdict(dict) # will map position to a list of fragment lengths with 3' ends at that position
        self.fragment_count = 0
        self.length_dist = defaultdict(int)
        self.assign_read_ends_from_sam(sam_file)

    def get_first_jxn_in_CDS(self):
        for exon_start in self.exon_starts:
            if exon_start> self.cds_start:
                return exon_start
        return None

    def assign_read_ends_from_sam(self, sam_file):
        all_mapping_reads = sam_file.fetch(reference = self.sequence_name)
        reads_reversed = self.experiment_settings.get_property('reads_reversed')
        for read in [r for r in all_mapping_reads if (not r.is_secondary) and (reads_reversed == r.is_reverse)]:
            if reads_reversed:
                fragment_start = read.reference_end
                fragment_end = read.reference_start
            else:
                fragment_start = read.reference_start #0-based start of fragment
                fragment_end = read.reference_end
            # get read length from sequence, or CIGAR string if unavailable
            fragment_length = read.infer_query_length(always=False)
            assert fragment_length != 0
            self.fragment_5p_ends_at_position[fragment_start] += 1
            self.fragment_3p_ends_at_position[fragment_end] += 1
            self.fragment_count += 1
            self.length_dist[fragment_length] += 1
            if fragment_length not in self.fragment_5p_lengths_at_position[fragment_start]:
                self.fragment_5p_lengths_at_position[fragment_start][fragment_length] = 0
            self.fragment_5p_lengths_at_position[fragment_start][fragment_length] += 1
            if fragment_length not in self.fragment_3p_lengths_at_position[fragment_end]:
                self.fragment_3p_lengths_at_position[fragment_end][fragment_length] = 0
            self.fragment_3p_lengths_at_position[fragment_end][fragment_length] += 1

    def contains_subsequence(self, subsequence):
        if subsequence in self.full_sequence:
            return True
        else:
            return False

    def positions_of_subsequence(self, subsequence):
        #this regex will NOT return overlapping sequences
        return [m.start() for m in re.finditer(subsequence, self.full_sequence)]

    def get_read_end_positions(self, read_end = '5p', read_lengths = 'all'):
        '''
        :param read_end:
        :param read_lengths:
        :return:
        '''
        if read_lengths == 'all':
            if read_end == '5p':
                read_dict = self.fragment_5p_ends_at_position
            elif read_end == '3p':
                read_dict = self.fragment_3p_ends_at_position
            else:
                print 'unidentified read end option', read_end
                sys.exit()
        else:
            assert isinstance(read_lengths, list)
            read_dict = {}
            if read_end == '5p':
                super_dict = self.fragment_5p_lengths_at_position
            elif read_end == '3p':
                super_dict = self.fragment_3p_lengths_at_position
            else:
                print 'unidentified read end option', read_end
                sys.exit()
            for position in super_dict:
                read_dict[position] = sum([super_dict[position][read_length] for read_length in read_lengths if read_length in super_dict[position]])
        return read_dict

    def get_positional_length_sums(self, read_end = '5p'):
        '''
        :param read_end:
        :param read_lengths:
        :return:
        '''
        lengths_dict = {}
        if read_end == '5p':
            super_dict = self.fragment_5p_lengths_at_position
        elif read_end == '3p':
            super_dict = self.fragment_3p_lengths_at_position
        else:
            print 'unidentified read end option', read_end
            sys.exit()
        for position in super_dict:
            lengths_dict[position] = sum([super_dict[position][read_length]*read_length for read_length in
                                          super_dict[position] if read_length in super_dict[position]])
        return lengths_dict

    def get_read_counts_array(self, anchor, left_offset, right_offset, read_end ='5p', read_lengths ='all'):
        '''

        :param anchor: the "reference" or "zero" position for the array. For example the CDS start or stop
        :param left_offset: position relative to anchor to start the array, negative numbers will be upstream
        :param right_offset: position relative to anchor to end the array, negative numbers will be upstream
        :param read_end: '5p' or 3p', which end of reads to count.
        :param read_lengths: read_lengths: read lengths to include in the count. must be 'all', or an array of ints.
        :return:
        counts_array: an array of the read counts over the window
        inclusion_array: an array of 1 and zero to indicate which positions are within the transcript boundaries.
        '''
        assert right_offset > left_offset
        read_dict = self.get_read_end_positions(read_end = read_end, read_lengths = read_lengths)
        counts_array = np.array([read_dict[position] if position in read_dict and position>=0 else 0
                        for position in range(anchor + left_offset, anchor + right_offset + 1)])
        inclusion_array = np.array([1 if position>=0 and position < self.tx_length else 0
                           for position in range(anchor + left_offset, anchor + right_offset + 1)])
        assert len(counts_array) == len (inclusion_array)
        return counts_array, inclusion_array

    def get_avg_read_lengths_array(self, anchor, left_offset, right_offset, read_end ='5p'):
        '''

        :param anchor: the "reference" or "zero" position for the array. For example the CDS start or stop
        :param left_offset: position relative to anchor to start the array, negative numbers will be upstream
        :param right_offset: position relative to anchor to end the array, negative numbers will be upstream
        :param read_end: '5p' or 3p', which end of reads to count.
        :param read_lengths: read_lengths: read lengths to include in the count. must be 'all', or an array of ints.
        :return:
        counts_array: an array of the read counts over the window
        length_sum_array: an array of the total length of reads at each position
        '''
        assert right_offset > left_offset
        collapsed_read_dict = self.get_read_end_positions(read_end = read_end, read_lengths = 'all')
        read_length_sums_dict = self.get_positional_length_sums(read_end = read_end)
        counts_array = np.array([collapsed_read_dict[position] if position in collapsed_read_dict and position>=0 else 0
                        for position in range(anchor + left_offset, anchor + right_offset + 1)])
        length_sum_array = np.array([read_length_sums_dict[position] if position in read_length_sums_dict and position>=0 else 0
                        for position in range(anchor + left_offset, anchor + right_offset + 1)])
        assert len(counts_array) == len (length_sum_array)
        return length_sum_array, counts_array

    def get_cds_read_count(self, start_offset, stop_offset, read_end='5p', read_lengths='all'):
        '''

        :param start_offset: position relative to CDS start to open the count window, negative numbers will be upstream
        :param stop_offset: position relative to CDS stop to close the count window, negative numbers will be upstream
        :param read_end: '5p' or 3p', which end of reads to count.
        :param read_lengths: read lengths to include in the count. must be 'all', or an array of ints.
        :return:
        '''
        read_dict = self.get_read_end_positions(read_end=read_end, read_lengths=read_lengths)
        return sum([read_dict[position] for position in read_dict
                    if position>=self.cds_start+start_offset and position<=self.cds_end+stop_offset])

    def get_tx_read_count(self, start_offset, stop_offset, read_end='5p', read_lengths='all'):
        '''

        :param start_offset: position relative to tx start to open the count window, negative numbers will be upstream
        :param stop_offset: position relative to tx end to close the count window, negative numbers will be upstream
        :param read_end: '5p' or 3p', which end of reads to count.
        :param read_lengths: read lengths to include in the count. must be 'all', or an array of ints.
        :return:
        '''
        read_dict = self.get_read_end_positions(read_end=read_end, read_lengths=read_lengths)
        return sum([read_dict[position] for position in read_dict
                    if position>=start_offset and position<=self.tx_length+stop_offset])

    def get_read_frame_counts(self, start, end, read_end='5p', read_lengths='all'):
        '''
        :return the number of reads in each CDS frame. the first nucleotide of the  start codon is frame zero.
        :param start: tx position (zero-indexed, inclusive) to start counting at
        :param stop: tx position (zero-indexed, inclusive) to stop counting at
        :param read_end: '5p' or 3p', which end of reads to count.
        :param read_lengths: read lengths to include in the count. must be 'all', or an array of ints.
        :return:
        '''
        read_dict = self.get_read_end_positions(read_end=read_end, read_lengths=read_lengths)
        frame_counts = np.zeros(3)
        for position in range(start, end+1):
            if position >= 0 and position <= self.tx_length and position in read_dict:
                cds_rel_position = position - self.cds_start
                frame_counts[cds_rel_position%3] += read_dict[position]

        return frame_counts

    def second_stop_position(self):
        #find the position of the first nt of the first in-frame stop after the canonical stop codon,
        # relative to the start of the tx

        stop_codon = self.full_sequence[self.cds_end-3: self.cds_end]

        if ribo_utils.GENETIC_CODE[stop_codon] == '_':
            for position in range(self.cds_end, self.tx_length, 3):
                codon = self.full_sequence[position: position+3]
                if codon in ribo_utils.GENETIC_CODE and ribo_utils.GENETIC_CODE[codon] == '_':
                    return position
        else:
            #pass
            print self.gene_id, self.strand, self.cds_end, self.full_sequence[self.cds_end-3: self.cds_end]
        return None

    def second_stop_codon(self):
        #return the triplet code of the second stop codon
        position = self.second_stop_position()
        if  position == None:
            return None
        else:
            return self.full_sequence[position: position+3]
    def readthrough_extension_length(self, pre_extension_stop_buffer=0, post_cds_stop_buffer=0):
        if self.second_stop_codon() == None:
            return None
        return float((self.second_stop_position()+3-pre_extension_stop_buffer)-(self.cds_end+post_cds_stop_buffer))

    def get_readthrough_counts(self, p_offset=0, read_end='5p', read_lengths='all', pre_extension_stop_buffer=0, post_cds_stop_buffer=0):
        read_dict = self.get_read_end_positions(read_end=read_end, read_lengths=read_lengths)
        second_stop = self.second_stop_position()+3
        return sum([read_dict[position] for position in read_dict if position>=self.cds_end+(post_cds_stop_buffer-p_offset) and position<=second_stop-pre_extension_stop_buffer-p_offset])



    def compute_readthrough_ratio(self, p_offset, read_end='5p', read_lengths='all', cds_cutoff=128, log=True,
                                  post_cds_start_buffer=12, pre_cds_stop_buffer = 15, pre_extension_stop_buffer=15,
                                  post_cds_stop_buffer=9):
        """

        :param p_offset: distance of p-site from designated read end, positive for 5p, negative for 3p
        :param read_end:
        :param read_lengths:
        :param cds_cutoff:
        :param log:
        :param post_cds_start_buffer: omit this many nucleotides at edge from counting
        :param pre_cds_stop_buffer: omit this many nucleotides at edge from counting
        :param pre_extension_stop_buffer: omit this many nucleotides at edge from counting
        :param post_cds_stop_buffer: omit this many nucleotides at edge from counting
        :return:
        """
        cds_counts = self.get_cds_read_count(post_cds_start_buffer-p_offset, (-1*pre_cds_stop_buffer)-p_offset, read_end=read_end,
                                                   read_lengths=read_lengths)
        cds_read_density =  cds_counts/float(self.cds_length)
        second_stop = self.second_stop_position()
        if not second_stop == None and cds_counts>=cds_cutoff:
            ex_length = self.readthrough_extension_length(pre_extension_stop_buffer=pre_extension_stop_buffer,
                                              post_cds_stop_buffer=post_cds_stop_buffer)
            if ex_length > 0:
                second_cds_density = self.get_readthrough_counts(p_offset=p_offset, read_end=read_end, read_lengths=read_lengths,
                                                                 pre_extension_stop_buffer=pre_extension_stop_buffer,
                                                                 post_cds_stop_buffer=post_cds_stop_buffer)/ex_length

                readthrough = second_cds_density/cds_read_density
                if log:
                    if readthrough>0:
                        return math.log(readthrough, 10)
                    else:
                        return None
                else:
                    return readthrough
        return None

    def stop_codon_context(self):
        #return the canonical stop codon sequence

        stop_codon = self.full_sequence[self.cds_end-3: self.cds_end]

        #assert ribo_utils.GENETIC_CODE[stop_codon] == '_'
        return self.full_sequence[self.cds_end-15: self.cds_end+12]

    def trailer_sequence(self):
        #return sequence of 3' trailer (3' UTR)
        return self.full_sequence[self.cds_end: self.tx_length]

    def trailer_monomer_fraction(self, monomer):
        assert monomer.upper() in 'ATCG'
        trailer_seq = self.trailer_sequence().upper()
        fraction = float(trailer_seq.count(monomer))/float(len(trailer_seq))
        return fraction

