from collections import defaultdict
import re

import pysam
import ribo_utils
import numpy as np


def initialize_pool_sequence_mappings(experiment_settings, lib_settings):
    lib_settings.write_to_log('counting reads or loading counts')
    if experiment_settings.get_property('force_recount') or not lib_settings.sequence_counts_exist():
        print "counting BAM reads"
        transcripts = {}
        tx_annotations = ribo_utils.tsv_to_dict(experiment_settings.get_property('canonical_tx_features'))
        tx_seqs = ribo_utils.convertFastaToDict(experiment_settings.get_property('canonical_tx_seqs'))
        samfile = pysam.AlignmentFile(lib_settings.get_transcript_mapped_reads(), "rb")
        for tx_id in tx_annotations:
            transcripts[tx_id] = transcript(tx_id, tx_annotations[tx_id], tx_seqs[tx_id])

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
        self.pool_sequence_mappings = ribo_utils.unPickle(self.lib_settings.get_transcript_counts())
        #this summing has to happen after library initialization and unpickling
        self.total_mapped_fragments = sum([mapping.fragment_count for mapping in self.pool_sequence_mappings.values()])
        self.log_mapping_tags()

    def name_sorted_counts(self):
        #returns counts for each sequence in pool, sorted by their sequence names, alphabetically
        return np.array([self.pool_sequence_mappings[sequence_name].fragment_count for sequence_name in
                         sorted(self.pool_sequence_mappings)])

    def name_sorted_rpms(self):
        # returns reads per million for each sequence in pool, sorted by their sequence names, alphabetically
        return (10**6)*self.name_sorted_counts()/float(self.total_mapped_fragments)

    def sorted_names(self):
        return sorted(self.pool_sequence_mappings.keys())

    def get_single_TL_mappings(self, names_only = False):
        single_TL_mappings = set()
        single_TL_names = set()
        for mapping_name in self.pool_sequence_mappings:
            if self.pool_sequence_mappings[mapping_name].is_only_tl:
                single_TL_mappings.add(self.pool_sequence_mappings[mapping_name])
                single_TL_names.add(mapping_name)

        if names_only:
            return single_TL_names
        else:
            return single_TL_mappings

    def compute_lib_fractions(self):
        total_passing_library_reads = float(sum([pool_sequence_mapping.total_passing_reads for pool_sequence_mapping in self.pool_sequence_mappings.values()]))
        for pool_sequence_mapping in self.pool_sequence_mappings.values():
            pool_sequence_mapping.lib_fraction = pool_sequence_mapping.total_passing_reads/total_passing_library_reads

    def calculate_enrichments(self, input_lib):
        for pool_sequence_mapping in self.pool_sequence_mappings.values():
            input_pool_sequence_mapping = input_lib.get_pool_sequence_mapping(pool_sequence_mapping.sequence_name)
            assert input_pool_sequence_mapping != None
            try:
                pool_sequence_mapping.enrichment = pool_sequence_mapping.lib_fraction/input_pool_sequence_mapping.lib_fraction
            except:
                pool_sequence_mapping.enrichment = 0

    def get_pool_sequence_mapping(self, sequence_name):
        if sequence_name in self.pool_sequence_mappings:
            return self.pool_sequence_mappings[sequence_name]
        else:
            return None

    def get_sample_name(self):
        return self.lib_settings.sample_name

    def get_counts(self, sequence_name):
        return self.pool_sequence_mappings[sequence_name].fragment_count

    def get_mappings_with_minimum_reads(self, minimum_reads, names_only = False):
        passing_mappings = set()
        for mapping in self.pool_sequence_mappings.values():
            if mapping.get_number_rt_stops() >= minimum_reads:
                passing_mappings.add(mapping)

        if names_only:
            return set([passing_mapping.sequence_name for passing_mapping in passing_mappings])
        else:
            return passing_mappings

    def get_all_fragment_length_counts(self):
        all_length_counts = defaultdict(int)
        for pool_sequence_mapping in self.pool_sequence_mappings.values():
            for fragment_length in pool_sequence_mapping.length_dist:
                all_length_counts[fragment_length] += pool_sequence_mapping.length_dist[fragment_length]
        return all_length_counts

    def get_fragment_count_distribution(self):
        return [self.pool_sequence_mappings[sequence].fragment_count for sequence in self.pool_sequence_mappings]

    def get_rpm(self, sequence_name):
        return (10**6)*self.get_pool_sequence_mapping(sequence_name).fragment_count/float(self.total_mapped_fragments)
    def log_mapping_tags(self):
        "writes aggregate paired-end mapping tag data to log file as a sanity check"
        aggregate_dict = defaultdict(int)
        for pool_seq in self.pool_sequence_mappings.values():
            for mapping_tag in pool_seq.paired_end_mapping_tags:
                aggregate_dict[mapping_tag] += pool_seq.paired_end_mapping_tags[mapping_tag]
        self.lib_settings.write_to_log('PE mapping tags- %s' % (', '.join(['%s:%d' % (tag, aggregate_dict[tag]) for
                                                                           tag in aggregate_dict])))

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
    def __init__(self, sequence_name, full_sequence, sam_file):
        self.sequence_name = sequence_name
        self.full_sequence = full_sequence
        self.fragment_5p_ends_at_position = defaultdict(int) #will map position to # of reads there
        self.fragment_3p_ends_at_position = defaultdict(int) #will map position to # of reads there
        self.fragment_lengths_at_position = defaultdict(list)  # will map position to a list of fragment lengths with 5' ends at that position
        #self.fragment_lengths = []
        self.paired_end_mapping_tags = defaultdict(int)
        self.fragment_count = 0
        self.length_dist = defaultdict(int)

        self.assign_read_ends_from_sam(sam_file)

    def assign_read_ends_from_sam(self, sam_file):
        all_mapping_reads = sam_file.fetch(reference = self.sequence_name)
        for read in all_mapping_reads:
            # only need to look at the forward mapping read in each pair, since the necessary mate info is there
            if read.is_read1:
                #read1 should be on the forawrd strans since --norc should be specified
                #this alignment should be the primary one. IF this throws erros, I will need to write more logic.
                assert (not read.is_reverse) and (not read.is_secondary)
                pe_mapping_tag = read.get_tag('YT') #this should always return 'CP' for a concordantly-mapped pair
                self.paired_end_mapping_tags[pe_mapping_tag]+=1
                fragment_start = read.reference_start #0-based start of fragment
                fragment_length = read.template_length
                fragment_end = fragment_start + fragment_length
                self.fragment_5p_ends_at_position[fragment_start] += 1
                self.fragment_3p_ends_at_position[fragment_end] += 1
                self.fragment_count += 1
                self.length_dist[fragment_length] += 1
                #self.fragment_lengths_at_position[fragment_start].append(fragment_length)
                #self.fragment_lengths.append(fragment_length)

    def contains_subsequence(self, subsequence):
        if subsequence in self.full_sequence:
            return True
        else:
            return False

    def positions_of_subsequence(self, subsequence):
        #this regex will NOT return overlapping sequences
        return [m.start() for m in re.finditer(subsequence, self.full_sequence)]

    def fraction_at_position(self, position):
        if position < 0 or position > len(self.full_sequence)-1:
            return None
        else:
            #return self.reads_at_position[position]/float(self.total_passing_reads)
            if self.get_number_rt_stops() == 0:
                return 0
            else:
                return self.fragment_5p_ends_at_position[position] / self.get_number_rt_stops()