from collections import defaultdict
import re

import pysam
import ribo_utils
import numpy as np
import sys

def initialize_pool_sequence_mappings(experiment_settings, lib_settings):
    lib_settings.write_to_log('counting reads or loading counts')
    if experiment_settings.get_property('force_recount') or not lib_settings.sequence_counts_exist():
        print "counting BAM reads"
        transcripts = {}
        tx_annotations = ribo_utils.tsv_to_dict(experiment_settings.get_property('canonical_tx_features'))
        tx_seqs = ribo_utils.convertFastaToDict(experiment_settings.get_property('canonical_tx_seqs'))
        samfile = pysam.AlignmentFile(lib_settings.get_transcript_mapped_reads(), "rb")
        for tx_id in tx_annotations:
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
        self.log_mapping_tags()

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

    def get_rpm(self, sequence_name):
        return (10**6)*self.transcripts(sequence_name).fragment_count/float(self.total_mapped_fragments)

    def get_rpkm(self, sequence_name):
        pass
        #TODO: implement RPKM computation
        #return (10**6)*self.transcripts(sequence_name).fragment_count/float(self.total_mapped_fragments)

    def log_mapping_tags(self):
        "writes aggregate paired-end mapping tag data to log file as a sanity check"
        aggregate_dict = defaultdict(int)
        for pool_seq in self.transcripts.values():
            for mapping_tag in pool_seq.paired_end_mapping_tags:
                aggregate_dict[mapping_tag] += pool_seq.paired_end_mapping_tags[mapping_tag]
        self.lib_settings.write_to_log('mapping tags- %s' % (', '.join(['%s:%d' % (tag, aggregate_dict[tag]) for
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
    def __init__(self, tx_id, tx_info, tx_seq, sam_file):
        #TODO: collect info from tx_info
        self.tx_length = None
        self.exon_starts = []
        self.exon_ends = []
        self.cds_start = None
        self.cds_end = None

        self.sequence_name = tx_id
        self.full_sequence = tx_seq
        self.fragment_5p_ends_at_position = defaultdict(int) #will map position to # of reads there
        self.fragment_3p_ends_at_position = defaultdict(int) #will map position to # of reads there
        self.fragment_5p_lengths_at_position = defaultdict(dict)  # will map position to a list of fragment lengths with 5' ends at that position
        self.fragment_3p_lengths_at_position = defaultdict(dict)
        self.fragment_count = 0
        self.length_dist = defaultdict(dict)

        self.assign_read_ends_from_sam(sam_file)

    def assign_read_ends_from_sam(self, sam_file):
        all_mapping_reads = sam_file.fetch(reference = self.sequence_name)
        for read in all_mapping_reads:
            if read.is_read1:
                #read1 should be on the forawrd strand since these are transcript-relative
                #this alignment should be the primary one. IF this throws erros, I will need to write more logic.
                assert (not read.is_reverse) and (not read.is_secondary)
                fragment_start = read.reference_start #0-based start of fragment
                fragment_length = read.template_length
                fragment_end = fragment_start + fragment_length
                self.fragment_5p_ends_at_position[fragment_start] += 1
                self.fragment_3p_ends_at_position[fragment_end] += 1
                self.fragment_count += 1
                self.length_dist[fragment_length] += 1
                if fragment_length not in self.fragment_5p_lengths_at_position[fragment_start]:
                    self.fragment_5p_lengths_at_position[fragment_start][fragment_length]+=1
                if fragment_length not in self.fragment_3p_lengths_at_position[fragment_start]:
                    self.fragment_3p_lengths_at_position[fragment_end][fragment_length]+=1
    def contains_subsequence(self, subsequence):
        if subsequence in self.full_sequence:
            return True
        else:
            return False

    def positions_of_subsequence(self, subsequence):
        #this regex will NOT return overlapping sequences
        return [m.start() for m in re.finditer(subsequence, self.full_sequence)]

    def get_read_counts_array(self, relative_start, upstream, downstream, read_end = '5p', read_lengths = 'all'):
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
                read_dict[position] = sum([super_dict[position][read_length] for read_length in read_lengths])
        counts_array = [read_dict[position] if position in read_dict and position>=0 else 0
                        for position in range(relative_start-upstream, relative_start+downstream+1)]
        inclusion_array = [1 if position>=0 and position < self.tx_length else 0
                        for position in range(relative_start-upstream, relative_start+downstream+1)]
        return  counts_array, inclusion_array