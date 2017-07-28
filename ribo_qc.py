from collections import defaultdict
import matplotlib.pyplot as plt

import uniform_colormaps

plt.rcParams['pdf.fonttype'] = 42
import scipy.stats as stats
import os
import ribo_utils
import numpy as np
import itertools
import math
import gzip
import pysam

class ribo_qc:
    def __init__(self, experiment, experiment_settings, threads):
        """
        Constructor for Library class
        """
        self.threads = threads
        self.experiment = experiment
        self.experiment_settings = experiment_settings
        self.experiment_settings.write_to_log('initializing QC engine')
        self.get_property = self.experiment_settings.get_property
        self.get_rdir = experiment_settings.get_rdir
        self.genome = ribo_utils.genome_sequence(self.experiment_settings.get_genome_sequence_files())
        self.GTF_annotations = ribo_utils.gtf_data(self.experiment_settings.get_annotation_GTF_file())
        self.lib_QCs = [single_lib_qc(self, lib_settings) for lib_settings in self.experiment_settings.iter_lib_settings()]
        ribo_utils.make_dir(self.experiment.rdir_path('QC'))

    def write_mapping_summary(self, output_file):

        f = open(output_file, 'w')

        f.write('sample name\ttotal reads\tread1 adaptors\tread2 adaptors\ttoo short\tpass trimming filter\tmapping input\tunaligned pairs\tuniquely aligned pairs\tmultiply aligned pairs\ttotal alignment %\n')
        for lib_settings in self.experiment_settings.iter_lib_settings():
            mapping_input, paired_reads, unaligned_pairs, uniquely_aligned_pairs, multiply_aligned_pairs,\
            overall_alignment_percent = self.parse_paired_end_mapping_stats(lib_settings.get_pool_mapping_stats())
            total_reads, read1_adaptors, read2_adaptors, too_short, passing_filter = self.get_adaptor_trimming_stats(lib_settings.get_log())
            f.write('%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n' % (lib_settings.sample_name, total_reads, read1_adaptors,
                                                      read2_adaptors, too_short, passing_filter, mapping_input,
                                                      unaligned_pairs, uniquely_aligned_pairs,
                                                      multiply_aligned_pairs, overall_alignment_percent))
            f.write('%s percents\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n' % (lib_settings.sample_name,
                                                                               100*total_reads/float(total_reads),
                                                                               100*read1_adaptors/float(total_reads),
                                                                               100*read2_adaptors/float(total_reads),
                                                                               100*too_short/float(total_reads),
                                                                               100*passing_filter/float(total_reads),
                                                                               100*mapping_input/float(total_reads),
                                                                               100*unaligned_pairs/float(total_reads),
                                                                               100*uniquely_aligned_pairs/float(total_reads),
                                                                               100*multiply_aligned_pairs/float(total_reads),
                                                                               100*(uniquely_aligned_pairs+multiply_aligned_pairs)/float(total_reads)))
        f.close()

    def write_trimming_stats_summary(self):
        out_file = open(self.experiment_settings.get_trimming_count_summary(), 'w')
        out_percent_file = open(self.experiment_settings.get_trimming_percent_summary(), 'w')

        header_list = ['trimming stat'] + [single_lib_qc.lib_settings.sample_name for single_lib_qc in self.lib_QCs]
        header_line = '%s\n' % ('\t'.join(header_list))
        out_file.write(header_line)
        out_percent_file.write(header_line)
        for stat_name in ['reads processed','reads filtered out by quality control',
                          'short reads filtered out after trimming by size control',
                          'empty reads filtered out after trimming by size control',
                          'reads available', 'trimmed reads available']:
            out_list = [stat_name] + [str(single_lib_qc.adaptor_stats[stat_name]) for single_lib_qc in self.lib_QCs]
            out_percent_list = [stat_name] +\
                               [str(100.*single_lib_qc.adaptor_stats[
                                        stat_name]/float(single_lib_qc.adaptor_stats['reads processed']))
                                for single_lib_qc in self.lib_QCs]
            out_line = '%s\n' % ('\t'.join(out_list))
            out_file.write(out_line)
            out_percent_line = '%s\n' % ('\t'.join(out_percent_list))
            out_percent_file.write(out_percent_line)

        out_file.close()
        out_percent_file.close()

class single_lib_qc():
    def __init__(self, parent_qc, lib_settings):
        """
        Constructor for Library class
        """
        self.parent_qc = parent_qc
        self.lib_settings = lib_settings
        self.get_adaptor_trimming_stats()
        #self.samfile = pysam.AlignmentFile(lib_settings.get_genome_mapped_reads(), "rb")
        #self.get_mapping_multiplicity_stats()
        #self.samfile.close()
    def get_adaptor_trimming_stats(self):
        """
        Parses the outpur from skewer to consolidate stats
        :param lib_settings: 
        :return: 
        """
        self.lib_settings.write_to_log('parsing adpator trimming stats in %s' % (self.lib_settings.get_adaptor_trimming_log()))
        self.adaptor_stats = {}
        trimming_log = open(self.lib_settings.get_adaptor_trimming_log())
        for line in trimming_log:
            if "reads processed; of these:" in line:
                self.adaptor_stats['reads processed'] = int(line.strip().split()[0])
                line = trimming_log.next()
                self.adaptor_stats['reads filtered out by quality control'] = int(line.strip().split()[0])
                line = trimming_log.next()
                self.adaptor_stats['short reads filtered out after trimming by size control'] = int(line.strip().split()[0])
                line = trimming_log.next()
                self.adaptor_stats['empty reads filtered out after trimming by size control'] = int(line.strip().split()[0])
                line = trimming_log.next()
                self.adaptor_stats['reads available'] = int(line.strip().split()[0])
                line = trimming_log.next()
                self.adaptor_stats['trimmed reads available'] = int(line.strip().split()[0])
            if "Length distribution of reads after trimming:" in line:
                self.adaptor_stats['post_trimming_lengths'] = {}
                line = trimming_log.next()
                while True:
                    try:
                        line = trimming_log.next()
                        ll = line.strip().split()
                        length = int(ll[0])
                        count = int(ll[1])
                        self.adaptor_stats['post_trimming_lengths'][length] = count
                    except:
                        break
        trimming_log.close()
        self.lib_settings.write_to_log('done parsing adpator trimming stats')

    def get_mapping_multiplicity_stats(self):
        #TODO: remove print statements and have this written to a summary file
        total_alignments = 0
        primary_alignments = 0
        secondary_alignments = 0
        unmapped_alignments = 0
        multiplicity = defaultdict(int)
        for alignment in self.samfile.fetch():
            if alignment.is_secondary:
                secondary_alignments += 1
            else:
                primary_alignments += 1
                multiplicity[int(alignment.get_tag('NH:i'))] += 1
            if alignment.is_unmapped:
                unmapped_alignments += 1
            total_alignments += 1
        print self.lib_settings.sample_name
        print 'total: ', total_alignments
        print 'primary: ', primary_alignments
        print 'unique mapping: ', multiplicity[1]
        multimapping = sum([multiplicity[mult] for mult in multiplicity.keys() if mult !=1])
        print 'multiply mapping: ', multimapping
        print 'secondary: ', secondary_alignments
        print 'unmapped: ', unmapped_alignments


        #for mult in sorted(multiplicity.keys()):
        #    print '%d: %d' % (mult, multiplicity[mult])

