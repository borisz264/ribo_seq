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

class ribo_qc:
    def __init__(self, experiment, experiment_settings, threads):
        """
        Constructor for Library class
        """
        self.threads = threads
        self.experiment = experiment
        self.experiment_settings = experiment_settings
        self.get_property = self.experiment_settings.get_property
        self.get_rdir = experiment_settings.get_rdir
        self.get_wdir = experiment_settings.get_wdir
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

    def get_adaptor_trimming_stats(self, log_file):

        return