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

    def parse_paired_end_mapping_stats(self, alignment_summary_file):
        '''
        example alignment summary:

        random stuff here

        100000 reads; of these:
          100000 (100.00%) were paired; of these:
            37476 (37.48%) aligned concordantly 0 times
            60871 (60.87%) aligned concordantly exactly 1 time
            1653 (1.65%) aligned concordantly >1 times
            ----
            37476 pairs aligned 0 times concordantly or discordantly; of these:
              74952 mates make up the pairs; of these:
                66796 (89.12%) aligned 0 times
                1106 (1.48%) aligned exactly 1 time
                7050 (9.41%) aligned >1 times
        66.60% overall alignment rate

        more stuff here
        '''
        f = open(alignment_summary_file)
        for line in f:
            if line.strip().endswith('reads; of these:'):
                total_reads = int(line.strip().split()[0])
                line=f.next()
                paired_reads = int(line.strip().split()[0])
                line =f.next()
                unaligned_pairs = int(line.strip().split()[0])
                line =f.next()
                uniquely_aligned_pairs = int(line.strip().split()[0])
                line =f.next()
                multiply_aligned_pairs = int(line.strip().split()[0])
                line =f.next()
                overall_alignment_percent = float(line.strip().split()[0][:-1])
        f.close()
        return total_reads, paired_reads, unaligned_pairs, uniquely_aligned_pairs, multiply_aligned_pairs, overall_alignment_percent

    def get_adaptor_trimming_stats(self, log_file):
        '''
        example:
        This is cutadapt 1.10 with Python 2.7.10
        Command line parameters: -a GCTGCACGGTGACGTCTCNNNTGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC --overlap 5 -u 7 -U 21 -q 10 --trim-n --minimum-length 30 --pair-filter=both -o /Users/boris/Gilbert_Lab/Book_4/4.146/80S_seq_pilot_pe/adaptor_removed/80S_1_1.fastq.gz -p /Users/boris/Gilbert_Lab/Book_4/4.146/80S_seq_pilot_pe/adaptor_removed/80S_1_2.fastq.gz /Users/boris/Gilbert_Lab/Book_4/4.146/truncated_pe_data/160519Gil_D16-4658_1_sequence.fastq.gz /Users/boris/Gilbert_Lab/Book_4/4.146/truncated_pe_data/160519Gil_D16-4658_2_sequence.fastq.gz
        Trimming 1 adapter with at most 10.0% errors in paired-end mode ...
        Finished in 10.02 s (40 us/read; 1.50 M reads/minute).

        === Summary ===

        Total read pairs processed:            250,000
          Read 1 with adapter:                  60,994 (24.4%)
          Read 2 with adapter:                       0 (0.0%)
        Pairs that were too short:              60,994 (24.4%)
        Pairs written (passing filters):       189,006 (75.6%)

        Total basepairs processed:    20,000,000 bp
          Read 1:    10,000,000 bp
          Read 2:    10,000,000 bp
        Quality-trimmed:                   3,840 bp (0.0%)
          Read 1:             7 bp
          Read 2:         3,833 bp
        Total written (filtered):      9,825,137 bp (49.1%)
          Read 1:     6,237,184 bp
          Read 2:     3,587,953 bp

        === First read: Adapter 1 ===

        Sequence: GCTGCACGGTGACGTCTCNNNTGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC; Type: regular 3'; Length: 54; Trimmed: 60994 times.

        No. of allowed errors:
        0-9 bp: 0; 10-19 bp: 1; 20-29 bp: 2; 30-39 bp: 3; 40-49 bp: 4; 50-54 bp: 5

        Bases preceding removed adapters:
          A: 0.1%
          C: 0.2%
          G: 0.1%
          T: 0.2%
          none/other: 99.3%

        Overview of removed sequences
        length	count	expect	max.err	error counts
        5	21	244.1	0	21
        6	6	61.0	0	6
        7	18	15.3	0	18
        8	10	3.8	0	10
        9	11	1.0	0	11
        10	5	0.2	1	5
        11	8	0.1	1	7 1
        12	12	0.0	1	12
        13	8	0.0	1	8
        14	9	0.0	1	9
        15	9	0.0	1	8 1
        16	6	0.0	1	6
        17	9	0.0	1	8 1
        18	11	0.0	1	9 1 1
        19	13	0.0	1	13
        20	7	0.0	2	7
        21	15	0.0	2	14 0 1
        22	19	0.0	2	16 3
        23	10	0.0	2	8 1 1
        24	18	0.0	2	15 2 1
        25	23	0.0	2	18 2 3
        26	23	0.0	2	22 1
        27	32	0.0	2	29 3
        28	14	0.0	2	14
        29	21	0.0	2	20 1
        30	14	0.0	3	9 3 2
        31	31	0.0	3	12 1 16 2
        32	23	0.0	3	21 1 1
        33	60588	0.0	3	102 52134 7310 1042
        '''
        f = open(log_file)
        for line in f:
            if line.strip().startswith('Total read pairs processed:'):
                total_reads = int(line.strip().split()[-1].replace(',', ''))
                line = f.next()
                read1_adaptors = int(line.strip().split()[-2].replace(',', ''))
                line = f.next()
                read2_adaptors = int(line.strip().split()[-2].replace(',', ''))
                line = f.next()
                too_short = int(line.strip().split()[-2].replace(',', ''))
                line = f.next()
                passing_filter = int(line.strip().split()[-2].replace(',', ''))
        f.close()
        return total_reads, read1_adaptors, read2_adaptors, too_short, passing_filter


    def get_collapsed_read_fractions(self, lib_settings):
        out_name =  os.path.join(
          self.experiment_settings.get_rdir(),
          'QC','collapsed_fracs',
          '%(sample_name)s.collapsed_read_fractions.pkl' % {'sample_name': lib_settings.sample_name})
        if not ribo_utils.file_exists(out_name) and not self.experiment_settings.get_property('force_recollapse'):
            collapsed_reads_file = lib_settings.get_collapsed_reads()
            read_counts = []
            f = gzip.open(collapsed_reads_file)
            for line in f:
                if not line.strip() == '' and not line.startswith('#'):#ignore empty lines and commented out lines
                    if line.startswith('>'):#> marks the start of a new sequence
                        num_reads = int(line[1:].strip().split('-')[1])
                        read_counts.append(num_reads)
                    else:
                        continue
            f.close()
            read_fractions = np.array(read_counts)/float(sum(read_counts))
            ribo_utils.makePickle(read_fractions, out_name)
        else:
            read_fractions = ribo_utils.unPickle(out_name)

        return (lib_settings.sample_name, read_fractions)

    def get_library_enrichment_correlation(self, lib1, lib2):
        lib1_enrichments = []
        lib2_enrichments = []
        for sequence in lib1.pool_sequence_mappings:
            lib1_enrichments.append(lib1.pool_sequence_mappings[sequence].enrichment)
            lib2_enrichments.append(lib2.pool_sequence_mappings[sequence].enrichment)
        spearmanR, spearmanP = stats.spearmanr(lib1_enrichments, lib2_enrichments)
        pearsonR, pearsonP = stats.pearsonr(lib1_enrichments, lib2_enrichments)
        return pearsonR, spearmanR, pearsonP, spearmanP

    def get_library_count_correlation(self, lib1, lib2):
        lib1_counts = []
        lib2_counts = []
        for sequence_name in lib1.pool_sequence_mappings:
            lib1_counts.append(lib1.pool_sequence_mappings[sequence_name].fragment_count)
            lib2_counts.append(lib2.pool_sequence_mappings[sequence_name].fragment_count)
        spearmanR, spearmanP = stats.spearmanr(lib1_counts, lib2_counts)
        pearsonR, pearsonP = stats.pearsonr(lib1_counts, lib2_counts)
        return pearsonR, spearmanR, pearsonP, spearmanP

    def print_library_count_concordances(self):
        out_name =  os.path.join(self.experiment_settings.get_rdir(), 'QC',
          'count_concordances.txt')
        f = open(out_name, 'w')
        header = 'sample1\tsample2\tpearson r\t pearson p\t spearman r\t spearman p\n'
        f.write(header)
        for libi, libj in itertools.combinations(self.experiment.libs, 2):
            pearsonR, spearmanR, pearsonP, spearmanP = self.get_library_count_correlation(libi, libj)
            line = '%s\t%s\t%f\t%f\t%f\t%f\n' % (libi.get_sample_name(), libj.get_sample_name(),
                                                         pearsonR, pearsonP, spearmanR, spearmanP)
            f.write(line)
        f.close()


    def plot_average_read_positions(self):
        for lib in self.experiment.libs:
            self.plot_average_read_positions_one_lib(lib)

    def plot_average_read_positions_one_lib(self, lib, min_x = 0, max_x = 150):
        positions = np.array(range(min_x, max_x+1))
        counts5p = np.array([sum([pool_sequence_mapping.fragment_5p_ends_at_position[position] for pool_sequence_mapping
                         in lib.pool_sequence_mappings.values()]) for position in positions])
        counts3p = np.array([sum([pool_sequence_mapping.fragment_3p_ends_at_position[position] for pool_sequence_mapping
                         in lib.pool_sequence_mappings.values()]) for position in positions])
        assert sum(counts3p) == sum(counts3p)
        averages5p = counts5p/float(sum(counts5p))
        averages3p = counts3p/float(sum(counts3p))
        fig = plt.figure(figsize=(8,8))
        plot = fig.add_subplot(111)
        plot.bar(positions, averages5p,color=ribo_utils.rainbow[0], lw=0, label = "5' end")
        plot.bar(positions, averages3p, color=ribo_utils.rainbow[1], lw=0, label = "3' end")
        plot.set_xticks(positions[::10]+0.5)
        plot.set_xticklabels(positions[::10])
        lg = plt.legend(loc=2, prop={'size': 12}, labelspacing=0.2)
        lg.draw_frame(False)
        plot.set_xlabel("position of fragment end from RNA 5' end")
        plot.set_ylabel("average read fraction")
        plot.set_title(lib.lib_settings.sample_name)
        plot.set_ylim(0, 1)
        out_name =  os.path.join(
          self.experiment_settings.get_rdir(),
          'QC',
          '%(sample_name)s.read_positions.pdf' % {'sample_name': lib.get_sample_name ()})
        plt.savefig(out_name, transparent='True', format='pdf')
        plt.clf()

    def plot_fragment_length_distributions(self):
        fig = plt.figure(figsize=(8, 8))
        plot = fig.add_subplot(111)
        colormap = uniform_colormaps.viridis
        color_index = 0
        for lib in self.experiment.libs:
            sample_name = lib.lib_settings.sample_name
            fragment_length_counts = lib.get_all_fragment_length_counts()
            bins = range(0, max(fragment_length_counts.keys()))
            frangment_length_fractions = np.array([fragment_length_counts[length] for length in bins])/float(lib.total_mapped_fragments)
            # note that all but the last bin exclude the right (larger) edge of the bin. So I add an extra bin.
            plot.plot(bins, frangment_length_fractions, color=colormap((color_index-1)/float(len(self.experiment.libs))), lw=2, label = sample_name)
            color_index += 1
            plot.set_xlabel("fragment length", fontsize=14)
            plot.set_ylabel("fraction of fragments", fontsize=14)
            #plot.set_ylim(0,0.5)
        lg = plt.legend(loc=2, prop={'size': 12}, labelspacing=0.2)
        lg.draw_frame(False)
        out_name = os.path.join(
            self.experiment_settings.get_rdir(),
            'QC',
            'fragment_length_distributions.pdf')
        plt.savefig(out_name, transparent='True', format='pdf')
        plt.clf()


    def plot_count_distributions(self):
        num_libs = len(self.experiment.libs)
        fig = plt.figure(figsize=(16,16))
        plot_index = 1
        cutoff = 64
        log_bins = [0]+[2**x for x in range(26)]
        #hbins = np.arange(0, 4000, 10)
        #hbins = np.append(hbins, 10000000)
        for lib in self.experiment.libs:
            plot = fig.add_subplot(math.sqrt(ribo_utils.next_square_number(num_libs)), math.sqrt(ribo_utils.next_square_number(num_libs)), plot_index)
            sample_name = lib.lib_settings.sample_name
            dist = lib.get_fragment_count_distribution()
            plot.hist(dist, bins = log_bins, color=ribo_utils.skyBlue, histtype='stepfilled', edgecolor = None, lw = 0)
            plot.set_xlabel("# reads", fontsize = 10)
            plot.set_ylabel("# sequences (%d/%d have >= %d reads)" % (ribo_utils.number_passing_cutoff(dist, cutoff), len(dist), cutoff), fontsize = 10)
            plot.set_xlim(0, 10**8.1)
            plot.set_xscale('symlog', linthreshx=1, basex=2)
            plot.set_xticks(log_bins[::2])
            plot.axvline(cutoff, ls = 'dashed')
            plot.set_title(sample_name, fontsize = 14)
            plot_index += 1
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.15, wspace=0.4, hspace=0.6)
        out_name =  os.path.join(
          self.experiment_settings.get_rdir(),
          'QC',
          'tl_count_distributions.pdf')
        plt.savefig(out_name, transparent='True', format='pdf')
        plt.clf()

    def read_cutoff_choice_plot(self):
        '''
        This should be similar to the supplementary figure from Ingolia 2009 to determine read cutoff.

        :return:
        '''
        output_file = os.path.join(
            self.experiment.settings.get_rdir(),
            'QC',
            'cutoff_choice.pdf')
        libs = self.experiment.monosome_libs
        num_libs = len(libs)
        num_plots_wide = num_libs-1
        num_plots_high = num_libs-1
        fig = plt.figure(figsize=(8,8))

        bins = []
        currentBinCeiling = 16
        while currentBinCeiling <= (2**12)+1:  # add 1 integer of wiggle room because of floating point error in sqrt(2) calculation
            bins.append(currentBinCeiling)
            currentBinCeiling = currentBinCeiling * math.sqrt(2)


        for i in range(len(libs)):
            for j in range(i+1, len(libs)):
                plot_index = (j-1)*(num_plots_wide)+(i+1)
                plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index)
                binnedGenes = {}
                for binCeiling in bins:
                    binnedGenes[binCeiling] = []
                for sequence_name in libs[i].sorted_names():
                    total = libs[i].get_counts(sequence_name) + libs[j].get_counts(sequence_name)
                    if total <= 0:
                        continue
                    else:
                        for bini in range(len(bins)):
                            if total <= bins[bini]:
                                binnedGenes[bins[bini]].append(sequence_name)
                                break
                                # if it's bigger than all the bins, add it to the last one
                        if total > bins[-1]:
                            binnedGenes[bins[-1]].append(sequence_name)
                binStDevs = {}
                for bin in binnedGenes:
                    values = []
                    for sequence_name in binnedGenes[bin]:
                        rep1Val = libs[i].get_counts(sequence_name)
                        total = rep1Val + libs[j].get_counts(sequence_name)
                        ratio = rep1Val / float(total)
                        values.append(ratio)
                    binStDevs[bin] = np.std(values)


                x = sorted(binStDevs.keys())
                y = np.array([binStDevs[bin] for bin in x])


                plot.set_xlabel("%s+%s counts" % (libs[i].lib_settings.sample_name, libs[j].lib_settings.sample_name))
                plot.set_ylabel("std dev %s/(%s+%s) counts" % (libs[i].lib_settings.sample_name,
                                                               libs[i].lib_settings.sample_name,
                                                               libs[j].lib_settings.sample_name))
                t = np.arange(2.0, 2.**15, 1.)
                p = float(libs[i].total_mapped_fragments) / float(libs[j].total_mapped_fragments + libs[i].total_mapped_fragments)
                expected = np.sqrt(p*(1.0-p)/t)  # will plot expected value from binomial counting


                plot.set_xscale('log', basex=2)
                #plot.set_yscale('symlog', linthreshy=0.01)
                plot.scatter(x, y, color=ribo_utils.blue, s=4)
                plot.plot(t, expected, color=ribo_utils.vermillion, lw = 1, linestyle='dashed')
                #plot.set_xlim(0, 1000)
                #plot.set_ylim(0, 1000)
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.4, hspace=0.4)
        plt.savefig(output_file, transparent='True', format='pdf')