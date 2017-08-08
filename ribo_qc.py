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
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

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
        ribo_utils.make_dir(self.experiment.rdir_path('QC'))
        self.genome = ribo_utils.genome_sequence(self.experiment_settings.get_genome_sequence_files())
        self.GTF_annotations = ribo_utils.gtf_data(self.experiment_settings.get_annotation_GTF_file())
        self.lib_QCs = [single_lib_qc(self, lib_settings) for lib_settings in self.experiment_settings.iter_lib_settings()]
        self.plot_rRNA_size_distributions()
        #self.write_read_annotations_table()

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

    def plot_rRNA_size_distributions(self):
        dfs = []
        for qc_lib in self.lib_QCs:
            frag_dict = qc_lib.rrna_fragment_lengths
            frag_lengths = sorted(frag_dict.keys())
            frag_length_counts = [frag_dict[length] for length in frag_lengths]
            d = {'fragment length': frag_lengths, '# reads': frag_length_counts,
                 '% reads': 100. * np.array(frag_length_counts) / sum(frag_length_counts),
                 'sample': [qc_lib.lib_settings.sample_name] * len(frag_length_counts)}
            temp_df = pd.DataFrame(data=d)
            dfs.append(temp_df)
        frag_length_df = pd.concat(dfs)
        outname = os.path.join(self.experiment_settings.get_rdir(), 'QC', 'rRNA_fragment_sizes.tsv')
        frag_length_df.to_csv(outname, sep='\t')
        sns.set(style="white", color_codes=True, font_scale=1)
        sns.set_style("ticks", {"xtick.major.size": 3, "xtick.minor.size": 2, "ytick.major.size": 2})
        g = sns.factorplot(x='fragment length', y='% reads', hue="sample", data=frag_length_df, legend_out=True,
                           kind="point", size=8)
        outname = os.path.join(self.experiment_settings.get_rdir(), 'QC', 'rRNA_fragment_sizes_single_plot.pdf')
        g.savefig(outname, transparent=True)
        g = sns.factorplot(x='fragment length', y='% reads', col="sample", col_wrap=3, data=frag_length_df, legend_out=True,
                           kind="point", size=8)
        outname = os.path.join(self.experiment_settings.get_rdir(), 'QC', 'rRNA_fragment_sizes_multi_plot.pdf')
        g.savefig(outname, transparent=True)

    def plot_read_annotations_summary(self):
        dfs = []
        for qc_lib in self.lib_QCs:
            dict_list = [] # a list pf tuples that I wil later cast to a dict
            dict_list.append(('sample',qc_lib.lib_settings.sample_name))
            dict_list.append(('short or low quality', qc_lib.adaptor_stats['reads processed'] - qc_lib.adaptor_stats['trimmed reads available'] ))
            for ncrna_reference in qc_lib.ncrna_reference_counts:
                dict_list.append((ncrna_reference, qc_lib.ncrna_reference_counts[ncrna_reference]))




            frag_dict = qc_lib.rrna_fragment_lengths
            frag_lengths = sorted(frag_dict.keys())
            frag_length_counts = [frag_dict[length] for length in frag_lengths]
            d = {'fragment length': frag_lengths, '# reads': frag_length_counts,
                 '% reads': 100. * np.array(frag_length_counts) / sum(frag_length_counts),
                 'sample': [qc_lib.lib_settings.sample_name] * len(frag_length_counts)}
            temp_df = pd.DataFrame(data=d)
            dfs.append(temp_df)
        frag_length_df = pd.concat(dfs)


class single_lib_qc():
    def __init__(self, parent_qc, lib_settings):
        """
        Constructor for Library class
        """
        self.parent_qc = parent_qc
        self.lib_settings = lib_settings
        self.get_adaptor_trimming_stats()
        self.ncrna_samfile = pysam.AlignmentFile(self.lib_settings.get_ncrna_mapped_reads(), "rb")
        self.count_ncrna_mapping_reads()
        self.ncrna_samfile.close()

        self.genome_samfile = pysam.AlignmentFile(self.lib_settings.get_genome_mapped_reads(), "rb")
        #self.get_mapping_multiplicity_stats()
        self.get_mapping_annotation_summary()
        self.genome_samfile.close()

    def count_ncrna_mapping_reads(self):
        self.ncrna_reference_counts = defaultdict(int)
        self.ncrna_sequence_multiplicities = {}
        self.rrna_fragment_lengths = defaultdict(int)
        self.total_rRNA_alignments = 0
        self.lib_settings.write_to_log('counting ncrna-mapped reads')
        for reference in self.ncrna_samfile.references:
                for alignment in self.ncrna_samfile.fetch(reference=reference):
                    if not alignment.is_secondary:
                        self.ncrna_reference_counts[reference] += 1
                        if 'rRNA' in reference:
                            self.total_rRNA_alignments += 1
                            self.rrna_fragment_lengths[alignment.query_length] += 1
                        if alignment.query_sequence in self.ncrna_sequence_multiplicities:
                            self.ncrna_sequence_multiplicities[alignment.query_sequence]['count'] += 1
                        else:
                            self.ncrna_sequence_multiplicities[alignment.query_sequence] = {}
                            self.ncrna_sequence_multiplicities[alignment.query_sequence]['reference'] = reference
                            self.ncrna_sequence_multiplicities[alignment.query_sequence]['count'] = 1
        out_file = open(self.lib_settings.get_ncrna_most_common_reads(), 'w')
        for sequence in sorted(self.ncrna_sequence_multiplicities.keys(), key=lambda x:self.ncrna_sequence_multiplicities[x]['count'], reverse=True):
            out_file.write('%s\t%s\t%d\n' % (self.ncrna_sequence_multiplicities[sequence]['reference'], sequence, self.ncrna_sequence_multiplicities[sequence]['count']))
        out_file.close()
        self.lib_settings.write_to_log('done counting rrna-mapped reads')


    def get_adaptor_trimming_stats(self):
        """
        Parses the outpur from skewer to consolidate stats
        :param lib_settings: 
        :return: 
        """
        self.lib_settings.write_to_log('parsing adaptor trimming stats in %s' % (self.lib_settings.get_adaptor_trimming_log()))
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
        for alignment in self.genome_samfile.fetch():
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

    def get_mapping_annotation_summary(self, reversed_reads=False):
        '''
        produce a summary of where all of the reads in a dataset map.
        :return: 
        '''
        self.lib_settings.write_to_log('making mapping annotation_summary for all mapping reads')
        self.annotation_mapping_counts = defaultdict(lambda : defaultdict(int))
        self.total_primary_alignments = 0
        print self.genome_samfile.references
        print sorted(self.parent_qc.GTF_annotations.chr_to_entry.keys())
        for chromosome in self.genome_samfile.references:
            self.lib_settings.write_to_log(chromosome)
            for alignment in self.genome_samfile.fetch(reference=chromosome):
                if not alignment.is_secondary:
                    self.total_primary_alignments += 1
                    multiplicity = int(alignment.get_tag('NH:i'))
                    assert multiplicity > 0
                    if multiplicity > 1:
                        uniqueness = 'multiple mapping'
                    else:
                        uniqueness = 'unique mapping'
                    if not alignment.is_reverse:  # alignment on + strand
                        strand = '+'
                    else:  # alignment on - strand
                        strand = '-'
                    annotation_entry = self.parent_qc.GTF_annotations.find_smallest_annotation_at_position(chromosome, strand, alignment.reference_start, alignment.reference_end)
                    if annotation_entry == None:
                        self.annotation_mapping_counts[uniqueness]['not annotated']+=1
                    elif annotation_entry.get_value('type') in ['transcript', 'exon']:
                        self.annotation_mapping_counts[uniqueness][annotation_entry.get_value('transcript_type')]+=1
                    elif annotation_entry.get_value('type') == 'UTR':
                        # need to differentiate 5' from 3' UTR
                        type = self.parent_qc.GTF_annotations.utr_type(annotation_entry)
                        assert type is not None
                        self.annotation_mapping_counts[uniqueness][type] += 1
                    else:
                        self.annotation_mapping_counts[uniqueness][annotation_entry.get_value('type')]+=1
        self.lib_settings.write_to_log('%s' % str(self.annotation_mapping_counts))
        self.lib_settings.write_to_log('total_aligned: %d' % self.total_primary_alignments)
        #self.lib_settings.write_to_log('total_annotated: %d' % sum(self.annotation_mapping_counts.values()))
