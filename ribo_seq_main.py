__author__ = 'boris zinshteyn'
"""
Intended for processing of ribosome footprint profiling data from mammalian cells
"""
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42 #leaves most text as actual text in PDFs, not outlines
import os
import argparse
import subprocess

import ribo_settings
import ribo_utils
import ribo_lib
import ribo_qc
import ribo_plotting
import ribo_tables
from collections import defaultdict
import numpy as np
import scipy.stats as stats

class experiment:
    def __init__(self, settings, threads):
        self.threads = threads
        self.settings = settings
        self.num_libs = len([x for x in settings.iter_lib_settings()])
        self.remove_adaptor()
        self.map_reads()
        self.initialize_libs()

    def initialize_libs(self):
        self.settings.write_to_log('initializing libraries, counting reads')
        ribo_utils.make_dir(self.rdir_path('transcript_counts'))
        self.libs = []

        ribo_utils.parmap(lambda lib_settings: ribo_lib.assign_tx_reads(self.settings, lib_settings),
                          self.settings.iter_lib_settings(), nprocs=self.threads)
        map(lambda lib_settings: self.initialize_lib(lib_settings), self.settings.iter_lib_settings())
        self.settings.write_to_log('initializing libraries, counting reads, done')


    def find_lib_by_sample_name(self, sample_name):
        for lib in self.libs:
            if lib.lib_settings.sample_name == sample_name:
                return lib
        assert False #if this triggers, your settings file is broken.


    def initialize_lib(self, lib_settings):
        lib = ribo_lib.ribo_lib(self.settings, lib_settings)
        self.libs.append(lib)

    def make_tables(self):
        ribo_utils.make_dir(self.rdir_path('tables'))
        ribo_tables.make_readthrough_table(self)

    def make_plots(self):
        ribo_utils.make_dir(self.rdir_path('plots'))

        ribo_plotting.plot_start_codon_average(self, up=100, down=100, min_cds_reads=self.settings.get_property('comparison_read_cutoff'),
                                               read_end='3p', read_lengths='all')
        ribo_plotting.plot_stop_codon_average(self, up=100, down=100, min_cds_reads=128,
                                               read_end='3p', read_lengths='all')
        ribo_plotting.plot_first_exon_average(self, up=100, down=100, min_cds_reads=128,
                                               read_end='3p', read_lengths='all')
        ribo_plotting.plot_start_codon_average(self, up=100, down=100,
                                               min_cds_reads=self.settings.get_property('comparison_read_cutoff'),
                                               read_end='3p', read_lengths=[29, 30])
        ribo_plotting.plot_stop_codon_average(self, up=100, down=100, min_cds_reads=128,
                                          read_end='3p', read_lengths=[29, 30])

        ribo_plotting.plot_start_codon_average(self, up=100, down=100,
                                               min_cds_reads=self.settings.get_property('comparison_read_cutoff'),
                                               read_end='5p', read_lengths='all')

        ribo_plotting.plot_first_exon_average(self, up=100, down=100,
                                               min_cds_reads=self.settings.get_property('comparison_read_cutoff'),
                                               read_end='5p', read_lengths='all')

        ribo_plotting.plot_stop_codon_average(self, up=100, down=100, min_cds_reads=128,
                                              read_end='5p', read_lengths='all')
        ribo_plotting.plot_start_codon_average(self, up=100, down=100,
                                               min_cds_reads=self.settings.get_property('comparison_read_cutoff'),
                                               read_end='5p', read_lengths=[29, 30])
        ribo_plotting.plot_stop_codon_average(self, up=100, down=100, min_cds_reads=128,
                                              read_end='5p', read_lengths=[29, 30])
        ribo_plotting.plot_fragment_length_distributions(self)
        ribo_plotting.plot_frame_distributions(self)

        ribo_plotting.plot_stop_positional_read_lengths(self, up=100, down=100, min_cds_reads=128, read_end='3p')
        ribo_plotting.plot_stop_positional_read_lengths(self, up=100, down=100, min_cds_reads=128, read_end='5p')
        ribo_plotting.plot_first_exon_positional_read_lengths(self, up=100, down=100, min_cds_reads=128, read_end='5p')
        ribo_plotting.plot_start_positional_read_lengths(self, up=100, down=100, min_cds_reads=128, read_end='3p')
        ribo_plotting.plot_start_positional_read_lengths(self, up=100, down=100, min_cds_reads=128, read_end='5p')
        ribo_plotting.plot_first_exon_positional_read_lengths(self, up=100, down=100, min_cds_reads=128, read_end='5p')
        ribo_plotting.plot_readthrough_box(self)



    def remove_adaptor(self):
        self.settings.write_to_log('trimming adaptors')
        if not self.settings.get_property('force_retrim'):
            for lib_settings in self.settings.iter_lib_settings():
                if not lib_settings.trimmed_reads_exist():
                    break
            else:
                return

        if self.settings.get_property('trim_adaptor'):
            ribo_utils.make_dir(self.rdir_path('trimmed_reads'))
            map(lambda lib_setting: self.remove_adaptor_one_lib(lib_setting, self.threads),
                              self.settings.iter_lib_settings())
        self.settings.write_to_log('done trimming adaptors')

    def remove_adaptor_one_lib(self, lib_settings, threads):
        lib_settings.write_to_log('adaptor trimming')
        """
        -x specifies the 3' adaptor to trim from the forward read
        -Q specifies the lowest acceptable mean read quality before trimming
        -l specifies the minimum post-trimming read length
        -L specifies the maximum post-trimming read length
        -o is the output prefix
        --threads specifies number of threads to use
        """
        command_to_run = 'skewer -x %s -Q %d  -l %d -L %d -o %s --quiet --threads %s %s 1>>%s 2>>%s' % (
            self.settings.get_property('adaptor_3p_sequence'),
            self.settings.get_property('quality_cutoff'), self.settings.get_property('min_insert_length'),
            self.settings.get_property('max_insert_length'),
            lib_settings.get_trimmed_reads(prefix_only=True),
            threads,
            lib_settings.get_fastq_gz_file(),
            lib_settings.get_log(), lib_settings.get_log())
        subprocess.Popen(command_to_run, shell=True).wait()
        subprocess.Popen('gzip %s-trimmed.fastq' % (lib_settings.get_trimmed_reads(prefix_only=True)), shell=True).wait()
        lib_settings.write_to_log('adaptor trimming done')


    def map_reads(self):
        """
        map all reads using bowtie
        :return:
        """
        self.settings.write_to_log('mapping reads')
        if not self.settings.get_property('force_remapping'):
            for lib_settings in self.settings.iter_lib_settings():
                if not lib_settings.mapped_reads_exist():
                    break
            else:
                return
        ribo_utils.make_dir(self.rdir_path('mapped_reads'))
        map(lambda lib_setting: self.map_one_library(lib_setting, self.threads), self.settings.iter_lib_settings())
        self.settings.write_to_log( 'finished mapping reads')

    def map_one_library(self, lib_settings, threads):
        lib_settings.write_to_log('mapping_reads')
        command_to_run = 'STAR --runThreadN %d --genomeDir %s --readFilesIn %s --readFilesCommand gunzip -c ' \
                         '--outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 1 --outWigType wiggle read1_5p --outFileNamePrefix %s' \
                         ' --quantMode TranscriptomeSAM --outReadsUnmapped FastX 1>>%s 2>>%s' %\
                         (threads, self.settings.get_star_genome_dir(), lib_settings.get_trimmed_reads(),
                          lib_settings.get_mapped_reads_prefix(), lib_settings.get_log(), lib_settings.get_log())

        subprocess.Popen(command_to_run, shell=True).wait()
        #sort transcript-mapped bam file

        subprocess.Popen('samtools sort %s %s.temp_sorted 1>>%s 2>>%s' % (lib_settings.get_transcript_mapped_reads(), lib_settings.get_transcript_mapped_reads(),
                                                                          lib_settings.get_log(), lib_settings.get_log()), shell=True).wait()
        subprocess.Popen('mv %s.temp_sorted.bam %s' % (lib_settings.get_transcript_mapped_reads(),
                                                                          lib_settings.get_transcript_mapped_reads()), shell = True).wait()

        subprocess.Popen('samtools index %s' % (lib_settings.get_transcript_mapped_reads()), shell = True).wait()
        subprocess.Popen('samtools index %s' % (lib_settings.get_genome_mapped_reads()), shell=True).wait()
        lib_settings.write_to_log('mapping_reads done')

    def rdir_path(self, *args):
        return os.path.join(self.settings.get_rdir(), *args)

    def get_rdir_fhandle(self, *args):
        """
        returns a filehandle to the fname in the rdir
        """
        out_path = self.rdir_path(*args)
        out_dir = os.path.dirname(out_path)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        return ribo_utils.aopen(out_path, 'w')

    def perform_qc(self):
        qc_engine = ribo_qc.ribo_qc(self, self.settings, self.threads)
        qc_engine.write_mapping_summary(self.settings.get_overall_mapping_summary())
        qc_engine.print_library_count_concordances()
        qc_engine.plot_average_read_positions()
        qc_engine.plot_fragment_length_distributions()
        qc_engine.plot_count_distributions()
        qc_engine.read_cutoff_choice_plot()

    def make_counts_table(self, fractional=False):
        """
        write out number of fragments mapping to each TL in each dataset
        :param fractional: if True, replace raw counts with library fraction in reads per million
        :return:
        """
        if fractional:
            summary_file = open(os.path.join(
                self.rdir_path('tables'),
                'rpm.txt'), 'w')
        else:
            summary_file = open(os.path.join(
                self.rdir_path('tables'),
                'raw_counts.txt'), 'w')

        header = 'sequence name\t' + '\t'.join([lib.lib_settings.sample_name for lib in self.libs]) + '\n'
        summary_file.write(header)
        if fractional:
            for sequence_name in self.libs[0].pool_sequence_mappings:
                out_line = '%s\t%s\n' % (sequence_name,
                                         '\t'.join(['%f' % ((10**6)*lib.pool_sequence_mappings[sequence_name].fragment_count/float(lib.total_mapped_fragments)) for lib in self.libs]))
                summary_file.write(out_line)
        else:
            for sequence_name in self.libs[0].pool_sequence_mappings:
                out_line = '%s\t%s\n' % (sequence_name,
                                         '\t'.join(['%f' %
                                                    lib.pool_sequence_mappings[sequence_name].fragment_count
                                                    for lib in self.libs]))
                summary_file.write(out_line)
        summary_file.close()



    def write_sequence_subset(self, recruitment_cutoff, read_cutoff=128, as_RNA=True):
        """
        write out fasta of all sequences that pass a certain recruitment cutoff in all libraries
        :return:
        """
        output_file = open(os.path.join(
            self.rdir_path('tables'),
            'recruitment_above_%f.fasta' % recruitment_cutoff), 'w')
        trimmed_sequences = ribo_utils.convertFastaToDict(self.settings.get_trimmed_pool_fasta())

        for sequence_name in self.monosome_libs[0].pool_sequence_mappings:
            rec_scores = [(self.monosome_libs[i].get_rpm(sequence_name) / (self.monosome_libs[i].get_rpm(sequence_name) +
                                                                           self.mrnp_libs[i].get_rpm(sequence_name)))
                          for i in range(len(self.monosome_libs)) if (self.monosome_libs[i].get_counts(sequence_name) +
                              self.mrnp_libs[i].get_counts(sequence_name)) >= read_cutoff]
            if (len(rec_scores) == len(self.monosome_libs)):
                average_score = np.average(rec_scores)
                if average_score >=recruitment_cutoff:
                    output_file.write('>%s_rec_%f\n' % (sequence_name, average_score))
                    seq = trimmed_sequences[sequence_name]
                    if as_RNA:
                        seq = ribo_utils.rna(seq)
                    output_file.write('%s\n' % (seq))
        output_file.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("settings_file")
    parser.add_argument("--make-tables",
                        help="Makes tables.",
                        action='store_true')
    parser.add_argument("--perform-qc",
                        help="performs quality control analysis.",
                        action='store_true')
    parser.add_argument("--make-plots",
                        help="Makes plots.",
                        action='store_true')
    parser.add_argument("--comparisons",
                        help="Does comparisons to other experiments",
                        action='store_true')
    parser.add_argument("--all-tasks",
                        help="Makes plots, tables, folding and comparisons",
                        action='store_true')
    parser.add_argument("--threads",
                        help="Max number of processes to use",
                        type = int, default = 8)
    args = parser.parse_args()

    return args

def main():
    """
    """
    args = parse_args()
    settings = ribo_settings.ribo_settings(args.settings_file)
    ribo_experiment = experiment(settings, args.threads)
    print 'experiment ready'
    if args.perform_qc or args.all_tasks:
        print 'QC'
        settings.write_to_log('performing QC')
        ribo_experiment.perform_qc()
        settings.write_to_log('done performing QC')
    if args.make_tables or args.all_tasks:
        print 'tables'
        settings.write_to_log('making tables')
        ribo_experiment.make_tables()
        settings.write_to_log('done making tables')

    if args.make_plots or args.all_tasks:
        print 'plots'
        settings.write_to_log('making plots')
        ribo_experiment.make_plots()
        settings.write_to_log('done making plots')
    '''
    if args.comparisons or args.all_tasks:
        settings.write_to_log('doing comparisons')
        ribo_experiment.compare_all_other_experiments()
    '''

main()