import os
import ConfigParser
import simplejson
import itertools
import shutil
import datetime

import ribo_utils

class ribo_settings:
    def __init__(self, settings_file):
        self.settings_file = settings_file
        self.process_settings(settings_file)
    def get_force_recount(self, count_type):
        return self.settings['force_%s_recount' % count_type]
    def get_settings_file(self):
        return self.settings_file
    def get_property(self, property, default=None):
        try:
            if not property in self.settings and default != None:
                return default
            return self.settings[property]
        except:
            print self.settings
            raise  ValueError('cannot find %s' % property)
    def get_rdir(self):
        ribo_utils.make_dir(self.rdir)
        return self.rdir
    def get_wdir(self):
        ribo_utils.make_dir(self.wdir)
        return self.wdir
    def get_input_barcode(self):
        return self.settings['library_seq_barcode']

    def iter_lib_settings(self):
        for i in range(len(self.sample_names)):
            yield ribo_lib_settings(self,
                                    self.sample_names[i],
                                    self.fastq_gz_file_handles[i])

    def process_settings(self, settings_file):
        """
        - reads the settings file and converts str to float, list, etc.
        - stores result in self.settings as a dict()
        """
        int_keys = ['comparison_read_cutoff', 'min_insert_length', 'max_insert_length',
                     'quality_cutoff']
        #float_keys = []
        str_keys = ['adaptor_3p_sequence', 'star_genome_dir', 'canonical_tx_features', 'canonical_tx_seqs']
        boolean_keys = ['force_remapping', 'force_recount', 'force_index_rebuild', 'force_retrim', 'trim_adaptor', 'make_interactive_plots']
        list_str_keys = ['fastq_gz_files', 'sample_names']
        #list_float_keys = ['concentrations', 'input_rna']
        extant_files = ['canonical_tx_features', 'canonical_tx_seqs']
        config = ConfigParser.ConfigParser()
        config.read(settings_file)
        settings = {}
        for section in config.sections():
            for option in config.options(section):
                settings[option] = config.get(section, option)
                settings[section] = True
        for k in int_keys:
            settings[k] = int(settings[k])
        for k in str_keys:
            settings[k] = settings[k]
        #for k in float_keys:
        #    settings[k] = float(settings[k])
        for k in boolean_keys:
            if not settings[k].lower() in ['true', 'false']:
                raise ValueError(
                  'Boolean value %s must be "true" or "false"' % k)
            settings[k] = settings[k].lower() == 'true'
        #for k in list_float_keys:
        #    settings[k] = map(float, simplejson.loads(settings[k]))
        #for k in list_int_keys:
        #    settings[k] = map(int, simplejson.loads(settings[k]))
        for k in list_str_keys:
            settings[k] = simplejson.loads(settings[k])
        self.fqdir = settings['fastq_dir']
        self.sample_names = settings['sample_names']
        #for paired end reads, there are now 2 fastq files per sample, at least that's how MIT provides the data.
        self.fastq_gz_files = settings['fastq_gz_files']
        self.fastq_gz_file_handles = [os.path.join(self.fqdir, fastq_gz_file) for fastq_gz_file in self.fastq_gz_files]
        for file_handle in self.fastq_gz_file_handles:
            assert ribo_utils.file_exists(file_handle)
        for k in extant_files:
            assert ribo_utils.file_exists(settings[k])
        self.settings = settings
        self.rdir = settings['results_dir']
        ribo_utils.make_dir(self.rdir)
        shutil.copy(settings_file, self.rdir)

    def check_barcode_lens(self):
        """
        verifies that all the barcodes are the same length
        """
        barcode_lens = set(map(len, self.settings['barcodes']))
        if 1 != len(barcode_lens):
            raise ValueError('all barcodes must be the same length')
        self.barcode_len = barcode_lens.pop()
        self.settings['barcode_len'] = self.barcode_len

    def check_barcodes_are_separated(self):
        """
        makes sure the barcodes are all totally distinguishable
        """
        for b1, b2 in itertools.combinations(self.settings['barcodes'], 2):
            hamming_dist = ribo_utils.hamming_distance(b1, b2)
            if hamming_dist < 2:
                raise ValueError('The barcodes supplied are not well '
                  'separated: %s-%s' % (b1, b2))

    def get_bowtie_index(self):
        index = os.path.join(
          self.get_rdir(),
          'bowtie_indices',
          'pool_index')
        return index

    def get_rRNA_bowtie_index(self):
        index = index = self.get_property('rrna_index')
        return index

    def get_star_genome_dir(self):
        index = self.get_property('star_genome_dir')
        return index

    def get_log(self):
        log = os.path.join(
          self.get_rdir(),
          'log.txt')
        return log

    def write_to_log(self, text, add_time = True):
        f = open(self.get_log(), 'a')
        now = datetime.datetime.now()
        time = now.strftime("%Y-%m-%d %H:%M")
        if add_time:
            f.write('[%s] %s\n' % (time, text))
        else:
            f.write(text)
        f.close()

    def get_overall_mapping_summary(self):
        summary_file = os.path.join(
          self.get_rdir(),
          'QC',
          'mapping_summary.txt')
        return summary_file

class ribo_lib_settings:
    def __init__(self, experiment_settings, sample_name, fastq_gz_filehandle):
        self.experiment_settings = experiment_settings
        self.sample_name = sample_name
        self.fastq_gz_filehandle = fastq_gz_filehandle

    def get_property(self, property):
        return self.experiment_settings.get_property(property)

    def get_log(self):
        ribo_utils.make_dir(os.path.join(self.experiment_settings.get_rdir(), 'logs'))
        log = os.path.join(
          self.experiment_settings.get_rdir(),
          'logs',
          '%(sample_name)s.log' %
           {'sample_name': self.sample_name})
        return log
    def write_to_log(self, text, add_time = True):
        f = open(self.get_log(), 'a')
        now = datetime.datetime.now()
        time = now.strftime("%Y-%m-%d %H:%M")
        if add_time:
            f.write('[%s] %s\n' % (time, text))
        else:
            f.write(text)
        f.close()

    def get_fastq_gz_file(self):
        return self.fastq_gz_filehandle

    def get_collapsed_reads(self):
        collapsed_reads = os.path.join(
          self.experiment_settings.get_rdir(),
          'collapsed_reads',
          '%(sample_name)s.fasta.gz' %
           {'sample_name': self.sample_name})
        return collapsed_reads

    def get_pool_mapping_stats(self):
        pool_mapping_stats = os.path.join(self.experiment_settings.get_rdir(), 'mapping_stats', '%(sample_name)s.pool.txt' % {'sample_name': self.sample_name})
        return pool_mapping_stats

    def get_mapped_reads_prefix(self):
        mapped_reads = os.path.join(self.experiment_settings.get_rdir(), 'mapped_reads', '%(sample_name)s' % {'sample_name': self.sample_name})
        return mapped_reads

    def get_genome_mapped_reads(self):
        mapped_reads = os.path.join(self.experiment_settings.get_rdir(), 'mapped_reads', '%(sample_name)sAligned.sortedByCoord.out.bam' % {'sample_name': self.sample_name})
        return mapped_reads

    def get_transcript_mapped_reads(self):
        mapped_reads = os.path.join(self.experiment_settings.get_rdir(), 'mapped_reads', '%(sample_name)sAligned.toTranscriptome.out.bam' % {'sample_name': self.sample_name})
        return mapped_reads

    def get_trimmed_reads(self, prefix_only = False):
        if prefix_only:
            trimmed_reads = os.path.join(
                self.experiment_settings.get_rdir(),
                'trimmed_reads',
                '%(sample_name)s' %
                {'sample_name': self.sample_name})
        else:
            trimmed_reads = os.path.join(
              self.experiment_settings.get_rdir(),
              'trimmed_reads',
              '%(sample_name)s-trimmed.fastq.gz' %
               {'sample_name': self.sample_name})
        return trimmed_reads

    def get_trimming_log(self):
        """
        :return: name of trimming log from Skewer
        """
        trimming_log = os.path.join(
          self.experiment_settings.get_rdir(),
          'trimmed_reads',
          '%(sample_name)s-trimmed.log' %
           {'sample_name': self.sample_name})
        return trimming_log

    def get_transcript_counts(self):
        sequence_counts = os.path.join(
          self.experiment_settings.get_rdir(),
          'transcript_counts',
          '%(sample_name)s.counts.pkl' %
           {'sample_name': self.sample_name})
        return sequence_counts

    def get_overall_contamination_summary(self):
        summary_file = os.path.join(
          self.experiment_settings.get_rdir(),
          'QC',
          '%(sample_name)s.contamination_summary.txt' %
           {'sample_name': self.sample_name})
        return summary_file

    def collapsed_reads_exist(self):
        collapsed_reads = self.get_collapsed_reads()
        return ribo_utils.file_exists(collapsed_reads)

    def trimmed_reads_exist(self):
        trimmed_reads = self.get_trimmed_reads()
        return ribo_utils.file_exists(trimmed_reads)

    def mapped_reads_exist(self):
        mapped_reads = self.get_genome_mapped_reads()
        return ribo_utils.file_exists(mapped_reads)

    def sequence_counts_exist(self):
        sequence_counts = self.get_transcript_counts()
        return ribo_utils.file_exists(sequence_counts)