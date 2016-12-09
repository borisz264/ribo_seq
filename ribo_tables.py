import ribo_utils
import numpy as np
import os
import scipy.stats as stats
import math

def make_readthrough_table(experiment):
    all_genes = set()
    sample_names = []
    for lib in experiment.libs:
        sample_name = lib.lib_settings.sample_name
        sample_names.append(sample_name)
        tx_w_data = set([tx.sequence_name for tx in lib.transcripts.values() if not tx.compute_readthrough_ratio(16, read_end='3p',
                                                                                           read_lengths='all',
                                                                                           cds_cutoff=128) == None ])
        all_genes = all_genes.union(tx_w_data)
    out_name = os.path.join(experiment.settings.get_rdir(), 'tables', 'readthrough_fractions.tsv')
    f = open(out_name, 'w')
    f.write('tx_id\t%s\n' % '\t'.join(sample_names))
    for tx_name in all_genes:
        values = [str(lib.transcripts[tx_name].compute_readthrough_ratio(16, read_end='3p',read_lengths='all', cds_cutoff=128, log = False)) if
                  tx_name in lib.transcripts and lib.transcripts[tx_name].compute_readthrough_ratio(16, read_end='3p',
                                                                                                    read_lengths='all', cds_cutoff=128, log = False) != None
                  else '' for lib in experiment.libs]
        f.write('%s\t%s\n' % (tx_name, '\t'.join(values)))
    f.close()

def transcriptome_features_table(experiment):
    all_transcripts = set()
    first_lib = experiment.libs[0]
    all_tx = set(first_lib.transcripts.values())
    out_name = os.path.join(experiment.settings.get_rdir(), 'tables', 'transcript_features.tsv')
    f = open(out_name, 'w')
    f.write('tx_id\tstop_codon_context\tsecond_stop_codon\tUTR_length\tTL_length\tCDS_length\t'
            'tx_length\textension_nt_length\tUTR_A_percent\tUTR_T_percent\tUTR_C_percent\tUTR_G_percent\n')
    for tx in all_tx:
        values = []
        values.append(tx.stop_codon_context())
        f.write('%s\t%s\n' % (tx.sequence_name, '\t'.join(values)))
    f.close()
