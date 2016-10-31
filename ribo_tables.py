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
    f = open(out_name)
    f.write('tx_id\t%s\n' % '\t'.join(sample_names))
    for tx_name in all_genes:
        values = [str(lib.transcripts[tx_name].compute_readthrough_ratio(16, read_end='3p',read_lengths='all', cds_cutoff=128)) if
                  tx_name in lib.transcripts and lib.transcripts[tx_name].compute_readthrough_ratio(16, read_end='3p',
                                                                                                    read_lengths='all', cds_cutoff=128) != None
                  else '' for lib in experiment.libs]
        f.write('%s\t%s\n' % (tx_name, '\t'.join(values)))
    f.close()
