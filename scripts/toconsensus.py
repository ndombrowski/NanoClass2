#!/usr/bin/env python3

__author__ = "Evelien Jongepier"
__copyright__ = "Copyright 2020, Evelien Jongepier"
__email__ = "e.jongepier@uva.nl"
__license__ = "MIT"

import sys
import re
import argparse
import collections


class ToConsensus(object):

    def __init__(self, taxlists, consensustax, pconsensus):
        self.taxlists = taxlists
        self.consensustax = consensustax
        self.pconsensus = pconsensus
        return


    def taxlist_per_sample_per_run_dict(self):
        '''
        Example outout: {run : {sample1 : [path1, path2, ...],
          sample2 : [path3, path4, ...], ...}
        Stores the path to each of the input taxlists for each
        of the samples for each of the runs.
        '''
        sampledict = collections.defaultdict(
          dict = collections.defaultdict(list))
        for idx, taxlist in enumerate(self.taxlists):
            run = taxlist.split('/')[1]
            sample = taxlist.split('/')[3].split('.')[0]
            if run not in sampledict:
                sampledict[run] = {sample:[taxlist]}
            elif sample not in sampledict[run]:
                sampledict[run][sample] = [taxlist]
            else:
                sampledict[run][sample].append(taxlist)
        return sampledict

    def run_all_runs_n_samples(self):
        '''
        '''
        sampledict = self.taxlist_per_sample_per_run_dict()
        for run in sampledict:
            for sample in sampledict[run]:
                ftaxlists = sampledict[run][sample]
                samplereaddict = self.get_read_dict_per_sample(ftaxlists)
                nlists = len(ftaxlists)
                consensusreaddict = self.comp_consensus_tax(samplereaddict, nlists)
                for ftaxlist in ftaxlists:
                    precisiondict = self.comp_precision_per_list(consensusreaddict, ftaxlist)
                    path = self.get_output_path(ftaxlist)
                    self.write_precision_tax(precisiondict, path)
                #self.write_consensus_tax(
        return


    def get_output_path(self, ftaxlist):
        dir = '/'.join(ftaxlist.split('/')[:3])
        method = ftaxlist.split('/')[2]
        sample = ftaxlist.split('/')[3].split('.')[0]
        outpath = dir + '/' + sample + '.' + method + '.precision'
        return outpath


    def get_read_dict_per_list(self, taxlist):
        '''
        Example output: readdict = {read_id: {Domain:tool1},
           {Phylum:tool1},{Class:tool1},...}
        Stores tax classification of a single taxlist (i.a the taxo-
        nomic classification based on a single tool) into a dict
        of dict. Outer dict is readid as key and dict of tax as
        value. Inner / nested dicts have tax level is a key with
        tax annotation at that level for that read in that taxlist
        as str.
        '''
        readdict = collections.defaultdict(dict)
        with open (taxlist) as f:
            for l in f.readlines():
                if l.startswith("#"):
                    headers = l.strip().split('\t')[1:]
                else:
                    readid = l.strip().split('\t')[0]
                    taxlist = l.strip().split('\t')[1:7]
                    readdict[readid] = dict(
                        (headers[idx], tax) for idx, tax in enumerate(taxlist))
        return readdict

    def get_read_dict_per_sample(self, ftaxlists):
        '''
        Example output: readdict = {read_id: {Domain:[tool1, tool2,
           tool3,...]},{Phylum:[tool1, tool2,tool3,...]},...}
        Stores tax classification of all taxlists (i.a the taxo-
        nomic classification based on a all tools) into a dict
        of dict of lists. Outer dict is readid as key and dict of tax
        as value. Inner / nested dicts have tax level is a key with
        tax annotation at that level for that read in each of the
        taxlist (thus tools) as list.
        '''
        samplereaddict = collections.defaultdict(
          dict = collections.defaultdict(list))
        for taxlist in ftaxlists:
            readdict = self.get_read_dict_per_list(taxlist)
            for readid in readdict:
                if readid not in samplereaddict:
                    samplereaddict[readid] = readdict[readid]
                    # change tax into list so all taxlists can be appended
                    for level in readdict[readid]:
                        samplereaddict[readid][level] = [samplereaddict[readid][level]]
                else:
                    for level in readdict[readid]:
                        tax = readdict[readid][level]
                        samplereaddict[readid][level].append(tax)
        return samplereaddict


    def comp_consensus_tax(self, samplereaddict, nlists):
        '''
        Example output: {readid = {Domain:cons},{Phylum:cons},...}
        Compute consensus for each read at each taxonomic level
        based on user defined level (0.5 = majority) and the taxo-
        nomic classificatons of each tool.
        '''
        consensusreaddict = collections.defaultdict(dict)

        for readid in samplereaddict:
            for level in samplereaddict[readid]:
                freqdict = collections.Counter(samplereaddict[readid][level])
                maxkey = max(freqdict, key=freqdict.get)
                if (freqdict[maxkey] / float(nlists)) > self.pconsensus:
                    consensusreaddict[readid][level] = maxkey
                else:
                    consensusreaddict[readid][level] = "NA"
        return consensusreaddict


    def comp_precision_per_list(self, consensusreaddict, ftaxlist):
        '''
        Example output: {readid = {Domain:bool}, {Phylum:bool}, ...}
        For each read and each taxomic level, a boolian whether
        the classification matches the consesnus clasification or
        not.
        '''
        readdict = self.get_read_dict_per_list(ftaxlist)
        precisiondict = collections.defaultdict(dict)
        for readid in readdict:
            for level in readdict[readid]:
                if readdict[readid][level] == 'NA' or consensusreaddict[readid][level] == 'NA':
                    precisiondict[readid][level] =  'NA'
                elif readdict[readid][level] == consensusreaddict[readid][level]:
                    precisiondict[readid][level] = '1'
                else:
                    precisiondict[readid][level] = '0'
        return precisiondict


    def write_consensus_tax(self):
        '''
        Example out: "readid \t Domain \t Phylum ..."
        Write consensus to file, incl. header line.
        '''
        consensusreaddict = self.comp_consensus_tax()
        levels = ['Domain','Phylum','Class','Order','Family','Genus']
        with open(self.consensustax, 'w') as f:
            f.write('#readid' + '\t' + '\t'.join(levels) + '\n')
            for readid in consensusreaddict:
                tax = []
                for level in levels:
                    tax.append(consensusreaddict[readid][level])
                f.write(readid + '\t' + '\t'.join(tax) + '\n')
        return


    def write_precision_tax(self, precisiondict, path):
        '''
        Example out: "readid \t Domain \t Phylum ..."
        Write precision (0/1) to file, incl. header line.
        '''
        levels = ['Domain','Phylum','Class','Order','Family','Genus']
        with open(path, 'w') as f:
            f.write('#readid' + '\t' + '\t'.join(levels) + '\n')
            for readid in precisiondict:
                tax = []
                for level in levels:
                    tax.append(precisiondict[readid][level])
                f.write(readid + '\t' + '\t'.join(tax) + '\n')
        return


def main(taxlists, consensustax, pconsensus):
    m = ToConsensus(
        taxlists = taxlists,
        consensustax = consensustax,
        pconsensus = pconsensus
    )
    m.run_all_runs_n_samples()
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Obtain consensus and '
                        'majority taxonomic classifications based on a '
                        'list of taxlists as generated by various classi'
                        'fication tools. Taxlists can be obtained with '
                        'tomat.py')
    parser.add_argument('-l', '--taxlists',  nargs='+', default=None,
                        help='List of paths to taxlists, which contain '
                        'an taxonomic identifier in the first column followed '
                        'by tab-delimited 6-level taxonomic classifications '
                        'sensu Silva.')
    parser.add_argument('-o', '--consensustax', default=None,
                        help='Paths to where output consensus taxonomy '
                        'should be written.')
    parser.add_argument('-c', '--pconsensus', type=float, default=0.5,
                        help='Consensus treshold ranging between 0.5 and 1. '
                        'Minimum fraction of studies that need to show '
                        'consensus about each taxonomic classification at '
                        'each taxonomic level. Majority obtained with 0.5. '
                        'Taxonomic annotations for which there is no consensus '
                        'are returned as NA.')
    args = parser.parse_args()

    assert (len(args.taxlists) > 2), '''
        \n Please provide at least three taxlists to determine
        consensus taxonomy. Read the help (--help) for further information.
        '''

    main(
        taxlists = args.taxlists,
        consensustax = args.consensustax,
        pconsensus = args.pconsensus
        )
