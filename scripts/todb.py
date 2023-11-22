#!/usr/bin/env python

__author__ = "Evelien Jongepier"
__copyright__ = "Copyright 2020, Evelien Jongepier"
__email__ = "e.jongepier@uva.nl"
__license__ = "MIT"

import sys
import re
import argparse
import collections
from Bio import SeqIO


class ToDb(object):

    def __init__(self, refseq, reftax, method, outseq, outtax):
        self.refseq = refseq
        self.reftax = reftax
        self.method = method
        self.outseq = outseq
        self.outtax = outtax
        return


    def get_ref_dict(self):
        '''
        Example out: {id1:"domain1;phylum1;...", id2:"domain2;phylum2;..."}
        Create a dictionary with sequence id as key and taxonomic string as
        value. The exact formatting of the taxonomic string differs between
        methods.
        '''
        refdict = {}
        with open(self.reftax) as t:
            for line in t.readlines():
                id, tax = line.strip().split('\t')[:2]
                tax = re.sub(r'D_[0-6]__', '', tax)
                tax = re.sub(r' ', '_', tax)
                if self.method == 'centrifuge':
                    refdict[id] = tax
                else:
                    tax = tax.split(';')[:6]
                    if self.method == 'mothur':
                        refdict[id] = ';'.join(tax)
                    if self.method == 'rdp':
                        refdict[id] = ';'.join(tax)
                    if self.method == 'spingo':
                        tax.reverse()
                        refdict[id] = '\t'.join(tax)
        return refdict


    def get_out_files(self):
        '''
        Reformat fasta header to make it compatible with the selected method.
        '''
        with open(self.outseq, "w") as outs, open(self.outtax, "w") as outt:
            refdict = self.get_ref_dict()
            for record in SeqIO.parse(open(self.refseq), 'fasta'):
                if record.id in refdict:
                    record.description = ''
                    outt.write(record.id + '\t' + refdict[record.id] + '\n')
                    if self.method == 'mothur':
                        record.id = record.id + '\tNA\t' + refdict[record.id]
                    if self.method == 'rdp':
                        record.id = refdict[record.id] + ';'
                    if self.method == 'spingo':
                        record.id = record.id + "]\t" + refdict[record.id]
                    if self.method == 'centrifuge':
                        record.id = record.id + ' ' + refdict[record.id]
                    SeqIO.write(record, outs, "fasta")



def main(refseq, reftax, method, outseq, outtax):
    m = ToDb(
        refseq = refseq,
        reftax = reftax,
        method = method,
        outseq = outseq,
        outtax = outtax
    )
    m.get_out_files()
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert qiime2-formatted '
                        'database, consisting of a sequence file and a 7-level '
                        'taxonomy file, into database file(s) compatible with '
                        'other taxonomic classification methods.')
    parser.add_argument('-s', '--refseq',  default=None,
                        help='Path to fasta file of reference nucleotide sequences. '
                        'Headers should contain only the sequence identifier.')
    parser.add_argument('-t', '--reftax', default=None,
                        help='Path to tab-separated table with reference taxonomy. '
                        'The first column should contain the sequence identifier '
                        'as given in "refseq", the second column a semi-colon-'
                        'separated 7-level taxonomy.')
    parser.add_argument('-m', '--method', default=None,
                        help='A string specifying for the which taxonomic classi'
                        'fication tool the output db files should be re-formatted. '
                        'Legit values are "mothur", "centrifuge", "rdp" or "spingo".')
    parser.add_argument('-S', '--outseq',  default=None,
                        help='Path to output fasta file with re-formatted sequences.')
    parser.add_argument('-T', '--outtax', default=None,
                        help='Path to tab separated output table with re-formatted '
                        'taxonomy.')


    args = parser.parse_args()

    main(
        refseq = args.refseq,
        reftax = args.reftax,
        method = args.method,
        outseq = args.outseq,
        outtax = args.outtax
        )


