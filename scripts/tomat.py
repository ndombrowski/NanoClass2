#!/usr/bin/env python3

__author__ = "Evelien Jongepier"
__copyright__ = "Copyright 2020, Evelien Jongepier"
__email__ = "e.jongepier@uva.nl"
__license__ = "MIT"

import sys
import re
import argparse
from collections import defaultdict


class ToMat(object):

    def __init__(self, taxlist, bblastlike, centrifugelike, krakenlike, dbtaxonomy, dbfasta, dbmap):
        self.taxlist = taxlist
        self.bblastlike = bblastlike
        self.centrifugelike = centrifugelike
        self.krakenlike = krakenlike
        self.dbtaxonomy = dbtaxonomy
        self.dbfasta = dbfasta
        self.dbmap = dbmap

        if taxlist is not None:
            self.pathfile = taxlist
        elif bblastlike is not None:
            self.pathfile = bblastlike
        elif centrifugelike is not None:
            self.pathfile = centrifugelike
        elif krakenlike is not None:
            self.pathfile = krakenlike
        self.dir = '/'.join(self.pathfile.split('/')[:3])
        self.run = self.pathfile.split('/')[1]
        self.method = self.pathfile.split('/')[2]
        self.sample = self.pathfile.split('/')[3].split('.')[0]
        self.outpath = self.dir + '/' + self.sample + '.' + self.method
        ## if input is bblastlike then create otumat and taxmat from taxlist created here
        if any([bblastlike, centrifugelike, krakenlike]):
            self.taxlist = self.outpath + ".taxlist"
        return


    def get_reftax_dict(self):
        '''
        Example output: {refid : [Domain, Phylum, ...]}
        Create a dict of taxdb identifyers and taxdb taxonomy.
        Input file can be taxonomy table or fasta with taxonomy
        in header, depending which one is provided.
        '''
        refdict = {}
        if self.dbtaxonomy is not None:
            with open(self.dbtaxonomy) as reftax:
                for l in reftax.readlines():
                    refid, taxstr = l.strip().split('\t')[:2]
                    taxstr = re.sub(r'D_[0-6]__', '', taxstr)
                    refdict[refid] = taxstr.split(';')[:6]
        elif self.dbfasta is not None:
            with open(self.dbfasta) as fasta:
                for l in fasta.readlines():
                    if l.startswith('>'):
                        line = l.strip().split('>')
                        refid = line[1].split(' ')[0]
                        taxstr = line[1].split(' ')[1]
                        refdict[refid] = taxstr.split(';')[:6]
        return refdict


    def get_krakentax_dict(self):
        '''
        Example output: {ncbiid : [silvaid1,silvaid2,...]}
        For krakenlike input the ncbi refids should first be
        translated into silva refids. This creates a dictionary
        of silva to ncbi id. Each silvaid has only one ncbiid
        but same ncbiid's can have multiple silvaid's.
        '''
        krakendict = defaultdict(list)
        with open (self.dbmap) as map:
            for l in map.readlines():
                line = l.strip().split('\t')
                ncbiid = line[1]
                krakendict[ncbiid].append(line[0])
        return krakendict


    def get_krakenlike_read_dict(self):
        '''
        Example output: {readid : [refid1, refid2, ...]}
        For krakenlike input the ncbi refids in the input table
        should first be replaces by the silva refids.
        There are multiple silva ref ids associated with each
        readid hense a dict of lists.
        '''
        with open(self.krakenlike) as f:
            readdict = defaultdict(list)
            krakendict = self.get_krakentax_dict()
            for l in f.readlines():
                if l.startswith("#") or not l.strip():
                    continue
                else:
                    readid, refrefid = l.strip().split('\t')[1:3]
                    if refrefid in krakendict:
                        readdict[readid].extend(krakendict[refrefid])
        return readdict


    def get_read_dict(self):
        '''
        Example output: {readid : [refid1, refid2, ...]}
        Read a bblastlike file and store reads as keys in a
        dictionary with tadid's as values list. Some classifiers
        return multiple refids per read, hense the list (e.g.
        centrifuge, but not bestblast).
        '''
        readdict = defaultdict(list)
        with open(self.pathfile) as f:
            for l in f.readlines():
                if l.startswith("#") or not l.strip():
                    continue
                else:
                    readid, refid = l.strip().split('\t')[:2]
                    readdict[readid].append(refid)
        return readdict


    def parse_read_dict(self):
        '''
        Parse correct read dict depending on input data type.
        '''
        if self.krakenlike is not None:
            readdict = self.get_krakenlike_read_dict()
        else:
            readdict = self.get_read_dict()
        return readdict


    def get_taxlist_dict(self):
        '''
        Example output: {readid : [Domain, Phylum, Class..]}
        Read a bblast like file and link taxdb identifyers
        to taxdb taxonomy in refdict to store in taxlist dictionary.
        '''
        readdict = self.parse_read_dict()
        refdict = self.get_reftax_dict()
        taxdict = defaultdict(list)
        for readid in readdict:
            refidlist = readdict[readid]
            for refid in refidlist:
                if refid in refdict and len(refdict[refid]) > 5:
                    ''' If already in dict, compare for each tax
                    level whether there is concensus. If not replace
                    with NA. '''
                    if readid in taxdict:
                        for level in range(5, -1, -1):
                            if taxdict[readid][level] == refdict[refid][level]:
                                break
                            else:
                                taxdict[readid][level] = "NA"
                    else:
                        taxdict[readid] = refdict[refid]
        return taxdict


    def write_taxlist(self):
        '''
        Example out: "readid \t Domain \t Phylum ..."
        Write taxlist to file, incl. header line.
        '''
        taxlist = self.get_taxlist_dict()
        with open(self.outpath + ".taxlist", 'w') as f:
            f.write('#readid' + '\t' + 'Domain' + '\t' +
                    'Phylum' + '\t' + 'Class' + '\t' +
                    'Order' + '\t' + 'Family' + '\t' +
                    'Genus' + '\n')
            for readid in taxlist:
                tax = taxlist[readid]
                f.write(readid + '\t' + '\t'.join(tax) + '\n')


    def get_taxmat_dict(self):
        '''
        Example out: {taxstr: [Domain,Phylum,...]}
        Read taxlist, concatenate taxonomy into taxid
        and create lib with key taxid and value taxonomy.
        '''
        with open(self.taxlist) as f:
            taxdict = {}
            for l in f.readlines():
                if l.startswith("#") or not l.strip():
                    continue
                else:
                    line = l.strip().split('\t')
                    tax = line[1:7]
                    taxid = ';'.join(tax)
                    taxdict[taxid] = tax
        return taxdict


    def write_taxmat(self):
        '''
        Example out: "taxstr \t Domain \t Phylum ..."
        Write taxmat to file, incl. header line.
        '''
        with open(self.outpath + ".taxmat", 'w') as f:
            f.write('#taxid' + '\t' + 'Domain' + '\t' +
                    'Phylum' + '\t' + 'Class' + '\t' +
                    'Order' + '\t' + 'Family' + '\t' +
                    'Genus' + '\n')
            for taxid in self.get_taxmat_dict():
                tax = self.get_taxmat_dict()[taxid]
                taxlab = re.sub(r'[_NA]*_NA$', '', taxid)
                taxstr = '\t'.join(tax)
                taxstr = re.sub(r'[\tNA]*\tNA$', '', taxstr)
                taxstr = re.sub(r'_', ' ', taxstr)
                f.write(taxlab + '\t' + taxstr + '\n')


    def get_otumat_dict(self):
        '''
        Example out: {taxstr: count}
        Read taxlist, concatenate taxonomy into taxid
        and create lib witk key taxid and value count
        where count is the frequency each taxid occured
        in taxlist
        '''
        with open(self.taxlist) as f:
            otudict = defaultdict(int)
            for l in f.readlines():
                if l.startswith("#") or not l.strip():
                    continue
                else:
                    line = l.strip().split('\t')
                    tax = line[1:7]
                    taxid = ';'.join(tax)
                    otudict[taxid] += 1
        return otudict


    def write_otumat(self):
        '''
        Example output: "taxstr \t count"
        Write otumat to file, incl. header line.
        '''
        with open(self.outpath + ".otumat", 'w') as f:
            f.write('#taxid' + '\t' + self.method + '_' + self.sample + '\n')
            for taxid in self.get_otumat_dict():
                count = self.get_otumat_dict()[taxid]
                taxlab = re.sub(r'[_NA]*_NA$', '', taxid)
                f.write(taxlab + '\t' + str(count) + '\n')


def main(taxlist, bblastlike, centrifugelike, krakenlike, dbtaxonomy, dbfasta, dbmap):
    m = ToMat(
        taxlist = taxlist,
        bblastlike = bblastlike,
        centrifugelike = centrifugelike,
        krakenlike = krakenlike,
        dbtaxonomy = dbtaxonomy,
        dbfasta = dbfasta,
        dbmap = dbmap
    )
    if any([bblastlike, centrifugelike, krakenlike]):
        m.write_taxlist()
    m.write_taxmat()
    m.write_otumat()
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert output tables '
                        'from classifyers in nanoclass snakemake pipeline '
                        'to a taxonomy and otu matrix. These serve to '
                        'generate taxonomic barplots. If not given, the '
                        'script also produces a taxlist which serves as input '
                        'to assign consensus and majority classifications.')
    parser.add_argument('-l', '--taxlist', default=None,
                        help='Path to a file with read IDs in the '
                        'first column, followed by at least 6 columns '
                        'with tab-delimited 7-level taxonomic classifications '
                        'sensu Silva.')
    parser.add_argument('-b', '--bblastlike', default=None,
                        help='Path to a best-blast-like output file '
                        'with the first column listing the readID and '
                        'the second column the taxonomic identifyer.')
    parser.add_argument('-c', '--centrifugelike', default=None,
                        help='Path to a centrifuge-like output file. '
                        'First column should contain read ID and second '
                        'the taxonomic identifyer.')
    parser.add_argument('-k', '--krakenlike', default=None,
                        help='Path to a kraken-like output file. '
                        'Second column should contain read ID and third '
                        'the ncbi taxonomic identifyer.')
    parser.add_argument('-t', '--dbtaxonomy', default=None,
                        help='Path to the taxonomy of the database used '
                        'to get taxonomic classifications. Only required '
                        'in combination with --bblastlike.')
    parser.add_argument('-f', '--dbfasta', default=None,
                        help='Path to the fasta of database that was used '
                        'to get taxonomic classifications. Fasta headers '
                        'should contain taxids and taxo string. Only required '
                        'in combination with --centrifugelike.')
    parser.add_argument('-m', '--dbmap', default=None,
                        help='Path to the map of database that was used '
                        'to get taxonomic classifications. Map contains '
                        'link between taxonomy in third column of '
                        '--krakenlike and taxids in --dbfasta. Only '
                        'required in combination with --krakenlike.')
    args = parser.parse_args()

    assert any([args.taxlist, args.bblastlike, args.centrifugelike,
                args.krakenlike]), '''
        \n Please provide the path to the taxlist or bblastlike.
        Read the help (--help) for further information.
        '''

    main(
        taxlist = args.taxlist,
        bblastlike = args.bblastlike,
        centrifugelike = args.centrifugelike,
        krakenlike = args.krakenlike,
        dbtaxonomy = args.dbtaxonomy,
        dbfasta = args.dbfasta,
        dbmap = args.dbmap
        )
