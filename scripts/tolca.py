#!/usr/bin/env python3

__author__ = "Evelien Jongepier"
__copyright__ = "Copyright 2021, Evelien Jongepier"
__email__ = "e.jongepier@uva.nl"
__license__ = "MIT"

import sys
import re
import argparse
import collections



class ToLCA(object):

    def __init__(self, blast, taxonomy, lca, pconsensus):
        self.blast = blast
        self.taxonomy = taxonomy
        self.lca = lca
        self.pconsensus = pconsensus
        self.levels = ["Domain","Phylum","Class","Order","Family","Genus"]
        return


    def get_blastdict(self):
        '''
        Example output: readdict = {read_id: [hitid1, hitid2, 
            hitid3],...}
        Stores list of blast hit id for each read in dict.
        '''
        readdict = collections.defaultdict(list)
        with open (self.blast) as f:
            for l in f.readlines():
                readid = l.strip().split('\t')[0]
                hitid = l.strip().split('\t')[1]
                if readid not in readdict:
                    readdict[readid] = [hitid]
                else:
                    readdict[readid].append(hitid)
        return readdict


    def get_dbdict(self):
        '''
        Example output: dbdict = {dbid: [domain, order,
            class, ...], ...}
        Stores list of taxonomic annotations (i.e. 7 levels)
        for each refernce in the db.
        '''
        dbdict = collections.defaultdict(list)
        with open (self.taxonomy) as f:
            for l in f.readlines():
                dbid, dbtaxstr = l.strip().split('\t')[:2]
                dbtaxstr = re.sub(r'D_[0-6]__', '', dbtaxstr)
                dbtax = dbtaxstr.split(';')[:6]
                dbdict[dbid] = dbtax
        return dbdict


    def get_classdict(self):
        '''
        Example output: classdict = {readid: {domain: [hit1, 
            hit2, hit3, ...], ...}, ...}
        Stores list of classifications for each of the 7
        levels of the taxonomy for each dbhit of each read.
        Redundant hits, not yet consensus.
        '''

        classdict = collections.defaultdict(
            dict = collections.defaultdict(list)
        )

        readdict = self.get_blastdict()
        dbdict = self.get_dbdict()

        for readid in readdict:
              for hitid in readdict[readid]:
                   dbtaxs = dbdict[hitid]

                   if readid not in classdict:
                       classdict[readid] = dict(
                           (self.levels[idx], [tax]) for idx, tax in enumerate(dbtaxs))
                   else:
                       for idx, tax in enumerate(dbtaxs):
                           classdict[readid][self.levels[idx]].append(tax)
        return classdict


    def comp_lca(self):
        '''
        Example output: {readid = {Domain:cons},{Phylum:cons},...}
        Compute consensus for each read at each taxonomic level
        based on user defined level (0.5 = majority) and the taxo-
        nomic classificatons of each hit.
        '''

        lcadict = collections.defaultdict(dict)
        classdict = self.get_classdict()

        for readid in classdict:
            for level in classdict[readid]:
                freqdict = collections.Counter(classdict[readid][level])
                totfreq = sum(freqdict.values())
                maxkey = max(freqdict, key=freqdict.get)
                if (freqdict[maxkey] / float(totfreq)) > self.pconsensus:
                    lcadict[readid][level] = maxkey
                else:
                    lcadict[readid][level] = "NA"
        return lcadict


    def write_lca(self):
        '''
        Example out: "readid \t Domain \t Phylum ..."
        Write consensus to file, incl. header line.
        '''
        lcadict = self.comp_lca()
        with open(self.lca, 'w') as f:
            f.write('#readid' + '\t' + '\t'.join(self.levels) + '\n')
            for readid in lcadict:
                tax = []
                for level in self.levels:
                    tax.append(lcadict[readid][level])
                f.write(readid + '\t' + '\t'.join(tax) + '\n')
        return


def main(blast, taxonomy, lca, pconsensus):
    m = ToLCA(
        blast = blast,
        taxonomy = taxonomy,
        lca = lca,
        pconsensus = pconsensus
    )
    m.write_lca()
    return      


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Obtain LCA based on blast '
                        'tabular and user defined consensus treshold.')
    parser.add_argument('-b', '--blast', default=None,
                        help='Path to blast tabular with read id and reference '
                        'hit id.')
    parser.add_argument('-t', '--taxonomy',  default=None,
                        help='Path to taxonomy with reference id and 7-level '
                        'taxonomy sensu Silva.')
    parser.add_argument('-l', '--lca', default=None,
                        help='Path to where output lca consensus taxonomy '
                        'should be written.')
    parser.add_argument('-c', '--pconsensus', type=float, default=0.5,
                        help='Consensus treshold ranging between 0.5 and 0.99. '
                        'Minimum fraction of studies that need to show '
                        'consensus about each taxonomic classification at '
                        'each taxonomic level. Majority obtained with 0.5. '
                        'Taxonomic annotations for which there is no consensus '
                        'are returned as NA.')
    args = parser.parse_args()

    main(
	blast = args.blast,
        taxonomy = args.taxonomy,
        lca = args.lca,
        pconsensus = args.pconsensus
        )
