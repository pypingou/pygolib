#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This project is licensed under the New BSD License:

Copyright (c) 2012, Pierre-Yves Chibon

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.
* Neither the name of the Wageningen University nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ''AS IS'' AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.
"""

"""
This script computes the distance between two GO terms.
This distance is computed as a combination of the number of edges
between the two terms of the tree weighted by one tens of the level
difference between the two terms.
"""


import datetime
import os
import urllib

GOURL = 'http://geneontology.org/ontology/obo_format_1_2/gene_ontology_ext.obo'


def write_down_ontology(ontodata, datafile):
    """ Writes to disk the provided ontology.
    :arg ontodata, the hash of hash containing all the ontology
    :arg datafile, the name of the file to which write the ontology.
    """
    stream = open(datafile, 'w')
    cnt = 0
    for key in ontodata.keys():
        stream.write('\n[Term] \n')
        cnt = cnt + 1
        info = ontodata[key]
        for infokey in info.keys():
            if isinstance(info[infokey], list):
                for entry in info[infokey]:
                    stream.write(infokey + ': ' + entry + '\n')
            else:
                stream.write(infokey + ': ' + info[infokey] + '\n')
    stream.close()
    print '%s terms written to the file %s' %(cnt, datafile)


class GoDistanceCounterException(Exception):
    """ This is the class used for all potential exception which can be
    generated while running this project.
    """
    pass


class GoDistanceCounter(object):
    """ This class is the main class of the project to compute the
    distance between two GO terms.
    """

    def __init__(self):
        """ Constructor. """
        self.godata = None
        self.goterms = {}
        self.biological_process = {}

    def __do_handle_parent(self, termid, level, pred, paths, verbose=False):
        """ Handle the output for one parent of a term.
        It will retrieve the GO term for the given ID, add it to the pred
        building up the path and keep building the upper part of the 
        tree.
        If verbose is True, it will print the line of the tree.

        :arg termid, GO term identifier (GO:XXXX)
        :arg level, the level in which we are while building the tree
        :arg pred, the precedant part of the paths browsed
        :arg paths, the list of paths already browsed
        :kwarg verbose, a boolean to actually print the tree or not.
        """
        if verbose: print " " * level, "\_", termid
        if '!' in termid:
            parentid = termid.split('!')[0].strip()
        else:
            parentid = termid
        pred = pred + "," + parentid
        parent = self.goterms[parentid]
        try:
            parent = self.get_path(parent, level=level + 1,
                pred=pred, paths=paths, verbose=verbose)
        except KeyError:
            paths.append(pred)

    def get_biological_process(self, verbose=False):
        """ From the list of GO terms, retrieve all the one which have
        for parent 'GO:0008150: biological_process'.
        """
        for key in self.goterms.keys():
            term = self.goterms[key]
            if verbose: print term['id']
            pathways = self.get_path(term, pred=term['id'], paths=[],
                verbose=verbose)
            for el in pathways:
                #print "  ", el
                if el.endswith('GO:0008150') and \
                        term['id'] not in self.biological_process.keys():
                    self.biological_process[term['id']] = term

    def get_go_data(self, go_file, force_dl=False):
        """ Retrieve the GO data from the specified file on the
        filesystem is provided or from the web or using the local
        version if dated from the day.

        :arg go_file, file on the filesystem containing the GO
        annotation.
        :kwarg force_dl, boolean to force the (re)download of the GO
        annotation file from the geneontology.org website. Defaults to
        False.
        """
        possible_go_file = 'geneontology-%s.obo' % \
            datetime.datetime.now().strftime('%Y%m%d')
        if not force_dl and not go_file \
                and not os.path.exists(possible_go_file):
            print "Retrieving GO from %s" % GOURL
            urllib.urlretrieve(GOURL, possible_go_file)
            go_file = possible_go_file

        if not go_file and os.path.exists(possible_go_file):
            go_file = possible_go_file
        
        if not os.path.exists(go_file):
            raise GoDistanceCounterException(
                'No file containing the GO ontology found or provided.')

        print 'Reading file: %s' % go_file
        gostream = open(go_file)
        self.godata = gostream.read()
        gostream.close()

    def get_go_terms(self):
        """ From the GO annotation file, extract all the GO terms and
        store them into dictionnary.
        """
        print 'Loading GO terms...'
        for entry in self.godata.split("\n\n"):
            if '[Term]' in entry:
                info = {}
                rows = entry.split('\n')
                for row in rows[1:]:
                    if row and ':' in row:
                        (key, value) = row.split(':', 1)
                        key = key.strip()
                        if key in info.keys():
                            if isinstance(info[key], str):
                                info[key] = [info[key], value.strip()]
                            elif isinstance(info[key], list):
                                info[key].append(value.strip())
                        else:
                            info[key] = value.strip()
                if info['id'] not in self.goterms.keys():
                    self.goterms[info['id']] = info

        print "%s GO terms retrieved" % len(self.goterms.keys())

    def get_path(self, term, level=0, pred="", paths=[], verbose=False):
        """ This is an iterative method which is used to retrieve the top
        parent of a given term.
        """
        #print term['id'], paths, pred
        if 'is_a' in term.keys():
            if isinstance(term['is_a'], list):
                before = pred
                for p in term['is_a']:
                    self.__do_handle_parent(p, level, before, paths,
                        verbose=verbose)
            else:
                self.__do_handle_parent(term['is_a'], level, pred, paths,
                    verbose=verbose)
            return paths
        else:
            paths.append(pred)
            return paths

    def main(self, go_file=None):
        """ Main method of the class.
        Download the latest GO annotation and save it to the current
        directory if no file was specified.
        Reads in the GO annotation.
        Builds the tree.
        """
        starttime = datetime.datetime.now()
        self.get_go_data(go_file)
        self.get_go_terms()
        self.get_biological_process()
        print "%s terms found in the biological process branch" % \
            len(self.biological_process.keys())
        outputfile =  'geneontology_biological_process-%s.obo' % \
            datetime.datetime.now().strftime('%Y%m%d')
        write_down_ontology(self.biological_process, outputfile)

        endtime  = datetime.datetime.now()
        print "Time spent: ", endtime - starttime, "minutes"
        self.godata=None


if __name__ == '__main__':
    gdc = GoDistanceCounter()
    #gdc.main()
    gdc.main(go_file='geneontology-20120117.all.obo')
    #gdc.main('multi_is_a.obo')
