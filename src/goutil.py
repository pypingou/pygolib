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
        if pred == "":
            pred = parentid
        else:
            pred = pred + "," + parentid
        parent = self.goterms[parentid]
        try:
            parent = self.get_path(parent, level=level + 1,
                pred=pred, paths=paths, verbose=verbose)
        except KeyError:
            paths.append(pred)

    def __get_ancester(self, path1, path2):
        """ For two given path, return the first common ancester.
        :arg path1, a list of nodes.
        :arg path2, a list of nodes.
        """
        for el in path1:
            if el in path2:
                return el
        return None

    def __score_cousins(self, goid1, goid2, path1=None, path2=None):
        """ For two given GO term ID and the list of their path, return
        the score between them.
        :arg goid1, GO term ID (ie: GO:XXXX).
        :arg goid2, GO term ID (ie: GO:XXXX).
        :kwarg path1, the list path from the first GO term to the top of
        the tree, as returned by get_path().
        :kwarg path2, the list path from the second GO term to the top
        of the tree, as returned by get_path().
        """
        if path1 is None:
            path1 = self.get_path(self.goterms[goid1])
        if path2 is None:
            path2 = self.get_path(self.goterms[goid2])

        print goid1, path1
        print goid2, path2
        ancester = None
        mindist = None
        deltalevel = None
        for path in path1:
            step1 = path.split(',')
            for opath in path2:
                step2 = opath.split(',')
                inter = self.__get_ancester(step1, step2)
                if inter:
                    index1 = step1.index(inter)
                    index2 = step2.index(inter)
                    dist = index1 + index2
                    deltaleveltmp = abs(index1 - index2)
                    if not mindist or dist < mindist:
                        mindist = dist
                        ancester = inter
                        deltalevel = deltaleveltmp
        print ancester, deltalevel, mindist
        return mindist + deltalevel / 10.0

    def __score_parents(self, goid1, goid2, path1=None, path2=None):
        """ For two given GO term ID and the list of their path, return
        the score between them if one is parent of the other.
        :arg goid1, GO term ID (ie: GO:XXXX).
        :arg goid2, GO term ID (ie: GO:XXXX).
        :kwarg path1, the list path from the first GO term to the top of
        the tree, as returned by get_path().
        :kwarg path2, the list path from the second GO term to the top
        of the tree, as returned by get_path().
        """
        #print goid1, path1
        #print goid2, path2
        if path1 is None:
            path1 = self.get_path(self.goterms[goid1])
        scores = []
        for paths in path1:
            el = paths.split(',')
            if goid2 in el:
                start = el.index(goid1)
                stop = el.index(goid2)
                score = abs(stop - start)
                score = score + score / 10.0
                scores.append(score)
        if path2 is None:
            path2 = self.get_path(self.goterms[goid2])
        for paths in path2:
            el = paths.split(',')
            if goid1 in el:
                start = el.index(goid1)
                stop = el.index(goid2)
                score = abs(stop - start)
                score = score + score / 10.0
                scores.append(score)
        return scores

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
        #self.get_biological_process()
        #print "%s terms found in the biological process branch" % \
            #len(self.biological_process.keys())
        #outputfile =  'geneontology_biological_process-%s.obo' % \
            #datetime.datetime.now().strftime('%Y%m%d')
        #write_down_ontology(self.biological_process, outputfile)

        #term = self.goterms['11']
        #term = self.goterms['4']
        #print term['id']
        #self.get_path(term, pred=term['id'], paths=[], verbose=True)

        self.scores('11', '0')

        endtime  = datetime.datetime.now()
        print "Time spent: ", endtime - starttime, "minutes"
        self.godata=None

    def scores(self, id1, id2):
        """Returns the score between two given GO terms.
        :arg id1, identifier of a GO term (ie: GO:XXX).
        :arg id2, identifier of a GO term (ie: GO:XXX).
        """
        goterm1 = self.goterms[id1]
        path1 = self.get_path(goterm1, pred=goterm1['id'], paths=[])
        goterm2 = self.goterms[id2]
        path2 = self.get_path(goterm2, pred=goterm2['id'], paths=[])
        scores = self.__score_parents(id1, id2, path1, path2)
        if scores:
            print "%s and %s are parents" % (id1, id2)
            print "The score between %s and %s is: %s" % (id1, id2, min(scores))
        else:
            scores = self.__score_cousins(id1, id2, path1, path2)
            if scores:
                print "The score between %s and %s is: %s" % (id1, id2, scores)
        

if __name__ == '__main__':
    gdc = GoDistanceCounter()
    #gdc.main()
    #gdc.main(go_file='geneontology-20120117.all.obo')
    gdc.main('test.obo')
    #gdc.main('multi_is_a.obo')