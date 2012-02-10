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


import os
import sys

try:
    from pygolib import get_logger, PyGoLib
except ImportError:
    sys.path.insert(0, os.path.abspath('../'))
    from src import get_logger, PyGoLib


class GsesameGO(object):
    """ This class re-implement in python the algorithm used in the
    g-sesame program to compare two GO term to each other.
    see: http://bioinformatics.clemson.edu/G-SESAME/
    """

    def __init__(self, data=None):
        """ Constructor.
        :arg data, the graph of ontologies
        """
        self.goterms = data
        if self.goterms is None:
            self.goterms = {}
        self.log = get_logger()

    def __get_all_ancesters(self, path):
        """ For a given path, return all the ancesters in it. """
        ancesters = []
        for el in path:
            for item in el.split(','):
                if item not in ['is_a', 'part_of'] \
                    and item not in ancesters:
                    ancesters.append(item)
        return ancesters

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
        if mindist != None and deltalevel != None:
            return mindist + deltalevel / 10.0
        else:
            return None

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

    def semantic_value(self, id1):
        """ Returns the semantic values of all the parents of a given
        term.
        :arg id1, identifier of a GO term (ie: GO:XXX).
        """
        sem_values = self.semantic_values(id1)
        return sum(sem_values.values())

    def semantic_values(self, id1):
        """ Returns the semantic values of all the parents of a given
        term.
        :arg id1, identifier of a GO term (ie: GO:XXX).
        """
        golib = PyGoLib(self.goterms)
        goterm1 = self.goterms[id1]
        path1 = golib.get_path(goterm1, pred=goterm1['id'], paths=[],
            details=True)

        semantic_values= {}
        for ancester in self.__get_all_ancesters(path1):
            if ancester == goterm1['id']:
                semantic_values[ancester] = 1
                continue
            tmp_score = []
            for item in path1:
                tmp_cnt = 1
                if ancester in item:
                    path_el = item.split(',')
                    ind = path_el.index(ancester)
                    for el in range(ind -1, 0, -2):
                        if path_el[el] == 'is_a':
                            tmp_cnt = tmp_cnt * 0.8
                        elif path_el[el] == 'part_of':
                            tmp_cnt = tmp_cnt * 0.6
                    tmp_score.append(tmp_cnt)
            semantic_values[ancester] = max(tmp_score)
        return semantic_values

    def scores(self, id1, id2):
        """Returns the score between two given GO terms.
        :arg id1, identifier of a GO term (ie: GO:XXX).
        :arg id2, identifier of a GO term (ie: GO:XXX).
        """
        golib = PyGoLib(self.goterms)
        #golib.fix_GO_graph()
        goterm1 = self.goterms[id1]
        path1 = golib.get_path(goterm1, pred=goterm1['id'], paths=[],
            details=True)
        #print goterm1['id'], len(path1)
        ancester1 = self.__get_all_ancesters(path1)
        semantic_values1 = self.semantic_values(goterm1['id'])

        goterm2 = self.goterms[id2]
        path2 = golib.get_path(goterm2, pred=goterm2['id'], paths=[],
            details=True)
        #print goterm2['id'], len(path2)
        ancester2 = self.__get_all_ancesters(path2)
        semantic_values2 = self.semantic_values(goterm2['id'])
        
        common_ancester = list(set(ancester1).intersection(set(ancester2)))
        sum_comm_anc = 0
        for ancester in common_ancester:
            sum_comm_anc = sum_comm_anc + semantic_values2[ancester] + \
                semantic_values1[ancester]

        score = sum_comm_anc / (sum(semantic_values1.values()) 
            + sum(semantic_values2.values()))
        return score

class GsesameGene(object):
    """ This class re-implement in python the algorithm used in the
    g-sesame program to compare two genes based on their GO annotation.
    see: http://bioinformatics.clemson.edu/G-SESAME/
    """

    def __init__(self, data=None):
        """ Constructor.
        :arg data, the graph of ontologies
        """
        self.goterms = data
        if self.goterms is None:
            self.goterms = {}
        self.log = get_logger()

    def __get_go_score(self, goid, golist):
        """ For a given GO term return the semantic similarity between
        this GO term and the given GO list.
        :arg goid, GO term identifier (ie: GO:XXXXX)
        :arg golist, list of GO terms
        """
        scores = []
        sesamego = GsesameGO(self.goterms)
        for go in golist:
            scores.append(sesamego.scores(goid, go))
        return max(scores)

    def scores(self, gene1, gene2):
        """ For two list of GO term associated with two genes, computes
        the semantic similarities of the genes.
        :arg gene1, list of GO term associated with the gene1
        :arg gene2, list of GO term associated with the gene2
        """
        sim1 = 0
        for go in gene1:
            sim1 = sim1 + self.__get_go_score(go, gene2)
        sim2 = 0
        for go in gene2:
            sim2 = sim2 + self.__get_go_score(go, gene1)
        score = (sim1 + sim2) / (len(gene1) + len(gene2))
        return score

if __name__ == '__main__':
    from oboio import OboIO
    from src import download_GO_graph
    from src import get_logger, PyGoLib

    obio = OboIO()
    ontology = download_GO_graph()
    terms = obio.get_graph(ontology)

    # G-Sesame GO example
    gdc = GsesameGO(terms)
    print 'GO:0043229 semantic value:', gdc.semantic_value('GO:0043229')
    print 'GO:0043231 semantic value:', gdc.semantic_value('GO:0043231')
    print 'Score for: GO:0043229 - GO:0043231:', \
        gdc.scores('GO:0043229','GO:0043231')

    # G-Sesame Gene example
    sesamegene = GsesameGene(terms)
    gene1 = ['GO:0004022', 'GO:0004024', 'GO:0004174', 'GO:0046872',
        'GO:0008270', 'GO:0004023']
    gene2 = ['GO:0009055', 'GO:0005515', 'GO:0046872', 'GO:0008270',
        'GO:0020037']
    print 'Semantic similarities between the two genes:', \
        sesamegene.scores(gene1, gene2)
