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

def _get_ancester(path1, path2):
    """ For two given path, return the first common ancester.
    :arg path1, a list of nodes.
    :arg path2, a list of nodes.
    """
    for step in path1:
        if step in path2:
            return step
    return None

def _get_all_ancesters(path):
    """ For a given path, return all the ancesters in it. """
    ancesters = []
    for step in path:
        for item in step.split(','):
            if item not in ['is_a', 'part_of'] \
                and item not in ancesters:
                ancesters.append(item)
    return ancesters


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
        self.pygo = PyGoLib(self.goterms)

    def semantic_value(self, id1):
        """ Returns the semantic values of all the parents of a given
        term.
        :arg id1, identifier of a GO term (ie: GO:0043231, or whatever
            identifier is in your ontology).
        """
        sem_values = self.semantic_values(id1)
        return sum(sem_values.values())

    def semantic_values(self, id1):
        """ Returns the semantic values of all the parents of a given
        term.
        :arg id1, identifier of a GO term (ie: GO:0043231, or whatever
            identifier is in your ontology).
        """
        golib = PyGoLib(self.goterms)
        goterm1 = self.goterms[id1]
        path1 = golib.get_path(goterm1, pred=goterm1['id'], paths=[],
            details=True)

        semantic_values = {}
        for ancester in _get_all_ancesters(path1):
            if ancester == goterm1['id']:
                semantic_values[ancester] = 1
                continue
            tmp_score = []
            for item in path1:
                tmp_cnt = 1
                if ancester in item:
                    path_el = item.split(',')
                    ind = path_el.index(ancester)
                    for step in range(ind -1, 0, -2):
                        if path_el[step] == 'is_a':
                            tmp_cnt = tmp_cnt * 0.8
                        elif path_el[step] == 'part_of':
                            tmp_cnt = tmp_cnt * 0.6
                    tmp_score.append(tmp_cnt)
            semantic_values[ancester] = max(tmp_score)
        return semantic_values

    def scores(self, id1, id2):
        """Returns the score between two given GO terms.
        :arg id1, identifier of a GO term (ie: GO:0043231, or whatever
            identifier is in your ontology).
        :arg id2, identifier of a GO term (ie: GO:0043229, or whatever
            identifier is in your ontology).
        """
        golib = PyGoLib(self.goterms)
        #golib.fix_go_graph()
        goterm1 = self.goterms[id1]
        path1 = golib.get_path(goterm1, pred=goterm1['id'], paths=[],
            details=True)
        #print goterm1['id'], len(path1)
        ancester1 = _get_all_ancesters(path1)
        semantic_values1 = self.semantic_values(goterm1['id'])

        goterm2 = self.goterms[id2]
        path2 = golib.get_path(goterm2, pred=goterm2['id'], paths=[],
            details=True)
        #print goterm2['id'], len(path2)
        ancester2 = _get_all_ancesters(path2)
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
        :arg goid, identifier of a GO term (ie: GO:0043229, or whatever
            identifier is in your ontology).
        :arg golist, list of GO terms
        """
        scores = []
        sesamego = GsesameGO(self.goterms)
        for goterm in golist:
            scores.append(sesamego.scores(goid, goterm))
        return max(scores)

    def scores(self, gene1, gene2):
        """ For two list of GO term associated with two genes, computes
        the semantic similarities of the genes.
        :arg gene1, list of GO term associated with the gene1
        :arg gene2, list of GO term associated with the gene2
        """
        sim1 = 0
        for goterm in gene1:
            sim1 = sim1 + self.__get_go_score(goterm, gene2)
        sim2 = 0
        for goterm in gene2:
            sim2 = sim2 + self.__get_go_score(goterm, gene1)
        score = (sim1 + sim2) / (len(gene1) + len(gene2))
        return score

if __name__ == '__main__':
    from oboio import OboIO
    from src import download_go_graph

    OBIO = OboIO()
    ONTO = download_go_graph()
    TERMS = OBIO.get_graph(ONTO)

    # G-Sesame GO example
    GSGO = GsesameGO(TERMS)
    print 'GO:0043229 semantic value:', GSGO.semantic_value('GO:0043229')
    print 'GO:0043231 semantic value:', GSGO.semantic_value('GO:0043231')
    print 'Score for: GO:0043229 - GO:0043231:', \
        GSGO.scores('GO:0043229','GO:0043231')

    # G-Sesame Gene example
    GSGENE = GsesameGene(TERMS)
    GENE1 = ['GO:0004022', 'GO:0004024', 'GO:0004174', 'GO:0046872',
        'GO:0008270', 'GO:0004023']
    GENE2 = ['GO:0009055', 'GO:0005515', 'GO:0046872', 'GO:0008270',
        'GO:0020037']
    print 'Semantic similarities between the two genes:', \
        GSGENE.scores(GENE1, GENE2)
