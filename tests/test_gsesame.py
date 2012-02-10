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
Unit-tests for the goutil library.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.abspath('../'))
from src import PyGoLib
from src.gsesame import GsesameGO, GsesameGene
from src.oboio import OboIO

if os.path.dirname(__file__):
    folder = os.path.dirname(__file__)
else:
    folder = '.'
# Ontology for the term GO:0043231 used in the G-Sesame article from 2007:
# James Z. Wang, Zhidian Du, Rapeeporn Payattakool, Philip S. Yu and
# Chin-Fu Chen, A New Method to Measure the Semantic Similarity of GO
# Terms, Bioinformatics, 2007, 23: 1274-1281;
GOFILE = '%s/test2.obo' % folder
# Ontology for the term GO:0043231 at 2012-02-10
GOFILE2 = '%s/test3.obo' % folder


class GsesameGOTests(unittest.TestCase):
    """ GsesameGO tests. """

    def __init__(self, methodName='runTest'):
        """ Constructor. """
        unittest.TestCase.__init__(self, methodName)

    def test_semantic_value(self):
        """ Test the semantic_value function. """
        obio = OboIO()
        terms = obio.get_graph(GOFILE2)
        gsgo = GsesameGO(terms)
        output = 4.7440000000000015
        self.assertEqual(output, gsgo.semantic_value('0043229'))
        output = 5.5952
        self.assertEqual(output, gsgo.semantic_value('0043231'))

    def test_semantic_value_article(self):
        """ Test the semantic_value function. """
        obio = OboIO()
        terms = obio.get_graph(GOFILE)
        gsgo = GsesameGO(terms)
        output = 3.4000000000000004
        self.assertEqual(output, gsgo.semantic_value('0043229'))
        output = 4.520000000000001
        self.assertEqual(output, gsgo.semantic_value('0043231'))

    def test_scores(self):
        """ Test the scores function. """
        obio = OboIO()
        terms = obio.get_graph(GOFILE2)
        gsgo = GsesameGO(terms)
        output = 0.8259052924791086
        self.assertEqual(output, gsgo.scores('0043229','0043231'))

    def test_scores_article(self):
        """ Test the scores function. """
        obio = OboIO()
        terms = obio.get_graph(GOFILE)
        gsgo = GsesameGO(terms)
        output = 0.7727272727272726
        self.assertEqual(output, gsgo.scores('0043229','0043231'))


class GsesameGeneTests(unittest.TestCase):
    """ GsesameGene tests. """

    def __init__(self, methodName='runTest'):
        """ Constructor. """
        unittest.TestCase.__init__(self, methodName)

    def test_scores(self):
        """ Test the scores function. """
        obio = OboIO()
        terms = obio.get_graph(GOFILE2)
        sesamegene = GsesameGene(terms)
        gene1 = ['0043229', '0044424']
        gene2 = ['0043231', '0043227']
        output = 0.6743128041470686
        self.assertEqual(output, sesamegene.scores(gene1, gene2))


suite = unittest.TestLoader().loadTestsFromTestCase(GsesameGOTests)
suite = unittest.TestLoader().loadTestsFromTestCase(GsesameGeneTests)
unittest.TextTestRunner(verbosity=2).run(suite)
