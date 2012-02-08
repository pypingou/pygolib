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
from src.godistance import GoDistanceCounter
from src.oboio import OboIO

GOFILE = '%s/test.obo' % os.path.dirname(__file__)

class GoDistanceCounterTests(unittest.TestCase):
    """ GoDistanceCounter tests. """

    def __init__(self, methodName='runTest'):
        """ Constructor. """
        unittest.TestCase.__init__(self, methodName)

    def test_get_go_terms(self):
        """ Test the get_go_terms function. """
        obio = OboIO()
        terms = obio.get_graph(GOFILE)

    def test_get_path(self):
        """ Test the get_path function. """
        obio = OboIO()
        terms = obio.get_graph(GOFILE)
        golib = PyGoLib(terms)
        term = terms['4']
        print term['id']
        output = ['4,2,1,0']
        self.assertEqual(output, 
            golib.get_path(term, pred=term['id'], paths=[], verbose=True))
        self.assertEqual(output, 
            golib.get_path(term, pred=term['id'], paths=[]))
        output = ['2,1,0']
        self.assertEqual(output, 
                golib.get_path(term, paths=[]))

    def test_scores(self):
        """ Test the scores function. """
        obio = OboIO()
        terms = obio.get_graph(GOFILE)
        gdc = GoDistanceCounter(terms)
        self.assertEqual(5.5, gdc.scores('11', '0'))
        self.assertEqual(4.4, gdc.scores('11', '1'))
        self.assertEqual(6.0, gdc.scores('9', '5'))
        self.assertEqual(5.1, gdc.scores('7', '5'))
        self.assertEqual(6.0, gdc.scores('8', '5'))

    def test_alt_id(self):
        """ Test that a GO term with an alt_id is added to the list. """
        obio = OboIO()
        terms = obio.get_graph(GOFILE)
        gdc = GoDistanceCounter(terms)
        term = terms['12']
        self.assertEqual('7', term['id'])
        term = terms['13']
        self.assertEqual('8', term['id'])
        gdc = GoDistanceCounter(terms)
        self.assertEqual(5.1, gdc.scores('12', '5'))
        self.assertEqual(6.0, gdc.scores('13', '5'))

suite = unittest.TestLoader().loadTestsFromTestCase(GoDistanceCounterTests)
unittest.TextTestRunner(verbosity=2).run(suite)
