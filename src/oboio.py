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

import sys
import os

try:
    from pygolib import get_logger
except ImportError:
    sys.path.insert(0, os.path.abspath('../'))
    from src import get_logger

class OboIO (object):
    """ This class handles the reading and writing of OBO files. """

    def __init__(self, graph=None):
        """ Constructor. """
        self.graph = graph
        if self.graph is None:
            self.graph = {}
        self.log = get_logger()

    def download_GO_graph(self, force_dl=False):
        """ Retrieve the GO data from the specified file on the
        filesystem is provided or from the web or using the local
        version if dated from the day.

        :kwarg force_dl, boolean to force the (re)download of the GO
        annotation file from the geneontology.org website. Defaults to
        False.
        """
        possible_go_file = 'geneontology-%s.obo' % \
            datetime.datetime.now().strftime('%Y%m%d')
        if not force_dl and not os.path.exists(possible_go_file):
            self.log.info("Retrieving GO from %s" % GOURL)
            urllib.urlretrieve(GOURL, possible_go_file)
        else:
            self.log.info(
                "%s already exists, no need to re-download it" % GOURL)
        return possible_go_file

    def get_go_terms(self, filename):
        """ From the GO annotation file, extract all the GO terms and
        store them into dictionnary.
        """
        stream = open(filename)
        data = stream.read()
        stream.close()
        self.log.info('Loading GO terms...')
        for entry in data.split("\n\n"):
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
                if info['id'] not in self.graph.keys():
                    self.graph[info['id']] = info

        self.log.info("%s GO terms retrieved" % len(self.graph.keys()))
        return self.graph

    def write_down_ontology(self, datafile):
        """ Writes graph to disk.
        :arg datafile, the name of the file to which write the ontology.
        """
        stream = open(datafile, 'w')
        cnt = 0
        for key in self.graph.keys():
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
        print '%s terms written to the file %s' % (cnt, datafile)
