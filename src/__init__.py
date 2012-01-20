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

import datetime
import logging
import os
import urllib

__version__ = '0.1.0'

GOURL = 'http://geneontology.org/ontology/obo_format_1_2/gene_ontology_ext.obo'

logging.basicConfig()
LOG = logging.getLogger('golib')
LOG.setLevel(logging.INFO)

class PyGoLibException(Exception):
    """ This is the class used for all potential exception which can be
    generated while running this project.
    """
    pass


class PyGoLib(object):
    """ Utility class with some functions to play with the graph."""

    def __init__(self, graph=None):
        """ Constructor. """
        self.graph = graph
        if not self.graph:
            self.graph = {}
        self.log = get_logger()
        self.subgraph = {}

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
        if verbose:
            print " " * level, "\_", termid
        if '!' in termid:
            parentid = termid.split('!')[0].strip()
        else:
            parentid = termid
        if pred == "":
            pred = parentid
        else:
            pred = pred + "," + parentid
        parent = self.graph[parentid]
        try:
            parent = self.get_path(parent, level=level + 1,
                pred=pred, paths=paths, verbose=verbose)
        except KeyError:
            paths.append(pred)

    def get_sub_graph(self, graph, termid, verbose=False):
        """ From the list of GO terms, retrieve all the one which have
        for parent the provided termid.
        """
        for key in graph.keys():
            term = graph[key]
            if verbose:
                print term['id']
            pathways = self.get_path(term, pred=term['id'], paths=[],
                verbose=verbose)
            for el in pathways:
                if termid in el:
                    for tid in el.split(termid)[0].split(','):
                        if tid:
                            term = self.graph[tid]
                            if term['id'] not in self.subgraph.keys():
                                self.subgraph[term['id']] = term
        return self.subgraph

    def get_path(self, term, level=0, pred="", paths=[], verbose=False):
        """ This is an iterative method which is used to retrieve the top
        parent of a given term.
        """
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


def get_logger():
    """ Return the logger. """
    return LOG


def set_logger(quiet=False, debug=False):
    """ Set the logger level. """
    if debug:
        LOG.setLevel(logging.DEBUG)
    elif quiet:
        LOG.setLevel(logging.WARNING)


def download_GO_graph(outputfile=None, force_dl=False):
    """ Retrieve the GO data from the specified file on the
    filesystem is provided or from the web or using the local
    version if dated from the day.

    :kwarg outputfile, name of the file in which the geneontology will
    be saved. If not provided it will be of the form:
    'geneontology-DATE.obo'
    :kwarg force_dl, boolean to force the (re)download of the GO
    annotation file from the geneontology.org website. Defaults to
    False.
    """
    LOG.debug('download_GO_graph: outputfile %s - force_dl %s' %
        (outputfile, force_dl))
    if not outputfile:
        go_file = 'geneontology-%s.obo' % \
            datetime.datetime.now().strftime('%Y%m%d')
    else:
        go_file = outputfile
    LOG.debug('Downloading GO term into: %s' % go_file)
    if force_dl or not os.path.exists(go_file):
        LOG.info('Retrieving GO from %s' % GOURL)
        urllib.urlretrieve(GOURL, go_file)
    else:
        LOG.info(
            '%s already exists, no need to re-download it' % go_file)
    return go_file
