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
