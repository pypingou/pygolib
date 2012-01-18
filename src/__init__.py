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

import logging

__version__ = '0.1.0'


logging.basicConfig()
LOG = logging.getLogger('golib')


def get_logger():
    """ Return the logger. """
    return LOG

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
