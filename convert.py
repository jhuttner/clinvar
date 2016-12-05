#!/usr/bin/env python2.7

import json
import Queue
import random
import sys
import threading

from bs4 import BeautifulSoup
from lxml import etree


PATH_DELIMITER = '/'
INPUT_FILENAME = '/tmp/ClinVarFullRelease_00-latest.xml'
OUTPUT_FILE = open('/tmp/out.txt', 'w')

# 1 for 100%, 0.1 for 10%, 0.01 for 1%, etc.
SAMPLE = 0.01 
QUEUE = Queue.Queue()
TITLE = None


class Extended(object):

    def __init__(self, soup):
        self.soup = soup

    def split_path(self, path):
        return path.split(PATH_DELIMITER)

    def find(self, path):
        if isinstance(path, str):
            path = self.split_path(path)

        result = self.soup

        for dir_ in path:
            result = result.find(dir_)

        result = Extended(result)
        return result

    def find_all(self, path):
        if isinstance(path, str):
            path = self.split_path(path)

        result = self.soup

        for index, dir_ in enumerate(path):
            # Last item - call findAll so we return a list
            if index == len(path) - 1:
                result = result.findAll(dir_)
            else:
                result = result.find(dir_)

            # Needed to address data inconsistency in the dataset for
            # MeasureSet/Measure/MeasureRelationship/XRef

            if result is None:
                return []

        result = [Extended(i) for i in result]
        return result
    
    def get_abs(self, path):
        pass

    def get_attr(self, attr):
        return self.soup[attr]

    def get_content(self):
        return self.soup.string


def assertion_to_json(assertion):
    # print 'assertion is ', assertion
    result = {}

    # rcvaccession
    result['rcvaccession'] = (assertion.find('ClinVarAccession')
                                       .get_attr('Acc'))

    # rcvaccession_version
    result['rcvaccession_version'] = (assertion.find('ClinVarAccession')
                                               .get_attr('Version'))

    # rcvaccession_full
    result['rcvaccession_full'] = '.'.join([result['rcvaccession'],
                                            result['rcvaccession_version']])

    # uuid
    result['uuid'] = result['rcvaccession_full']

    # title
    result['title'] = TITLE

    # preferred_name
    result['preferred_name'] = (assertion.find('MeasureSet/Measure/Name'
                                               '/ElementValue')
                                         .get_content())

    # hgvs
    hgvs_values = []
    all_attr_sets = assertion.find_all('MeasureSet/Measure/AttributeSet')

    for attr_set in all_attr_sets:
        all_attrs = attr_set.find_all('Attribute')
        for attr in all_attrs:
            types = [i.strip() for i in attr.get_attr('Type').split(',')]
            if 'HGVS' in types:
                hgvs_values.append(attr.get_content())

    result['hgvs'] = hgvs_values

    # type
    result['type'] = (assertion.find('MeasureSet/Measure').get_attr('Type'))

    # clinical_significance
    result['clinical_significance'] = (assertion.find('ClinicalSignificance'
                                                      '/Description')
                                                .get_content())

    # entrez_gene_id
    xrefs = assertion.find_all('MeasureSet/Measure/MeasureRelationship/XRef')
    entrez_gene_id_values = []

    for xref in xrefs:
        if xref.get_attr('DB') == 'Gene':
            entrez_gene_id_values.append(xref.get_attr('ID'))

    result['entrez_gene_id'] = entrez_gene_id_values

    # rs_id
    xrefs = assertion.find_all('MeasureSet/Measure/XRef')
    rs_id_values = []

    for xref in xrefs:
        if xref.get_attr('DB') == 'dbSNP':
            rs_id_values.append('rs_' + xref.get_attr('ID'))

    result['rs_id'] = rs_id_values

    return json.dumps(result)


class TitleTarget(object):
    def __init__(self):
        self.text = []

    def start(self, tag, attrib):
        self.is_assertion = (tag == 'ReferenceClinVarAssertion')

    def end(self, tag):
        pass

    def data(self, data):
        if self.is_assertion:
            print 'appending', data
            self.text.append(data.encode('utf-8'))

    def close(self):
        return self.text


# Copied from: http://stackoverflow.com/questions/4695826
def fast_iter(context, func, *args, **kwargs):
    """
    fast_iter is useful if you need to free memory while iterating through a
    very large XML file.

    http://lxml.de/parsing.html#modifying-the-tree
    Based on Liza Daly's fast_iter
    http://www.ibm.com/developerworks/xml/library/x-hiperfparse/
    See also http://effbot.org/zone/element-iterparse.htm
    """
    for event, elem in context:
        func(elem, *args, **kwargs)
        # It's safe to call clear() here because no descendants will be
        # accessed
        elem.clear()
        # Also eliminate now-empty references from the root node to elem
        for ancestor in elem.xpath('ancestor-or-self::*'):
            while ancestor.getprevious() is not None:
                del ancestor.getparent()[0]
    del context


def enqueue_element(elt):
    if random.random() <= SAMPLE:
        text = etree.tostring(elt)
        QUEUE.put(text)


def process_element(text):
    try:
        soup = Extended(BeautifulSoup(text, 'xml'))
        OUTPUT_FILE.write(assertion_to_json(soup) + '\n')
    except:
        raise


def worker():
    while True:
        item = QUEUE.get()
        if item is None:
            break
        process_element(item)
        QUEUE.task_done()


def set_title():
    '''
    This is a hack to get the <Title> out.  

    There is some better way to do this using lxml, but I gave up trying
    to figure it out.

    TODO:Remove this hack
    '''
    global TITLE
    with open(INPUT_FILENAME) as f:
        lines = [f.readline() for _ in range(100)]
        soup = Extended(BeautifulSoup('\n'.join(lines), 'xml'))
        element = soup.find('ReleaseSet/ClinVarSet/Title')
        if element is None:
            raise Exception('Could not find <Title> in first 100 lines')
        else:
            TITLE = element.get_content()


def main(args):
    set_title()

    threads = []

    # Increasing this value beyond 1 doesn't seem to make a difference for this
    # workload, likely because this work is CPU-bound.  I've left the code in
    # just to demonstrate what I attempted.
    num_worker_threads = 2 

    for i in range(num_worker_threads):
        t = threading.Thread(target=worker)
        t.thread_index = 1
        t.setDaemon(True)
        t.start()
        threads.append(t)

    context = etree.iterparse(INPUT_FILENAME, events=('end',), 
                              tag='ReferenceClinVarAssertion')
    fast_iter(context, enqueue_element)

    # block until all tasks are done
    QUEUE.join()

    # stop workers
    for i in range(num_worker_threads):
        QUEUE.put(None)
    for t in threads:
        t.join()

    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
