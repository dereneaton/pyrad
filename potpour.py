#!/usr/bin/env python2

""" extra class objects """

import multiprocessing
#import cPickle as pickle

 
class Worker(multiprocessing.Process):
    """ multiprocessing object """
 
    def __init__(self, work_queue, result_queue, func):
 
        # base class initialization
        multiprocessing.Process.__init__(self)
 
        # job management stuff
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False
        self.func = func
 
    def run(self):
        while not self.kill_received:
            # get a task
            if self.work_queue.empty():
                break
            else:
                #job = self.work_queue.get_nowait()
                job = self.work_queue.get()
 
            # the actual processing
            res = self.func(*job)
 
            # store the result
            self.result_queue.put(res)



class Locusobj():
    """ store all the depth data """
    def __init__(self, init_id):
        ## init assignments
        self.init_id = init_id
        self.final_id = ""

        ## data
        self.consobjs = []

        ## post align
        self.trimr = None
        self.triml = None
        self.indels = []
        self.filter = []

        ## clustering records
        self.seed = ""
        self.hits = []
        self.orients = []

    # def alignment():
    #     """returns the .loci alignment"""
    #     ## create locus

    #     ## trim edges

    #     ## make snpstring
    #     pass

class Consobj():
    """ store all the depth data """
    def __init__(self):
        ## init assignments
        self.seq = ""
        self.Cs = []
        self.As = []
        self.Ts = []
        self.Gs = []

