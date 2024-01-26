"""
Created on Oct 10, 2012

@author: smirarab
"""
from sepp.algorithm import AbstractAlgorithm
from sepp.config import options
from sepp.tree import PhylogeneticTree
from sepp.alignment import MutableAlignment, ExtendedAlignment,\
    hamming_distance
from sepp.problem import SeppProblem, RootProblem
from dendropy.datamodel.treemodel import Tree
from sepp.jobs import HMMBuildJob, HMMSearchJob, HMMAlignJob, PplacerJob,\
    MergeJsonJob
from sepp.scheduler import JobPool, Join, Job
from sepp import get_logger
from sepp.math_utils import lcm
import copy
import random 
import time
import numpy as np
_LOG = get_logger(__name__)


def get_placement_job_name(chunk_number):
    return "pplacer_%d" % chunk_number


class QuerySetJob(Job):
    def __init__(self, fc):
        Job.__init__(self)
        self.fc = fc
        self.qseqs = None
        self.query_seqs = None
        self.kmers_to_check = None
        self.set_to_seqs = None
        self.top_sets_to_check = None
        self.rounds = None
        self.int_score = None
        
        self.qseq_to_K = None
        self.K_to_set_kmers = None
        self.K_to_set_num_kmers = None
        
    def setup_problem(self, qseqs, query_seqs, qseq_to_K, K_to_set_kmers, K_to_set_num_kmers, set_to_seqs, use_fraction_batch, batch_size_fraction, min_batch_size, kmers_to_check, top_sets_to_check, rounds, int_score, exact):
        self.qseqs = qseqs    
        self.query_seqs = query_seqs
        self.qseq_to_K = qseq_to_K
        self.K_to_set_kmers = K_to_set_kmers
        self.K_to_set_num_kmers = K_to_set_num_kmers
        self.set_to_seqs = set_to_seqs
        self.use_fraction_batch = use_fraction_batch 
        self.batch_size_fraction = batch_size_fraction 
        self.min_batch_size = min_batch_size
        self.kmers_to_check = kmers_to_check
        self.top_sets_to_check = top_sets_to_check   
        self.rounds = rounds
        self.int_score = int_score
        self.exact = exact   

    def run(self):
        qseqs = self.qseqs
        query_seqs = self.query_seqs
        set_to_seqs = self.set_to_seqs
        kmers_to_check = self.kmers_to_check
        set_to_seqs = self.set_to_seqs
        top_sets_to_check = self.top_sets_to_check
        rounds = self.rounds
        int_score = self.int_score
        exact = self.exact
        
        qseq_to_K = self.qseq_to_K
        K_to_set_kmers = self.K_to_set_kmers 
        K_to_set_num_kmers = self.K_to_set_num_kmers 
        
        
        qseq_too_short = {}        

        
        num_sets = len(K_to_set_kmers[list(K_to_set_kmers.keys())[0]])
        qseq_to_kmers = {}
        if len(qseq_to_K) == 0:
            K = list(K_to_set_kmers.keys())[0]
            for qseq in qseqs:
                if len(query_seqs[qseq]) < K:
                    qseq_too_short[qseq] = True
                    continue
                else:
                    qseq_to_kmers[qseq] = {}
                    for i in range(0, len(query_seqs[qseq]) - K + 1):
                        kmer = query_seqs[qseq][i:i+K] 
                        if kmer in qseq_to_kmers[qseq]:
                            qseq_to_kmers[qseq][kmer] += 1
                        else:
                            qseq_to_kmers[qseq][kmer] = 1   
        else:
            for qseq in qseqs:
                if len(query_seqs[qseq]) < K:
                    qseq_too_short[qseq] = True
                else:            
                    qseq_to_kmers[qseq] = {}
                    K = qseq_to_K[qseq]
                    for i in range(0, len(query_seqs[qseq]) - K + 1):
                        kmer = query_seqs[qseq][i:i+K] 
                        if kmer in qseq_to_kmers[qseq]:
                            qseq_to_kmers[qseq][kmer] += 1
                        else:
                            qseq_to_kmers[qseq][kmer] = 1               
        

        
        #print(K_to_set_num_kmers)
        #print(qseq_to_K)
        query_set_comparison = {}
        #print(qseq)
        #query_set_comparison_data[qseq] = {}
        #print("num_sets: ", num_sets)

        K = list(K_to_set_kmers.keys())[0]
        for qseq in qseq_to_kmers:
            if len(qseq_to_K) != 0:
            
                K = qseq_to_K[qseq]
            
            
            if qseq in qseq_too_short:
                query_set_comparison[qseq] = {}    
                for set_num in range(num_sets):
                    # the first entry of this vector should not be 0, but it doesn't matter because it's not used anyway.
                    query_set_comparison[qseq][set_num] = [0, 0, 0]
                
            else:    
                set_to_kmers = K_to_set_kmers[K]
                set_to_num_kmers = K_to_set_num_kmers[K]
                
                if self.use_fraction_batch:
                    kmers_to_check = max(self.min_batch_size, int(np.ceil(self.batch_size_fraction*len(query_seqs[qseq]))))
                    #print(len(query_seqs[qseq]), kmers_to_check)
                
                #print(qseq)
                #print(K)
                #print(" ")
                #print(set_to_kmers)
                #print(" ")
                #print(qseq_to_kmers[qseq])
                if (not exact) and (num_sets > top_sets_to_check):
                    cands = {}
                    for set_num in range(num_sets):
                        cands[set_num]  = 0
                    
                    qseq_cands = {}
                    cands_sorted = {}
                    #qseq_to_kmers_sorted = sorted(qseq_to_kmers[qseq].items(), key = lambda kv: kv[1], reverse=True)
                    num_kmers_in_q = len(qseq_to_kmers[qseq])
                    
                    for rnd in range(rounds):
                        #print(qseq, len(cands), top_sets_to_check)
                        if len(cands) > top_sets_to_check:
                            
                            # sample without replacement
                            kmers_to_check2 = min(kmers_to_check, len(qseq_to_kmers[qseq]))
                            qseq_to_kmers_sample = random.sample(qseq_to_kmers[qseq].keys(), kmers_to_check2)
                            
                            # sample with replacement
                            #kmers_to_check2 = kmers_to_check
                            #qseq_to_kmers_sample = random.choices([i for i in range(len(query_seqs[qseq]) - K + 1)], k = kmers_to_check2)
                            for set_num in cands:
                                intersection1 = 0
                                for i1 in range(kmers_to_check2):
                                    #kmer1 = query_seqs[qseq][qseq_to_kmers_sample[i1]: qseq_to_kmers_sample[i1] + K]
                                    kmer1 = qseq_to_kmers_sample[i1]
                                    if kmer1 in set_to_kmers[set_num]:
                                        intersection1 += min(set_to_kmers[set_num][kmer1] / len(set_to_seqs[set_num]), qseq_to_kmers[qseq][kmer1]) 
                                if int_score: 
                                    cands[set_num] += intersection1 * (num_kmers_in_q / kmers_to_check2)
                                else:
                                    cands[set_num] += intersection1 * (num_kmers_in_q / kmers_to_check2)  /  ((len(query_seqs[qseq]) - K + 1)  + set_to_num_kmers[set_num]/ len(set_to_seqs[set_num]))
                                    #print(qseq, set_num, (len(query_seqs[qseq]) - K + 1), set_to_num_kmers[set_num], len(set_to_seqs[set_num]),  (len(query_seqs[qseq]) - K + 1) / kmers_to_check2,  ((len(query_seqs[qseq]) - K + 1)  + set_to_num_kmers[set_num]/ len(set_to_seqs[set_num])), ((len(query_seqs[qseq]) - K + 1) / kmers_to_check2)  /  ((len(query_seqs[qseq]) - K + 1)  + set_to_num_kmers[set_num]/ len(set_to_seqs[set_num])))
                            cands_sorted = sorted(cands.items(), key = lambda kv: kv[1], reverse=True)   
                            new_cands = {}
                            for i1 in range(max(int(len(cands_sorted)/2), top_sets_to_check)):
                                new_cands[cands_sorted[i1][0]] = cands_sorted[i1][1]
                            #print(qseq, rnd, cands)
                            #print(" ")
                            cands = new_cands
                
                    # for rnd in range(rounds):
                    #     if len(cands) > top_sets_to_check:
                    #         kmers_to_check2 = min(kmers_to_check, len(qseq_to_kmers[qseq]))
                    #         qseq_to_kmers_sample = random.sample(qseq_to_kmers[qseq].keys(), kmers_to_check2)
                    #         for set_num in cands:
                    #             intersection1 = 0
                    #             for i1 in range(kmers_to_check2):
                    #                 #kmer1 = qseq_to_kmers_sorted[i1][0]
                    #                 kmer1 = qseq_to_kmers_sample[i1]
                    #                 if kmer1 in set_to_kmers[set_num]:
                    #                     intersection1 += min(set_to_kmers[set_num][kmer1] / len(set_to_seqs[set_num]), qseq_to_kmers[qseq][kmer1])
                    #             cands[set_num] += intersection1
                    #         cands_sorted = sorted(cands.items(), key = lambda kv: kv[1], reverse=True)   
                    #         new_cands = {}
                    #         for i1 in range(max(int(len(cands_sorted)/2), top_sets_to_check)):
                    #             new_cands[cands_sorted[i1][0]] = cands_sorted[i1][1]
                    #         cands = new_cands
                            
                        
                    qseq_cands[qseq] = cands_sorted            
                    
                    
                    query_set_comparison[qseq] = {}
                    #for set_num in range(num_sets):
                    for ts in range(top_sets_to_check):  
                        set_num = qseq_cands[qseq][ts][0]  
                        
                        intersection = {}
                        qseq_minus_set = {}
                        set_minus_qseq = {}
                        total_in_qseq = 0
                        
                        total_in_set = 0
                        for kmer in set(set_to_kmers[set_num].keys()).union(set(qseq_to_kmers[qseq].keys())):
                            if kmer in set_to_kmers[set_num] and kmer in qseq_to_kmers[qseq]:
                                s_norm = set_to_kmers[set_num][kmer]/len(set_to_seqs[set_num])
                                q_norm = qseq_to_kmers[qseq][kmer]
                                intersection[kmer] = min(s_norm, q_norm)
                                if s_norm > q_norm:
                                    set_minus_qseq[kmer] = s_norm - q_norm
                                elif q_norm > s_norm:
                                    qseq_minus_set[kmer] = q_norm - s_norm
                                total_in_qseq += q_norm
                                total_in_set += s_norm
                            
                            elif kmer in set_to_kmers[set_num]:
                                s_norm = set_to_kmers[set_num][kmer]/len(set_to_seqs[set_num])
                                set_minus_qseq[kmer] = s_norm 
                                total_in_set += s_norm
                                
                            elif kmer in qseq_to_kmers[qseq]:
                                q_norm = qseq_to_kmers[qseq][kmer]
                                qseq_minus_set[kmer] = q_norm
                                total_in_qseq += q_norm
                                    
                        total_in_union = total_in_qseq + total_in_set
                        total_in_intersection =  sum([intersection[kmer] for kmer in intersection])
                        
                        #query_set_comparison_data[qseq][set_num] = [intersection, qseq_minus_set, set_minus_qseq] 
                        
                        query_set_comparison[qseq][set_num] = [total_in_union, total_in_intersection, total_in_intersection/total_in_union]
               
                
                else:
                    query_set_comparison[qseq] = {}
                    
                    for set_num in range(num_sets):
                        
                        intersection = {}
                        qseq_minus_set = {}
                        set_minus_qseq = {}
                        total_in_qseq = 0
                        
                        total_in_set = 0
                        for kmer in set(set_to_kmers[set_num].keys()).union(set(qseq_to_kmers[qseq].keys())):
                            if kmer in set_to_kmers[set_num] and kmer in qseq_to_kmers[qseq]:
                                s_norm = set_to_kmers[set_num][kmer]/len(set_to_seqs[set_num])
                                q_norm = qseq_to_kmers[qseq][kmer]
                                intersection[kmer] = min(s_norm, q_norm)
                                if s_norm > q_norm:
                                    set_minus_qseq[kmer] = s_norm - q_norm
                                elif q_norm > s_norm:
                                    qseq_minus_set[kmer] = q_norm - s_norm
                                total_in_qseq += q_norm
                                total_in_set += s_norm
                            
                            elif kmer in set_to_kmers[set_num]:
                                s_norm = set_to_kmers[set_num][kmer]/len(set_to_seqs[set_num])
                                set_minus_qseq[kmer] = s_norm 
                                total_in_set += s_norm
                                
                            elif kmer in qseq_to_kmers[qseq]:
                                q_norm = qseq_to_kmers[qseq][kmer]
                                qseq_minus_set[kmer] = q_norm
                                total_in_qseq += q_norm
                                    
                        total_in_union = total_in_qseq + total_in_set
                        total_in_intersection =  sum([intersection[kmer] for kmer in intersection])
                        
                        #query_set_comparison_data[qseq][set_num] = [intersection, qseq_minus_set, set_minus_qseq] 
                        
                        query_set_comparison[qseq][set_num] = [total_in_union, total_in_intersection, total_in_intersection/total_in_union]

                
        query_assignment = {}
        if int_score:
            for qseq in query_set_comparison:
                best_seen = None
                best_val = 0
                for set_num in query_set_comparison[qseq]:
                    if query_set_comparison[qseq][set_num][1] > best_val:
                        best_seen = set_num
                        best_val = query_set_comparison[qseq][set_num][1]
                        #print(best_seen, best_val)
                query_assignment[qseq] = best_seen
        else:
            for qseq in query_set_comparison:
                best_seen = None
                best_val = 0
                for set_num in query_set_comparison[qseq]:
                    if query_set_comparison[qseq][set_num][2] > best_val:
                        best_seen = set_num
                        best_val = query_set_comparison[qseq][set_num][2]
                        #print(best_seen, best_val)
                query_assignment[qseq] = best_seen
        _LOG.info("Finished Fast Search Job with input: " + str(self.fc) + ",     time: " + str(time.time()))
        return query_assignment



# class TestJob(Job):
#     def run(self):
#         print("hi")
#         return 1 + 1
    
# def factorial(num):
#     print("hi2")
#     return num + 2

        
# class JoinFastSearchJobs(Join):  
#     def __init__(self):
#         Join.__init__(self)
#         self.done1 = False
#     def perform(self):
#         self._lock.acquire()
#         self.done1 = True
#         print("Fast query searches are done")
#         print(self.done1)
#         self._lock.release()

class JoinSearchJobs(Join):
    """
    After all search jobs have finished on tips, we need to figure out which
    fragment goes  to which subset and start aligning fragments.
    This join takes care of that step.
    """
    def __init__(self):
        Join.__init__(self)
        self.root_problem = None
        self.fc_jobs = None

    def setup_with_root_problem(self, root_problem, fc_jobs):
        self.root_problem = root_problem
        self.fc_jobs = fc_jobs 
        for p in root_problem.iter_leaves():
            self.add_job(p.jobs["hmmsearch"])

    def figureout_fragment_subset_og(self):
        """ Figure out which fragment should go to which subproblem"""
        if "fragments.distribution.done" in self.root_problem.annotations:
            return
        max_evalues = dict(
            [(name, (None, None)) for name in
             self.root_problem.fragments.keys()])
        for fragment_chunk_problem in self.root_problem.iter_leaves():
            align_problem = fragment_chunk_problem.get_parent()
            assert isinstance(align_problem, SeppProblem)
            '''For each subproblem start with an empty set of fragments,
            and add to them as we encounter new best hits for that
            subproblem'''
            if align_problem.fragments is None:
                align_problem.fragments = \
                    self.root_problem.fragments.get_soft_sub_alignment([])
            search_res = fragment_chunk_problem.get_job_result_by_name(
                "hmmsearch")
            for key, val in search_res.items():
                (best_value, prev_align_problem) = max_evalues[key]
                ''' If this is better than previous best hit, remove this
                fragment from the previous hit, and add it to this subproblem
                '''
                if best_value is None or (best_value < val[1]):
                    max_evalues[key] = (val[1], align_problem)

        # TODO: is the following efficient enough? Do we need to make lists
        # and then turn them to sets?
        not_scored = []
        for key, v in max_evalues.items():
            if v[1] is None:
                not_scored.append(key)
            else:
                v[1].fragments.seq_names.add(key)

        self.root_problem.annotations["fragments.distribution.done"] = 1

        ''' Make sure all fragments are in at least one subproblem.
        TODO: what to do with those that are not?  For now, only output
        warning message'''
        # notScored = [k for k, v in max_evalues.iteritems() if v[1] is None]
        _LOG.warning(
            "Fragments %s are not scored against any subset" % str(not_scored))
        # assert len(notScored) == 0, "Fragments %s are not scored against
        # any subset" %str(notScored)
                


    def perform(self):
        """
        Distributes fragments to alignments subsets with best score,
        and runs align jobs on those. Also, creates new chunks of fragments
        for better parallelism.
        """

        ''' Figure out which fragment should go to which subproblem'''
        """ Figure out which fragment should go to which subproblem"""
        if "fragments.distribution.done" in self.root_problem.annotations:
            return

        _LOG.info("Beginning Query Set Assignments, time: " + str(time.time()))
        
        fragments = self.root_problem.fragments     
        fragment_keys = fragments.keys()
        
        backbone_seqs = copy.deepcopy(self.root_problem.subalignment)
        backbone_seqs.degap()
        
        seqs = {}
        avg_frag_len = 0
        for seq in fragments:
            seqs[seq] = fragments[seq]
            avg_frag_len += len(fragments[seq])
        for seq in backbone_seqs:
            seqs[seq] = backbone_seqs[seq]
        avg_frag_len = avg_frag_len / len(fragments)
        
        set_to_seqs = {}
        align_problems = []
        #K = 10  #4
        K = options().kmer_size
        
        use_fraction_batch = options().use_fraction_batch
        batch_size_fraction = options().batch_size_fraction
        min_batch_size = options().min_batch_size
        
        kmers_to_check = options().kmers_to_check
        top_sets_to_check = options().top_sets_to_check
        rounds = options().sample_rounds
        int_score = options().int_score
        exact = options().exact
        choose_K = options().choose_K
        Ks_to_check = [int(i) for i in options().Ks_to_check.split(",")]
        

        # if use_fraction_batch:
        #     kmers_to_check = max(min_batch_size, np.ceil(batch_size_fraction*avg_frag_len))
            
        
        
        i = 0
        for fragment_chunk_problem in self.root_problem.iter_leaves():
            align_problem = fragment_chunk_problem.get_parent()
            if align_problem not in align_problems:
                assert isinstance(align_problem, SeppProblem)
                '''For each subproblem start with an empty set of fragments'''
                if align_problem.fragments is None:
                    align_problem.fragments = \
                        self.root_problem.fragments.get_soft_sub_alignment([])
                align_problems.append(align_problem)
                set_to_seqs[i] = {}
                for seq in align_problem.taxa:
                    set_to_seqs[i][seq] = True
                i += 1
            
        num_sets = len(align_problems)
        
        query_seqs = {} 
        for frag in fragment_keys:
            query_seqs[frag] = seqs[frag]       
        
        
        
                
        # if choose_K:
        #     bestK = None
        #     Ks_to_check.sort(reverse = True)
        #     for i in range(len(Ks_to_check)):
        #         K1 = Ks_to_check[i]
        #         if i < len(Ks_to_check) - 1:
        #             backbone_to_kmers = {}
        #             for set_num in range(num_sets):
        #                 for seq_ID in set_to_seqs[set_num]:
        #                     seq = seqs[seq_ID]
        #                     for j in range(0, len(seq) - K1 + 1):
        #                         backbone_to_kmers[seq[j:j+K1]] = True 
    
        #             num_unmatched = 0
        #             for qseq in fragment_keys:
        #                 qseq_seq = query_seqs[qseq]
        #                 match1 = False
        #                 for j in range(0, len(qseq_seq) - K1 + 1):
        #                     if qseq_seq[j:j+K1] in backbone_to_kmers:
        #                         match1 = True
        #                         break
        #                 if match1 == False:
        #                     num_unmatched += 1
        #             print("K = ", K1, num_unmatched / len(query_seqs) )
        #             if num_unmatched / len(query_seqs) < seqs_unmatched_thresh:
        #                 bestK = K1
        #                 break
        #         else:
        #             bestK = K1
        #     K = bestK
                    
                            
        K_to_set_kmers = {}
        K_to_set_num_kmers = {}
        Ks_used = {}
        qseq_to_K = {}
        K_to_qseq = {}
        if choose_K:
            Ks_to_check.sort(reverse = True)
            K_to_kmers = {}
            for K1 in Ks_to_check:
                K_to_qseq[K1] = []
                backbone_to_kmers = {}
                for set_num in range(num_sets):
                    for seq_ID in set_to_seqs[set_num]:
                        seq = seqs[seq_ID]
                        for j in range(0, len(seq) - K1 + 1):
                            backbone_to_kmers[seq[j:j+K1]] = True 
                K_to_kmers[K1] = backbone_to_kmers
            for qseq in fragment_keys:
                qseq_seq = query_seqs[qseq]
                for i in range(len(Ks_to_check)):
                    K1 = Ks_to_check[i]
                    if i < len(Ks_to_check) - 1:
                        match1 = False
                        for j in range(0, len(qseq_seq) - K1 + 1):
                            if qseq_seq[j:j+K1] in K_to_kmers[K1]:
                                match1 = True
                                break        
                        if match1 == True:
                            qseq_to_K[qseq] = K1
                            Ks_used[K1] = True
                            K_to_qseq[K1].append(qseq)
                            break
                    else:
                        qseq_to_K[qseq] = K1
                        Ks_used[K1] = True
                        K_to_qseq[K1].append(qseq)
            print("")
            for K1 in Ks_to_check:
                print("K = " + str(K1) + ": " + str(len(K_to_qseq[K1])))
            print("")
        else:
            Ks_used[K] = True
                        
                        
                        
                        
        for K1 in Ks_used:
            set_to_num_kmers = {} 
            for set_num in range(num_sets):
                set_to_num_kmers[set_num] = 0
                for seq in  set_to_seqs[set_num]:
                    set_to_num_kmers[set_num] += len(seqs[seq]) - K1 + 1
            
            set_to_kmers = {}
            for set_num in range(num_sets):
                set_to_kmers[set_num] = {}
                for seq_ID in set_to_seqs[set_num]:
                    seq = seqs[seq_ID]
                    for i in range(0, len(seq) - K1 + 1):
                        if seq[i:i+K1-1] in set_to_kmers[set_num]:
                            set_to_kmers[set_num][seq[i:i+K1]] += 1
                        else:
                            set_to_kmers[set_num][seq[i:i+K1]] = 1
            
            K_to_set_kmers[K1] = set_to_kmers
            K_to_set_num_kmers[K1] = set_to_num_kmers
                
            
            
             
        
        for set1 in set_to_seqs:
            print(set1, len(set_to_seqs[set1].keys()))
        

        #print(K_to_set_kmers)
            
    
        fc_fragments = {}    
        i = 0
        for fragment_chunk_problem in self.root_problem.get_children()[0].get_children()[0].get_children():
            fc_fragments[i] = list(MutableAlignment().read_filepath(fragment_chunk_problem.fragments).keys())
            i +=1
            
        for fc in self.fc_jobs:
            fcjob = self.fc_jobs[fc]
            fc_seqs = fc_fragments[fc]
            # print(fc)
            # print(fcjob)
            # print(fc_seqs)
            # print(" ")
            fcjob.setup_problem(fc_seqs, query_seqs, qseq_to_K, K_to_set_kmers, K_to_set_num_kmers, set_to_seqs, use_fraction_batch, batch_size_fraction, min_batch_size, kmers_to_check, top_sets_to_check, rounds, int_score, exact)
            #fcjob.setup_problem(fc_seqs, query_seqs, set_to_seqs, set_to_kmers, set_to_num_kmers, K, kmers_to_check, top_sets_to_check, rounds, int_score, exact)
        
        
        # while True:
        #     time.sleep(1)
        
        for fc in self.fc_jobs:
            JobPool().enqueue_job(self.fc_jobs[fc])

    def __str__(self):
        return "join search jobs for all tips of ", self.root_problem

class JoinFastSearchJobs(Join):
    """
    After all search jobs have finished on tips, we need to figure out which
    fragment goes  to which subset and start aligning fragments.
    This join takes care of that step.
    """
    def __init__(self):
        Join.__init__(self)
        self.root_problem = None
        

    def setup_with_root_problem(self, root_problem, fc_jobs):
        self.root_problem = root_problem
        self.fc_jobs = fc_jobs
        for fc in fc_jobs:
            self.add_job(fc_jobs[fc])

        
    def perform(self):
        """
        Distributes fragments to alignments subsets with best score,
        and runs align jobs on those. Also, creates new chunks of fragments
        for better parallelism.
        """

        ''' Figure out which fragment should go to which subproblem'''
        fragments = self.root_problem.fragments     
        fragment_keys = fragments.keys()
        
        backbone_seqs = copy.deepcopy(self.root_problem.subalignment)
        backbone_seqs.degap()
        
        seqs = {}
        for seq in fragments:
            seqs[seq] = fragments[seq]
        for seq in backbone_seqs:
            seqs[seq] = backbone_seqs[seq]

        
        #set_to_seqs = {}
        align_problems = []
        #K = options().kmer_size
        #kmers_to_check = options().kmers_to_check
        #top_sets_to_check = options().top_sets_to_check
             
        
        i = 0
        for fragment_chunk_problem in self.root_problem.iter_leaves():
            align_problem = fragment_chunk_problem.get_parent()
            if align_problem not in align_problems:
                assert isinstance(align_problem, SeppProblem)
                '''For each subproblem start with an empty set of fragments'''
                if align_problem.fragments is None:
                    align_problem.fragments = \
                        self.root_problem.fragments.get_soft_sub_alignment([])
                align_problems.append(align_problem)
                # set_to_seqs[i] = {}
                # for seq in align_problem.taxa:
                #     set_to_seqs[i][seq] = True
                # i += 1
            
        #num_sets = len(align_problems)

        query_seqs = {} 
        for frag in fragment_keys:
            query_seqs[frag] = seqs[frag]


        fc_jobs = self.fc_jobs
        
        query_assignment = {}
        set_to_queries = {}
        for fc in fc_jobs:
            fc_job = fc_jobs[fc]
            fc_job_res = fc_job.result
            for qseq in fc_job_res:
                best_seen = fc_job_res[qseq]
                query_assignment[qseq] = best_seen
                if best_seen in set_to_queries:
                    set_to_queries[best_seen].append(qseq)
                else:
                    set_to_queries[best_seen] = [qseq]
        
        print("")
        print(set_to_queries)
        print("")
        #print(align_problems)
        
        not_scored = []
        for i in range(len(align_problems)):
            if i in set_to_queries:
                for qseq in set_to_queries[i]:
                    align_problems[i].fragments.seq_names.add(qseq)

        if None in set_to_queries:        
            for qseq in set_to_queries[None]:
                not_scored.append(qseq)
                align_problems[0].fragments.seq_names.add(qseq)
                
        # for qseq in query_assignment:
        #     if query_assignment[qseq] == None:
        #         not_scored.append(qseq)
         
        self.root_problem.annotations["fragments.distribution.done"] = 1
        
        for set_num in set_to_queries:
            print(set_num, ":", len(set_to_queries[set_num]))
        print("")

        ''' Make sure all fragments are in at least one subproblem.
        TODO: what to do with those that are not?  For now, only output
        warning message'''
        # notScored = [k for k, v in max_evalues.iteritems() if v[1] is None]
        
        _LOG.info("Finished Query Set Assignments, time: " + str(time.time()))
        _LOG.warning(
            "Fragments %s are not scored against any subset" % str(not_scored))
        # assert len(notScored) == 0, "Fragments %s are not scored against
        # any subset" %str(notScored)      

        ''' For each alignment subproblem,
        1) make sure its fragments are evenly distributed to fragment chunks.
        2) Setup alignment jobs for its children and enqueue them'''
        alg_problems = [alg for p in self.root_problem.children
                        for alg in p.children]
        for alg_problem in alg_problems:
            assert isinstance(alg_problem, SeppProblem)
            chunks = len(alg_problem.get_children())
            fragment_chunks = alg_problem.fragments.divide_to_equal_chunks(
                chunks)

            ''' Now setup alignment jobs and enqueue them'''
            for (i, fragment_chunk_problem) in enumerate(alg_problem.children):
                fragment_chunk_problem.fragments = fragment_chunks[i]
                aj = fragment_chunk_problem.jobs['hmmalign']
                assert isinstance(aj, HMMAlignJob)
                ''' First Complete setting up alignments'''
                aj.hmmmodel = alg_problem.get_job_result_by_name('hmmbuild')
                aj.base_alignment = alg_problem.jobs["hmmbuild"].infile

                if (
                    fragment_chunk_problem.fragments is None or
                    fragment_chunk_problem.fragments.is_empty()
                   ):
                    aj.fake_run = True
                else:
                    fragment_chunk_problem.fragments.write_to_path(
                        aj.fragments)
                ''' Now the align job can be put on the queue '''
                JobPool().enqueue_job(aj)

    def __str__(self):
        return "join search jobs for all tips of ", self.root_problem

    


class JoinAlignJobs(Join):
    """
    After all alignments jobs for a placement subset have finished,
    we need to build those extended alignments and start placing fragments.
    This join takes care of that step.
    """
    def __init__(self):
        Join.__init__(self)

    def setup_with_placement_problem(self, placement_problem):
        self.placement_problem = placement_problem
        self.root_problem = placement_problem.parent
        for p in placement_problem.iter_leaves():
            self.add_job(p.jobs["hmmalign"])

    def merge_subalignments(self):
        """
        Merge alignment subset extended alignments to get one extended
        alignment for current placement subset.
        """
        pp = self.placement_problem
        _LOG.info("Merging sub-alignments for placement problem : %s." %
                  (pp.label))
        ''' First find fragments assigned to this placement problem'''
        pp.fragments = pp.parent.fragments.get_soft_sub_alignment([])
        for ap in pp.get_children():
            pp.fragments.seq_names |= set(ap.fragments)

        ''' Then, gather a list of all alignments relevant to this placement
        subset'''
        fragfilesperap = dict()
        for ap in pp.children:
            assert isinstance(ap, SeppProblem)
            ''' Get all fragment chunk alignments for this alignment subset'''
            aligned_files = [fp.get_job_result_by_name('hmmalign') for
                             fp in ap.children]
            fragfilesperap[ap] = aligned_files

        ''' Now, build an extended alignment *per each fragment chunk*.
            Simply merge all hmmalign results for fragment chunk numbered i'''
        extendedAlignments = []
        for i in range(0, self.root_problem.fragment_chunks):
            extendedAlignment = ExtendedAlignment(pp.fragments.seq_names)
            for ap in pp.children:
                # _LOG.debug("Merging fragment chunks for subalignment : %s."
                # %(ap.label))
                if fragfilesperap[ap][i]:
                    ap_alg = ap.read_extendend_alignment_and_relabel_columns(
                        ap.jobs["hmmbuild"].infile, [fragfilesperap[ap][i]])
                else:
                    ap_alg = ap.read_extendend_alignment_and_relabel_columns(
                        ap.jobs["hmmbuild"].infile, [])
                _LOG.debug(
                    ("Merging alignment subset into placement subset for "
                     "chunk %d: %s.") % (i, ap.label))
                extendedAlignment.merge_in(ap_alg, convert_to_string=False)
            '''Extended alignmnts have all fragments. remove the ones that
               don't belong to thsi chunk'''
            extendedAlignment.remove_missing_fragments()
            extendedAlignment.from_bytearray_to_string()
            extendedAlignments.append(extendedAlignment)
        return extendedAlignments

    def perform(self):
        pp = self.placement_problem
        fullExtendedAlignments = self.merge_subalignments()

        for i in range(0, self.root_problem.fragment_chunks):
            fullExtendedAlignment = fullExtendedAlignments[i]
            # Split the backbone alignment and query sequences
            # into separate files
            queryExtendedAlignment = \
                fullExtendedAlignment.get_fragments_readonly_alignment()
            baseAlignment = fullExtendedAlignment.get_base_readonly_alignment()
            pj = pp.jobs[get_placement_job_name(i)]
            assert isinstance(pj, PplacerJob)
            if (queryExtendedAlignment.is_empty()):
                pj.fake_run = True

            # Write out the extended alignments, split into query and full-
            # length for pplacer
            queryExtendedAlignment.write_to_path(pj.extended_alignment_file)
            baseAlignment.write_to_path(pj.backbone_alignment_file)

            # But keep the extended alignment on everything
            pj.set_attribute(
                "full_extended_alignment_object", fullExtendedAlignment)

            JobPool().enqueue_job(pj)

    def __str__(self):
        return "join align jobs for tips of ", self.placement_problem


class ExhaustiveAlgorithm(AbstractAlgorithm):
    """
    This implements the exhaustive algorithm where all alignments subsets
    are searched for every fragment.
    """
    def __init__(self):
        AbstractAlgorithm.__init__(self)
        self.place_nomatch_fragments = False
        ''' Hardcoded E-Lim for hmmsearch '''  # TODO: what to do with this
        self.elim = 99999999
        self.filters = False
        self.symfrac = True
        self.strategy = options().exhaustive.strategy
        self.decomp_strategy = options().decomp_strategy
        self.minsubsetsize = int(options().exhaustive.minsubsetsize)
        # Temp fix for now,
        self.molecule = self.options.molecule
        self.distances = dict()

    def compute_distances(self, sequences):
        for seq1, val1 in sequences.items():
            for seq2, val2 in sequences.items():
                if ("".join([seq1, seq2]) not in self.distances):
                    self.distances["".join([seq1, seq2])] = \
                        hamming_distance(val1, val2)
                    self.distances["".join([seq2, seq1])] = \
                        self.distances["".join([seq1, seq2])]

    def merge_results(self):
        assert isinstance(self.root_problem, RootProblem)

        '''Generate single extended alignment'''
        fullExtendedAlignment = ExtendedAlignment(
            self.root_problem.fragments.keys())
        # self.root_problem.get_children()[0].jobs[get_placement_job_name(0)]\
        # .get_attribute("full_extended_alignment_object")
        for pp in self.root_problem.get_children():
            for i in range(0, self.root_problem.fragment_chunks):
                extended_alignment = pp.jobs[
                    get_placement_job_name(i)].get_attribute(
                        "full_extended_alignment_object")
                fullExtendedAlignment.merge_in(
                    extended_alignment, convert_to_string=True)
        self.results = fullExtendedAlignment

        mergeinput = []
        '''Append main tree to merge input'''
        mergeinput.append(
            "%s;" % (self.root_problem.subtree.compose_newick(labels=True)))
        for pp in self.root_problem.get_children():
            assert isinstance(pp, SeppProblem)
            for i in range(0, self.root_problem.fragment_chunks):
                if (pp.get_job_result_by_name(
                        get_placement_job_name(i)) is None):
                    continue
                '''Append subset trees and json locations to merge input'''
                mergeinput.append(
                    "%s;\n%s" % (
                        pp.subtree.compose_newick(labels=True),
                        pp.get_job_result_by_name(get_placement_job_name(i))))
        mergeinput.append("")
        mergeinput.append("")
        meregeinputstring = "\n".join(mergeinput)
        mergeJsonJob = MergeJsonJob()
        mergeJsonJob.setup(meregeinputstring,
                           self.get_output_filename("placement.json"))
        mergeJsonJob.run()

    def output_results(self):
        """ Merged json file is already saved in merge_results function and
            full extended alignment already created in merge_results function
        """
        outfilename = self.get_output_filename("alignment.fasta")
        self.results .write_to_path(outfilename)
        self.results.remove_insertion_columns()
        outfilename = self.get_output_filename("alignment_masked.fasta")
        self.results.write_to_path(outfilename)
        namerev_script = self.root_problem.subtree.rename_script()
        if namerev_script:
            outfilename = self.get_output_filename("rename-json.py")
            with open(outfilename, 'w') as s:
                s.write(namerev_script)

    def check_options(self, supply=[]):
        if (options().info_file is None):
            supply = supply + ["raxml file"]
        AbstractAlgorithm.check_options(self, supply)

    def modify_tree(self, a_tree):
        pass

    def build_subproblems(self):
        (alignment, tree) = self.read_alignment_and_tree()

        if options().distance != 1:
            self.compute_distances(alignment)

        assert isinstance(tree, PhylogeneticTree)
        assert isinstance(alignment, MutableAlignment)

        tree.get_tree().resolve_polytomies()
        # Label edges with numbers so that we could assemble things back
        # at the end
        tree.lable_edges()

        ''' Make sure size values are set, and are meaningful. '''
        self.check_and_set_sizes(alignment.get_num_taxa())

        self._create_root_problem(tree, alignment)

        ''' Decompose the tree based on placement subsets'''
        placement_tree_map = PhylogeneticTree(
            Tree(tree.den_tree)).decompose_tree(
                self.options.placement_size,
                strategy=self.strategy,
                minSize=self.options.placement_size / int(
                    self.options.exhaustive.placementminsubsetsizefacotr),
                tree_map={}, pdistance=1,
                decomp_strategy=self.decomp_strategy,
                distances=self.distances,
                maxDiam=None)
        assert len(placement_tree_map) > 0, (
            "Tree could not be decomposed"
            " given the following settings; strategy:%s minsubsetsize:%s"
            " placement_size:%s"
            % (self.strategy, self.minsubsetsize, self.options.placement_size))
        _LOG.info("Breaking into %d placement subsets." % len(
            placement_tree_map))

        ''' For placement subsets create a placement subproblem,
            and decompose further'''
        for (p_key, p_tree) in placement_tree_map.items():
            assert isinstance(p_tree, PhylogeneticTree)
            placement_problem = SeppProblem(
                p_tree.leaf_node_names(), self.root_problem)
            placement_problem.subtree = p_tree
            placement_problem.label = "P_%s" % str(p_key)
            _LOG.debug(
                "Placement subset %s has %d nodes" %
                (placement_problem.label, len(p_tree.leaf_node_names())))
            ''' Further decompose to alignment subsets '''
            alignment_tree_map = PhylogeneticTree(
                Tree(p_tree.den_tree)).decompose_tree(
                    self.options.alignment_size,
                    strategy=self.strategy,
                    minSize=self.minsubsetsize,
                    tree_map={},
                    decomp_strategy=self.options.decomp_strategy,
                    pdistance=options().distance,
                    distances=self.distances,
                    maxDiam=self.options.maxDiam)
            assert len(alignment_tree_map) > 0, (
                "Tree could not be decomposed"
                " given the following settings; strategy:%s"
                " minsubsetsize:%s alignmet_size:%s" %
                (self.strategy, self.minsubsetsize,
                 self.options.alignment_size))

            _LOG.debug("Placement subset %s has %d alignment subsets: %s" %
                       (placement_problem.label, len(alignment_tree_map),
                        str(sorted(alignment_tree_map.keys()))))
            _LOG.debug("Placement subset %s has %d taxa:" %
                       (placement_problem.label,
                        sum([len(a_tree.leaf_node_names())
                             for a_tree in alignment_tree_map.values()])))
            for (a_key, a_tree) in alignment_tree_map.items():
                assert isinstance(a_tree, PhylogeneticTree)
                self.modify_tree(a_tree)
                alignment_problem = SeppProblem(a_tree.leaf_node_names(),
                                                placement_problem)
                alignment_problem.subtree = a_tree
                alignment_problem.label = "A_%s_%s" % (str(p_key), str(a_key))

        _LOG.info("Breaking into %d alignment subsets." %
                  (len(list(self.root_problem.iter_leaves()))))

        ''' Divide fragments into chunks, to help achieve better parallelism'''
        fragment_chunk_files = self.create_fragment_files()
        self.root_problem.fragment_chunks = len(fragment_chunk_files)
        for alignment_problem in self.root_problem.iter_leaves():
            for afc in range(0, self.root_problem.fragment_chunks):
                frag_chunk_problem = SeppProblem(alignment_problem.taxa,
                                                 alignment_problem)
                frag_chunk_problem.subtree = alignment_problem.subtree
                frag_chunk_problem.label = alignment_problem.label.replace(
                    "A_", "FC_") + "_" + str(afc)
                frag_chunk_problem.fragments = fragment_chunk_files[afc]

        _LOG.info("Breaking each alignment subset into %d fragment chunks." %
                  self.root_problem.fragment_chunks)
        _LOG.debug("Subproblem structure: %s" % str(self.root_problem))
        return self.root_problem

    def create_fragment_files(self):
        alg_subset_count = len(list(self.root_problem.iter_leaves()))
        frag_chunk_count = lcm(
            alg_subset_count, self.options.cpu) // alg_subset_count
        return self.read_and_divide_fragments(frag_chunk_count)

    def _get_new_Join_Align_Job(self):
        return JoinAlignJobs()

    def _log_pipe(self):
        if hasattr(self.options.hmmsearch, "piped"):
            pipe = self.options.hmmsearch.piped.\
                strip().lower() == "true"
        else:
            pipe = True
        results_on_temp = not pipe
        _LOG.debug("HmmSearch: Piped?: %s and keep on temp?: %s" % (
            str(pipe), str(results_on_temp)))

    def build_jobs(self):
        assert isinstance(self.root_problem, SeppProblem)

        self._log_pipe()

        for placement_problem in self.root_problem.get_children():
            ''' Create placer jobs'''
            for i in range(0, self.root_problem.fragment_chunks):
                pj = PplacerJob()
                pj.partial_setup_for_subproblem(
                    placement_problem, self.options.info_file, i)
                placement_problem.add_job(get_placement_job_name(i), pj)

            '''For each alignment subproblem, ...'''
            for alg_problem in placement_problem.children:
                assert isinstance(alg_problem, SeppProblem)
                ''' create the build model job'''
                bj = HMMBuildJob()
                _LOG.debug("Alignment subproblem: %s" % str(alg_problem))
                bj.setup_for_subproblem(
                    alg_problem, symfrac=self.symfrac,
                    molecule=self.molecule,
                    **vars(self.options.hmmbuild))
                alg_problem.add_job(bj.job_type, bj)
                ''' create the search jobs'''
                for fc_problem in alg_problem.get_children():
                    sj = HMMSearchJob()
                    sj.partial_setup_for_subproblem(
                        fc_problem.fragments,
                        fc_problem, self.elim, self.filters)
                    fc_problem.add_job(sj.job_type, sj)
                    ''' create the align job'''
                    aj = HMMAlignJob()
                    fc_problem.add_job(aj.job_type, aj)
                    aj.partial_setup_for_subproblem(
                        fc_problem, molecule=self.molecule)

    def connect_jobs_og(self):
        """ a callback function called after hmmbuild jobs are finished"""
        def enq_job_searchfragment(result, search_job):
            search_job.hmmmodel = result
            JobPool().enqueue_job(search_job)
        assert isinstance(self.root_problem, SeppProblem)
        for placement_problem in self.root_problem.get_children():
            '''For each alignment subproblem, ...'''
            for alg_problem in placement_problem.children:
                assert isinstance(alg_problem, SeppProblem)
                ''' create the build model job'''
                bj = alg_problem.jobs["hmmbuild"]
                ''' create the search jobs'''
                for fc_problem in alg_problem.get_children():
                    sj = fc_problem.jobs["hmmsearch"]
                    ''' connect bulid and search jobs'''
                    bj.add_call_Back(
                        lambda result, next_job=sj: enq_job_searchfragment(
                            result, next_job))
            '''Join all align jobs of a placement subset (enqueues
                placement job)'''
            jaj = self._get_new_Join_Align_Job()
            jaj.setup_with_placement_problem(placement_problem)
        ''' Join all search jobs together (enqueues align jobs)'''
        jsj = JoinSearchJobs()
        jsj.setup_with_root_problem(self.root_problem)
        _LOG.debug("Jobs joined successfully")

    def connect_jobs(self):
        """ a callback function called after hmmbuild jobs are finished"""
        def enq_job_searchfragment(result, search_job):
            search_job.hmmmodel = result
            JobPool().enqueue_job(search_job)
        assert isinstance(self.root_problem, SeppProblem)
        for placement_problem in self.root_problem.get_children():
            '''For each alignment subproblem, ...'''
            for alg_problem in placement_problem.children:
                assert isinstance(alg_problem, SeppProblem)
                ''' create the build model job'''
                bj = alg_problem.jobs["hmmbuild"]
                ''' create the search jobs'''
                for fc_problem in alg_problem.get_children():
                    sj = fc_problem.jobs["hmmsearch"]

                    # print(fc_problem)
                    # print(fc_problem.fragments)
                    # print(" ")

                    ''' connect bulid and search jobs'''
                    bj.add_call_Back(
                        lambda result, next_job=sj: enq_job_searchfragment(
                            result, next_job))
            '''Join all align jobs of a placement subset (enqueues
                placement job)'''
            jaj = self._get_new_Join_Align_Job()
            jaj.setup_with_placement_problem(placement_problem)
        ''' Join all search jobs together (enqueues align jobs)'''
        
        fc_jobs = {}
        i = 0
        #  iterate over all the fragment chunks, not all the fragment chunk problems
        #  the fragments variable is a fasta file, need to open that using alignment to get qseqs
        for fragment_chunk_problem in self.root_problem.get_children()[0].get_children()[0].get_children():
            fc_jobs[i] = QuerySetJob(i)
            i += 1
        # for qseq in self.root_problem.fragments:
        #     qseqjob = QuerySetJob(qseq)
        #     qseq_jobs[qseq] = qseqjob
            
        jsj = JoinSearchJobs()
        jsj.setup_with_root_problem(self.root_problem, fc_jobs)
        
        jfsj = JoinFastSearchJobs()
        jfsj.setup_with_root_problem(self.root_problem, fc_jobs)
        
        _LOG.debug("Jobs joined successfully")

    def enqueue_firstlevel_job(self):
        for p in self.root_problem.children:
            for ap in p.children:
                JobPool().enqueue_job(ap.jobs["hmmbuild"])


if __name__ == '__main__':
    ExhaustiveAlgorithm().run()
