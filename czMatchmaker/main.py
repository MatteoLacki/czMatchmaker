import pandas as pd
from dplython import DplyFrame, X, sift, select, mutate, rename, summarize, group_by, DelayFunction
from fragmentator import collect_fragments
from pprint import pprint
from collections import Counter
from misc import crossprod, Round

D = DplyFrame(pd.read_csv('data/synapt_substanceP_wh_wv_parsed.csv'))
D = D >> sift(X.estimates > 10) >> mutate( estimates = Round(X.estimates) )
precursor_aa_seq = 'RPKPQQFFGLM'
Q = 3

data = D >> sift( X.wH == 80, X.wV == 300 )

def analyze(data):
	precursors = data >> \
		sift( X.tag == 'precursor' ) >> \
		select( X.active, X.neutral, X.estimates, X.wH, X.wV, X.fixed )

	fragments = data >> sift( X.tag != 'precursor' ) >> \
		group_by( X.wH, X.wV, X.fixed, X.tag, X.fragType, X.active, X.aa_break_from_N_term ) >> \
		summarize( estimates = X.estimates.sum() )

	I_on_fragments 	= {}
	optiminfos 		= {}

	for break_point, data in fragments.groupby('aa_break_from_N_term'):
		pairing, optiminfo 			= collect_fragments(data, Q)
		I_on_fragments[break_point] = pairing
		optiminfos[break_point] 	= optiminfo

	cations_fragmented_I = sum(
		sum( I_on_fragments[bP][p] for p in I_on_fragments[bP])
		for bP in I_on_fragments )

	try:
		I_no_reactions = precursors >> \
			sift( X.active==Q, X.neutral == 0) >> \
			select( X.estimates )
	except:
		pprint(precursors)

	I_no_reactions = I_no_reactions.values.flatten()[0]

		# This looks very very fishy.... 
		# ETnoD should be equal to this result.
	prec_ETnoD_PTR_I = precursors >> \
		sift( X.active != Q ) >> \
		rename( PTR = X.neutral, I = X.estimates ) >> \
		mutate( ETnoD = Q - X.PTR - X.active ) >> \
		select( X.ETnoD, X.PTR, X.I )

	I_prec_no_frag = prec_ETnoD_PTR_I >> \
		summarize( I = X.I.sum() )

	I_prec_no_frag = I_prec_no_frag.values.flatten()[0]

	precursorNoReactions = precursors >> \
		sift( X.active == Q ) >> \
		select( X.estimates )

	prec_ETnoD_PTR_I = prec_ETnoD_PTR_I >> mutate(
			I_PTR 	= crossprod(X.PTR, X.I), \
			I_ETnoD = crossprod(X.ETnoD, X.I) ) >> \
		summarize( I_PTR = X.I_PTR.sum(), I_ETnoD = X.I_ETnoD.sum() )

	I_PTR_no_frag, I_ETnoD_no_frag = prec_ETnoD_PTR_I.values.flatten()

	prob_PTR   = I_PTR_no_frag/(I_PTR_no_frag + I_ETnoD_no_frag)
	prob_ETnoD = 1. - prob_PTR

	I_frags = dict(
		( bP, sum( I_on_fragments[bP][pairing] for pairing in I_on_fragments[bP] ) )
		for bP in I_on_fragments 	)

	I_frag_total = sum( I_frags[bP] for bP in I_frags )

	prob_frag = Counter(dict( (int(bP), I_frags[bP]/I_frag_total) for bP in I_frags))
	prob_frag = [ prob_frag[i] for i in range(len(precursor_aa_seq)) ]

	I_frags_PTRETnoD_total = sum(
		( Q-1-sum( q for cz, q in pairing ) ) * I_on_fragments[bP][pairing] for bP in I_on_fragments for pairing in I_on_fragments[bP] 	)

	anion_meets_cation = I_frags_PTRETnoD_total + I_PTR_no_frag + I_ETnoD_no_frag
	prob_fragmentation = I_frags_PTRETnoD_total / anion_meets_cation
	prob_no_fragmentation = 1 - prob_fragmentation

	prob_no_reaction = I_no_reactions/( I_no_reactions + I_frag_total + I_prec_no_frag )
	prob_reaction = 1. - prob_no_reaction

	res = [ prob_reaction, prob_no_reaction, prob_fragmentation, prob_no_fragmentation, prob_PTR, prob_ETnoD ]

	res.extend( prob_frag )
	return res

results = []
for (wH, wV), data in D.groupby(['wH', 'wV']):
	res = [wH, wV]
	res.extend( analyze(data) )
	results.append(res)

results = pd.DataFrame(results)
varNames = ['wH','wV','reaction','no_reaction','frag','no_frag','PTR','ETnoD']
varNames.extend( [ 'AA_'+str(i) for i in range(len(precursor_aa_seq)) ] )
results.columns = varNames
results.to_csv('results/analysis_results.csv', na_rep='NA', index=False)
