# Independent functional-graph terminal oracle

All 3413 functional graphs with one through five states are checked against
an independent per-start trajectory oracle. Fixed endpoint identities and
cycle labels agree for every start. Invalid successor indices, booleans,
nonintegral entries and empty graphs now raise ValueError rather than being
silently indexed or mislabelled. NumPy integer-array inputs remain supported.
Existing 58 rule-counterfactual basin results are unchanged. This tests an
aggregation algorithm, not empirical biological validity or independent
network-model validation. Full suite: 72 passed, one optional mpbn-dependent
test skipped, one existing empty-slice warning.
