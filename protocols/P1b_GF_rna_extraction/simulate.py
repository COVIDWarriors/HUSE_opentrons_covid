import opentrons.simulate
protocol_file = open('p1b_GF_rna_extraction.py')
opentrons.simulate.simulate(protocol_file)
