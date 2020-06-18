import opentrons.simulate
protocol_file = open('P1_GF_rna_extraction.py')
opentrons.simulate.simulate(protocol_file)
