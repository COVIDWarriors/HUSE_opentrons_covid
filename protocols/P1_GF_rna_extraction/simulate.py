import opentrons.simulate
protocol_file = open('rna_extraction.py')
opentrons.simulate.simulate(protocol_file)
