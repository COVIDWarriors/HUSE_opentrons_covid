import opentrons.simulate
protocol_file = open('p1a_GF_prekingfisher.py')
opentrons.simulate.simulate(protocol_file)
