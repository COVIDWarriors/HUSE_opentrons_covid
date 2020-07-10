import opentrons.simulate
protocol_file = open('p1b_KF_multi_prekingfisher_biel.py')
opentrons.simulate.simulate(protocol_file)
