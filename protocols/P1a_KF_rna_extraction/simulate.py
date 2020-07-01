import opentrons.simulate
protocol_file = open('p1a_KF_prekingfisher.py')
opentrons.simulate.simulate(protocol_file)