import opentrons.simulate
protocol_file = open('p1_KF_prekingfisher.py')
opentrons.simulate.simulate(protocol_file)
