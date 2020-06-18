import opentrons.simulate
protocol_file = open('p1kf_prekingfisher.py')
opentrons.simulate.simulate(protocol_file)
