import opentrons.simulate
protocol_file = open('p2b_mmix.py')
opentrons.simulate.simulate(protocol_file)
