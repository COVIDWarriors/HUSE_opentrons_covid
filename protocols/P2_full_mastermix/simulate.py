import opentrons.simulate
protocol_file = open('p2_mmix.py')
opentrons.simulate.simulate(protocol_file)
