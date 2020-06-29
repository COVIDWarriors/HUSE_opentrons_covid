import opentrons.simulate
protocol_file = open('p2a_mmix.py')
opentrons.simulate.simulate(protocol_file)
