import opentrons.simulate
protocol_file = open('mmix_new.py')
opentrons.simulate.simulate(protocol_file)
