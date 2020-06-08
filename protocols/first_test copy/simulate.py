import opentrons.simulate
protocol_file = open('thermocycler_csv.py')
opentrons.simulate.simulate(protocol_file)
