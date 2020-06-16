from opentrons import protocol_api

metadata = {'apiLevel': '2.2'}

def run(protocol: protocol_api.ProtocolContext):
    plate = protocol.load_labware('corning_96_wellplate_360ul_flat', 1)
    tiprack_1 = protocol.load_labware('opentrons_96_tiprack_300ul', 2)
    tiprack_2 = protocol.load_labware('opentrons_96_tiprack_300ul', 3)
    reservoir = protocol.load_labware('usascientific_12_reservoir_22ml', 4)
    p20 = protocol.load_instrument('p20_single_gen2', 'right', tip_racks=[tiprack_1, tiprack_2])
    p20.distribute(50, reservoir['A12'], plate.wells())  # dilutent

    # loop through each row
    for i in range(8):

        # save the source well and destination column to variables
        source = reservoir.wells()[i]
        row = plate.rows()[i]

        # transfer 30uL of source to first well in column
        p20.transfer(10, source, row[0], mix_after=(3, 5))

        # dilute the sample down the column
        p20.transfer(
            10, row[:11], row[1:],
            mix_after=(3, 5))

    