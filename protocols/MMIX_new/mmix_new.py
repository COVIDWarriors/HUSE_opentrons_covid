import math
from opentrons.types import Point
from opentrons import protocol_api
import time
import os
import numpy as np
from timeit import default_timer as timer
import json
from datetime import datetime
import csv

# metadata
metadata = {
    'protocolName': 'Per Version 2',
    'author': 'Matias Bonet Fullana & Antoni Morla. based on: Malen Aguirregabiria,Aitor Gastaminza & José Luis Villanueva (jlvillanueva@clinic.cat)',
    'source': 'Hospital Son Espases Palma',
    'apiLevel': '2.2',
    'description': 'Protocol for Marter mix'

}
'''
'technician': '$technician',
'date': '$date'
'''
# Defined variables
##################
NUM_SAMPLES = 8
# NUM_SAMPLES = NUM_SAMPLES -1 #Remove last sample (PC), done manually

air_gap_vol = 0
air_gap_mmix = 0
air_gap_sample = 0
run_id = '$run_id'

# Tune variables
size_transfer = 4  # Number of wells the distribute function will fill
volume_sample = 20  # Volume of the sample
# (NUM_SAMPLES * 1.5 * volume_mmix)  # Total volume of first screwcap
extra_dispensal = 0  # Extra volume for master mix in each distribute transfer
diameter_screwcap = 8.25  # Diameter of the screwcap
temperature = 10  # Temperature of temp module
volume_cone = 50  # Volume in ul that fit in the screwcap cone
x_offset = [0, 0]
pipette_allowed_capacity = 18

#############################################################
# Available master mastermixes
#############################################################

MMIX_available = {1: 'SonEspases'}

mmix_selection = 1  # select the mastermix to be used

# volume of mastermixes per sample and number of wells in which is distributed
MMIX_vol = {1: [20, 1]}
MMIX_recipe = {1: [6.25, 1.25, 8.25]}  # Reactive volumes for the mmix

MMIX_make_location = 4  # Cell D1 in which the first tube for the MMIX will be placed

# Volume of transfered master mix per well
volume_mmix = MMIX_vol[mmix_selection][0]

MMIX_make = {}
for mmix_type in MMIX_recipe.keys():
    MMIX_make[mmix_type] = []
    for needed_vol in MMIX_recipe[mmix_type]:
        MMIX_make[mmix_type].append(needed_vol * NUM_SAMPLES * 1.1)

# Total volume of mastermix that will be prepared
volume_mmix_available = (NUM_SAMPLES * 1.1 * MMIX_vol[mmix_selection][0])

#############################################
# Calculated variables
area_section_screwcap = (np.pi * diameter_screwcap**2) / 4
h_cone = (volume_cone * 3 / area_section_screwcap)
num_cols = math.ceil(NUM_SAMPLES / 8)  # Columns we are working on


def run(ctx: protocol_api.ProtocolContext):
    import os
    # from opentrons.drivers.rpi_drivers import gpio
    # gpio.set_rail_lights(False) #Turn off lights (termosensible reagents)
    comment(ctx, 'Actual used columns: ' + str(num_cols))

    # Define the STEPS of the protocol
    STEP = 0
    STEPS = {  # Dictionary with STEP activation, description, and times
        1: {'Execute': False, 'description': 'Make MMIX'},
        2: {'Execute': True, 'description': 'Transfer MMIX'},
        3: {'Execute': True, 'description': 'Transfer elution'}
    }

    for s in STEPS:  # Create an empty wait_time
        if 'wait_time' not in STEPS[s]:
            STEPS[s]['wait_time'] = 0

    # Folder and file_path for log time
    folder_path = '/var/lib/jupyter/notebooks/'+run_id
    if not ctx.is_simulating():
        if not os.path.isdir(folder_path):
            os.mkdir(folder_path)
        file_path = folder_path + '/KC_qPCR_time_log.txt'

    # Reagents and their characteristics
    taq_path = Reagent(name='R1_MM_TaqPath',
                       rinse=False,
                       flow_rate_aspirate=1,
                       flow_rate_dispense=1,
                       reagent_reservoir_volume=1000,
                       num_wells=1,  # change with num samples
                       delay=0,
                       h_cono=h_cone,
                       v_fondo=volume_cone  # V cono
                       )

    covid_assay = Reagent(name='R2_AM_TaqPath_Covid19_assay',
                          rinse=False,
                          flow_rate_aspirate=1,
                          flow_rate_dispense=1,
                          reagent_reservoir_volume=150,
                          num_wells=1,  # change with num samples
                          delay=0,
                          h_cono=h_cone,
                          v_fondo=volume_cone  # V cono
                          )

    mmix_water = Reagent(name='R3_Water',
                         rinse=False,
                         flow_rate_aspirate=1,
                         flow_rate_dispense=1,
                         reagent_reservoir_volume=1000,
                         num_wells=1,  # change with num samples
                         delay=0,
                         h_cono=h_cone,
                         v_fondo=volume_cone  # V cono
                         )

    MMIX = Reagent(name='Master Mix',
                   rinse=False,
                   flow_rate_aspirate=1,
                   flow_rate_dispense=1,
                   reagent_reservoir_volume=volume_mmix_available,
                   num_wells=1,  # change with num samples
                   delay=0,
                   h_cono=h_cone,
                   v_fondo=volume_cone  # V cono
                   )

    Samples = Reagent(name='Samples',
                      rinse=False,
                      flow_rate_aspirate=1,
                      flow_rate_dispense=1,
                      reagent_reservoir_volume=0,
                      delay=0,
                      num_wells=num_cols,  # num_cols comes from available columns
                      h_cono=0,
                      v_fondo=0
                      )
    Elution = Reagent(name='Elution',
                      rinse=False,
                      flow_rate_aspirate=1,
                      flow_rate_dispense=1,
                      reagent_reservoir_volume=50,
                      delay=0,
                      num_wells=num_cols,  # num_cols comes from available columns
                      h_cono=0,
                      v_fondo=0
                      )

    # Assign class type reactives to a summary
    MMIX_components = [taq_path, covid_assay, mmix_water]

    ####################################
    # load labware and modules
    ############################################
    # tempdeck
    tempdeck = ctx.load_module('tempdeck', '4')
    tempdeck.set_temperature(temperature)

    ##################################
    # 24 well rack
    tuberack = tempdeck.load_labware(
        'opentrons_24_aluminumblock_generic_2ml_screwcap')

    ##################################
    # Sample plate - comes from King fisher/ Manual / Other
    source_plate = ctx.load_labware(
        'biorad_96_wellplate_200ul_pcr', '2')

    # Elution
    # Final result
    destination_plate = ctx.load_labware(
        'biorad_96_wellplate_200ul_pcr', '3')

    ##################################
    # Load Tipracks
    tips20 = [
        ctx.load_labware('opentrons_96_filtertiprack_20ul', slot)
        for slot in ['5', '7']
    ]

    ################################################################################
    # Declare which reagents are in each reservoir as well as deepwell and elution plate
    # 1 row, 2 columns (first ones)
    MMIX.reagent_reservoir = tuberack.rows()[0][:MMIX.num_wells]
    MMIX_components_location = tuberack.wells()[MMIX_make_location:(
        MMIX_make_location + len(MMIX_make[mmix_selection]))]

    comment(ctx, 'Wells in: ' + str(tuberack.rows()
                                    [0][:MMIX.num_wells]) + ' element: '+str(MMIX.reagent_reservoir[MMIX.col]))

    # setup up sample sources and destinations
    samples = source_plate.wells()[:NUM_SAMPLES]
    samples_multi = source_plate.rows()[0][:num_cols]
    pcr_wells = destination_plate.wells()[:NUM_SAMPLES]
    pcr_wells_multi = destination_plate.rows()[0][:num_cols]

    # pipettes
    m20 = ctx.load_instrument(
        'p20_multi_gen2', mount='left', tip_racks=tips20)

    p20 = ctx.load_instrument(
        'p20_single_gen2', mount='right', tip_racks=tips20)

    # used tips counter
    tip_track = {
        'counts': {p20: 0, m20: 0},
        'maxes': {p20: len(tips20) * 96, m20: len(tips20)*96}
    }

    ##########
    ############################################################################
    # STEP 1: Make Master MIX
    ############################################################################

    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        # Check if among the pipettes, p20_single is installed
        used_vol = []
        comment(ctx, 'Selected MMIX: ' +
                MMIX_available[mmix_selection], add_hash=True)

        for i, [source, vol] in enumerate(zip(MMIX_components_location, MMIX_make[mmix_selection])):

            pick_up(p20, tip_track)
            comment(ctx, 'Add component: ' +
                    MMIX_components[i].name, add_hash=True)

            # because 20ul is the maximum volume of the tip we will choose 17
            if (vol + air_gap_vol) > pipette_allowed_capacity:
                # calculate what volume should be transferred in each step
                vol_list = divide_volume(vol, pipette_allowed_capacity)
                comment(ctx, vol_list)
                for vol in vol_list:
                    move_vol_multichannel(ctx, p20, reagent=MMIX_components[i], source=source, dest=MMIX.reagent_reservoir[0],
                                          # should be changed with picku_up_height calculation
                                          vol=vol, air_gap_vol=air_gap_vol, x_offset=x_offset, pickup_height=1,
                                          rinse=False, disp_height=-10, blow_out=True, touch_tip=True)
            else:
                move_vol_multichannel(ctx, p20, reagent=MMIX_components[i], source=source, dest=MMIX.reagent_reservoir[0],
                                      # should be changed with picku_up_height calculation
                                      vol=vol, air_gap_vol=air_gap_vol, x_offset=x_offset, pickup_height=1,
                                      rinse=False, disp_height=-10, blow_out=True, touch_tip=True)

            if i+1 < len(MMIX_components):
                p20.drop_tip()
            else:
                comment(ctx, 'Final mix', add_hash=True)

                custom_mix(p20, reagent=MMIX, location=MMIX.reagent_reservoir[0], vol=18, rounds=5,
                           blow_out=True, mix_height=2, x_offset=x_offset)

                p20.drop_tip()

            tip_track['counts'][p20] += 1

        end = datetime.now()
        time_taken = (end - start)
        comment(ctx, 'Step ' + str(STEP) + ': ' +
                STEPS[STEP]['description'] + ' took ' + str(time_taken), add_hash=True)
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 2: Transfer Master MIX
    ############################################################################
    ctx._hw_manager.hardware.set_lights(
        rails=True)  # set lights off when using MMIX
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        p20.pick_up_tip()

        for dest in pcr_wells:
            [pickup_height, col_change] = calc_height(ctx,
                                                      MMIX, area_section_screwcap, volume_mmix)

            move_vol_multichannel(ctx, p20, reagent=MMIX, source=MMIX.reagent_reservoir[0],
                                  dest=dest, vol=volume_mmix, air_gap_vol=air_gap_mmix, x_offset=x_offset,
                                  pickup_height=pickup_height, disp_height=-10, rinse=False,
                                  blow_out=True, touch_tip=True)

            # used_vol.append(used_vol_temp)

        p20.drop_tip()
        tip_track['counts'][p20] += 1
        #MMIX.unused_two = MMIX.vol_well

        end = datetime.now()
        time_taken = (end - start)
        comment(ctx, 'Step ' + str(STEP) + ': ' +
                STEPS[STEP]['description'] + ' took ' + str(time_taken))
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 3: TRANSFER Samples
    ############################################################################
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        comment(ctx, 'pcr_wells')
        # Loop over defined wells
        for s, d in zip(samples_multi, pcr_wells_multi):
            comment(ctx, "%s %s" % (s, d))
            m20.pick_up_tip()
            # Source samples
            move_vol_multichannel(ctx, m20, reagent=Elution, source=s, dest=d,
                                  vol=volume_sample, air_gap_vol=air_gap_sample, x_offset=x_offset,
                                  pickup_height=0.5, disp_height=-10, rinse=False,
                                  blow_out=True, touch_tip=False, post_airgap=True)

            # ADD Custom mix
            m20.drop_tip()
            tip_track['counts'][m20] += 8

        end = datetime.now()
        time_taken = (end - start)

        comment(ctx, 'Step ' + str(STEP) + ': ' +
                STEPS[STEP]['description'] + ' took ' + str(time_taken))

        STEPS[STEP]['Time:'] = str(time_taken)

    # Export the time log to a tsv file
    if not ctx.is_simulating():
        with open(file_path, 'w') as f:
            f.write('STEP\texecution\tdescription\twait_time\texecution_time\n')
            for key in STEPS.keys():
                row = str(key)
                for key2 in STEPS[key].keys():
                    row += '\t' + format(STEPS[key][key2])
                f.write(row + '\n')
        f.close()

    ############################################################################
    # Light flash end of program

    time.sleep(2)
    import os
    #os.system('mpg123 -f -8000 /etc/audio/speaker-test.mp3 &')

    '''if STEPS[1]['Execute'] == True:
        total_used_vol = np.sum(used_vol)
        total_needed_volume = total_used_vol
        comment(ctx,'Total Master Mix used volume is: ' + str(total_used_vol) + '\u03BCl.')
        comment(ctx,'Needed Master Mix volume is ' +
                    str(total_needed_volume + extra_dispensal*len(dests)) +'\u03BCl')
        comment(ctx,'Used Master Mix volumes per run are: ' + str(used_vol) + '\u03BCl.')
        comment(ctx,'Master Mix Volume remaining in tubes is: ' +
                    format(np.sum(MMIX.unused)+extra_dispensal*len(dests)+MMIX.vol_well) + '\u03BCl.')
        comment(ctx,'200 ul Used tips in total: ' + str(tip_track['counts'][p20]))
        comment(ctx,'200 ul Used racks in total: ' + str(tip_track['counts'][p20] / 96))'''

    if STEPS[2]['Execute'] == True:
        comment(ctx, '20 ul Used tips in total: ' +
                str(tip_track['counts'][m20]))
        comment(ctx, '20 ul Used racks in total: ' +
                str(tip_track['counts'][m20] / 96))

    blink(ctx)
    comment(ctx, 'Finished! \nMove plate to PCR')


##################
# Custom functions
##################
# Define Reagents as objects with their properties
class Reagent:
    def __init__(self, name, flow_rate_aspirate, flow_rate_dispense, rinse,
                 reagent_reservoir_volume, delay, num_wells, h_cono, v_fondo,
                 tip_recycling='none'):
        self.name = name
        self.flow_rate_aspirate = flow_rate_aspirate
        self.flow_rate_dispense = flow_rate_dispense
        self.rinse = bool(rinse)
        self.reagent_reservoir_volume = reagent_reservoir_volume
        self.delay = delay
        self.num_wells = num_wells
        self.col = 0
        self.h_cono = h_cono
        self.v_cono = v_fondo
        self.unused = []
        self.tip_recycling = tip_recycling
        self.vol_well_original = reagent_reservoir_volume / num_wells
        self.vol_well = self.vol_well_original


def divide_volume(volume, max_vol):
    num_transfers = math.ceil(volume/max_vol)
    vol_roundup = math.ceil(volume/num_transfers)
    last_vol = volume - vol_roundup*(num_transfers-1)
    vol_list = [vol_roundup for v in range(1, num_transfers)]
    vol_list.append(last_vol)
    return vol_list


def divide_destinations(l, n):
    a = []
    # Divide the list of destinations in size n lists.
    for i in range(0, len(l), n):
        a.append(l[i:i + n])

    return a


def distribute_custom(ctx, pipette, reagent, volume, src, dest, waste_pool, pickup_height, extra_dispensal, disp_height=0):
    # Custom distribute function that allows for blow_out in different location and adjustement of touch_tip
    air_gap = 10
    pipette.aspirate((len(dest) * volume) +
                     extra_dispensal, src.bottom(pickup_height), rate=reagent.flow_rate_aspirate)
    pipette.touch_tip(speed=20, v_offset=-5)
    pipette.move_to(src.top(z=5))
    pipette.aspirate(air_gap, rate=reagent.flow_rate_aspirate)  # air gap

    for d in dest:
        pipette.dispense(air_gap, d.top(), rate=reagent.flow_rate_dispense)
        drop = d.top(z=disp_height)
        pipette.dispense(volume, drop, rate=reagent.flow_rate_dispense)
        # pause for x seconds depending on reagent
        ctx.delay(seconds=reagent.delay)
        pipette.move_to(d.top(z=5))
        pipette.aspirate(air_gap, d.top(
            z=5), rate=reagent.flow_rate_aspirate)  # air gap

    try:
        pipette.blow_out(waste_pool.wells()[0].bottom(pickup_height + 3))
    except:
        pipette.blow_out(waste_pool.bottom(pickup_height + 3))
    return (len(dest) * volume)

##########
# pick up tip and if there is none left, prompt user for a new rack


def pick_up(pip, tip_track):
    if not ctx.is_simulating():
        if tip_track['counts'][pip] == tip_track['maxes'][pip]:
            ctx.pause('Replace ' + str(pip.max_volume) + 'µl tipracks before \
            resuming.')
            pip.reset_tipracks()
            tip_track['counts'][pip] = 0

    if not pip.hw_pipette['has_tip']:
        pip.pick_up_tip()


def move_vol_multichannel(ctx, pipette, reagent, source, dest, vol, air_gap_vol, x_offset,
                          pickup_height, rinse, disp_height, blow_out, touch_tip,
                          post_dispense=False, post_dispense_vol=20,
                          post_airgap=False, post_airgap_vol=10):
    '''
    x_offset: list with two values. x_offset in source and x_offset in destination i.e. [-1,1]
    pickup_height: height from bottom where volume
    rinse: if True it will do 2 rounds of aspirate and dispense before the tranfer
    disp_height: dispense height; by default it's close to the top (z=-2), but in case it is needed it can be lowered
    blow_out, touch_tip: if True they will be done after dispensing
    '''
    # Rinse before aspirating
    if rinse == True:
        custom_mix(pipette, reagent, location=source, vol=vol,
                   rounds=2, blow_out=True, mix_height=0,
                   x_offset=x_offset)
    # SOURCE
    s = source.bottom(pickup_height).move(Point(x=x_offset[0]))
    # aspirate liquid
    pipette.aspirate(vol, s, rate=reagent.flow_rate_aspirate)
    if air_gap_vol != 0:  # If there is air_gap_vol, switch pipette to slow speed
        pipette.aspirate(air_gap_vol, source.top(z=-2),
                         rate=reagent.flow_rate_aspirate)  # air gap
    # GO TO DESTINATION
    drop = dest.top(z=disp_height).move(Point(x=x_offset[1]))
    pipette.dispense(vol + air_gap_vol, drop,
                     rate=reagent.flow_rate_dispense)  # dispense all
    # pause for x seconds depending on reagent
    ctx.delay(seconds=reagent.delay)
    if blow_out == True:
        pipette.blow_out(dest.top(z=-2))
    if post_dispense == True:
        pipette.dispense(post_dispense_vol, dest.top(z=-2))
    if touch_tip == True:
        pipette.touch_tip(speed=20, v_offset=-5, radius=0.9)
    if post_airgap == True:
        pipette.dispense(post_airgap_vol, dest.top(z=5), rate=2)


def custom_mix(pipet, reagent, location, vol, rounds, blow_out, mix_height,
               x_offset, source_height=3, post_airgap=False, post_airgap_vol=10,
               post_dispense=False, post_dispense_vol=20,):
    '''
    Function for mixing a given [vol] in the same [location] a x number of [rounds].
    blow_out: Blow out optional [True,False]
    x_offset = [source, destination]
    source_height: height from bottom to aspirate
    mix_height: height from bottom to dispense
    '''
    if mix_height == 0:
        mix_height = 3
    pipet.aspirate(1, location=location.bottom(
        z=source_height).move(Point(x=x_offset[0])), rate=reagent.flow_rate_aspirate)
    for _ in range(rounds):
        pipet.aspirate(vol, location=location.bottom(
            z=source_height).move(Point(x=x_offset[0])), rate=reagent.flow_rate_aspirate)
        pipet.dispense(vol, location=location.bottom(
            z=mix_height).move(Point(x=x_offset[1])), rate=reagent.flow_rate_dispense)
    pipet.dispense(1, location=location.bottom(
        z=mix_height).move(Point(x=x_offset[1])), rate=reagent.flow_rate_dispense)
    if blow_out == True:
        pipet.blow_out(location.top(z=-2))  # Blow out
    if post_dispense == True:
        pipet.dispense(post_dispense_vol, location.top(z=-2))
    if post_airgap == True:
        pipet.dispense(post_airgap_vol, location.top(z=5))


def calc_height(ctx, reagent, cross_section_area, aspirate_volume, min_height=0.5, extra_volume=30):

    comment(ctx, 'Remaining volume ' + str(reagent.vol_well) +
            '< needed volume ' + str(aspirate_volume) + '?')
    if reagent.vol_well < aspirate_volume + extra_volume:
        reagent.unused.append(reagent.vol_well)
        comment(ctx, 'Next column should be picked')
        comment(ctx, 'Previous to change: ' + str(reagent.col))
        # column selector position; intialize to required number
        reagent.col = reagent.col + 1
        comment(ctx, str('After change: ' + str(reagent.col)))
        reagent.vol_well = reagent.vol_well_original
        comment(ctx, 'New volume:' + str(reagent.vol_well))
        height = (reagent.vol_well - aspirate_volume -
                  reagent.v_cono) / cross_section_area
        # - reagent.h_cono
        reagent.vol_well = reagent.vol_well - aspirate_volume
        comment(ctx, 'Remaining volume:' + str(reagent.vol_well))
        if height < min_height:
            height = min_height
        col_change = True
    else:
        height = (reagent.vol_well - aspirate_volume -
                  reagent.v_cono) / cross_section_area  # - reagent.h_cono
        reagent.vol_well = reagent.vol_well - aspirate_volume
        comment(ctx, 'Calculated height is ' + str(height))
        if height < min_height:
            height = min_height
        comment(ctx, 'Used height is ' + str(height))
        col_change = False
    return height, col_change


def comment(ctx, comment, add_hash=False):
    hash_string = '#######################################################'
    if not ctx.is_simulating():
        if (add_hash):
            ctx.comment(hash_string)
        ctx.comment(comment)
        if (add_hash):
            ctx.comment(hash_string)
    else:
        if (add_hash):
            print(hash_string)
        print(comment)
        if (add_hash):
            print(hash_string)


def blink(ctx):
    for i in range(3):
        ctx._hw_manager.hardware.set_lights(rails=False)
        # ctx._hw_manager.hardware.set_button_light(1,0,0)
        time.sleep(0.3)
        ctx._hw_manager.hardware.set_lights(rails=True)
        # ctx._hw_manager.hardware.set_button_light(0,0,1)
        time.sleep(0.3)
        ctx._hw_manager.hardware.set_lights(rails=False)
