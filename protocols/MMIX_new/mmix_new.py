import math
from opentrons.types import Point
from opentrons import protocol_api
from opentrons import labware
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
NUM_SAMPLES = 96

# NUM_SAMPLES = NUM_SAMPLES -1 #Remove last sample (PC), done manually

air_gap_vol = 10
air_gap_mmix = 0
air_gap_sample = 0
run_id = '$run_id'

# Tune variables
size_transfer = 4  # Number of wells the distribute function will fill
volume_sample = 10  # Volume of the sample
# (NUM_SAMPLES * 1.5 * volume_mmix)  # Total volume of first screwcap
extra_dispensal = 0  # Extra volume for master mix in each distribute transfer
diameter_screwcap = 8.1  # Diameter of the screwcap
temperature = 10  # Temperature of temp module
volume_cone = 50  # Volume in ul that fit in the screwcap cone
x_offset = [0, 0]

#############################################
# Calculated variables
area_section_screwcap = (np.pi * diameter_screwcap**2) / 4
h_cone = (volume_cone * 3 / area_section_screwcap)
num_cols = math.ceil(NUM_SAMPLES / 8)  # Columns we are working on


#############################################################
# Available master mastermixes
#############################################################
MMIX_available = {1: 'SonEspases'}
mmix_selection = 1  # select the mastermix to be used
# volume of mastermixes per sample and number of wells in which is distributed
MMIX_recipe = {1: [8.25, 6.25, 1.25]}  # Reactive volumes for the mmix
MMIX_vol = {1: [20, 1]}
MMIX_make_location = 4  # Cell A4 in which the first tube for the MMIX will be placed

# Volume of transfered master mix per well
volume_mmix = MMIX_vol[mmix_selection][0]

MMIX_make = {}
for mmix_type in MMIX_recipe.keys():
    MMIX_make[mmix_type] = []
    for needed_vol in MMIX_recipe[mmix_type]:
        MMIX_make[mmix_type].append(needed_vol * NUM_SAMPLES * 1.1)

# Total volume of mastermix that will be prepared
volume_mmix_available = (NUM_SAMPLES * 1.1 * MMIX_vol[mmix_selection][0])


def run(ctx: protocol_api.ProtocolContext):

    # Init protocol run
    run = ProtocolRun(ctx)

    run.addStep(execute=True, description="Make MMIX")
    run.addStep(execute=True, description="Transfer MMIX")
    run.addStep(execute=True, description="Make MMIX")

    ##################################
    # Define desk
    tempdeck = ctx.load_module('tempdeck', '4')
    tuberack = tempdeck.load_labware(
        'opentrons_24_aluminumblock_generic_2ml_screwcap')

    # tempdeck.set_temperature(temperature)

    # PCR
    pcr_plate = ctx.load_labware(
        'opentrons_96_aluminumblock_generic_pcr_strip_200ul', '1')

    # Eluted from King fisher/ Manual / Other
    elution_plate = ctx.load_labware(
        'biorad_96_wellplate_200ul_pcr', '2')

    # Tipracks20_multi
    tips20 = ctx.load_labware('opentrons_96_tiprack_20ul', 3)
    tips300 = ctx.load_labware('opentrons_96_filtertiprack_200ul', 7)

    # Mount pippets and set racks
    run.mount_right_pip('p20_single_gen2', tip_racks=[tips20], capacity=20)
    run.mount_left_pip('p300_single_gen2', tip_racks=[tips300], capacity=300)

    # Define wells interaction
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

    pcr_well = Reagent(name='Samples',
                       rinse=False,
                       flow_rate_aspirate=1,
                       flow_rate_dispense=1,
                       reagent_reservoir_volume=0,
                       delay=0,
                       num_wells=num_cols,  # num_cols comes from available columns
                       h_cono=0,
                       v_fondo=0
                       )

    elution_well = Reagent(name='Elution',
                           rinse=False,
                           flow_rate_aspirate=1,
                           flow_rate_dispense=1,
                           reagent_reservoir_volume=50,
                           delay=0,
                           num_wells=num_cols,  # num_cols comes from available columns
                           h_cono=0,
                           v_fondo=0
                           )

    MMIX_components = [mmix_water, taq_path, covid_assay]

    ################################################################################
    # Declare which reagents are in each reservoir as well as deepwell and elution plate
    # 1 row, 2 columns (first ones)
    MMIX.reagent_reservoir = tuberack.rows()[0][:MMIX.num_wells]
    MMIX_components_location = tuberack.wells()[MMIX_make_location:(
        MMIX_make_location + len(MMIX_make[mmix_selection]))]

    # setup up sample sources and destinations
    pcr_wells = pcr_plate.wells()[:NUM_SAMPLES]
    elution_wells = elution_plate.wells()[:NUM_SAMPLES]

    ############################################################################
    # STEP 1: Make Master MIX
    ############################################################################
    if (run.next_step()):
        # Check if among the pipettes, p20_single is installed
        used_vol = []
        run.comment('Selected MMIX: ' +
                    MMIX_available[mmix_selection], add_hash=True)
        run.set_pip("left")
        run.pick_up()
        drop = False
        for i, [source, vol] in enumerate(zip(MMIX_components_location, MMIX_make[mmix_selection])):

            run.comment('Add component: ' +
                        MMIX_components[i].name, add_hash=True)

            # because 20ul is the maximum volume of the tip we will choose 17
            if (vol + air_gap_vol) > run.get_pip_capacity():
                # calculate what volume should be transferred in each step
                vol_list = run.divide_volume(vol, run.get_pip_capacity())
                for vol in vol_list:
                    # If not in first step we need to change everytime
                    if(i > 0):
                        run.pick_up()

                    run.move_vol_multichannel(reagent=MMIX_components[i], source=source, dest=MMIX.reagent_reservoir[0],
                                              # should be changed with picku_up_height calculation
                                              vol=vol, air_gap_vol=air_gap_vol, x_offset=x_offset, pickup_height=1,
                                              rinse=False, disp_height=-10, blow_out=True, touch_tip=False)

                    # If not in first step we need to change everytime
                    if(i > 0):
                        run.drop_tip()
                        drop = True

            else:
                if(i > 0):
                    run.pick_up()
                run.move_vol_multichannel(reagent=MMIX_components[i], source=source, dest=MMIX.reagent_reservoir[0],
                                          # should be changed with picku_up_height calculation
                                          vol=vol, air_gap_vol=air_gap_vol, x_offset=x_offset, pickup_height=1,
                                          rinse=False, disp_height=-10, blow_out=True, touch_tip=False)
                if(i > 0):
                    run.drop_tip()
                    drop = True

            if i+1 < len(MMIX_components):
                if(not drop):
                    run.drop_tip()

            else:
                run.pick_up()
                run.comment('Final mix', add_hash=True)

                run.custom_mix(reagent=MMIX, location=MMIX.reagent_reservoir[0], vol=50, rounds=5,
                               blow_out=True, mix_height=2, x_offset=x_offset)
                run.drop_tip()

        run.finish_step()

    ############################################################################
    # STEP 2: Transfer Master MIX
    ############################################################################
    run.start_lights()
    if (run.next_step()):
        run.set_pip("rigth")
        run.pick_up_tip()
        for dest in pcr_wells:
            [pickup_height, col_change] = run.calc_height(
                MMIX, area_section_screwcap, volume_mmix)

            run.move_vol_multichannel(reagent=MMIX, source=MMIX.reagent_reservoir[0],
                                      dest=dest, vol=volume_mmix, air_gap_vol=air_gap_mmix, x_offset=x_offset,
                                      pickup_height=pickup_height, disp_height=-3, rinse=False,
                                      blow_out=True, touch_tip=True)

            # used_vol.append(used_vol_temp)

        run.drop_tip()
        run.finish_step()

    ############################################################################
    # STEP 3: TRANSFER Samples
    ############################################################################
    if(run.next_steps()):
        run.comment('pcr_wells')
        run.set_pip("right")
        # Loop over defined wells
        for s, d in zip(elution_wells, pcr_wells):
            run.comment("%s %s" % (s, d))
            run.pick_up_tip()
            # Source samples
            run.move_vol_multichannel(reagent=elution_well, source=s, dest=d,
                                      vol=volume_sample, air_gap_vol=air_gap_sample, x_offset=x_offset,
                                      pickup_height=3, disp_height=-0, rinse=False,
                                      blow_out=True, touch_tip=True, post_airgap=True,)

            # ADD Custom mix
            run.drop_tip()

        run.finish_step()

    if STEPS[2]['Execute'] == True:
        comment(ctx, '20 ul Used tips in total: ' +
                str(tip_track['counts'][p20]))
        comment(ctx, '20 ul Used racks in total: ' +
                str(tip_track['counts'][p20] / 96))

    ############################################################################
    # Light flash end of program
    run.log_steps_time()
    run.blink()
    run.comment('Finished! \nMove plate to PCR')

##################
# Custom functionss
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


class ProtocolRun:
    def __init__(self, ctx):
        self.ctx = ctx
        self.step_list = []
        self.step = 0

        # Folder and file_path for log time
        folder_path = '/var/lib/jupyter/notebooks/'+run_id
        if not self.ctx.is_simulating():
            if not os.path.isdir(folder_path):
                os.mkdir(folder_path)
            self.file_path = folder_path + '/KC_qPCR_time_log.txt'
        self.selected_pip = "right"
        self.pips = {"right": {}, "left": {}}

    def addStep(self, execute, description, wait_time=0):
        self.step_list.append(
            {'Execute': execute, 'description': description, 'wait_time': wait_time})

    def next_step(self):
        if self.step_list[self.step]['Execute'] == False:
            self.step += 1
            return False
        self.start = datetime.now()
        return True

    def finish_step(self):
        end = datetime.now()
        time_taken = (end - self.start)
        self.comment('Step ' + str(STEP) + ': ' +
                     self.step_list[self.step]['description'] + ' took ' + str(time_taken), add_hash=True)

        self.step_list[self.step]['Time'] = str(time_taken)
        self.step += 1

    def mount_right_pip(self, type, tip_racks, capacity):
        self.pips["right"]["pip"] = self.ctx.load_instrument(
            type, mount='right', tip_racks=tip_racks)
        self.pips["right"]["capacity"] = capacity
        self.pips["right"]["counts"] = 0
        self.pips["right"]["maxes"] = len(tip_racks)

    def mount_left_pip(self, type, tip_racks, capacity):
        self.pips["left"]["pip"] = self.ctx.load_instrument(
            type, mount='left', tip_racks=tip_racks)
        self.pips["right"]["capacity"] = capacity
        self.pips["left"]["counts"] = 0
        self.pips["left"]["maxes"] = len(tip_racks)

    def get_current_pip(self):
        return self.pips[self.selected_pip]["pip"]

    def get_pip_count(self):
        return self.pips[self.selected_pip]["count"]

    def reset_pip_count(self):
        self.pips[self.selected_pippet]["count"] = 0

    def get_pip_maxes(self):
        return self.pips[self.selected_pip]["maxes"]

    def get_pip_capacity(self):
        print(self.selected_pip)
        return self.pips[self.selected_pip]["capacity"]

    def set_pip(self, position):
        self.selected_pippet = position

    def custom_mix(self, reagent, location, vol, rounds, blow_out, mix_height,
                   x_offset, source_height=3, post_airgap=False, post_airgap_vol=10,
                   post_dispense=False, post_dispense_vol=20,):
        '''
        Function for mixing a given [vol] in the same [location] a x number of [rounds].
        blow_out: Blow out optional [True,False]
        x_offset = [source, destination]
        source_height: height from bottom to aspirate
        mix_height: height from bottom to dispense
        '''
        pip = self.get_current_pip()
        if mix_height == 0:
            mix_height = 3
        pip.aspirate(1, location=location.bottom(
            z=source_height).move(Point(x=x_offset[0])), rate=reagent.flow_rate_aspirate)
        for _ in range(rounds):
            pip.aspirate(vol, location=location.bottom(
                z=source_height).move(Point(x=x_offset[0])), rate=reagent.flow_rate_aspirate)
            pip.dispense(vol, location=location.bottom(
                z=mix_height).move(Point(x=x_offset[1])), rate=reagent.flow_rate_dispense)
        pip.dispense(1, location=location.bottom(
            z=mix_height).move(Point(x=x_offset[1])), rate=reagent.flow_rate_dispense)
        if blow_out == True:
            pip.blow_out(location.top(z=-2))  # Blow out
        if post_dispense == True:
            pip.dispense(post_dispense_vol, location.top(z=-2))
        if post_airgap == True:
            pip.dispense(post_airgap_vol, location.top(z=5))

    def pick_up(self):
        pip = self.get_current_pip()
        if not self.ctx.is_simulating():
            if self.get_pip_count() == self.get_pip_maxes():
                self.ctx.pause('Replace ' + str(pip.max_volume) + 'µl tipracks before \
                resuming.')
                pip.reset_tipracks()
                self.reset_pip_count()

        if not pip.hw_pipette['has_tip']:
            pip.pick_up_tip()

    def drop_tip(self):
        pip = self.get_current_pip()
        pip.drop_tip()
        pip['counts'] += 1

    def comment(self, comment, add_hash=False):
        hash_string = '#######################################################'
        if not self.ctx.is_simulating():
            if (add_hash):
                self.ctx.comment(hash_string)
            self.ctx.comment(comment)
            if (add_hash):
                self.ctx.comment(hash_string)
        else:
            if (add_hash):
                print(hash_string)
            print(comment)
            if (add_hash):
                print(hash_string)

    def move_vol_multichannel(self, reagent, source, dest, vol, air_gap_vol, x_offset,
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
        pip = self.get_current_pip()

        # Rinse before aspirating
        if rinse == True:
            run.custom_mix(reagent, location=source, vol=vol,
                           rounds=2, blow_out=True, mix_height=0,
                           x_offset=x_offset)
        # SOURCE
        s = source.bottom(pickup_height).move(Point(x=x_offset[0]))
        # aspirate liquid
        pip.aspirate(vol, s, rate=reagent.flow_rate_aspirate)
        if air_gap_vol != 0:  # If there is air_gap_vol, switch pipette to slow speed
            pip.aspirate(air_gap_vol, source.top(z=-2),
                         rate=reagent.flow_rate_aspirate)  # air gap
        # GO TO DESTINATION
        drop = dest.top(z=disp_height).move(Point(x=x_offset[1]))
        pip.dispense(vol + air_gap_vol, drop,
                     rate=reagent.flow_rate_dispense)  # dispense all
        # pause for x seconds depending on reagent
        ctx.delay(seconds=reagent.delay)
        if blow_out == True:
            pip.blow_out(dest.top(z=-2))
        if post_dispense == True:
            pip.dispense(post_dispense_vol, dest.top(z=-2))
        if touch_tip == True:
            pip.touch_tip(speed=20, v_offset=-5, radius=0.9)
        if post_airgap == True:
            pip.dispense(post_airgap_vol, dest.top(z=5), rate=2)

    def divide_volume(self, volume, max_vol):
        num_transfers = math.ceil(volume/max_vol)
        vol_roundup = math.ceil(volume/num_transfers)
        last_vol = volume - vol_roundup*(num_transfers-1)
        vol_list = [vol_roundup for v in range(1, num_transfers)]
        vol_list.append(last_vol)
        return vol_list

    def divide_destinations(self, l, n):
        a = []
        # Divide the list of destinations in size n lists.
        for i in range(0, len(l), n):
            a.append(l[i:i + n])

        return a

    def distribute_custom(self, reagent, volume, src, dest, waste_pool, pickup_height, extra_dispensal, disp_height=0):
        # Custom distribute function that allows for blow_out in different location and adjustement of touch_tip
        pip = self.get_current_pip()
        air_gap = 10
        pip.aspirate((len(dest) * volume) +
                     extra_dispensal, src.bottom(pickup_height), rate=reagent.flow_rate_aspirate)
        pip.touch_tip(speed=20, v_offset=-5)
        pip.move_to(src.top(z=5))
        pip.aspirate(air_gap, rate=reagent.flow_rate_aspirate)  # air gap

        for d in dest:
            pip.dispense(air_gap, d.top(), rate=reagent.flow_rate_dispense)
            drop = d.top(z=disp_height)
            pip.dispense(volume, drop, rate=reagent.flow_rate_dispense)
            # pause for x seconds depending on reagent
            self.ctx.delay(seconds=reagent.delay)
            pip.move_to(d.top(z=5))
            pip.aspirate(air_gap, d.top(
                z=5), rate=reagent.flow_rate_aspirate)  # air gap

        try:
            pip.blow_out(waste_pool.wells()[0].bottom(pickup_height + 3))
        except:
            pip.blow_out(waste_pool.bottom(pickup_height + 3))
        return (len(dest) * volume)

    def custom_mix(self, reagent, location, vol, rounds, blow_out, mix_height,
                   x_offset, source_height=3, post_airgap=False, post_airgap_vol=10,
                   post_dispense=False, post_dispense_vol=20,):
        '''
        Function for mixing a given [vol] in the same [location] a x number of [rounds].
        blow_out: Blow out optional [True,False]
        x_offset = [source, destination]
        source_height: height from bottom to aspirate
        mix_height: height from bottom to dispense
        '''
        pip = self.get_current_pip()

        if mix_height == 0:
            mix_height = 3
        pip.aspirate(1, location=location.bottom(
            z=source_height).move(Point(x=x_offset[0])), rate=reagent.flow_rate_aspirate)
        for _ in range(rounds):
            pip.aspirate(vol, location=location.bottom(
                z=source_height).move(Point(x=x_offset[0])), rate=reagent.flow_rate_aspirate)
            pip.dispense(vol, location=location.bottom(
                z=mix_height).move(Point(x=x_offset[1])), rate=reagent.flow_rate_dispense)
        pip.dispense(1, location=location.bottom(
            z=mix_height).move(Point(x=x_offset[1])), rate=reagent.flow_rate_dispense)
        if blow_out == True:
            pip.blow_out(location.top(z=-2))  # Blow out
        if post_dispense == True:
            pip.dispense(post_dispense_vol, location.top(z=-2))
        if post_airgap == True:
            pip.dispense(post_airgap_vol, location.top(z=5))

    def calc_height(self, reagent, cross_section_area, aspirate_volume, min_height=0.5, extra_volume=30):

        self.comment('Remaining volume ' + str(reagent.vol_well) +
                     '< needed volume ' + str(aspirate_volume) + '?')
        if reagent.vol_well < aspirate_volume + extra_volume:
            reagent.unused.append(reagent.vol_well)
            self.comment('Next column should be picked')
            self.comment('Previous to change: ' + str(reagent.col))
            # column selector position; intialize to required number
            reagent.col = reagent.col + 1
            self.comment(str('After change: ' + str(reagent.col)))
            reagent.vol_well = reagent.vol_well_original
            self.comment('New volume:' + str(reagent.vol_well))
            height = (reagent.vol_well - aspirate_volume -
                      reagent.v_cono) / cross_section_area
            # - reagent.h_cono
            reagent.vol_well = reagent.vol_well - aspirate_volume
            self.comment('Remaining volume:' + str(reagent.vol_well))
            if height < min_height:
                height = min_height
            col_change = True
        else:
            height = (reagent.vol_well - aspirate_volume -
                      reagent.v_cono) / cross_section_area  # - reagent.h_cono
            reagent.vol_well = reagent.vol_well - aspirate_volume
            self.comment('Calculated height is ' + str(height))
            if height < min_height:
                height = min_height
            self.comment('Used height is ' + str(height))
            col_change = False
        return height, col_change

    def start_lights(self):
        ctx._hw_manager.hardware.set_lights(
            rails=True)  # set lights off when using MMIX

    def stop_lights(self):
        ctx._hw_manager.hardware.set_lights(
            rails=False)  # set lights off when using MMIX

    def blink(self):
        for i in range(3):
            self.stop_lights()
            # ctx._hw_manager.hardware.set_button_light(1,0,0)
            time.sleep(0.3)
            self.start_lights()
            # ctx._hw_manager.hardware.set_button_light(0,0,1)
            time.sleep(0.3)
            self.stop_lights()

    def log_steps_time(self):
        # Export the time log to a tsv file
        if not self.ctx.is_simulating():
            with open(self.file_path, 'w') as f:
                f.write('STEP\texecution\tdescription\twait_time\texecution_time\n')
                for row in self.step_list:
                    for key in self.step_list[row].keys():
                        row += '\t' + format(self.step_list[row][key])
                    f.write(row + '\n')
            f.close()
