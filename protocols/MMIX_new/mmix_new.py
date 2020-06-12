import math
from opentrons.types import Point
from opentrons import robot
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
NUM_SAMPLES = 8
steps = []  # Steps you want to execute


air_gap_vol = 10
air_gap_mmix = 0
air_gap_sample = 0
run_id = '$run_id'

# Tune variables
select_mmix = "SonEspases1"  # Now only one recipe available
temperature = 10  # Temperature of temp module
volume_elution = 10  # Volume of the sample
extra_dispensal = 0  # Extra volume for master mix in each distribute transfer
diameter_screwcap = 8.1  # Diameter of the screwcap
elution_initial_volume = 50  # True
volume_cone = 57  # Volume in ul that fit in the screwcap cone
area_section_screwcap = (np.pi * diameter_screwcap**2) / 4
h_cone = (volume_cone * 3 / area_section_screwcap)
num_cols = math.ceil(NUM_SAMPLES/8)

#############################################################
# Available master mastermixes
#############################################################
MMIX_available = {'SonEspases1':
                  {
                      "recipe": [8.25, 6.25, 1.25],
                      "sources": ["D3", "C3", "B3"],
                      "dest": "D6",
                      "volume_mmix": 15,

                  }
                  }


MMIX_make = MMIX_available[select_mmix]
MMIX_make["volumes"] = []
for needed_vol in MMIX_make["recipe"]:
    MMIX_make["volumes"].append(needed_vol * NUM_SAMPLES * 1.1)
# Total volume of mastermix that will be prepared
MMIX_make["volume_available"] = sum(MMIX_make["volumes"])


def run(ctx: protocol_api.ProtocolContext):

    # Init protocol run
    run = ProtocolRun(ctx)

    run.addStep(description="Make MMIX")
    run.addStep(description="Transfer MMIX")
    run.addStep(description="Make MMIX")
    run.addStep(description="Set up positive control")

    # execute avaliaible steps
    if(len(steps) > 0):
        for index in steps:
            if(index <= len(run.step_list)):
                run.setExecutionStep(index-1, True)
            else:
                print("Step index out of range")
    else:
        # print(enumerate(run.step_list))
        for step in run.step_list:
            step['Execute'] = True
            # run.setExecutionStep(index['Execute'],True)

    ##################################
    # Define desk
    tempdeck = ctx.load_module('tempdeck', '10')
    tuberack = tempdeck.load_labware(
        'opentrons_24_aluminumblock_generic_2ml_screwcap')

    # tempdeck.set_temperature(temperature)

    # PCR
    pcr_plate = ctx.load_labware(
        'opentrons_96_aluminumblock_generic_pcr_strip_200ul', '11')

    # Eluted from King fisher/ Manual / Other
    elution_plate = ctx.load_labware(
        'biorad_96_wellplate_200ul_pcr', '8')

    # Tipracks20_multi
    tips20 = ctx.load_labware('opentrons_96_tiprack_20ul', 9)
    tips300 = ctx.load_labware('opentrons_96_filtertiprack_200ul', 7)

    # Mount pippets and set racks
    run.mount_right_pip('p20_single_gen2', tip_racks=[tips20], capacity=20)
    run.mount_left_pip('p300_single_gen2', tip_racks=[tips300], capacity=300)

    # Define wells interaction
    # Reagents and their characteristics

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
    positive_control = Reagent(name='Positive control',
                               rinse=False,
                               flow_rate_aspirate=1,
                               flow_rate_dispense=1,
                               reagent_reservoir_volume=50,
                               num_wells=1,  # change with num samples
                               delay=0,
                               h_cono=h_cone,
                               v_fondo=volume_cone  # V cono
                               )

    MMIX = Reagent(name='Master Mix',
                   rinse=False,
                   flow_rate_aspirate=1,
                   flow_rate_dispense=1,
                   reagent_reservoir_volume=MMIX_make["volume_available"],
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
                           reagent_reservoir_volume=elution_initial_volume,
                           delay=0,
                           num_wells=num_cols,  # num_cols comes from available columns
                           h_cono=0,
                           v_fondo=0
                           )

    MMIX_components = [mmix_water, taq_path, covid_assay]

    ################################################################################
    # Declare which reagents are in each reservoir as well as deepwell and elution plate
    # 1 row, 2 columns (first ones)
    MMIX_destination = tuberack.wells(MMIX_make["dest"])
    MMIX_components_location = []
    for source in MMIX_make["sources"]:
        MMIX_components_location.append(
            tuberack.wells(source))

    # setup up sample sources and destinations
    pcr_wells = pcr_plate.wells()[:NUM_SAMPLES]
    elution_wells = elution_plate.wells()[:NUM_SAMPLES]

    ############################################################################
    # STEP 1: Make Master MIX
    ############################################################################
    if (run.next_step()):
        run.stop_lights()
        run.comment('Selected MMIX: ' +
                    select_mmix, add_hash=True)

        run.set_pip("left")
        run.pick_up()
        drop = False
        for i, [source] in enumerate(MMIX_components_location):
            run.comment('Add component: ' +
                        MMIX_components[i].name, add_hash=True)

            # Get volumen calculated
            vol = MMIX_make["volumes"][i]
            # because 20ul is the maximum volume of the tip we will choose 17
            if (vol + air_gap_vol) > run.get_pip_capacity():
                # calculate what volume should be transferred in each step
                vol_list = run.divide_volume(vol, run.get_pip_capacity())
                for vol in vol_list:
                    # If not in first step we need to change everytime
                    if(i > 0):
                        run.pick_up()

                    run.move_vol_multichannel(reagent=MMIX_components[i], source=source, dest=MMIX_destination[0],
                                              vol=vol, air_gap_vol=air_gap_vol, pickup_height=0, disp_height=-10,
                                              blow_out=True)

                    # If not in first step we need to change everytime
                    if(i > 0):
                        run.drop_tip()
                        drop = True

            else:
                if(i > 0):
                    run.pick_up()
                run.move_vol_multichannel(reagent=MMIX_components[i], source=source, dest=MMIX_destination[0],
                                          vol=vol, air_gap_vol=air_gap_vol, pickup_height=-10,
                                          disp_height=-10, blow_out=True)
                if(i > 0):
                    run.drop_tip()
                    drop = True

            if i+1 < len(MMIX_components):
                if(not drop):
                    run.drop_tip()

            else:
                run.pick_up()
                run.comment('Final mix', add_hash=True)

                run.custom_mix(reagent=MMIX, location=MMIX_destination[0], vol=50, rounds=5,
                               blow_out=True, mix_height=2)
                run.drop_tip()

        run.finish_step()

    ############################################################################
    # STEP 2: Transfer Master MIX
    ############################################################################
    # run.start_lights()
    if (run.next_step()):
        run.set_pip("right")
        run.pick_up()
        volumen_mmix = MMIX_make["volume_available"]
        for dest in pcr_wells:
            [pickup_height, col_change] = run.calc_height(
                MMIX, area_section_screwcap, MMIX_make["volume_mmix"])

            run.move_vol_multichannel(reagent=MMIX, source=MMIX_destination[0],
                                      dest=dest, vol=MMIX_make["volume_mmix"], air_gap_vol=air_gap_mmix,
                                      pickup_height=pickup_height, disp_height=-10,
                                      blow_out=True, touch_tip=True)
            # change
            # mmix to positive and negativo control
        #    -> Positive
        run.move_vol_multichannel(reagent=positive_control, source=tuberack.wells('D6')[0],
                                  dest=pcr_plate.wells('H12')[0],
                                  vol=volume_elution, air_gap_vol=air_gap_sample,
                                  pickup_height=3, disp_height=-10,
                                  blow_out=True, touch_tip=True, post_airgap=True,)

        #    -> Negative
        run.move_vol_multichannel(reagent=positive_control, source=tuberack.wells('D6')[0],
                                  dest=pcr_plate.wells('G12')[0],
                                  vol=volume_elution, air_gap_vol=air_gap_sample,
                                  pickup_height=3, disp_height=-10,
                                  blow_out=True, touch_tip=True, post_airgap=True,)

        run.drop_tip()
        run.finish_step()
        tempdeck.deactivate()

    ############################################################################
    # STEP 3: TRANSFER Samples
    ############################################################################
    if(run.next_step()):
        run.comment('pcr_wells')
        run.set_pip("right")
        # Loop over defined wells
        for s, d in zip(elution_wells, pcr_wells):
            run.comment("%s %s" % (s, d))
            run.pick_up()
            # Source samples
            run.move_vol_multichannel(reagent=elution_well, source=s, dest=d,
                                      vol=volume_elution, air_gap_vol=air_gap_sample,
                                      pickup_height=3, disp_height=-10,
                                      blow_out=False, touch_tip=True, post_airgap=True,)
            run.custom_mix(reagent=elution_well, location=d, vol=8, rounds=3,
                           blow_out=False, mix_height=2)

            # ADD Custom mix
            run.drop_tip()

        run.finish_step()

    ############################################################################
    # STEP 4: Set up positive control
    ############################################################################
    if(run.next_step()):
        run.comment('pcr_wells')
        run.set_pip("right")
        run.pick_up()

        # Negative control
        run.move_vol_multichannel(reagent=positive_control, source=elution_plate.wells('G12')[0],
                                  dest=pcr_plate.wells('G12')[0],
                                  vol=volume_elution, air_gap_vol=air_gap_sample,
                                  pickup_height=3, disp_height=-10,
                                  blow_out=True, touch_tip=True, post_airgap=True)
        run.custom_mix(reagent=positive_control, location=pcr_plate.wells('G12')[0], vol=8, rounds=3,
                       blow_out=False, mix_height=2)

        run.change_tip()

        # Positive Control

        run.move_vol_multichannel(reagent=positive_control, source=tuberack.wells('A6')[0],
                                  dest=pcr_plate.wells('H12')[0],
                                  vol=volume_elution, air_gap_vol=air_gap_sample,
                                  pickup_height=3, disp_height=-10,
                                  blow_out=True, touch_tip=True, post_airgap=True)
        run.custom_mix(reagent=positive_control, location=pcr_plate.wells('H12')[0], vol=8, rounds=3,
                       blow_out=False, mix_height=2)

        run.drop_tip()
        run.finish_step()

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

    def addStep(self, description, execute=False, wait_time=0):
        self.step_list.append(
            {'Execute': execute, 'description': description, 'wait_time': wait_time})

    def setExecutionStep(self, index, value):
        self.step_list[index]["Execute"] = value

    def next_step(self):
        robot.clear_commands()
        # print(self.step_list[self.step]['Execute'])
        if self.step_list[self.step]['Execute'] == False:
            self.step += 1
            return False
        self.start = datetime.now()
        return True

    def finish_step(self):
        for c in robot.commands():
            print(c)
        end = datetime.now()
        time_taken = (end - self.start)
        self.comment('Step ' + str(self.step + 1) + ': ' +
                     self.step_list[self.step]['description'] + ' took ' + str(time_taken), add_hash=True)

        self.step_list[self.step]['Time'] = str(time_taken)
        self.step += 1

    def mount_pip(self, position, type, tip_racks, capacity):
        self.pips[position]["pip"] = self.ctx.load_instrument(
            type, mount=position, tip_racks=tip_racks)
        self.pips[position]["capacity"] = capacity
        self.pips[position]["count"] = 0
        self.pips[position]["maxes"] = len(tip_racks)

    def mount_right_pip(self, type, tip_racks, capacity):
        self.mount_pip("right", type, tip_racks, capacity)

    def mount_left_pip(self, type, tip_racks, capacity):
        self.mount_pip("left", type, tip_racks, capacity)

    def get_current_pip(self):
        return self.pips[self.selected_pip]["pip"]

    def get_pip_count(self):
        return self.pips[self.selected_pip]["count"]

    def reset_pip_count(self):
        self.pips[self.selected_pip]["count"] = 0

    def add_pip_count(self):
        self.pips[self.selected_pip]["count"] + 1

    def get_pip_maxes(self):
        return self.pips[self.selected_pip]["maxes"]

    def get_pip_capacity(self):
        return self.pips[self.selected_pip]["capacity"]

    def set_pip(self, position):
        self.selected_pip = position

    def custom_mix(self, reagent, location, vol, rounds, mix_height, blow_out=False,
                   source_height=3, post_airgap=False, post_airgap_vol=10,
                   post_dispense=False, post_dispense_vol=20, x_offset=[0, 0]):
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
        self.add_pip_count()

    def change_tip(self):
        self.drop_tip()
        self.pick_up()

    def comment(self, comment, add_hash=False):
        hash_string = '#######################################################'
        if not self.ctx.is_simulating():
            if (add_hash):
                robot.comment(hash_string)
            robot.comment(('{}').format(comment))
            if (add_hash):
                robot.comment(hash_string)
        else:
            if (add_hash):
                print(hash_string)
            print(comment)
            if (add_hash):
                print(hash_string)

    def move_vol_multichannel(self, reagent, source, dest, vol, air_gap_vol,
                              pickup_height, disp_height, x_offset=[0, 0],
                              rinse=False,  blow_out=False, touch_tip=False,
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
        if (s.point.z < source.bottom().point.z):
            self.comment("Pickup height too low you will hit the bottom")
            self.comment(s.point.z)
            self.comment(source.bottom().point.z)
            return False

        if (s.point.z > source.top().point.z):
            self.comment("Pickup too high you will not get any liquid")
            self.comment(s.point.z)
            self.comment(source.top().point.z)
            return False

        # aspirate liquid
        pip.aspirate(vol, s, rate=reagent.flow_rate_aspirate)
        if air_gap_vol != 0:  # If there is air_gap_vol, switch pipette to slow speed
            pip.aspirate(air_gap_vol, source.top(z=-2),
                         rate=reagent.flow_rate_aspirate)  # air gap
        # GO TO DESTINATION
        drop = dest.top(z=disp_height).move(Point(x=x_offset[1]))
        if (drop.point.z < dest.bottom().point.z):
            self.comment("Dispense height too low you will hit the bottom")
            self.comment(drop.point.z)
            self.comment(dest.bottom().point.z)
            return False

        pip.dispense(vol + air_gap_vol, drop,
                     rate=reagent.flow_rate_dispense)  # dispense all
        # pause for x seconds depending on reagent
        self.ctx.delay(seconds=reagent.delay)
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
                   x_offset=[0, 0], source_height=3, post_airgap=False, post_airgap_vol=10,
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
        # if support_selected == pcr_support.index[1] : --> refdefine height (calculate_heigh(self))
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
        self.ctx._hw_manager.hardware.set_lights(
            rails=True)  # set lights off when using ººººººº

    def stop_lights(self):
        self.ctx._hw_manager.hardware.set_lights(
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
                row = ""
                '''for step in self.step_list:
                    row = ('{}\t{}\t{}\t{}\t{}').format(
                        row, step["Execution"], step["description"], step["wait_time"], step["execution_time"])
                    total_time += data["execution_time"]
                    f.write(row + '\n')'''
            f.close()
