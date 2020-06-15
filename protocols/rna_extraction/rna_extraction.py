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
    'protocolName': 'RNA Extraction Version 2',
    'author': 'Matias Bonet & Antoni Morla. based on: Malen Aguirregabiria,Aitor Gastaminza & José Luis Villanueva (jlvillanueva@clinic.cat)',
    'source': 'Hospital Son Espases Palma',
    'apiLevel': '2.3',
    'description': 'Protocol for rna extraction'
}

'''
'technician': 'Toni',
'date': '$date'
'''
# Defined variables
##################
NUM_SAMPLES = 8
steps = []  # Steps you want to execute
set_temp_on = False  # Do you want to start temperature module?
temperature = 65  # Set temperature. It will be uesed if set_temp_on is set to True
set_mag_on = False  # Do you want to start magnetic module?
mag_height = 14  # Height needed for NEST deepwell in magnetic deck

use_waits = True

num_cols = math.ceil(NUM_SAMPLES/8)

diameter_screwcap = 8.1  # Diameter of the screwcap
volume_cone = 57  # Volume in ul that fit in the screwcap cone
area_section_screwcap = (np.pi * diameter_screwcap**2) / 4
h_cone = (volume_cone * 3 / area_section_screwcap)


air_gap_vol = 10
air_gap_r1 = 0
air_gap_sample = 0
run_id = 'testing'


##################
# Custom function
##################
# Define Reagents as objects with their properties
class Reagent:
    def __init__(self, name, flow_rate_aspirate, reagent_volume, reagent_reservoir_volume,
                 flow_rate_dispense, flow_rate_aspirate_mix, flow_rate_dispense_mix,
                 air_gap_vol_bottom, air_gap_vol_top, disposal_volume,  max_volume_allowed,
                 num_wells, h_cono, v_fondo, area_section_screwcap = 8.254,rinse=False, delay=0, tip_recycling='none'):
        self.name = name
        self.flow_rate_aspirate = flow_rate_aspirate
        self.flow_rate_dispense = flow_rate_dispense
        self.flow_rate_aspirate_mix = flow_rate_aspirate_mix
        self.flow_rate_dispense_mix = flow_rate_dispense_mix
        self.air_gap_vol_bottom = air_gap_vol_bottom
        self.air_gap_vol_top = air_gap_vol_top
        self.disposal_volume = disposal_volume
        self.rinse = bool(rinse)
        self.max_volume_allowed = max_volume_allowed
        self.reagent_volume = reagent_volume
        self.reagent_reservoir_volume = reagent_reservoir_volume
        self.num_wells = num_wells
        self.col = 0
        self.h_cono = h_cono
        self.v_cono = v_fondo
        self.tip_recycling = tip_recycling
        self.vol_well_original = reagent_reservoir_volume / num_wells
        self.vol_well = self.vol_well_original
        self.unused = []
        self.delay = delay
        self.area_section_screwcap = area_section_screwcap


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

    def add_step(self, description, execute=False, wait_time=0):
        self.step_list.append(
            {'Execute': execute, 'description': description, 'wait_time': wait_time, 'execution_time': 0})

    def init_steps(self, steps):
        if(len(steps) > 0):
            for index in steps:
                if(index <= len(self.step_list)):
                    self.set_execution_step(index-1, True)
                else:
                    print("Step index out of range")
        else:
            for index, step in enumerate(self.step_list):
                self.set_execution_step(index, True)

    def set_execution_step(self, index, value):
        self.step_list[index]["Execute"] = value

    def get_current_step(self):
        return self.step_list[self.step]

    def next_step(self):
        if self.step_list[self.step]['Execute'] == False:
            self.step += 1
            return False
        self.start = datetime.now()
        return True

    def finish_step(self):
        if (self.get_current_step()["wait_time"] > 0 and use_waits):
            self.ctx.delay(seconds=int(self.get_current_step()[
                           "wait_time"]), msg=self.get_current_step()["description"])
        if (self.get_current_step()["wait_time"] > 0 and not use_waits):
            self.comment("We simulate a wait of:%s seconds" %
                         self.get_current_step()["wait_time"])
        end = datetime.now()
        time_taken = (end - self.start)
        self.comment('Step ' + str(self.step + 1) + ': ' +
                     self.step_list[self.step]['description'] + ' took ' + str(time_taken), add_hash=True)

        self.step_list[self.step]['execution_time'] = str(time_taken)
        self.step += 1
        self.log_steps_time()

    def log_steps_time(self):
        # Export the time log to a tsv file
        if not self.ctx.is_simulating():
            with open(self.file_path, 'w') as f:
                f.write('STEP\texecution\tdescription\twait_time\texecution_time\n')
                row = ""
                for step in self.step_list:
                    row = ('{}\t{}\t{}\t{}\t{}').format(
                        row, step["Execute"], step["description"], step["wait_time"], step["execution_time"])
                    f.write(row + '\n')
            f.close()

    def mount_pip(self, position, type, tip_racks, capacity, multi=False, size_tipracks=96):
        self.pips[position]["pip"] = self.ctx.load_instrument(
            type, mount=position, tip_racks=tip_racks)
        self.pips[position]["capacity"] = capacity
        self.pips[position]["count"] = 0
        self.pips[position]["maxes"] = len(tip_racks)*size_tipracks
        if(multi):
            self.pips[position]["increment_tips"] = 8
        else:
            self.pips[position]["increment_tips"] = 1

    def mount_right_pip(self, type, tip_racks, capacity, multi=False):
        self.mount_pip("right", type, tip_racks, capacity)

    def mount_left_pip(self, type, tip_racks, capacity, multi=False):
        self.mount_pip("left", type, tip_racks, capacity)

    def get_current_pip(self):
        return self.pips[self.selected_pip]["pip"]

    def get_pip_count(self):
        return self.pips[self.selected_pip]["count"]

    def reset_pip_count(self):
        self.pips[self.selected_pip]["count"] = 0

    def add_pip_count(self):
        self.pips[self.selected_pip]["count"] + \
            self.pips[self.selected_pip]["increment_tips"]

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
                self.ctx.comment(hash_string)
            self.ctx.comment(('{}').format(comment))
            if (add_hash):
                self.ctx.comment(hash_string)
        else:
            if (add_hash):
                print(hash_string)
            print(comment)
            if (add_hash):
                print(hash_string)

    def pause(self, comment):
        if not self.ctx.is_simulating():
            self.ctx.pause(comment)
        else:
            print("%s\n Press any key to continue " % comment)

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
            self.custom_mix(reagent, location=source, vol=vol,
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

    def calc_height(self, reagent, min_height=0.5, extra_volume=30):
        # if support_selected == pcr_support.index[1] : --> refdefine height (calculate_heigh(self))
        debug = False
        if (debug):
            self.comment('Remaining volume ' + str(reagent.vol_well) +
                         '< needed volume ' + str(aspirate_volume) + '?')
        
        cross_section_area= reagent.area_section_screwcap
        aspirate_volume = reagent.disposal_volume        

        if reagent.vol_well < aspirate_volume + extra_volume:
            reagent.unused.append(reagent.vol_well)
            if (debug):
                self.comment('Next column should be picked')
                self.comment('Previous to change: ' + str(reagent.col))
            # column selector position; intialize to required number
            reagent.col = reagent.col + 1
            if (debug):
                self.comment(str('After change: ' + str(reagent.col)))
            reagent.vol_well = reagent.vol_well_original
            if (debug):
                self.comment('New volume:' + str(reagent.vol_well))
            height = (reagent.vol_well - aspirate_volume -
                      reagent.v_cono) / cross_section_area
            # - reagent.h_cono
            reagent.vol_well = reagent.vol_well - aspirate_volume
            if (debug):
                self.comment('Remaining volume:' + str(reagent.vol_well))
            if height < min_height:
                height = min_height
            col_change = True
        else:
            height = (reagent.vol_well - aspirate_volume -
                      reagent.v_cono) / cross_section_area  # - reagent.h_cono
            reagent.vol_well = reagent.vol_well - aspirate_volume
            if (debug):
                self.comment('Calculated height is ' + str(height))
            if height < min_height:
                height = min_height
            if (debug):
                self.comment('Used height is ' + str(height))
            col_change = False
        return height, col_change

    def start_lights(self):
        self.ctx._hw_manager.hardware.set_lights(
            rails=True)  # set lights off when using MMIX

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


def run(ctx: protocol_api.ProtocolContext):

    # Init protocol run
    run = ProtocolRun(ctx)
    # yo creo que este tiene que ser manual o sacarlo a otro robot
    run.add_step(description="Transfer A6 - To AW_PLATE Single Slot1 -> Slot2")
    run.add_step(description="Wait until bell is done")  # INTERACTION
    run.add_step(description="Transfer BBUIX 3 - 2 Multi")
    run.add_step(description="Transfer Beats - To AW_PLATE Multi")
    run.add_step(
        description="Replace tips, empty trash, move Slot2 -> Slot 10")  # INTERACTION

    run.add_step(description="65C Incubation", wait_time=5 * 60)  # 5 minutos
    run.add_step(description="Transfer From temperature to magnet 485ul")
    run.add_step(description="Magnetic on: 10 minutes", wait_time=10) #10*60
    run.add_step(
        description="Extraer liquido sin tocar los beats. Slot 7 - Piscina Slot 3")
    run.add_step(description="Magnetic off")

    run.add_step(
        description="Replace tips, add WB, add ETOH, vaciar piscina y trash. Cambiar nuevo DW SLOT 10")  # INTERACTION

    # Add WB
    run.add_step(description="Add 500ul de WB  a los beats Slot 4 - 7 ")
    run.add_step(description="Magnetic on: 2 minutes", wait_time=2 ) # 2*60
    run.add_step(
        description="Extraer liquido sin tocar los beats. Slot 7 - Piscina Slot 3")
    run.add_step(description="Magnetic off")

    # Add ETOH First step
    run.add_step(description="Add 500ul de etoh a los beats Slot 8 - 7 ")
    run.add_step(description="Magnetic on: 2 minutes", wait_time=2) #2*60
    run.add_step(
        description="Extraer liquido sin tocar los beats. Slot 7 - Piscina Slot 3")
    run.add_step(description="Magnetic off")

    # Add ETOH First step
    run.add_step(description="Add 500ul de etoh a los beats Slot 8 - 7 ")
    run.add_step(description="Magnetic on: 2 minutes", wait_time=2) #2*60
    run.add_step(
        description="Extraer liquido sin tocar los beats. Slot 7 - Piscina Slot 3")
    run.add_step(description="Magnetic off")

    # Add ETOH Second step
    run.add_step(description="Add 250ul de etoh a los beats Slot 8 - 7 ")
    run.add_step(description="Magnetic on: 2 minutes", wait_time=10)# 2*60
    run.add_step(
        description="Extraer liquido sin tocar los beats. Slot 7 - Piscina Slot 3")
    run.add_step(description="Secar durante 10 minutos", wait_time=10) #10 * 60
    run.add_step(description="Magnetic off")
    run.add_step(
        description="Add elution move to temperature same tip 4 -> 7 -> 10")
    run.add_step(description="65C Incubation 10'", wait_time=10) # 10 * 60
    run.add_step(description="Move 50ul from temp to magnet 10-7")
    run.add_step(description="Magnetic on: 3 minutes", wait_time=3 * 60) # 3 * 60
    #run.addStep(description="Move 50ul Magnet Final destination 7-> 2")

    # execute avaliaible steps
    run.init_steps(steps)

    ##################################
    # Define desk
    moving_type = "biorad_96_wellplate_200ul_pcr"

    # Tube rack
    tube_rack = ctx.load_labware(
        'opentrons_24_tuberack_nest_1.5ml_screwcap', 1)

    # Destination plate SLOT 2
    aw_slot = ctx.load_labware(moving_type, 2)
    aw_wells = aw_slot.wells()[:NUM_SAMPLES]
    aw_wells_multi = aw_slot.rows()[0][:num_cols]

    # Magnetic Beads Pool
    beads_slot = ctx.load_labware(
        'nest_12_reservoir_15ml', 3)
    beads_wells_multi = beads_slot.rows()[0][:num_cols]

    # setup up sample sources and destinations
    # Wash Buffer Pool-
    wb_slot = ctx.load_labware(
        'nest_12_reservoir_15mL', 4)
    wb_wells_multi = wb_slot.rows()[0][:num_cols]

    # # Magnetic module plus NEST_Deep_well_reservoire
    magdeck = ctx.load_module('magnetic module gen2', 7)
    magdeck.disengage()
    mag_slot = magdeck.load_labware(moving_type)
    mag_wells_multi = mag_slot.rows()[0][:num_cols]

    # Ethanol Pool
    etoh_slot = ctx.load_labware(
        'nest_12_reservoir_15ml', 8)
    etoh_wells_multi = etoh_slot.rows()[0][:num_cols]

    # Temperature module plus NEST_Deep_well_reservoire
    tempdeck = ctx.load_module('tempdeck', 10)
    temp_slot = tempdeck.load_labware(moving_type)
    temp_wells_multi = temp_slot.rows()[0][:num_cols]

    # Mount pippets and set racks
    # Tipracks20_multi
    tips20 = ctx.load_labware('opentrons_96_tiprack_20ul', 11)
    tips300_9 = ctx.load_labware('opentrons_96_filtertiprack_200ul', "9")
    tips300_6 = ctx.load_labware('opentrons_96_filtertiprack_200ul', "6")
    tips300_5 = ctx.load_labware('opentrons_96_filtertiprack_200ul', "5")

    run.mount_right_pip('p20_single_gen2', tip_racks=[tips20], capacity=20)
    run.mount_left_pip('p300_multi_gen2', tip_racks=[
                       tips300_9, tips300_6,tips300_5], capacity=200, multi=True)

    # Reagents and their characteristics
    WB = Reagent(name='WB washing buffer',
                 flow_rate_aspirate=3,
                 flow_rate_dispense=3,
                 flow_rate_aspirate_mix=15,
                 flow_rate_dispense_mix=25,
                 air_gap_vol_bottom=5,
                 air_gap_vol_top=0,
                 disposal_volume=1,
                 max_volume_allowed=180,
                 reagent_volume=500,
                 reagent_reservoir_volume=(
                      NUM_SAMPLES + 5) * 500,  # 60000, #38400
                 # num_Wells max is 4
                 num_wells=math.ceil((NUM_SAMPLES + 5) * 500 / 13000),
                 h_cono=1.95,
                 v_fondo=750,  # 1.95 * multi_well_rack_area / 2, #Prismatic
                 tip_recycling='A1')


    aw_well = Reagent(name='dw_plate well',
                      num_wells=1,  # change with num samples
                      delay=0,
                      flow_rate_aspirate=3,  # Original 0.5
                      flow_rate_dispense=3,  # Original 1
                      flow_rate_aspirate_mix=15,
                      flow_rate_dispense_mix=25,
                      air_gap_vol_bottom=5,
                      air_gap_vol_top=0,
                      disposal_volume=1,
                      max_volume_allowed=150,
                      reagent_volume=50,
                      reagent_reservoir_volume=150,
                      h_cono=4,
                      v_fondo=4 * math.pi * 4 ** 3 / 3
                      )

    ############################################################################
    # STEP 1: Transfer A6 - To AW_PLATE
    ############################################################################
    if (run.next_step()):
        run.set_pip("right")  # single 20
        volumen_move = 5
        source = tube_rack.wells("A6")[0]
        liquid = Reagent(name='Proteinasa K',
                         num_wells=1,  # change with num samples
                         delay=0,
                         flow_rate_aspirate=3,  # Original 0.5
                         flow_rate_dispense=3,  # Original 1
                         flow_rate_aspirate_mix=15,
                         flow_rate_dispense_mix=25,
                         air_gap_vol_bottom=5,
                         air_gap_vol_top=0,
                         disposal_volume=volumen_move,
                         reagent_volume=volumen_move*NUM_SAMPLES,
                         max_volume_allowed=150,
                         reagent_reservoir_volume=volumen_move*NUM_SAMPLES,  # 14800,
                         h_cono=4,
                         v_fondo=4 * math.pi * 4 ** 3 / 3,
                         area_section_screwcap = (np.pi * 8.25**2) / 4
                         )


        
        run.pick_up()
        for dest in aw_wells:
            [pickup_height, col_change] = run.calc_height(liquid)

            run.move_vol_multichannel(reagent=liquid, source=source,
                                      dest=dest, vol=volumen_move, air_gap_vol=air_gap_r1,
                                      pickup_height=pickup_height, disp_height=-10,
                                      blow_out=True, touch_tip=True)

        
        run.drop_tip()
        run.finish_step()

    ############################################################################
    # STEP 2: Pause until the hood is done
    ############################################################################
    if (run.next_step()):
        run.blink()
        ctx.pause('Go to the hood to disable sample')

        run.finish_step()

    ############################################################################
    # STEP 3: Slot 3 -2 beats_PK AW
    ############################################################################
    if (run.next_step()):
        ############################################################################
        # Light flash end of program
        run.set_pip("left")  # p300 multi
         
        liquid = Reagent(name='Magnetic beads',
                       flow_rate_aspirate=1,
                       flow_rate_dispense=1.5,
                       flow_rate_aspirate_mix=1.5,
                       flow_rate_dispense_mix=5,
                       air_gap_vol_bottom=5,
                       air_gap_vol_top=0,
                       disposal_volume=275,
                       rinse=True,
                       max_volume_allowed=180,
                       reagent_volume=250,
                       reagent_reservoir_volume=NUM_SAMPLES * 250 * 1.1,  
                       num_wells=num_cols,
                       h_cono=1.95,
                       v_fondo=750)
        air_gap_vol = 3
        disposal_height = -5
        pickup_height = 1

        for source, destination in zip(aw_wells_multi, beads_wells_multi):
            run.pick_up()
            run.move_vol_multichannel(reagent=liquid, source=source,
                                      dest=destination, vol=150, air_gap_vol=air_gap_vol,
                                      pickup_height=pickup_height, disp_height=disposal_height,
                                      rinse=True)
            run.change_tip()
            run.move_vol_multichannel(reagent=liquid, source=source,
                                      dest=destination, vol=125, air_gap_vol=air_gap_vol,
                                      pickup_height=pickup_height, disp_height=disposal_height,rinse=True)
            run.drop_tip()

        run.finish_step()

    ############################################################################
    # STEP 4: Transfer Beats - To AW_PLATE
    ############################################################################
    if (run.next_step()):
        
        run.set_pip("right")  # single 20
        volumen_move = 5
        source = tube_rack.wells("B6")[0]
        liquid = Reagent(name='MS2',
                         num_wells=1,  # change with num samples
                         delay=0,
                         flow_rate_aspirate=3,  # Original 0.5
                         flow_rate_dispense=3,  # Original 1
                         flow_rate_aspirate_mix=15,
                         flow_rate_dispense_mix=25,
                         air_gap_vol_bottom=5,
                         air_gap_vol_top=0,
                         disposal_volume=volumen_move,
                         reagent_volume=volumen_move*NUM_SAMPLES,
                         max_volume_allowed=150,
                         reagent_reservoir_volume=volumen_move*NUM_SAMPLES,  # 14800,
                         h_cono=4,
                         v_fondo=4 * math.pi * 4 ** 3 / 3,
                         area_section_screwcap = (np.pi * 8.25**2) / 4
                         )
       
        
        for dest in aw_wells:
            run.pick_up()
            [pickup_height, col_change] = run.calc_height(liquid)

            run.move_vol_multichannel(reagent=liquid, source=source,
                                      dest=dest, vol=volumen_move, air_gap_vol=air_gap_r1,
                                      pickup_height=pickup_height, disp_height=-10,
                                      blow_out=True, touch_tip=True)

            run.drop_tip()

        run.finish_step()

    ############################################################################
    # STEP 5: Mix and Pause to replace
    ############################################################################
    if (run.next_step()):
        run.set_pip("left")
        liquid = Reagent(name='PK+Beads+MS2 Mix',
                         num_wells=1,  # change with num samples
                         delay=0,
                         flow_rate_aspirate=3,  # Original 0.5
                         flow_rate_dispense=3,  # Original 1
                         flow_rate_aspirate_mix=15,
                         flow_rate_dispense_mix=25,
                         air_gap_vol_bottom=5,
                         air_gap_vol_top=0,
                         disposal_volume=10,
                         reagent_volume=10,
                         max_volume_allowed=500,
                         rinse=False,
                         reagent_reservoir_volume=10*NUM_SAMPLES, 
                         h_cono=4,
                         v_fondo=4 * math.pi * 4 ** 3 / 3,
                         area_section_screwcap = (np.pi * 8.25**2) / 4
                         )

        for source in aw_wells_multi:
            run.pick_up()
            run.custom_mix(liquid, location=source, vol=100,
                           rounds=10, blow_out=True, mix_height=0)
            run.drop_tip()

        run.blink()
        ctx.pause('Replace tips, empty trash, move Slot2 -> Slot 10')
        run.finish_step()

    ############################################################################
    # STEP 6: Incubation at 65ºC
    ############################################################################
    if (run.next_step()):
        if (set_temp_on):
            tempdeck.set_temperature(temperature)
        run.finish_step()
        tempdeck.deactivate()

    ############################################################################
    # STEP 7: Transfer From temperature to magnet 485ul
    ############################################################################
    if (run.next_step()):
        
        run.set_pip("left")  # p300 multi 
        liquid = Reagent(name='MIX_HOT',
                       flow_rate_aspirate=1,
                       flow_rate_dispense=1.5,
                       flow_rate_aspirate_mix=1.5,
                       flow_rate_dispense_mix=5,
                       air_gap_vol_bottom=5,
                       air_gap_vol_top=0,
                       disposal_volume=485,
                       rinse=True,
                       max_volume_allowed=500,#No aplica
                       reagent_volume=0, #No aplica
                       reagent_reservoir_volume=NUM_SAMPLES * 250 * 1.1,  #No aplica
                       num_wells=num_cols, # multi
                       h_cono=1.95,
                       v_fondo=750)

        air_gap_vol = 3
        disposal_height = -5
        pickup_height = 1

        for source, destination in zip(temp_wells_multi, mag_wells_multi):
            run.pick_up()
            run.move_vol_multichannel(reagent=liquid, source=source,
                                      dest=destination, vol=175, air_gap_vol=air_gap_vol,
                                      pickup_height=pickup_height, disp_height=disposal_height,
                                      rinse=True)
            run.move_vol_multichannel(reagent=liquid, source=source,
                                      dest=destination, vol=175, air_gap_vol=air_gap_vol,
                                      pickup_height=pickup_height, disp_height=disposal_height,
                                      rinse=True)
            run.move_vol_multichannel(reagent=liquid, source=source,
                                      dest=destination, vol=135, air_gap_vol=air_gap_vol,
                                      pickup_height=pickup_height, disp_height=disposal_height,rinse=True)
            run.drop_tip()

        run.finish_step()

    # Extraer liquido sin tocar los beats. Slot 7 - Piscina Slot 3
    def move_magnet_to_trash(move_vol_steps = 3):
        run.set_pip("left")  # p300 multi
        # Sobre nadante primer paso
        liquid = Reagent(name='Sobrenadante',
                       flow_rate_aspirate=1,
                       flow_rate_dispense=1.5,
                       flow_rate_aspirate_mix=1.5,
                       flow_rate_dispense_mix=5,
                       air_gap_vol_bottom=5,
                       air_gap_vol_top=0,
                       disposal_volume=485,
                       rinse=True,
                       max_volume_allowed=500,#No aplica
                       reagent_volume=0, #No aplica
                       reagent_reservoir_volume=NUM_SAMPLES * 250 * 1.1,  #No aplica
                       num_wells=num_cols, # multi
                       h_cono=1.95,
                       v_fondo=750)

        air_gap_vol = 3
        pickup_height = 1
        disposal_height = 0
        # Hay que revisar los offsets para el movimiento este
        for source, destination in zip(mag_wells_multi, beads_wells_multi):
            # Replace this
            run.pick_up()
            run.move_vol_multichannel(reagent=liquid, source=source,
                                      dest=destination, vol=175, air_gap_vol=air_gap_vol,
                                      pickup_height=pickup_height, disp_height=disposal_height,
                                      rinse=True)
            
            # Patch for last step of etho 250ul instead of 500
            if(move_vol_steps==3):
                run.move_vol_multichannel(reagent=liquid, source=source,
                                        dest=destination, vol=175, air_gap_vol=air_gap_vol,
                                        pickup_height=pickup_height, disp_height=disposal_height,
                                        rinse=True)

            # We want to empty does not matter if we aspirate more                         
            run.move_vol_multichannel(reagent=liquid, source=source,
                                      dest=destination, vol=175, air_gap_vol=air_gap_vol,
                                      pickup_height=pickup_height, disp_height=disposal_height,rinse=True)
            run.drop_tip()

    ############################################################################
    # STEP 8: Magnet on 10 minutos
    ############################################################################
    if (run.next_step()):
        if (set_mag_on):
            magdeck.engage(height=mag_height)
        run.finish_step()

    ############################################################################
    # STEP 9: Extract liquid from magnet to liquid trash
    ############################################################################
    if (run.next_step()):
        move_magnet_to_trash()
        run.finish_step()

    ############################################################################
    # STEP 10: Magnet off
    ############################################################################
    if (run.next_step()):
        magdeck.disengage()
        run.finish_step()

    ############################################################################
    # STEP 11: Pause to replace
    ############################################################################
    if (run.next_step()):
        run.blink()
        ctx.pause(
            'Replace tips, add WB, add ETOH, vaciar piscina y trash. Cambiar nuevo DW SLOT 10')
        run.finish_step()

    ############################################################################
    # STEP 12: Add 500ul de WB a los bits 4 - 7
    ############################################################################
    if (run.next_step()):
        run.set_pip("left")  # p300 multi 
        liquid = Reagent(name='WB',
                       flow_rate_aspirate=1,
                       flow_rate_dispense=1.5,
                       flow_rate_aspirate_mix=1.5,
                       flow_rate_dispense_mix=5,
                       air_gap_vol_bottom=5,
                       air_gap_vol_top=0,
                       disposal_volume=485,
                       rinse=True,
                       max_volume_allowed=500,#No aplica
                       reagent_volume=0, #No aplica
                       reagent_reservoir_volume=NUM_SAMPLES * 250 * 1.1,  #No aplica
                       num_wells=num_cols, # multi
                       h_cono=1.95,
                       v_fondo=750)

        air_gap_vol = 3
        disposal_height = -1 # Arriba y el último paso lo hacemos dentro
        pickup_height = 1 

        for source, destination in zip(wb_wells_multi, mag_wells_multi):
            run.pick_up()
            run.move_vol_multichannel(reagent=liquid, source=source,
                                      dest=destination, vol=175, air_gap_vol=air_gap_vol,
                                      pickup_height=pickup_height, disp_height=disposal_height
                                     )
            run.move_vol_multichannel(reagent=liquid, source=source,
                                      dest=destination, vol=175, air_gap_vol=air_gap_vol,
                                      pickup_height=pickup_height, disp_height=disposal_height)

            # This will be drop inside
            [disposal_height,column_change] = run.calc_height(liquid)
            run.move_vol_multichannel(reagent=liquid, source=source,
                                      dest=destination, vol=135, air_gap_vol=air_gap_vol,
                                      pickup_height=pickup_height, disp_height=disposal_height-3)
            
            run.custom_mix(liquid, location=source, vol=150,
                           rounds=10, blow_out=True, mix_height=0)
            run.drop_tip()

        run.finish_step()


    ############################################################################
    # STEP 13: Magnet on 10 minutos
    ############################################################################
    if (run.next_step()):
        if (set_mag_on):
            magdeck.engage(height=mag_height)
        run.finish_step()
    ############################################################################
    # STEP 14: Extract liquid from magnet to liquid trash
    ############################################################################
    if (run.next_step()):
        move_magnet_to_trash()
        run.finish_step()

    ############################################################################
    # STEP 15: Magnet off
    ############################################################################
    if (run.next_step()):
        magdeck.disengage()
        run.finish_step()

    # Used twice in the next steps
    etoh = Reagent(name='Etoh',
                       flow_rate_aspirate=1,
                       flow_rate_dispense=1.5,
                       flow_rate_aspirate_mix=1.5,
                       flow_rate_dispense_mix=5,
                       air_gap_vol_bottom=5,
                       air_gap_vol_top=0,
                       disposal_volume=485,
                       rinse=True,
                       max_volume_allowed=500,#No aplica
                       reagent_volume=0, #No aplica
                       reagent_reservoir_volume=NUM_SAMPLES * 250 * 1.1,  #No aplica
                       num_wells=num_cols, # multi
                       h_cono=1.95,
                       v_fondo=750)
    ############################################################################
    # STEP 16: Add 500ul de etoh a los beats Slot 8 - 7
    ############################################################################
    if (run.next_step()):

        run.set_pip("left")  # p300 multi 
        liquid = etoh
        air_gap_vol = 3
        disposal_height = -1 # Arriba y el último paso lo hacemos dentro
        pickup_height = 1 

        for source, destination in zip(etoh_wells_multi, mag_wells_multi):
            run.pick_up()
            run.move_vol_multichannel(reagent=liquid, source=source,
                                      dest=destination, vol=175, air_gap_vol=air_gap_vol,
                                      pickup_height=pickup_height, disp_height=disposal_height
                                     )
            run.move_vol_multichannel(reagent=liquid, source=source,
                                      dest=destination, vol=175, air_gap_vol=air_gap_vol,
                                      pickup_height=pickup_height, disp_height=disposal_height)

            # This will be drop inside
            [disposal_height,column_change] = run.calc_height(liquid)
            run.move_vol_multichannel(reagent=liquid, source=source,
                                      dest=destination, vol=135, air_gap_vol=air_gap_vol,
                                      pickup_height=pickup_height, disp_height=disposal_height-3)
            
            run.custom_mix(liquid, location=source, vol=150,
                           rounds=10, blow_out=True, mix_height=0)

            run.drop_tip()

    ############################################################################
    # STEP 17: Magnet on 10 minutos
    ############################################################################
    if (run.next_step()):
        if (set_mag_on):
            magdeck.engage(height=mag_height)
        run.finish_step()
    ############################################################################
    # STEP 18: Extract liquid from magnet to liquid trash
    ############################################################################
    if (run.next_step()):
        move_magnet_to_trash()
        run.finish_step()

    ############################################################################
    # STEP 19: Magnet off
    ############################################################################
    if (run.next_step()):
        magdeck.disengage()
        run.finish_step()

    ############################################################################
    # STEP 20: Add 250 de etoh a los beats Slot 8 - 7
    ############################################################################
    if (run.next_step()):

        run.set_pip("left")  # p300 multi 
        liquid = etoh
        air_gap_vol = 3
        disposal_height = -1 # Arriba y el último paso lo hacemos dentro
        pickup_height = 1 

        for source, destination in zip(etoh_wells_multi, mag_wells_multi):
            run.pick_up()
            run.move_vol_multichannel(reagent=liquid, source=source,
                                      dest=destination, vol=175, air_gap_vol=air_gap_vol,
                                      pickup_height=pickup_height, disp_height=disposal_height
                                     )
            # This will be drop inside
            [disposal_height,column_change] = run.calc_height(liquid)
            run.move_vol_multichannel(reagent=liquid, source=source,
                                      dest=destination, vol=125, air_gap_vol=air_gap_vol,
                                      pickup_height=pickup_height, disp_height=disposal_height-3)
            
            run.drop_tip()

    ############################################################################
    # STEP 21: Magnet on 10 minutos
    ############################################################################
    if (run.next_step()):
        if (set_mag_on):
            magdeck.engage(height=mag_height)
        run.finish_step()

    ############################################################################
    # STEP 22: Extract liquid from magnet to liquid trash
    ############################################################################
    if (run.next_step()):
        move_magnet_to_trash(move_vol_steps=2)
        run.finish_step()

    ############################################################################
    # STEP 23: Magnet off
    ############################################################################
    if (run.next_step()):
        magdeck.disengage()
        run.finish_step()

    ############################################################################
    # STEP 24: Add elution move to temperature same tip 4 -> 7 -> 10
    ############################################################################
    #Used to move from temp to magnet and from magnet to destionation
    elu_beads = Reagent(name='BitsToHot',
                       flow_rate_aspirate=1,
                       flow_rate_dispense=1.5,
                       flow_rate_aspirate_mix=1.5,
                       flow_rate_dispense_mix=5,
                       air_gap_vol_bottom=5,
                       air_gap_vol_top=0,
                       disposal_volume=485,
                       rinse=True,
                       max_volume_allowed=500,#No aplica
                       reagent_volume=0, #No aplica
                       reagent_reservoir_volume=NUM_SAMPLES * 250 * 1.1,  #No aplica
                       num_wells=num_cols, # multi
                       h_cono=1.95,
                       v_fondo=750)

    if (run.next_step()):
        # to liquid types
        elution = Reagent(name='Elution',
                       flow_rate_aspirate=1,
                       flow_rate_dispense=1.5,
                       flow_rate_aspirate_mix=1.5,
                       flow_rate_dispense_mix=5,
                       air_gap_vol_bottom=5,
                       air_gap_vol_top=0,
                       disposal_volume=485,
                       rinse=True,
                       max_volume_allowed=500,#No aplica
                       reagent_volume=0, #No aplica
                       reagent_reservoir_volume=NUM_SAMPLES * 250 * 1.1,  #No aplica
                       num_wells=num_cols, # multi
                       h_cono=1.95,
                       v_fondo=750)

        for source,dest_source,destination in zip(wb_wells_multi,mag_wells_multi,temp_wells_multi):
            run.pick_up()
            run.move_vol_multichannel(reagent=elution, source=source,
                                      dest=dest_source, vol=50, air_gap_vol=air_gap_vol,
                                      pickup_height=pickup_height, disp_height=disposal_height
                                     )
            
            run.custom_mix(elu_beads, location=dest_source, vol=100,
                           rounds=10, blow_out=True, mix_height=0)

            # This will be drop inside
            [disposal_height,column_change] = run.calc_height(elu_beads)
            run.move_vol_multichannel(reagent=elu_beads, source=dest_source,
                                      dest=destination, vol=50, air_gap_vol=air_gap_vol,
                                      pickup_height=pickup_height, disp_height=disposal_height-3)            

            run.drop_tip()

        run.finish_step()

    ############################################################################
    # STEP 25: Incubation at 65ºC
    ############################################################################
    if (run.next_step()):
        if (set_temp_on):
            tempdeck.set_temperature(temperature)
        run.finish_step()

    ############################################################################
    # STEP 26: Move from temp to magnet
    ############################################################################
    if (run.next_step()):
        result = Reagent(name='Elution',
                       flow_rate_aspirate=1,
                       flow_rate_dispense=1.5,
                       flow_rate_aspirate_mix=1.5,
                       flow_rate_dispense_mix=5,
                       air_gap_vol_bottom=5,
                       air_gap_vol_top=0,
                       disposal_volume=485,
                       rinse=True,
                       max_volume_allowed=500,#No aplica
                       reagent_volume=0, #No aplica
                       reagent_reservoir_volume=NUM_SAMPLES * 250 * 1.1,  #No aplica
                       num_wells=num_cols, # multi
                       h_cono=1.95,
                       v_fondo=750)

        for source,destination in zip(mag_wells_multi,aw_wells_multi):
            run.pick_up()
            run.move_vol_multichannel(reagent=result, source=source,
                                      dest=dest_source, vol=50, air_gap_vol=air_gap_vol,
                                      pickup_height=pickup_height, disp_height=disposal_height
                                     )
            run.drop_tip()

        run.finish_step()

    ############################################################################
    # STEP 26: Magnet on 3 minutos
    ############################################################################
    if (run.next_step()):
        if (set_mag_on):
            magdeck.engage(height=mag_height)
        run.finish_step()

    ############################################################################
    # STEP 27: Move from magnet to final output slot 2
    ############################################################################
    if (run.next_step()):

        for source,destination in zip(temp_wells_multi,mag_wells_multi):
            run.pick_up()
            run.move_vol_multichannel(reagent=elu_beads, source=source,
                                      dest=dest_source, vol=50, air_gap_vol=air_gap_vol,
                                      pickup_height=pickup_height, disp_height=disposal_height
                                     )
            run.drop_tip()

        run.finish_step()

    run.log_steps_time()
    run.blink()
    ctx.comment('Finished! \nMove plate to PCR')
