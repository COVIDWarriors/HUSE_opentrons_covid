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
import subprocess

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
NUM_SAMPLES = 96
steps = []  # Steps you want to execut
set_temp_on = False  # Do you want to start temperature module?
temperature = 65  # Set temperature. It will be uesed if set_temp_on is set to True
set_mag_on = False  # Do you want to start magnetic module?
mag_height = 14  # Height needed for NEST deepwell in magnetic deck

if(NUM_SAMPLES>94):
    NUM_SAMPLES=94;

use_waits = True
num_cols = math.ceil(NUM_SAMPLES/8)

diameter_screwcap = 8.1  # Diameter of the screwcap
volume_cone = 57  # Volume in ul that fit in the screwcap cone
area_section_screwcap = (np.pi * diameter_screwcap**2) / 4
h_cone = (volume_cone * 3 / area_section_screwcap)

air_gap_vol = 10
air_gap_r1 = 0
air_gap_sample = 0

#Folder for the log files
run_id = 'prekingfisher'

##################
# Custom function
##################
# Define Reagents as objects with their properties
class Reagent:
    def __init__(self, name, flow_rate_aspirate, reagent_volume, reagent_reservoir_volume,
                 flow_rate_dispense, flow_rate_aspirate_mix, flow_rate_dispense_mix,
                 air_gap_vol_bottom, air_gap_vol_top, disposal_volume,  max_volume_allowed,
                 num_wells, h_cono, v_fondo, area_section_screwcap=8.254, rinse=False, delay=0, tip_recycling='none'):
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

ctx = None

class ProtocolRun:
    def __init__(self):
        global ctx
        print(ctx)
        self.step_list = []
        self.step = 0

        # Folder and file_path for log time
        folder_path = '/var/lib/jupyter/notebooks/'+run_id
        if not ctx.is_simulating():
            if not os.path.isdir(folder_path):
                os.mkdir(folder_path)
            self.file_path = folder_path + '/KC_qPCR_time_log.txt'
        self.selected_pip = "right"
        self.pips = {"right": {}, "left": {}}

        # To delete if not working
        self.VOICE_FILES_DICT = {
            'start': '/var/lib/jupyter/notebooks/data/sounds/started_process.mp3',
            'finish': '/var/lib/jupyter/notebooks/data/sounds/finished_process.mp3',
            'close_door': '/var/lib/jupyter/notebooks/data/sounds/close_door.mp3',
            'replace_tipracks': '/var/lib/jupyter/notebooks/data/sounds/replace_tipracks.mp3',
            'empty_trash': '/var/lib/jupyter/notebooks/data/sounds/empty_trash.mp3'
        }

    def voice_notification(self, action):
        if not ctx.is_simulating():
            fname = self.VOICE_FILES_DICT[action]
            if os.path.isfile(fname) is True:
                subprocess.run(
                    ['mpg123', fname],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE
                )
            else:
                ctx.comment(
                    f"Sound file does not exist. Call the technician")

    def add_step(self, description, execute=False, wait_time=0):
        self.step_list.append(
            {'execute': execute, 'description': description, 'wait_time': wait_time, 'execution_time': 0})

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
        self.step_list[index]["execute"] = value

    def get_current_step(self):
        return self.step_list[self.step]

    def next_step(self):
        if self.step_list[self.step]['execute'] == False:
            self.step += 1
            return False
        ctx.comment(self.step_list[self.step]['description'])
        self.comment(self.step_list[self.step]['description'],add_hash=True)
        self.start = datetime.now()
        return True

    def finish_step(self):
        if (self.get_current_step()["wait_time"] > 0 and use_waits):
            ctx.delay(seconds=int(self.get_current_step()[
                "wait_time"]), msg=self.get_current_step()["description"])
        if (self.get_current_step()["wait_time"] > 0 and not use_waits):
            ctx.comment("We simulate a wait of:%s seconds" %
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
        if not ctx.is_simulating():
            with open(self.file_path, 'w') as f:
                f.write('STEP\texecution\tdescription\twait_time\texecution_time\n')
                row = ""
                for step in self.step_list:
                    row = ('{}\t{}\t{}\t{}\t{}').format(
                        row, step["execute"], step["description"], step["wait_time"], step["execution_time"])
                    f.write(row + '\n')
            f.close()

    def mount_pip(self, position, type, tip_racks, capacity, multi=False, size_tipracks=96):
        self.pips[position]["pip"] = ctx.load_instrument(
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
        if not ctx.is_simulating():
            if self.get_pip_count() == self.get_pip_maxes():
                ctx.pause('Replace ' + str(pip.max_volume) + 'µl tipracks before \
                    resuming.')
                pip.reset_tipracks()
                self.reset_pip_count()

        if not pip.hw_pipette['has_tip']:
            pip.pick_up_tip()

    def drop_tip(self):
        pip = self.get_current_pip()
        pip.drop_tip(home_after=False)
        self.add_pip_count()

    def change_tip(self):
        self.drop_tip()
        self.pick_up()

    def comment(self, comment, add_hash=False):
        hash_string = '#######################################################'
        if (add_hash):
            ctx.comment(hash_string)
        ctx.comment(('{}').format(comment))
        if (add_hash):
            ctx.comment(hash_string)

        if ctx.is_simulating():
            if (add_hash):
                print(hash_string)
            print(comment)
            if (add_hash):
                print(hash_string)

    def pause(self, comment):
        if not ctx.is_simulating():
            ctx.pause(comment)
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
            ctx.delay(seconds=reagent.delay)
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

        cross_section_area = reagent.area_section_screwcap
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


def run(context: protocol_api.ProtocolContext):
    global ctx
    ctx = context
    # Init protocol run
    run = ProtocolRun()
    # yo creo que este tiene que ser manual o sacarlo a otro robot
    run.add_step(
        description="Transfer PK A6 - To AW_PLATE Single Slot1 -> Slot2")  # 1
    run.add_step(description="Wait until bell is done")  # INTERACTION 2
    run.add_step(description="Transfer Beats 3 - 2 Multi")  # 3
    run.add_step(description="Transfer MS2 B6 - To AW_PLATE Single 1->2")  # 4

    # execute avaliaible steps
    run.init_steps(steps)
    run.voice_notification('start')
    ##################################

    # Tube rack
    tube_rack = ctx.load_labware(
        'opentrons_24_tuberack_nest_1.5ml_screwcap', 1)

    # Destination plate SLOT 2
    aw_slot = ctx.load_labware(
            'axygen_96_wellplate_2000ul', 2)

    aw_wells = aw_slot.wells()[:NUM_SAMPLES]
    aw_wells_multi = aw_slot.rows()[0][:num_cols]

    # Magnetic Beads Pool
    beads_slot = ctx.load_labware(
        'nest_12_reservoir_15ml', 3)
    beads_wells_multi = beads_slot.rows()[0][:num_cols]

    # Mount pippets and set racks
    # Tipracks20_multi
    tips20 = ctx.load_labware('opentrons_96_tiprack_20ul', 4)
    tips300 =  [ctx.load_labware('opentrons_96_filtertiprack_200ul', slot) for slot in ['5','6','8']]
    
    run.mount_right_pip('p20_single_gen2', tip_racks=[tips20], capacity=20)
    run.mount_left_pip('p300_multi_gen2', tip_racks=
                       tips300, capacity=200, multi=True)

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
                         reagent_volume=volumen_move*NUM_SAMPLES*1.1,
                         max_volume_allowed=1500,
                         reagent_reservoir_volume=volumen_move*NUM_SAMPLES*1.1,  
                         h_cono=4,
                         v_fondo=4 * math.pi * 4 ** 3 / 3,
                         area_section_screwcap=(np.pi * 8.25**2) / 4
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
    # STEP 3: Slot 3 -2 beats_PK AW
    ############################################################################
    if (run.next_step()):
        ############################################################################
        # Light flash end of program
        run.set_pip("left")  # p300 multi
        volume = 275
        liquid = Reagent(name='Magnetic beads',
                         flow_rate_aspirate=1,
                         flow_rate_dispense=1.5,
                         flow_rate_aspirate_mix=1.5,
                         flow_rate_dispense_mix=5,
                         air_gap_vol_bottom=5,
                         air_gap_vol_top=0,
                         disposal_volume=volume,
                         rinse=True,
                         max_volume_allowed=15000,
                         reagent_volume=5000,
                         reagent_reservoir_volume=NUM_SAMPLES * volume * 1.1,
                         num_wells=num_cols,
                         h_cono=1.95,
                         v_fondo=750)

        air_gap_vol = 3
        disposal_height = -5
        pickup_height = 1

        for source, destination in zip(beads_wells_multi, aw_wells_multi):
            run.pick_up()
            run.move_vol_multichannel(reagent=liquid, source=source,
                                      dest=destination, vol=150, air_gap_vol=air_gap_vol,
                                      pickup_height=pickup_height, disp_height=disposal_height,
                                      rinse=True)
            run.change_tip()
            run.move_vol_multichannel(reagent=liquid, source=source,
                                      dest=destination, vol=125, air_gap_vol=air_gap_vol,
                                      pickup_height=pickup_height, disp_height=disposal_height, rinse=True)
            run.drop_tip()

        run.finish_step()

    ############################################################################
    # STEP 4: Transfer B6 MS2 - To AW_PLATE
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
                         reagent_volume=volumen_move*NUM_SAMPLES*1.1,
                         max_volume_allowed=150,
                         reagent_reservoir_volume=volumen_move*NUM_SAMPLES*1.1,  # 14800,
                         h_cono=4,
                         v_fondo=4 * math.pi * 4 ** 3 / 3,
                         area_section_screwcap=(np.pi * 8.25**2) / 4
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
    # STEP 5: Mix
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
                         area_section_screwcap=(np.pi * 8.25**2) / 4
                         )

        for source in aw_wells_multi:
            run.pick_up()
            run.custom_mix(liquid, location=source, vol=10,
                           rounds=2, blow_out=True, mix_height=0)
            run.drop_tip()

        run.blink()
        ctx.pause('Replace tips, empty trash, move Slot2 -> Slot 10')
        run.finish_step()

    run.log_steps_time()
    run.blink()
    for c in robot.commands():
        ctx.comment(c)
    ctx.comment('Finished! \nMove plate to PCR')