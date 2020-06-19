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
    'protocolName': 'RNA Extraction PreKingFisher Version 2',
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

# No quitar es seguridad por control + o -
if(NUM_SAMPLES > 94):
    NUM_SAMPLES = 94

num_cols = math.ceil(NUM_SAMPLES/8)

# Usar control general para las esperas para debug, siempre True
use_waits = True

diameter_screwcap = 8.1  # Diameter of the screwcap
volume_cone = 57  # Volume in ul that fit in the screwcap cone
area_section_screwcap = (np.pi * diameter_screwcap**2) / 4
h_cone = (volume_cone * 3 / area_section_screwcap)


# Folder for the log files
log_folder = 'prekingfisher'

##################
# Custom function
##################
class Reagent:
    def __init__(self, name, flow_rate_aspirate, flow_rate_dispense,
                 reagent_reservoir_volume, num_wells, h_cono, v_fondo, rinse=False, delay=0,
                 tip_recycling='none', rinse_loops=3, flow_rate_dispense_mix=2, flow_rate_aspirate_mix=2):

        self.name = name
        self.flow_rate_aspirate = flow_rate_aspirate
        self.flow_rate_dispense = flow_rate_dispense
        self.flow_rate_aspirate_mix = flow_rate_aspirate_mix
        self.flow_rate_dispense_mix = flow_rate_dispense_mix
        self.rinse = bool(rinse)
        self.reagent_reservoir_volume = reagent_reservoir_volume
        self.delay = delay  # Delay of reagent in dispense
        self.num_wells = num_wells
        self.col = 0
        self.h_cono = h_cono
        self.v_cono = v_fondo
        self.unused = []
        self.tip_recycling = tip_recycling
        self.vol_well_original = reagent_reservoir_volume / num_wells
        self.vol_well = reagent_reservoir_volume / num_wells
        self.rinse_loops = rinse_loops
    
    def get_current_position(self):
        return self.reagent_reservoir[self.col]
    
    def set_positions(self,labware_address):
        self.reagent_reservoir = labware_address

    def next_column(self):
        # Move to next position inside reagent
        self.vol_well = self.vol_well_original
        self.col =self.col+1

    def calc_height(self, cross_section_area, aspirate_volume,
                    min_height=0.3, extra_volume=50):

        if self.vol_well < aspirate_volume + extra_volume:
            self.unused.append(self.vol_well)
            # column selector position; intialize to required number
            self.next_column() 
        
        height = (self.vol_well - aspirate_volume -
                    self.v_cono) / cross_section_area  # - reagent.h_cono
        self.vol_well = self.vol_well - aspirate_volume
        if height < min_height:
            height = min_height
        return height

    def divide_volume(self, volume, max_vol):

        num_transfers = math.ceil(volume/max_vol)
        vol_roundup = math.ceil(volume/num_transfers)
        last_vol = volume - vol_roundup*(num_transfers-1)
        vol_list = [vol_roundup for v in range(1, num_transfers)]
        vol_list.append(last_vol)
        return vol_list

class ProtocolRun:
    def __init__(self, ctx):
        self.ctx = ctx
        self.step_list = []
        self.step = 0

        # Folder and file_path for log time
        folder_path = '/var/lib/jupyter/notebooks/'+log_folder
        if not self.ctx.is_simulating():
            if not os.path.isdir(folder_path):
                os.mkdir(folder_path)
            self.file_path = folder_path + \
                '/rna_extraction_%s.tsv' % datetime.now().strftime("%d_%m_%Y_%H_%M_%S")

        self.selected_pip = "right"
        self.pips = {"right": {}, "left": {}}

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


        self.comment("###############################################")
        self.comment("You are about to run %s samples" % (NUM_SAMPLES))
        for step in self.step_list:
            if(step['execute']):
                self.comment(step["description"])
        self.blink(5)
        self.pause("Are you sure the set up is correct? \n Check the desk before continue\n press resume")
        self.comment("###############################################")

    def set_execution_step(self, index, value):
        self.step_list[index]["execute"] = value

    def get_current_step(self):
        return self.step_list[self.step]

    def next_step(self):
        if self.step_list[self.step]['execute'] == False:
            self.step += 1
            return False

        self.comment(self.step_list[self.step]['description'], add_hash=True)
        self.start = datetime.now()
        return True

    def finish_step(self):
        if (self.get_current_step()["wait_time"] > 0 and use_waits):
            self.cdelay(seconds=int(self.get_current_step()[
                "wait_time"]), msg=self.get_current_step()["description"])
        if (self.get_current_step()["wait_time"] > 0 and not use_waits):
            self.ccomment("We simulate a wait of:%s seconds" %
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
                    row = ('{}\t{}\t{}\t{}').format(
                        step["execute"], step["description"], step["wait_time"], step["execution_time"])
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

    def reset_pip_count(self,pip):
        
        pip.reset_tipracks()
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
                   post_dispense=False, post_dispense_vol=10, x_offset=[0, 0]):
        '''
        Function for mixing a given [vol] in the same [location] a x number of [rounds].
        blow_out: Blow out optional [True,False]
        x_offset = [source, destination]
        source_height: height from bottom to aspirate
        mix_height: height from bottom to dispense
        '''
        pip = self.get_current_pip()
        vol = vol-1
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

    def pick_up(self, multi=None):
        pip = self.get_current_pip()
        self.multi = multi
        if not self.ctx.is_simulating():
            if self.get_pip_count() == self.get_pip_maxes():
                self.ctx.pause('Replace ' + str(pip.max_volume) + 'µl tipracks before \
                resuming.')
                self.reset_pip_count(pip)

        if multi != None:
            pip.pick_up_tip(self.multi)
        else:
            if not pip.hw_pipette['has_tip']:
                self.add_pip_count()
                pip.pick_up_tip()

    def drop_tip(self):
        pip = self.get_current_pip()
        pip.drop_tip(home_after=False)
        self.add_pip_count()

    def change_tip(self):
        self.drop_tip()
        self.pick_up()

    def comment(self, comment, add_hash=False):
        hash_string = "#######################################################"
        if (add_hash):
            self.ctx.comment(hash_string)
        self.ctx.comment(('{}').format(comment))
        if (add_hash):
            self.ctx.comment(hash_string)

        if self.ctx.is_simulating():
            if (add_hash):
                print(hash_string)
            print(comment)
            if (add_hash):
                print(hash_string)

    def pause(self, comment):
        self.ctx.pause(comment)
        if self.ctx.is_simulating():
            print("%s\n Press any key to continue " % comment)

    def move_volume(self, reagent, source, dest, vol, air_gap_vol,
                    pickup_height, disp_height, blow_out, touch_tip=False, rinse=False,
                    post_dispense=False, post_dispense_vol=20,
                    post_airgap=True, post_airgap_vol=10, x_offset=[0, 0]):
        # x_offset: list with two values. x_offset in source and x_offset in destination i.e. [-1,1]
        # pickup_height: height from bottom where volume
        # rinse: if True it will do 2 rounds of aspirate and dispense before the tranfer
        # disp_height: dispense height; by default it's close to the top (z=-2), but in case it is needed it can be lowered
        # blow_out, touch_tip: if True they will be done after dispensing

        # Rinse before aspirating
        pipet = self.get_current_pip()
        if rinse == True:
            self.custom_mix(reagent, location=source, vol=vol,
                            rounds=reagent.rinse_loops, blow_out=True, mix_height=0,
                            x_offset=x_offset)
        # SOURCE
        s = source.bottom(pickup_height).move(Point(x=x_offset[0]))
        # aspirate liquid
        pipet.aspirate(vol, s, rate=reagent.flow_rate_aspirate)
        if air_gap_vol != 0:  # If there is air_gap_vol, switch pipette to slow speed
            pipet.aspirate(air_gap_vol, source.top(z=-2),
                           rate=reagent.flow_rate_aspirate)  # air gap
        # GO TO DESTINATION
        drop = dest.top(z=disp_height).move(Point(x=x_offset[1]))
        pipet.dispense(vol + air_gap_vol, drop,
                       rate=reagent.flow_rate_dispense)  # dispense all
        # pause for x seconds depending on reagent
        self.ctx.delay(seconds=reagent.delay)
        if blow_out == True:
            pipet.blow_out(dest.top(z=-2))
        if post_airgap == True:
            pipet.dispense(post_airgap_vol, dest.top(z=-2))
        if post_dispense == True:
            pipet.dispense(post_dispense_vol, dest.top(z=-2))
        if touch_tip == True:
            pipet.touch_tip(speed=20, v_offset=-5, radius=0.9)

    def start_lights(self):
        self.ctx._hw_manager.hardware.set_lights(
            rails=True)  # set lights off when using MMIX

    def stop_lights(self):
        self.ctx._hw_manager.hardware.set_lights(
            rails=False)  # set lights off when using MMIX

    def blink(self, blink_number=3):
        for i in range(blink_number):
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
    

    # Define stesp
    run.add_step(
        description="Transfer PK A6 - To AW_PLATE Single Slot1 -> Slot2")  # 1
    run.add_step(description="Transfer MS2 B6 - To AW_PLATE Single 1->2")  # 2
    run.add_step(description="Transfer Beads 3 - 2 Multi and mix")  # 3
    run.add_step(description="Pause to replace tip racks")  # 4
    run.add_step(description="Slot 2 -> 1 Washing buffer to plate")  # 5
    run.add_step(description="Slot 2 -> 3 elution buffer to plate")  # 6
    run.add_step(description="Slot 10 -> 11 etoh to plate")  # 7

    # execute avaliaible steps
    run.init_steps(steps)

    ##################################

    # Tube rack
    tube_rack = ctx.load_labware(
        'opentrons_24_tuberack_nest_1.5ml_screwcap', 4)
        
    # Destination plate SLOT 2
    if(ctx.is_simulating()):
        labware_type = 'opentrons_96_aluminumblock_generic_pcr_strip_200ul'    
    else:
        labware_type = 'axygen_96_wellplate_2000ul'

    aw_slot = ctx.load_labware(labware_type, 5)
    aw_wells = aw_slot.wells()[:NUM_SAMPLES]
    aw_wells_multi = aw_slot.rows()[0][:num_cols]

    # Magnetic Beads Pool
    beads_slot = ctx.load_labware(
        'nest_12_reservoir_15ml', 6)
    beads_wells_multi = beads_slot.rows()[0][:num_cols]

    # Magnetic Beads Pool
    wbeb_slot = ctx.load_labware(
        'nest_12_reservoir_15ml', 2)

    etoh_pool_slot = ctx.load_labware(
        'nest_1_reservoir_195ml', 10)
    etoh_pool_wells_multi = etoh_pool_slot.rows()[0][:num_cols]

    wb_slot = ctx.load_labware(labware_type, 1)
    wb_wells_multi = wb_slot.rows()[0][:num_cols]

    eb_slot = ctx.load_labware(labware_type, 3)
    eb_wells_multi = eb_slot.rows()[0][:num_cols]

    etoh_slot = ctx.load_labware(labware_type, 11)
    etoh_wells_multi = etoh_slot.rows()[0][:num_cols]

    # Mount pippets and set racks
    # Tipracks20_multi
    tips20 = ctx.load_labware('opentrons_96_tiprack_20ul', 7)
    tips300_1 = ctx.load_labware('opentrons_96_filtertiprack_200ul', 8)
    tips300_2 = ctx.load_labware('opentrons_96_filtertiprack_200ul', 9)

    run.mount_right_pip('p20_single_gen2', tip_racks=[tips20], capacity=20)
    run.mount_left_pip('p300_multi_gen2', tip_racks=[
                       tips300_1, tips300_2], capacity=200, multi=True)

    ############################################################################
    # STEP 1: Transfer A6 - To AW_PLATE
    ############################################################################
    if (run.next_step()):
        run.set_pip("right")  # single 20
        volumen_move = 5
        source = tube_rack.wells("A6")[0]
        liquid = Reagent(name='Proteinasa K',
                         num_wells=1,  # change with num samples
                         flow_rate_aspirate=0.75,  # Original 0.5
                         flow_rate_dispense=3,  # Original 1
                         reagent_reservoir_volume=528,
                         h_cono=4,
                         v_fondo=4 * math.pi * 4 ** 3 / 3
                         )

        run.pick_up()
        for dest in aw_wells:
            pickup_height = liquid.calc_height(
                4.12*4.12*math.pi, volumen_move)
            run.move_volume(reagent=liquid, source=source,
                            dest=dest, vol=volumen_move, air_gap_vol=1,
                            pickup_height=pickup_height, disp_height=-10,
                            blow_out=True, post_dispense=True, post_dispense_vol=5)

        run.drop_tip()
        run.finish_step()

    ############################################################################
    # STEP 2: Transfer B6 MS2 - To AW_PLATE
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
                         reagent_reservoir_volume=528,
                         h_cono=4,
                         v_fondo=4 * math.pi * 4 ** 3 / 3
                         )
        run.pick_up()
        for dest in aw_wells:

            pickup_height = liquid.calc_height(
                4.12*4.12*math.pi, volumen_move)
            run.move_volume(reagent=liquid, source=source,
                            dest=dest, vol=volumen_move, air_gap_vol=1,
                            pickup_height=pickup_height, disp_height=-10,
                            blow_out=True, post_dispense=True, post_dispense_vol=5)

        run.drop_tip()

        run.finish_step()

    ############################################################################
    # STEP 3: Slot 3 -2 Beads_PK AW
    ############################################################################
    if (run.next_step()):
        ############################################################################
        # Light flash end of program
        run.set_pip("left")  # p300 multi
        volume = 275
        beads = Reagent(name='Magnetic beads',
                        flow_rate_aspirate=0.5,
                        flow_rate_dispense=0.5,
                        flow_rate_dispense_mix=4,
                        flow_rate_aspirate_mix=4,
                        rinse=True,
                        delay=2,
                        reagent_reservoir_volume=30000,
                        num_wells=3,
                        h_cono=1.95,
                        v_fondo=695,
                        rinse_loops=3)

        #First 3 rows in this case
        beads.set_positions(beads_slot.rows()[0][0:3])

        air_gap_vol = 5
        disposal_height = -5
        pickup_height = 1
        
        pool_area = 8.3*71.1

        for destination in aw_wells_multi:
            run.pick_up()
            vol = 150
            vol_min = 1000

            pickup_height = beads.calc_height(pool_area, vol*8, extra_volume=vol_min)
            run.move_volume(reagent=beads, source=beads.get_current_position(),
                            dest=destination, vol=vol, air_gap_vol=air_gap_vol,
                            pickup_height=pickup_height, disp_height=disposal_height,
                            rinse=True, blow_out=True)
            run.change_tip()
            vol = 125
            pickup_height = beads.calc_height(pool_area, vol*8, extra_volume=vol_min)
            run.move_volume(reagent=beads, source=beads.get_current_position(),
                            dest=destination, vol=vol, air_gap_vol=air_gap_vol,
                            pickup_height=pickup_height, disp_height=disposal_height,
                            rinse=True, blow_out=True)

            run.custom_mix(beads, location=destination, vol=150,
                           rounds=3, blow_out=True, mix_height=0)

            run.drop_tip()

        run.finish_step()

    ############################################################################
    # STEP 4: Pause, get DW replace tip racks
    ############################################################################
    if (run.next_step()):
        ############################################################################
        # Light flash end of program
        run.blink(5)
        run.pause("Get deepwell PK+MS2+Beads and replace tiprack")
        run.set_pip("left")
        run.reset_pip_count(run.get_current_pip())
        

        
    ############################################################################
    # STEP 5: Slot 2 -> 1 Washing buffer to plate
    ############################################################################
    if (run.next_step()):
        ############################################################################
        # Light flash end of program
        run.set_pip("left")  # p300 multi
        volume = 500
        wb = Reagent(name='WB Wash buffer',
                        flow_rate_aspirate=0.5,
                        flow_rate_dispense=0.5,
                        flow_rate_dispense_mix=4,
                        flow_rate_aspirate_mix=4,
                        rinse=True,
                        delay=2,
                        reagent_reservoir_volume=60000,
                        num_wells=5,
                        h_cono=1.95,
                        v_fondo=695,
                        rinse_loops=3)

        air_gap_vol = 5
        disposal_height = -5
        wb.set_positions(wbeb_slot.rows()[0][0:5])
        pool_area = 8.3*71.1

        run.pick_up()
        for destination in wb_wells_multi:            
            vol_min = 10
            for vol in wb.divide_volume(volume,180):
                pickup_height= wb.calc_height(
                    pool_area, vol*8, extra_volume=vol_min)

                run.move_volume(reagent=wb, source=wb.get_current_position(),
                                dest=destination, vol=vol, air_gap_vol=air_gap_vol,
                                pickup_height=pickup_height, disp_height=disposal_height,
                                rinse=True, blow_out=True)
            
        run.drop_tip()

        run.finish_step()
    
    ############################################################################
    # STEP 6: Slot 2 -> 3 elution buffer to plate
    ############################################################################
    if (run.next_step()):
        ############################################################################
        # Light flash end of program
        run.set_pip("left")  # p300 multi
        volume = 50
        elution = Reagent(name='Elution Buffer',
                        flow_rate_aspirate=0.5,
                        flow_rate_dispense=0.5,
                        flow_rate_dispense_mix=4,
                        flow_rate_aspirate_mix=4,
                        reagent_reservoir_volume=5000,
                        num_wells=1,
                        h_cono=1.95,
                        v_fondo=695,
                        rinse_loops=3)

        air_gap_vol = 1
        disposal_height = -5
        elution.set_positions(wbeb_slot.rows()[0][11:12])
        pool_area = 8.3*71.1

        run.pick_up()
        for destination in wb_wells_multi:            
            run.move_volume(reagent=elution, source=elution.get_current_position(),
                                dest=destination, vol=vol, air_gap_vol=air_gap_vol,
                                pickup_height=pickup_height, disp_height=disposal_height,
                                rinse=True, blow_out=True)
            
        run.drop_tip()

        run.finish_step()
    
    # ############################################################################
    # # STEP 7: Slot 10 -> 11 etoh to plate
    # ############################################################################
    if (run.next_step()):
        ############################################################################
        # Light flash end of program
        run.set_pip("left")  # p300 multi
        volume = 1000
        etoh = Reagent(name='ETOH',
                        flow_rate_aspirate=0.5,
                        flow_rate_dispense=0.5,
                        flow_rate_dispense_mix=4,
                        flow_rate_aspirate_mix=4,
                        reagent_reservoir_volume=30000,
                        num_wells=12,
                        h_cono=1.95,
                        v_fondo=695)

        air_gap_vol = 5
        disposal_height = -5
        
        pool_area = 8.3*71.1
        pickup_height = 1

        run.pick_up()
        for source,destination in zip(etoh_pool_wells_multi,etoh_wells_multi):            
            vol_min = 1000
            volume_list = etoh.divide_volume(volume,175)
            for vol in volume_list:
                run.move_volume(reagent=etoh, source=source,
                                dest=destination, vol=vol, air_gap_vol=air_gap_vol,
                                pickup_height=1, disp_height=disposal_height,
                                rinse=True, blow_out=True)
            
        run.drop_tip()

        run.finish_step()
    
    run.log_steps_time()
    run.blink()
    for c in robot.commands():
        ctx.comment(c)
    ctx.comment('Finished! \nMove plate to PCR')
