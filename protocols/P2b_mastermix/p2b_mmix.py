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
# Define boolean variable to ask if want deactivate termoblock after step 3
##################
remove_termoblock = False
stop_termoblock = True

# Check stop termoblock when remove termoblock
if remove_termoblock == True:
    stop_termoblock == True

# Defined variables
##################
NUM_SAMPLES = 94
VOL_SAMPLE = 200
steps = []  # Steps you want to execute

# Tune variables
volume_elution = 10  # Volume of the sample
elution_initial_volume = 50  # True


temp = 10  # Define termoblock temperature
num_blinks = 3  # Define number of advisor temperature blinks

log_folder = 'p2b_mmix'

# Correct num samples
if NUM_SAMPLES >= 95:
    NUM_SAMPLES = 94


air_gap_vol = 10
air_gap_mmix = 0
air_gap_sample = 0

extra_dispensal = 0  # Extra volume for master mix in each distribute transfer
diameter_screwcap = 8.1  # Diameter of the screwcap
volume_cone = 57  # Volume in ul that fit in the screwcap cone
area_section_screwcap = (np.pi * diameter_screwcap**2) / 4
h_cone = (volume_cone * 3 / area_section_screwcap)

num_cols = math.ceil(NUM_SAMPLES/8)

def run(ctx: protocol_api.ProtocolContext):

    # Init protocol run
    run = ProtocolRun(ctx)
    run.add_step(description="TRANSFER Samples")
    run.init_steps(steps)

    ##################################
    # Define desk
    tempdeck = ctx.load_module('tempdeck', '7')

    # PCR
    pcr_plate = tempdeck.load_labware(
        'opentrons_96_aluminumblock_generic_pcr_strip_200ul')
    # setup up sample sources and destinations
    pcr_wells_multi = pcr_plate.rows()[0][:num_cols]
    
    # Eluted from King fisher/ Manual / Other
    try:
        elution_plate = ctx.load_labware(
            'axygen_96_wellplate_2000ul', '5')
    except:
        elution_plate = ctx.load_labware(
            'opentrons_96_aluminumblock_generic_pcr_strip_200ul', '5')

    elution_wells_multi = elution_plate.rows()[0][:num_cols]

    # Tipracks20_multi
    tips20 = ctx.load_labware('opentrons_96_tiprack_20ul', 9)
    
    # Mount pippets and set racks
    run.mount_left_pip('p20_multi_gen2', tip_racks=[tips20], capacity=20)
    run.mount_right_pip('p20_single_gen2', tip_racks=[tips20], capacity=20)
   

    # check temperature to know if the protocol can start
    tempdeck.set_temperature(temp)
    if tempdeck.temperature == temp: run.blink(blink_number=num_blinks)


    ############################################################################
    # STEP 1: TRANSFER Samples
    ############################################################################
    if(run.next_step()):
        run.comment('pcr_wells')
        run.set_pip("left")

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


        for s, d in zip(elution_wells_multi, pcr_wells_multi):
            run.comment("%s %s" % (s, d))
            run.pick_up()
            # Source samples
            run.move_volume(reagent=elution_well, source=s, dest=d,
                            vol=volume_elution, air_gap_vol=air_gap_sample,
                            pickup_height=0, disp_height=-10,
                            blow_out=False)
            run.custom_mix(reagent=elution_well, location=d, vol=8, rounds=3,
                           blow_out=False, mix_height=2, touch_tip=True)
            # ADD Custom mix
            run.drop_tip()

        if NUM_SAMPLES <= 88:
            run.set_pip("right")
            run.pick_up(tips20['H12'])
            source = elution_plate.wells("G12")[0]
            destination = pcr_plate.wells("G12")[0]

            run.move_volume(reagent=elution_well, source=source, dest=destination,
                            vol=volume_elution, air_gap_vol=air_gap_sample,
                            pickup_height=0, disp_height=-10,
                            blow_out=False)
            run.custom_mix(reagent=elution_well, location=destination, vol=8, rounds=3,
                           blow_out=False, mix_height=2, touch_tip=True)

            run.drop_tip()

        run.finish_step()
        tempdeck.deactivate()

    ############################################################################
    # Light flash end of program
    run.log_steps_time()
    run.blink(blink_number=num_blinks)
    run.comment('Finished! \nMove plate to PCR')




##################
# Custom function
##################
class Reagent:
    def __init__(self, name, flow_rate_aspirate, flow_rate_dispense,
                 reagent_reservoir_volume,  h_cono, v_fondo, vol_well_max=12000, num_wells=-1,rinse=False, delay=0,
                 tip_recycling='none', rinse_loops=3, flow_rate_dispense_mix=2, flow_rate_aspirate_mix=2):

        self.name = name
        self.flow_rate_aspirate = flow_rate_aspirate
        self.flow_rate_dispense = flow_rate_dispense
        self.flow_rate_aspirate_mix = flow_rate_aspirate_mix
        self.flow_rate_dispense_mix = flow_rate_dispense_mix
        self.rinse = bool(rinse)
        self.reagent_reservoir_volume = reagent_reservoir_volume
        self.delay = delay  # Delay of reagent in dispense
            
        self.col = 0
        self.h_cono = h_cono
        self.v_cono = v_fondo
        self.tip_recycling = tip_recycling
        self.rinse_loops = rinse_loops

        
        if(num_wells!=-1):
            if(num_wells==1):
                self.num_wells = num_wells
                self.vol_well = self.reagent_reservoir_volume
                self.vol_last_well = self.vol_well
                self.vol_well_max = vol_well_max

            else:   
                self.num_wells = num_wells
                #IF manually made we set up all to have the same
                self.vol_well_max = self.reagent_reservoir_volume/self.num_wells
                self.vol_last_well = self.vol_well_max
                self.vol_well = self.vol_last_well
        else:
            self.vol_well_max = vol_well_max-self.v_cono
            num_wells = math.floor(self.reagent_reservoir_volume/self.vol_well_max)
            self.vol_last_well = math.ceil(self.reagent_reservoir_volume-num_wells*self.vol_well_max)
            if(self.vol_last_well>0):
                self.num_wells = num_wells+1
            else:
                self.num_wells = num_wells
            
            if(self.num_wells==1):
                self.vol_well = self.vol_last_well
                self.vol_well_max = self.vol_last_well
            else:
                self.vol_well = math.ceil(self.vol_well_max)
    
    def get_current_position(self):
        
        return self.reagent_reservoir[self.col]
    
    def set_positions(self,labware_address):
        self.reagent_reservoir = labware_address

    def get_volumes_fill_print(self):
        if(self.num_wells==1):
            return "===> '%s' has %s wells with %s Volume"%(
                                            self.name,
                                            self.num_wells,
                                            self.vol_last_well+self.v_cono
                                            )

        else:
            return "===> '%s' has %s wells with %s Volume. Volumen_last_well: %s"%(self.name,
                                            self.num_wells-1,
                                            self.vol_well_max+self.v_cono,
                                            self.vol_last_well+self.v_cono)

    def next_column(self):
        # Move to next position inside reagent
        self.col =self.col+1
        if(self.col<self.num_wells):
            self.vol_well = self.vol_well_max
        else:
            self.vol_well = self.vol_last_well


    def calc_height(self, cross_section_area, aspirate_volume,
                    min_height=0.3):

        if self.vol_well < aspirate_volume:
            # column selector position; intialize to required number
            self.next_column() 
        
        height = (self.vol_well - aspirate_volume) / cross_section_area - 5
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
                   source_height=3, post_dispense=0, x_offset=[0, 0],touch_tip=False):
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
            z=source_height).move(Point(x=x_offset[0])), rate=reagent.flow_rate_aspirate_mix)
        for _ in range(rounds):
            pip.aspirate(vol, location=location.bottom(
                z=source_height).move(Point(x=x_offset[0])), rate=reagent.flow_rate_aspirate_mix)
            pip.dispense(vol, location=location.bottom(
                z=mix_height).move(Point(x=x_offset[1])), rate=reagent.flow_rate_dispense_mix)
        pip.dispense(1, location=location.bottom(
            z=mix_height).move(Point(x=x_offset[1])), rate=reagent.flow_rate_dispense_mix)
        if blow_out == True:
            pip.blow_out(location.top(z=-2))  # Blow out
        if post_dispense > 0:
            pip.dispense(post_dispense, location.top(z=-2))
        
        if touch_tip == True:
            pip = self.get_current_pip()
            pip.touch_tip(speed=20, v_offset=-5, radius=0.9)

    def pick_up(self, position=None):
        pip = self.get_current_pip()
        
        if not self.ctx.is_simulating():
            if self.get_pip_count() == self.get_pip_maxes():
                self.ctx.pause('Replace ' + str(pip.max_volume) + 'µl tipracks before \
                resuming.')
                self.reset_pip_count(pip)
        if position != None:
            pip.pick_up_tip(position)
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
        self.blink(3)
        if self.ctx.is_simulating():
            print("%s\n Press any key to continue " % comment)

    def move_volume(self, reagent, source, dest, vol, 
                    pickup_height, disp_height, air_gap_vol = 0,blow_out=False, touch_tip=False, rinse=False,
                    post_dispense=0,x_offset=[0, 0]):
        # x_offset: list with two values. x_offset in source and x_offset in destination i.e. [-1,1]
        # pickup_height: height from bottom where volume
        # rinse: if True it will do 2 rounds of aspirate and dispense before the tranfer
        # disp_height: dispense height; by default it's close to the top (z=-2), but in case it is needed it can be lowered
        # blow_out, touch_tip: if True they will be done after dispensing

        # Rinse before aspirating
        pipet = self.get_current_pip()
        if rinse == True:
            self.custom_mix(reagent, location=source, vol=vol,
                            rounds=reagent.rinse_loops, blow_out=True, mix_height=1, source_height=pickup_height,
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
        if post_dispense >0:
            pipet.dispense(post_dispense, dest.top(z=-2))
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

