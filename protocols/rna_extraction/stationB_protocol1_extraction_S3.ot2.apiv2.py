from opentrons import protocol_api
from opentrons.types import Point
from opentrons.drivers.rpi_drivers import gpio
import time
import math
import os
import sys
import subprocess
import json
from datetime import datetime
custom_modules_path = "/var/user-packages/usr/lib/python3.7/site-packages"
if custom_modules_path not in sys.path:
    sys.path.append(custom_modules_path)
import requests

# metadata
metadata = {
    'protocolName': 'S3 Station B Protocol 1 extraction Version 2',
    'author': 'Nick <protocols@opentrons.com> Sara <smonzon@isciii.es> Miguel <mjuliam@isciii.es>',
    'source': 'Custom Protocol Request',
    'apiLevel': '2.3'
}

"""
REAGENT SETUP:
- slot 7 12-channel reservoir:
    - elution buffer: channel 1
    - magnetic beads: channel 2
    - bead buffer: channels 3-5
    - wash 1: channels 6-7
    - wash 2: channels 8-9
    - wash 3: channels 10-11
- slot 11 single-channel reservoir:
    - empty reservoir for liquid waste (supernatant removals)
"""

# Parameters to adapt the protocol
# Warning writing any Parameters below this line.
# It will be deleted if opentronsWeb is used.

NUM_SAMPLES = 96
REAGENT_LABWARE = 'nest 12 reservoir plate'
MAGPLATE_LABWARE = 'nest deep generic well plate'
WASTE_LABWARE = 'nest 1 reservoir plate'
ELUTION_LABWARE = 'opentrons aluminum nest plate'
DISPENSE_BEADS = False
REUSE_TIPS = True
LANGUAGE = 'esp'
RESET_TIPCOUNT = False
PROTOCOL_ID = "0000-AA"
URL = 'localhost'
# End Parameters to adapt the protocol

## global vars
## initialize robot object
robot = None
# default var for drop tip switching
switch = True
# initialize tip_log dictionary
tip_log = {}
tip_log['count'] = {}
tip_log['tips'] = {}
tip_log['max'] = {}

"""
NUM_SAMPLES is the number of samples, must be an integer number

REAGENT_LABWARE must be one of the following:
    nest 12 reservoir plate

MAGPLATE_LABWARE must be one of the following:
    opentrons deep generic well plate
    nest deep generic well plate
    vwr deep generic well plate
    ecogen deep generic well plate

WASTE labware
    nest 1 reservoir plate

ELUTION_LABWARE
    opentrons aluminum biorad plate
    opentrons aluminum nest plate
"""

# Calculated variables
if MAGPLATE_LABWARE == 'nest deep generic well plate':
    MAGNET_HEIGHT = 22
elif MAGPLATE_LABWARE == 'vwr deep generic well plate':
    MAGNET_HEIGHT = 22
elif MAGPLATE_LABWARE == 'ecogen deep generic well plate':
    MAGNET_HEIGHT = 21
else:
    MAGNET_HEIGHT = 22

# Config variables
ACTION = "StationB-protocol1-extraction"

# Constants
REAGENT_LW_DICT = {
    'nest 12 reservoir plate': 'nest_12_reservoir_15ml'
}

MAGPLATE_LW_DICT = {
    'opentrons deep generic well plate': 'usascientific_96_wellplate_2.4ml_deep',
    'nest deep generic well plate': 'nest_96_deepwellplate_2000ul',
    'ecogen deep generic well plate': 'ecogen_96_deepwellplate_2000ul',
    'vwr deep generic well plate': 'vwr_96_deepwellplate_2000ul'
}

WASTE_LW_DICT = {
    # Radius of each possible tube
    'nest 1 reservoir plate': 'nest_1_reservoir_195ml'
}

ELUTION_LW_DICT = {
    'opentrons aluminum biorad plate': 'opentrons_96_aluminumblock_biorad_wellplate_200ul',
    'opentrons aluminum nest plate': 'opentrons_96_aluminumblock_nest_wellplate_100ul'

}

LANGUAGE_DICT = {
    'esp': 'esp',
    'eng': 'eng'
}

if LANGUAGE_DICT[LANGUAGE] == 'eng':
    VOICE_FILES_DICT = {
        'start': './data/sounds/started_process.mp3',
        'finish': './data/sounds/finished_process.mp3',
        'close_door': './data/sounds/close_door.mp3',
        'replace_tipracks': './data/sounds/replace_tipracks.mp3',
        'empty_trash': './data/sounds/empty_trash.mp3'
    }
elif LANGUAGE_DICT[LANGUAGE] == 'esp':
    VOICE_FILES_DICT = {
        'start': './data/sounds/started_process_esp.mp3',
        'finish': './data/sounds/finished_process_esp.mp3',
        'close_door': './data/sounds/close_door_esp.mp3',
        'replace_tipracks': './data/sounds/replace_tipracks_esp.mp3',
        'empty_trash': './data/sounds/empty_trash_esp.mp3'
    }

# Function definitions

def write_to_error_log (info, reason):
    date = datetime.now().strftime("%Y_%m_%d")
    folder_date = os.path.join('/data/logs', date)
    time_now = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    json_file = time_now + '.json'
    folder_file_name = os.path.join(folder_date, json_file)
    folder_error_log = os.path.join(folder_date,'error.log')
    if not os.path.exists(folder_date):
        try:
            os.makedirs(folder_date)
        except:
            return
    try:
        # Create a new file for dumping json data
        with open (folder_file_name , 'w') as fh:
            json.dump(info, fh, indent=4)
        # Append status reason code to the log
        with open(folder_error_log, 'a') as fh:
            fh.write( time_now +  '  Unable to accept the requests get error : '+ reason + '\n')
    except:
        return

def run_info(start, end, parameters = dict()):
    info = {}
    hostname = subprocess.run(
        ['hostname'],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    ).stdout.decode('utf-8')

    info["RobotID"] = hostname
    info["executedAction"] = ACTION
    info["ProtocolID"] = PROTOCOL_ID
    info["StartRunTime"] = start
    info["FinishRunTime"] = end
    info["parameters"] = parameters

    headers = {'Content-type': 'application/json'}
    url_https = 'https://' + URL
    url_http = 'http://' + URL

    if not robot.is_simulating():
        try:
            r = requests.post(url_https, data=json.dumps(info), headers=headers)
        except:
            try:
                r = requests.post(url_http, data=json.dumps(info), headers=headers)
            except:
                write_to_error_log(info, 'Server communication error')
                return
        if r.status_code > 201 :
            write_to_error_log(info, str(r.status_code))

def check_door():
    return gpio.read_window_switches()

def confirm_door_is_closed():
    if not robot.is_simulating():
        #Check if door is opened
        if check_door() == False:
            #Set light color to red and pause
            gpio.set_button_light(1,0,0)
            robot.pause()
            voice_notification('close_door')
            time.sleep(5)
            confirm_door_is_closed()
        else:
            #Set light color to green
            gpio.set_button_light(0,1,0)

def start_run():
    voice_notification('start')
    gpio.set_button_light(0,1,0)
    now = datetime.now()
    # dd/mm/YY H:M:S
    start_time = now.strftime("%Y/%m/%d %H:%M:%S")
    return start_time

def finish_run():
    voice_notification('finish')
    #Set light color to blue
    gpio.set_button_light(0,0,1)
    now = datetime.now()
    # dd/mm/YY H:M:S
    finish_time = now.strftime("%Y/%m/%d %H:%M:%S")
    return finish_time

def voice_notification(action):
    if not robot.is_simulating():
        fname = VOICE_FILES_DICT[action]
        if os.path.isfile(fname) is True:
                subprocess.run(
                ['mpg123', fname],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
                )
        else:
            robot.comment(f"Sound file does not exist. Call the technician")

def reset_tipcount(file_path = '/data/B/tip_log.json'):
    if os.path.isfile(file_path):
        os.remove(file_path)

def retrieve_tip_info(pip,tipracks,file_path = '/data/B/tip_log.json'):
    global tip_log
    if not tip_log['count'] or pip not in tip_log['count']:
        tip_log['count'][pip] = 0
        if not robot.is_simulating():
            if os.path.isfile(file_path):
                with open(file_path) as json_file:
                    data = json.load(json_file)
                    if 'P1000' in str(pip):
                        tip_log['count'][pip] = data['tips1000']
                    elif 'P300' in str(pip) and 'Single-Channel' in str(pip):
                        tip_log['count'][pip] = data['tips300']
                    elif 'P300' in str(pip) and '8-Channel' in str(pip):
                        tip_log['count'][pip] = data['tipsm300']
                    elif 'P20' in str(pip) and 'Single-Channel' in str(pip):
                        tip_log['count'][pip] = data['tips20']
                    elif 'P20' in str(pip) and '8-Channel' in str(pip):
                        tip_log['count'][pip] = data['tipsm20']

        if "8-Channel" in str(pip):
            tip_log['tips'][pip] =  [tip for rack in tipracks for tip in rack.rows()[0]]
        else:
            tip_log['tips'][pip] = [tip for rack in tipracks for tip in rack.wells()]

        tip_log['max'][pip] = len(tip_log['tips'][pip])

    return tip_log

def save_tip_info(file_path = '/data/B/tip_log.json'):
    data = {}
    if not robot.is_simulating():
        if os.path.isfile(file_path):
            os.rename(file_path,file_path + ".bak")
        for pip in tip_log['count']:
            if "P1000" in str(pip):
                data['tips1000'] = tip_log['count'][pip]
            elif 'P300' in str(pip) and 'Single-Channel' in str(pip):
                data['tips300'] = tip_log['count'][pip]
            elif 'P300' in str(pip) and '8-Channel' in str(pip):
                data['tipsm300'] = tip_log['count'][pip]
            elif 'P20' in str(pip) and 'Single-Channel' in str(pip):
                data['tips20'] = tip_log['count'][pip]
            elif 'P20' in str(pip) and '8-Channel' in str(pip):
                data['tipsm20'] = tip_log['count'][pip]

        with open(file_path, 'a+') as outfile:
            json.dump(data, outfile)

def pick_up(pip,tiprack):
    if tip_log['count'][pip] == tip_log['max'][pip]:
        voice_notification('replace_tipracks')
        robot.pause('Replace ' + str(pip.max_volume) + 'µl tipracks before \
resuming.')
        confirm_door_is_closed()
        pip.reset_tipracks()
        tip_log['count'][pip] = 0
    pip.pick_up_tip(tip_log['tips'][pip][tip_log['count'][pip]])
    tip_log['count'][pip] += 1

def drop(pip):
    global switch
    if "8-Channel" not in str(pip):
        side = 1 if switch else -1
        drop_loc = robot.loaded_labwares[12].wells()[0].top().move(Point(x=side*20))
        pip.drop_tip(drop_loc,home_after=False)
        switch = not switch
    else:
        drop_loc = robot.loaded_labwares[12].wells()[0].top().move(Point(x=20))
        pip.drop_tip(drop_loc,home_after=False)

def mix_beads(reps, dests, pip, tiprack):
    ## Dispense beads to deep well plate.
    for i, m in enumerate(dests):
        if not pip.hw_pipette['has_tip']:
            pick_up(pip,tiprack)
        dispense_default_speed = pip.flow_rate.dispense
        pip.flow_rate.dispense = 600
        pip.mix(reps, 200, m.bottom(5))
        pip.flow_rate.dispense = dispense_default_speed
        # PENDING TO FIX THIS blow_out
        # pip.blow_out(m.top(-2))
        pip.aspirate(20, m.top(-2))
        drop(pip)

def dispense_beads(reps,sources,dests,pip,tiprack):
    ## Mix beads prior to dispensing.
    pick_up(pip,tiprack)
    for s in sources:
        pip.mix(reps, 200, s.bottom(20))

    ## Dispense beads to deep well plate.
    for i, m in enumerate(dests):
        if not pip.hw_pipette['has_tip']:
            pick_up(pip,tiprack)
        pip.transfer(200, dests[i//3], m.bottom(5), new_tip='never', air_gap=20)
        pip.blow_out(m.top(-2))
        drop(pip)
        pick_up(pip,tiprack)
        pip.transfer(200, dests[i//3], m.bottom(5), new_tip='never', air_gap=20)
        mix_beads(reps, dests, pip, tiprack)

def remove_supernatant(sources,waste,pip,tiprack):
    for i, m in enumerate(sources):
        loc = m.bottom(1.5)
        pick_up(pip,tiprack)
        pip.transfer(850, loc, waste, air_gap=100, new_tip='never')
        pip.blow_out(waste)
        drop(pip)

def wash_reuse(wash_sets,dests,waste,magdeck,pip,tiprack,tipreuse):
    wash_num = 0
    for wash_set in wash_sets:
        # transfer wash
        pick_up(pip,tiprack)
        for i, m in enumerate(dests):
            magdeck.disengage()
            wash_chan = wash_set[i//6]
            dispense_default_speed = pip.flow_rate.dispense
            pip.flow_rate.dispense = 200
            pip.transfer(200, wash_chan.bottom(2), m.top(), new_tip='never', air_gap=10)
            pip.flow_rate.dispense = dispense_default_speed

        # mix beads with wash
        tips_loc = 0
        for i, m in enumerate(dests):

            if tips_loc == 0:
                drop(pip)

            if  wash_num != 0:
                pip.pick_up_tip(tipreuse[0].rows()[0][tips_loc])
            else:
                pick_up(pip,tiprack)

            dispense_default_speed = pip.flow_rate.dispense
            pip.flow_rate.dispense = 1500
            pip.mix(7, 200, m.bottom(2))
            pip.flow_rate.dispense = dispense_default_speed
            if  wash_num != 0:
                pip.return_tip(home_after=False)
            else:
                pip.drop_tip(tipreuse[0].rows()[0][tips_loc], home_after=False)
            tips_loc += 1

        magdeck.engage(height_from_base=MAGNET_HEIGHT)
        robot.delay(seconds=75, msg='Incubating on magnet for 75 seconds.')

        wash_num += 1

        # remove supernatant
        tips_loc = 0
        for i, m in enumerate(dests):
            pip.pick_up_tip(tipreuse[0].rows()[0][tips_loc])
            aspire_default_speed = pip.flow_rate.aspirate
            pip.flow_rate.aspirate = 75
            asp_loc = m.bottom(1.5)
            pip.transfer(220, asp_loc, waste, new_tip='never', air_gap=20)
            pip.flow_rate.aspirate = aspire_default_speed
            pip.blow_out(waste)
            if  wash_num != 3:
                pip.return_tip(home_after=False)
            else:
                pip.drop_tip(home_after=False)
            tips_loc += 1

def wash(wash_sets,dests,waste,magdeck,pip,tiprack):
    for wash_set in wash_sets:
        for i, m in enumerate(dests):
            # transfer and mix wash with beads
            magdeck.disengage()
            wash_chan = wash_set[i//6]
            pick_up(pip,tiprack)
            pip.transfer(
                200, wash_chan.bottom(2), m.center(), new_tip='never', air_gap=20)
            # Mix heigh has to be really close to bottom, it was 5 now reduced to 2, maybe should be 1?
            dispense_default_speed = pip.flow_rate.dispense
            pip.flow_rate.dispense = 1500
            pip.mix(7, 200, m.bottom(2))
            pip.flow_rate.dispense = dispense_default_speed

            magdeck.engage(height_from_base=MAGNET_HEIGHT)
            robot.delay(seconds=75, msg='Incubating on magnet for 75 seconds.')

            # remove supernatant
            aspire_default_speed = pip.flow_rate.aspirate
            pip.flow_rate.aspirate = 75
            asp_loc = m.bottom(1.5)
            pip.transfer(220, asp_loc, waste, new_tip='never', air_gap=10)
            pip.flow_rate.aspirate = aspire_default_speed
            pip.blow_out(waste)
            drop(pip)

def elute_samples_reuse(sources,dests,buffer,magdeck,pip,tipracks,tipreuse):
    ## dispense buffer
    tips_loc = 0
    for i, m in enumerate(sources):
        pick_up(pip,tipracks)
        dispense_default_speed = pip.flow_rate.dispense
        pip.flow_rate.dispense = 1500
        pip.transfer(
            50, buffer.bottom(2), m.bottom(1), new_tip='never', air_gap=10)
        pip.mix(20, 200, m.bottom(1))
        pip.flow_rate.dispense = dispense_default_speed
        pip.drop_tip(tipreuse[0].rows()[0][tips_loc], home_after=False)
        tips_loc += 1

    ## Incubation steps
    robot.delay(minutes=5, msg='Incubating off magnet for 5 minutes.')
    magdeck.engage(height_from_base=MAGNET_HEIGHT)
    robot.delay(seconds=120, msg='Incubating on magnet for 120 seconds.')

    aspire_default_speed = pip.flow_rate.aspirate
    pip.flow_rate.aspirate = 50
    ## Dispense elutes in pcr plate.
    tips_loc = 0
    for i, (m, e) in enumerate(zip(sources, dests)):
        # tranfser and mix elution buffer with beads
        asp_loc = m.bottom(1.5)
        pip.pick_up_tip(tipreuse[0].rows()[0][tips_loc])
        # transfer elution to new plate
        pip.transfer(50, asp_loc, e, new_tip='never', air_gap=10)
        pip.blow_out(e.top(-2))
        drop(pip)
        tips_loc += 1
    pip.flow_rate.aspirate = aspire_default_speed

def elute_samples(sources,dests,buffer,magdeck,pip,tipracks):
    ## dispense buffer
    for i, m in enumerate(sources):
        pick_up(pip,tipracks)
        dispense_default_speed = pip.flow_rate.dispense
        pip.flow_rate.dispense = 1500
        pip.transfer(
            50, buffer.bottom(2), m.bottom(1), new_tip='never', air_gap=10)
        pip.mix(20, 200, m.bottom(1))
        pip.flow_rate.dispense = dispense_default_speed
        drop(pip)

    ## Incubation steps
    robot.delay(minutes=5, msg='Incubating off magnet for 5 minutes.')
    magdeck.engage(height_from_base=MAGNET_HEIGHT)
    robot.delay(seconds=120, msg='Incubating on magnet for 120 seconds.')

    aspire_default_speed = pip.flow_rate.aspirate
    pip.flow_rate.aspirate = 50
    ## Dispense elutes in pcr plate.
    for i, (m, e) in enumerate(zip(sources, dests)):
        # tranfser and mix elution buffer with beads
        asp_loc = m.bottom(1.5)
        pick_up(pip,tipracks)
        # transfer elution to new plate
        pip.transfer(50, asp_loc, e, new_tip='never', air_gap=10)
        pip.blow_out(e.top(-2))
        drop(pip)
    pip.flow_rate.aspirate = aspire_default_speed

def run(ctx: protocol_api.ProtocolContext):
    global robot
    robot = ctx

    # check if tipcount is being reset
    if RESET_TIPCOUNT:
        reset_tipcount()

    # confirm door is close
    robot.comment(f"Please, close the door")
    confirm_door_is_closed()

    # Begin run
    start_time = start_run()

    # load labware and modules
    ## ELUTION LABWARE
    if ELUTION_LABWARE not in ELUTION_LW_DICT:
        raise Exception('Invalid ELUTION_LABWARE. Must be one of the \
    following:\nopentrons aluminum biorad plate\nopentrons aluminum nest plate')

    elution_plate = robot.load_labware(
        ELUTION_LW_DICT[ELUTION_LABWARE], '1',
        'elution plate')

    ## MAGNETIC PLATE LABWARE
    magdeck = robot.load_module('magdeck', '10')
    magdeck.disengage()

    if MAGPLATE_LABWARE not in MAGPLATE_LW_DICT:
        raise Exception('Invalid MAGPLATE_LABWARE. Must be one of the \
following:\nopentrons deep generic well plate\nnest deep generic well plate\nvwr deep generic well plate')

    magplate = magdeck.load_labware(MAGPLATE_LW_DICT[MAGPLATE_LABWARE])

    ## WASTE LABWARE
    if WASTE_LABWARE not in WASTE_LW_DICT:
        raise Exception('Invalid WASTE_LABWARE. Must be one of the \
    following:\nnest 1 reservoir plate')

    waste = robot.load_labware(
        WASTE_LW_DICT[WASTE_LABWARE], '11', 'waste reservoir').wells()[0].top(-10)

    ## REAGENT RESERVOIR
    if REAGENT_LABWARE not in REAGENT_LW_DICT:
        raise Exception('Invalid REAGENT_LABWARE. Must be one of the \
    following:\nnest 12 reservoir plate')

    if REUSE_TIPS == True:
        reagent_res = robot.load_labware(
            REAGENT_LW_DICT[REAGENT_LABWARE], '4', 'reagent reservoir')
    else:
        reagent_res = robot.load_labware(
            REAGENT_LW_DICT[REAGENT_LABWARE], '7', 'reagent reservoir')

    ## TIPS
    # using standard tip definition despite actually using filter tips
    # so that the tips can accommodate ~220µl per transfer for efficiency
    if REUSE_TIPS == True:
        tips300 = [
            robot.load_labware(
                'opentrons_96_tiprack_300ul', slot, '200µl filter tiprack')
            for slot in ['8', '6', '2', '3']
        ]
        tipsreuse = [
            robot.load_labware(
                'opentrons_96_tiprack_300ul', slot, '200µl filter tiprack')
            for slot in ['7']
        ]
        tips1000 = [
            robot.load_labware('opentrons_96_filtertiprack_1000ul', slot,
                             '1000µl filter tiprack')
            for slot in ['5']
        ]
    else:
        tips300 = [
        robot.load_labware(
            'opentrons_96_tiprack_300ul', slot, '200µl filter tiprack')
            for slot in ['2', '3', '5', '6', '9','4']
        ]
        tips1000 = [
            robot.load_labware('opentrons_96_filtertiprack_1000ul', slot,
                             '1000µl filter tiprack')
            for slot in ['8']
        ]

    # reagents and samples
    num_cols = math.ceil(NUM_SAMPLES/8)
    mag_samples_m = magplate.rows()[0][:num_cols]
    mag_samples_s = magplate.wells()[:NUM_SAMPLES]
    elution_samples_m = elution_plate.rows()[0][:num_cols]
    elution_buffer = reagent_res.wells()[0]
    bead_buffer = reagent_res.wells()[1:5]
    wash_sets = [reagent_res.wells()[i:i+2] for i in [5, 7, 9]]

    # pipettes
    m300 = robot.load_instrument('p300_multi_gen2', 'left', tip_racks=tips300)
    p1000 = robot.load_instrument('p1000_single_gen2', 'right',
                                tip_racks=tips1000)

    ## retrieve tip_log
    retrieve_tip_info(p1000,tips1000)
    retrieve_tip_info(m300,tips300)

    m300.flow_rate.aspirate = 150
    m300.flow_rate.dispense = 300
    m300.flow_rate.blow_out = 300
    p1000.flow_rate.aspirate = 100
    p1000.flow_rate.dispense = 1000
    p1000.flow_rate.blow_out = 1000

    if(DISPENSE_BEADS):
        # premix, transfer, and mix magnetic beads with sample
        ## bead dests depending on number of samples
        bead_dests = bead_buffer[:math.ceil(num_cols/4)]
        dispense_beads(7,bead_dests,mag_samples_m,m300,tips300)
    else:
        # Mix bead
        mix_beads(7, mag_samples_m,m300,tips300)

    # incubate off the magnet
    robot.delay(minutes=10, msg='Incubating off magnet for 10 minutes.')

    ## First incubate on magnet.
    magdeck.engage(height_from_base=MAGNET_HEIGHT)
    robot.delay(minutes=7, msg='Incubating on magnet for 7 minutes.')

    # empty trash
    if NUM_SAMPLES >= 48:
        voice_notification('empty_trash')
        robot.pause(f"Please, empty trash")
        confirm_door_is_closed()

    # remove supernatant with P1000
    remove_supernatant(mag_samples_s,waste,p1000,tips1000)

    # empty trash
    if NUM_SAMPLES >= 48:
        voice_notification('empty_trash')
        robot.pause(f"Please, empty trash")
        confirm_door_is_closed()

    # 3x washes
    if REUSE_TIPS == True:
        wash_reuse(wash_sets,mag_samples_m,waste,magdeck,m300,tips300,tipsreuse)
    else:
        wash(wash_sets,mag_samples_m,waste,magdeck,m300,tips300)

    # elute samples
    magdeck.disengage()
    if REUSE_TIPS == True:
        elute_samples_reuse(mag_samples_m,elution_samples_m,elution_buffer,magdeck,m300,tips300,tipsreuse)
    else:
        elute_samples(mag_samples_m,elution_samples_m,elution_buffer,magdeck,m300,tips300)

    # track final used tip
    save_tip_info()

    magdeck.disengage()

    finish_time = finish_run()

    par = {
        "NUM_SAMPLES" : NUM_SAMPLES,
        "REAGENT_LABWARE" : REAGENT_LABWARE,
        "MAGPLATE_LABWARE" : MAGPLATE_LABWARE,
        "WASTE_LABWARE" : WASTE_LABWARE,
        "ELUTION_LABWARE" : ELUTION_LABWARE,
        "DISPENSE_BEADS" : DISPENSE_BEADS,
        "LANGUAGE" : LANGUAGE,
        "RESET_TIPCOUNT" : RESET_TIPCOUNT
    }

    run_info(start_time, finish_time, par)
