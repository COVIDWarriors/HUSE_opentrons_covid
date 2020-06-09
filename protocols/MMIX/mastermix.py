#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from opentrons import protocol_api,robot,labware,instruments

metadata = {'apiLevel': '2.2'}

def run(protocol: protocol_api.ProtocolContext):

    muestras = 96
    reac1 = 6.875
    reac2 = 1.375
    reac3 = 8.25
    r1 = (reac1 * muestras) + 2
    r2 = (reac3 * muestras) + 2
    r3 = (reac3 * muestras) + 2
    
    # Cargando tipracks de puntas en los slots definidos como segundo parámetro
    tiprack1 = protocol.load_labware('opentrons_96_tiprack_20ul', 3)
    tiprack2 = protocol.load_labware('opentrons_96_tiprack_20ul', 7)
    # cargando pipetas
    left = protocol.load_instrument('p20_multi_gen2', 'left',tip_racks=[tiprack1])
    right = protocol.load_instrument('p20_single_gen2', 'right',tip_racks=[tiprack2])

    # cargando placas
    plate_96_1 = protocol.load_labware('biorad_96_wellplate_200ul_pcr', '1')

    plate_96_2 = protocol.load_labware('opentrons_96_aluminumblock_generic_pcr_strip_200ul','2')
    #plate_96_2 = protocol.load_labware('nest_96_wellplate_2ml_deep', '2')
    
     
    
    # Cargando módulo de temperatura
    temp_mod = protocol.load_module('tempdeck', '4')
    plate_temp = temp_mod.load_labware('opentrons_24_aluminumblock_generic_2ml_screwcap')
    # Establecer temperatura a n grados -> temp_mod.set_temperature(4)
    # Desactivar modulo de temperatura -> temp_mod.deactivate()
    
    
    # # mezclando los reactivos del modulo de temperatura
    right.transfer(r1, plate_temp.wells('A1'), plate_temp.wells('D1'))
    right.transfer(r2, plate_temp.wells('A2'), plate_temp.wells('D1'),mix_before=(2,20))
    right.transfer(r3, plate_temp.wells('A3'), plate_temp.wells('D1'),mix_before=(2,20))
    
    # # Coger mezcla del pozillo de mezcla y dispensarlo en los 96 espacios del lab
    # right.transfer(15, plate_temp.wells('D1'), plate_96_1.wells()[:])
    
    # coger del plato2 con la multi i dispensar en el plato1
    #left.transfer(10, plate_96_2.wells()[:], plate_96_1.wells()[:],new_tip='once',mix_before=(2,20),mix_after=(2,10))
    
    # Coger control positivo(slot4) con p20Single y mandarlo a pos 96 de slot 1
    # right.transfer(100, plate_temp.wells('A4'), plate_96_1.wells('H12'))