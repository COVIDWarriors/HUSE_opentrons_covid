#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 21:48:34 2020

@author: amorvan23a
"""

from opentrons import protocol_api,robot,labware,instruments

metadata = {'apiLevel': '2.2'}

def run(protocol: protocol_api.ProtocolContext):
    
    # Cargando tipracks de puntas en los slots definidos como segundo parámetro
    tiprack1 = protocol.load_labware('opentrons_96_tiprack_1000ul', 7)
    tiprack2 = protocol.load_labware('opentrons_96_tiprack_1000ul', 5)
    
    # cargando placas
    plate_96_1 = labware.load('96-flat', '1')
    plate_96_2 = labware.load('96-flat', '2')
    
     # cargando pipetas
    left = protocol.load_instrument('p20_multi_gen2', 'left',tip_racks=[tiprack1])
    right = protocol.load_instrument('p20_single_gen2', 'right',tip_racks=[tiprack2])
    
    # Cargando módulo de temperatura
    temp_mod = protocol.load_module('temperature module Gen2', '4')
    plate_temp = temp_mod.load_labware('opentrons_24_aluminumblock_generic_2ml_screwcap')
    # Establecer temperatura a n grados -> temp_mod.set_temperature(4)
    # Desactivar modulo de temperatura -> temp_mod.deactivate()
    
    
    # mezclando los reactivos del modulo de temperatura
    right.transfer(6875, plate_temp.wells('A1'), plate_temp.wells('D1'))
    right.transfer(1375, plate_temp.wells('A2'), plate_temp.wells('D1'))
    right.transfer(8.25, plate_temp.wells('A3'), plate_temp.wells('D1'))
    
    # Coger mezcla del pozillo de mezcla y dispensarlo en los 96 espacios del lab
    right.transfer(15, plate_temp.wells('D1'), plate_96_1.columns('1', to='12'))
    
    # coger del plato2 con la multi i dispensar en el plato1
    left.transfer(50, plate_96_2.columns('1', to='12'), plate_96_1.columns('1', to='12'))
    
    # Coger control positivo(slot4) con p20Single y mandarlo a pos 96 de slot 1
    right.transfer(100, plate_temp.wells('A6'), plate_96_1.wells('H2'))