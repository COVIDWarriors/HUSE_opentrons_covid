import opentrons

metadata = {
    'protocolName': 'Thermocycler Example Protocol',
    'author': 'Opentrons <protocols@opentrons.com>',
    'source': 'Protocol Library',
    'apiLevel': '2.0'
}

def run(protocol):
	
	[STATIC_IP] = get_values( # noqa: F821
				'static_ip')

	opentrons.comment("Are you sure: %s "%STATIC_IP)