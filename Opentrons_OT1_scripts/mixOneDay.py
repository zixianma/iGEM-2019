from opentrons import robot, containers, instruments

# containers
plate = containers.load('trough-12row-short', 'B1')
tiprack = containers.load('tiprack-1000ul', 'E1')

# pipettes
pipette = instruments.Pipette(axis='b', max_volume=1000, tip_racks=[tiprack], aspirate_speed=600,
    dispense_speed=600)

# commands
pipette.pick_up_tip(tiprack.wells('A1'))
for i in range (0, 1000):
    pipette.mix(10, plate.wells('A1'))
    pipette.delay(minutes=1)
