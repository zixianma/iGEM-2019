from opentrons import robot, containers, instruments

# containers
plate = containers.load('trough-12row-short', 'B1')
regPlate = containers.load('384-plate', 'B1')
tiprack = containers.load('tiprack-1000ul', 'E1')
tiprack2 = containers.load('tiprack-200ul', 'E2')

# pipettes
pipette1000 = instruments.Pipette(axis='b', max_volume=1000, tip_racks=[tiprack], aspirate_speed=600,
    dispense_speed=600)
pipette200 = instruments.Pipette(axis='a', max_volume=200, tip_racks=[tiprack2])

# commands
pipette1000.pick_up_tip(tiprack.wells('A1'))
for i in range (0, 500):
    pipette1000.mix(10, plate.wells('A1'))
    pipette1000.delay(minutes=1)

pipette200.transfer(15, regPlate.wells('A1'), plate.wells('A1'))

for i in range (0, 200):
    pipette1000.mix(10, plate.wells('A1'))
    pipette1000.delay(minutes=1)

pipette1000.transfer(1000, plate.wells('A1'), plate.wells('A6'))

for i in range (0, 1000):
    pipette1000.mix(10, plate.wells('A6'))
    pipette1000.delay(minutes=1)
