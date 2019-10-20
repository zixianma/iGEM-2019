from opentrons import robot, containers, instruments

# containers
plate = containers.load('trough-12row-short', 'B1')
regPlate = containers.load('384-plate', 'D2')
tiprack = containers.load('tiprack-1000ul', 'D1')
tiprack2 = containers.load('tiprack-200ul', 'E2')

# pipettes
pipette1000 = instruments.Pipette(axis='b', max_volume=1000, tip_racks=[tiprack], aspirate_speed=600,
    dispense_speed=600)
pipette200 = instruments.Pipette(axis='a', max_volume=200, tip_racks=[tiprack2])

# commands
pipette1000.pick_up_tip(tiprack.wells('A2'))

pipette1000.mix(10, plate.wells('A1'))

#pipette200.pick_up_tip(tiprack2.wells('A1'))
pipette1000.aspirate(200, plate.wells('A1'))
pipette1000.move_to(plate.wells('A1').top(100))
pipette1000.move_to(plate.wells('A6'),strategy='direct')
pipette1000.dispense(200)
pipette1000.return_tip()

pipette1000.pick_up_tip(tiprack.wells('A1'))

for i in range (0, 1000):
    pipette1000.mix(10, plate.wells('A6'))
    pipette1000.delay(minutes=1)
