from opentrons import robot, containers, instruments

#containers.create(
#    'cuvette_rack',                    # name of you container
#    grid=(7, 1),                    # specify amount of (columns, rows)
#    spacing=(13, 1),               # distances (mm) between each (column, row)
#    diameter=6.5,                     # diameter (mm) of each well on the plate
#    depth=43)

cuvette_rack = containers.load('point', 'A2')

tiprack = containers.load('tiprack-1000ul', 'E1')

pipette = instruments.Pipette(axis='a', max_volume=1000, tip_racks=[tiprack], aspirate_speed=600,
        dispense_speed=600)

pipette.move_to(cuvette_rack.wells('A1').top(130))
#pipette.pick_up_tip(tiprack.wells('A1'))
for i in range (0, 3000):
    pipette.mix(10, cuvette_rack.wells('A1'))
    pipette.delay(minutes=1)
