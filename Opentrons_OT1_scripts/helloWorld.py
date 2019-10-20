from opentrons import robot, containers, instruments

# containers
plate = containers.load('384-plate', 'B1')
tiprack = containers.load('tiprack-200ul', 'E2')

# pipettes
pipette = instruments.Pipette(axis='a', max_volume=200, tip_racks=[tiprack])

# commands
pipette.transfer(50, plate.wells('A1'), plate.wells('B1'))
