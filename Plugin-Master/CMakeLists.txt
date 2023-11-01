include (plugin)

openflipper_plugin(
        DIRS ./eigen3/
        DEPS OpenVolumeMesh
        DEPS CGAL
        TYPES TRIANGLEMESH TETRAHEDRALMESH)
