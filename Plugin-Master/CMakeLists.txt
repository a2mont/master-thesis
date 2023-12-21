include (plugin)

openflipper_plugin(
        DIRS ./eigen3/
        DEPS OpenVolumeMesh
        TYPES TRIANGLEMESH TETRAHEDRALMESH)

string (REGEX MATCH "Plugin-.+[/\\]?$" _dir ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(${_dir} ortools)
