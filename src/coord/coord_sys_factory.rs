
use super::coord_sys_base::CoordSysBase;
use super::coord_sys_geographic::CoordSysGeographic;
use super::coord_sys_utm::CoordSysUTM;
use super::coord_sys_def::CoordSysDefinition;
use super::coordinate_type::CoordinateType;

pub struct CoordSysFactory {
}


impl CoordSysFactory {

    // This class takes a Coordinate-System Definition, then constructs a CoordSysBase-defined trait struct from it
    //
    // cs_def Coordinate System Definition with all required information to construct the appropriate CS type.
    pub fn create(cs_def: CoordSysDefinition) -> Option<Box<dyn CoordSysBase>> {

        match cs_def.coord_type {

            // Attempt to construct a Geographic coordinate system from the definition
            CoordinateType::GEOGRAPHIC => {

                return CoordSysGeographic::create( cs_def );

            },

            // Attempt to construct a UTM coordinate system from the definition
            CoordinateType::UTM => {

                return CoordSysUTM::create( cs_def );

            }
        }
    }
}