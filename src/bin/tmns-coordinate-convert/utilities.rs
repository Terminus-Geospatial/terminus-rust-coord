use tmns_coord::coord::coord_sys_base::CoordSysBase;
use tmns_coord::coord::coord_sys_def::CoordSysDefinition;
use tmns_coord::coord::coord_sys_factory::CoordSysFactory;

pub fn load_coordinate_system(epsg_code : i32 ) -> Box<dyn CoordSysBase>
{
    // Step 1:  Convert an EPSG code into a Coordinate System Definition.
    let cs_def: Option<CoordSysDefinition> = CoordSysDefinition::from_epsg( epsg_code );

    //     Log information about the definition
    match &cs_def {
        Some(x) => log::debug!("{}", x.to_log_string()),
        None => panic!("CS Was UNDEFINED!"),
    }

    // Step 2:  Convert the Coordinate System Definition into a real Coordinate-System object
    let cs_out : Option<Box<dyn CoordSysBase>> = match cs_def {
        None => panic!( "Destination Coordinate System was unable to be defined. Aborting!" ),
        Some(x) => CoordSysFactory::create( x ),
    };

    match cs_out {
        Some(x) => return x,
        None => panic!( "CS Out was undefined" ),
    };
}