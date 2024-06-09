/**
 * @file    tmns-coordinate-convert.rs
 * @author  Marvin Smith
 * @date    5/24/2024
 */

mod config;
mod utilities;

use tmns_coord::coord::coord_sys_base::CoordSysBase;
use tmns_coord::coord::transformer::Transformer;

fn main() {

    env_logger::init();

    // Parse Command-Line Options
    let _config = config::parse_command_line();

    // Create Source Coordinate System
    let src_cs : Box<dyn CoordSysBase> = utilities::load_coordinate_system( 32613 );

    // Create Destination Coordinate System
    let dst_cs : Box<dyn CoordSysBase> = utilities::load_coordinate_system( 4326 );


    // Build transformer to execute conversion
    let xform = Transformer::create( src_cs, dst_cs, );

    let result = xform.apply( vec![500000.00, 4400343.22, 1581.0 ] );

    let val = result.unwrap();
    log::info!( "{}", format!( "Result: {}, X: {}, Y: {}, Z: {}", val.len(), val[0] * 180.0 / std::f64::consts::PI, val[1]* 180.0 / std::f64::consts::PI, val[2] ) );
}