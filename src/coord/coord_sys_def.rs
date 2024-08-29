/**************************** INTELLECTUAL PROPERTY RIGHTS ****************************/
/*                                                                                    */
/*                           Copyright (c) 2024 Terminus LLC                          */
/*                                                                                    */
/*                                All Rights Reserved.                                */
/*                                                                                    */
/*          Use of this source code is governed by LICENSE in the repo root.          */
/*                                                                                    */
/***************************# INTELLECTUAL PROPERTY RIGHTS ****************************/

/**
 * @file    CoordSys.rs
 * @author  Marvin Smith
 * @date    5/24/2024
 */

use crate::coord::coordinate_type::epsg_to_coordinate_type;
use super::coordinate_type::CoordinateType;

pub struct CoordSysDefinition {

    /// For much easier processing, it is recommended to stick with EPSG codes, as that is far less
    /// error-prone then using things like WKT strings.
    pub epsg_code:  Option<i32>,

    /// Underlying structural type used internally by the API
    pub coord_type: CoordinateType,
}

impl CoordSysDefinition {

    pub fn to_log_string(&self) -> String {
        let mut output = String::from( "CoordSysDefinition:\n" );

        // Log the EPSG Code
        output = output + &format!( "    EPSG Code: defined: {}", self.epsg_code.is_some() );
        if self.epsg_code.is_some() {
            output = output + &format!( ", value: {}", self.epsg_code.unwrap() );
        }
        output.push_str( "\n" );

        // Log the Coordinate Type
        output = output + &format!( "    CoordinateType: {}", self.coord_type );

        return output;
    }

    pub fn from_epsg( epsg_code: i32 ) -> Option<CoordSysDefinition> {

        // Check if EPSG code is set
        let coordinate_type: Option<CoordinateType> = epsg_to_coordinate_type(epsg_code);

        match coordinate_type {
            None => return None,
            _ => () // do nothing,
        }

        // Extract final coordinate type
        let res_ctype: CoordinateType = match coordinate_type {
            Some(coordinate_type) => coordinate_type,
            None => return None,
        };

        // Determine type from EPSG Code
        // This library uses a provided configuration-file, which
        // if not available, will fallback to a large table. 
        let cs_def_out: CoordSysDefinition = CoordSysDefinition {
            epsg_code: Some(epsg_code),
            coord_type: res_ctype,
        };

        Some(cs_def_out)
    }
}