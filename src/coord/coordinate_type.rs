/**************************** INTELLECTUAL PROPERTY RIGHTS ****************************/
/*                                                                                    */
/*                           Copyright (c) 2024 Terminus LLC                          */
/*                                                                                    */
/*                                All Rights Reserved.                                */
/*                                                                                    */
/*          Use of this source code is governed by LICENSE in the repo root.          */
/*                                                                                    */
/***************************# INTELLECTUAL PROPERTY RIGHTS ****************************/

use std::fmt;

pub enum CoordinateType {
    GEOGRAPHIC,
    UTM,
}

// This takes an EPSG code and creates a CoordinateType enum from it.
pub fn epsg_to_coordinate_type(ctype: i32) -> Option<CoordinateType> {

    // UTM (32600 series is WGS84)
    match ctype {

        // Standard WGS84 Geographic
        4326 => return Some(CoordinateType::GEOGRAPHIC),

        // UTM (WGS84 Northern Hemisphere)
        32600..=32700 => return Some(CoordinateType::UTM),

        // This is bad
        _ => return None,
    };
}

impl fmt::Display for CoordinateType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            CoordinateType::GEOGRAPHIC => write!(f, "GEOGRAPHIC"),
            CoordinateType::UTM => write!(f, "UTM"),
        }
    }
}