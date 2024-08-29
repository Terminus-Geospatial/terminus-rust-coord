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

#[derive(Debug)]
pub enum ErrorCode {
    NotFound,
    LatitudeOutOfRange,
    LongitudeOutOfRange,
    TransverseMercatorEastingOutOfRange,
    TransverseMercatorNorthingOutOfRange,
    UtmZoneOverrideFailure,
}

impl fmt::Display for ErrorCode {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ErrorCode::NotFound                             => write!(f, "NotFound"),
            ErrorCode::LatitudeOutOfRange                   => write!(f, "LatitudeOutOfRange"),
            ErrorCode::LongitudeOutOfRange                  => write!(f, "LongitudeOutOfRange"),
            ErrorCode::TransverseMercatorEastingOutOfRange  => write!(f, "TransverseMercatorEastingOutOfRange"),
            ErrorCode::TransverseMercatorNorthingOutOfRange => write!(f, "TransverseMercatorNorthingOutOfRange"),
            ErrorCode::UtmZoneOverrideFailure               => write!(f, "UtmZoneOverrideFailure" ),
        }
    }
}