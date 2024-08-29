/**************************** INTELLECTUAL PROPERTY RIGHTS ****************************/
/*                                                                                    */
/*                           Copyright (c) 2024 Terminus LLC                          */
/*                                                                                    */
/*                                All Rights Reserved.                                */
/*                                                                                    */
/*          Use of this source code is governed by LICENSE in the repo root.          */
/*                                                                                    */
/***************************# INTELLECTUAL PROPERTY RIGHTS ****************************/

use super::error;

pub trait CoordSysBase {

    fn to_log_string(&self) -> String;

    fn to_geographic(&self, input_coord: Vec<f64>) -> Result<Vec<f64>,error::ErrorCode>;

    fn from_geographic(&self, input_coord: Vec<f64>) -> Result<Vec<f64>,error::ErrorCode>;

    fn is_geographic(&self) -> bool;

    fn is_projected(&self) -> bool;
}