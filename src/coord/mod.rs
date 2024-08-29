/**************************** INTELLECTUAL PROPERTY RIGHTS ****************************/
/*                                                                                    */
/*                           Copyright (c) 2024 Terminus LLC                          */
/*                                                                                    */
/*                                All Rights Reserved.                                */
/*                                                                                    */
/*          Use of this source code is governed by LICENSE in the repo root.          */
/*                                                                                    */
/***************************# INTELLECTUAL PROPERTY RIGHTS ****************************/

pub mod configuration;
pub mod coord_sys_base;
pub mod coord_sys_def;
pub mod coord_sys_factory;
pub mod coord_sys_geographic;

pub mod coord_sys_utm;

pub mod coordinate_type;

pub mod datum;

pub mod ellipsoid;

pub mod transformer;
mod error;
mod coord_sys_transverse_mercator;
mod trans_merc_coeffs;